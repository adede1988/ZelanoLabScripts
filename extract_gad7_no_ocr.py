# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 11:22:23 2026

@author: dtf8829
"""

import argparse, json, zipfile, shutil
from pathlib import Path

import numpy as np
import pandas as pd
import cv2
import fitz  # PyMuPDF


ITEM_COLS = [
    "Feeling nervous, anxious, or on edge",
    "Not being able to stop or control worrying",
    "Worrying too much about different things",
    "Trouble relaxing",
    "Being so restless that it's hard to sit still",
    "Becoming easily annoyed or irritable",
    "Feeling afraid as if something awful might happen",
]


def render_page_bgr(doc, pidx, zoom=1.0):
    page = doc[pidx]
    pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)
    img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
    if img.shape[2] == 4:
        img = img[:, :, :3]
    return cv2.cvtColor(img, cv2.COLOR_RGB2BGR)


def rotate_bgr(img, rot):
    if rot == 0:
        return img
    if rot == 90:
        return cv2.rotate(img, cv2.ROTATE_90_CLOCKWISE)
    if rot == 180:
        return cv2.rotate(img, cv2.ROTATE_180)
    if rot == 270:
        return cv2.rotate(img, cv2.ROTATE_90_COUNTERCLOCKWISE)
    raise ValueError(rot)


def binarize_adaptive(gray):
    gray = cv2.equalizeHist(gray)
    return cv2.adaptiveThreshold(
        gray, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY_INV, 41, 13
    )


def detect_horiz_lines(bin_img, min_frac=0.35, kernel_frac=0.14, close_frac=0.05):
    h, w = bin_img.shape
    kernel_w = max(50, int(w * kernel_frac))
    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (kernel_w, 1))
    horiz = cv2.morphologyEx(bin_img, cv2.MORPH_OPEN, kernel, iterations=1)

    close_k = cv2.getStructuringElement(cv2.MORPH_RECT, (max(10, int(w * close_frac)), 1))
    horiz = cv2.morphologyEx(horiz, cv2.MORPH_CLOSE, close_k, iterations=2)

    cnts, _ = cv2.findContours(horiz, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    lines = []
    for c in cnts:
        x, y, ww, hh = cv2.boundingRect(c)
        if ww < min_frac * w:
            continue
        if hh > 0.05 * h:
            continue
        yc = y + hh / 2
        if yc < 0.06 * h or yc > 0.97 * h:
            continue
        lines.append((x, x + ww, yc, ww))

    lines.sort(key=lambda t: t[2])

    # merge near-duplicate y lines
    merged = []
    ytol = max(5, int(0.008 * h))
    for x0, x1, yc, ww in lines:
        if not merged or abs(yc - merged[-1][2]) > ytol:
            merged.append([x0, x1, yc, ww])
        else:
            m = merged[-1]
            m[0] = min(m[0], x0)
            m[1] = max(m[1], x1)
            m[3] = m[1] - m[0]
            m[2] = (m[2] + yc) / 2

    return merged


def best_group(lines, shape):
    h, w = shape
    best = None
    for n in (8, 7):
        if len(lines) < n:
            continue
        for i in range(len(lines) - n + 1):
            grp = lines[i : i + n]
            ys = np.array([g[2] for g in grp], float)
            gaps = np.diff(ys)
            if gaps.size == 0:
                continue
            mg = float(gaps.mean())
            sg = float(gaps.std())
            span = float(np.mean([g[1] - g[0] for g in grp]))

            if mg < 0.02 * h or mg > 0.30 * h:
                continue
            if span < 0.60 * w:
                continue

            score = span - 2500 * sg
            if best is None or score > best["score"]:
                best = {"n": n, "grp": grp, "mg": mg, "ys": sorted(ys.tolist()), "score": score}
    return best


def col_peak_score(bin_img, bg):
    h, w = bin_img.shape
    ys = bg["ys"]
    x0 = int(min(g[0] for g in bg["grp"]))
    x1 = int(max(g[1] for g in bg["grp"]))
    y_top = int(max(0, ys[0]))
    y_bot = int(min(h, ys[-1]))

    rx0 = int(max(0.55 * w, x0 + 0.55 * (x1 - x0)))
    rx1 = int(min(w, x1))
    if rx1 - rx0 < 40 or y_bot - y_top < 40:
        return 0.0

    roi = bin_img[y_top:y_bot, rx0:rx1]
    xsum = roi.sum(axis=0).astype(np.float32)
    if xsum.size < 10:
        return 0.0

    xsum_s = cv2.GaussianBlur(xsum.reshape(1, -1), (1, 81), 0).ravel()
    peaks = []
    for i in range(1, len(xsum_s) - 1):
        if xsum_s[i] > xsum_s[i - 1] and xsum_s[i] > xsum_s[i + 1]:
            peaks.append(xsum_s[i])

    if len(peaks) < 4:
        return 0.0
    peaks.sort(reverse=True)
    base = float(np.median(xsum_s) + 1e-6)
    return float(np.mean(peaks[:4]) / base)


def edge_from_bin(bin_img):
    return cv2.Canny(bin_img, 50, 150)


def corrcoef(a, b):
    a = a.astype(np.float32).ravel()
    b = b.astype(np.float32).ravel()
    a -= a.mean()
    b -= b.mean()
    denom = (np.linalg.norm(a) * np.linalg.norm(b)) + 1e-9
    return float((a @ b) / denom)


def make_template_edge(pdf_path):
    doc = fitz.open(pdf_path)
    best_score = None
    best_edge = None

    for pidx in range(doc.page_count):
        bgr0 = render_page_bgr(doc, pidx, zoom=1.1)
        for rot in (0, 90, 180, 270):
            bgr = rotate_bgr(bgr0, rot)
            gray = cv2.cvtColor(bgr, cv2.COLOR_BGR2GRAY)
            bin_img = binarize_adaptive(gray)
            bg = best_group(detect_horiz_lines(bin_img), bin_img.shape)
            if bg is None:
                continue
            pk = col_peak_score(bin_img, bg)
            score = bg["score"] + pk * 1000

            if best_score is None or score > best_score:
                ys = bg["ys"]
                mg = bg["mg"]
                x0 = int(min(g[0] for g in bg["grp"]))
                x1 = int(max(g[1] for g in bg["grp"]))
                h, w = gray.shape
                hy0 = int(max(0, ys[0] - 1.3 * mg))
                hy1 = int(max(0, ys[0] - 0.1 * mg))
                hx0 = int(max(0, x0))
                hx1 = int(min(w, x0 + 0.92 * (x1 - x0)))
                if hy1 > hy0 and hx1 > hx0:
                    crop = gray[hy0:hy1, hx0:hx1]
                    best_edge = edge_from_bin(binarize_adaptive(crop))
                    best_score = score

    doc.close()
    if best_edge is None:
        return None

    target_w = 520
    scale = target_w / best_edge.shape[1]
    target_h = max(40, int(best_edge.shape[0] * scale))
    return cv2.resize(best_edge, (target_w, target_h), interpolation=cv2.INTER_AREA)


def header_match(gray, bg, template_edge):
    if template_edge is None:
        return 0.0

    ys = bg["ys"]
    mg = bg["mg"]
    x0 = int(min(g[0] for g in bg["grp"]))
    x1 = int(max(g[1] for g in bg["grp"]))
    h, w = gray.shape

    hy0 = int(max(0, ys[0] - 1.3 * mg))
    hy1 = int(max(0, ys[0] - 0.1 * mg))
    hx0 = int(max(0, x0))
    hx1 = int(min(w, x0 + 0.92 * (x1 - x0)))
    if hy1 <= hy0 or hx1 <= hx0:
        return 0.0

    crop = gray[hy0:hy1, hx0:hx1]
    ce = edge_from_bin(binarize_adaptive(crop))
    ce = cv2.resize(ce, (template_edge.shape[1], template_edge.shape[0]), interpolation=cv2.INTER_AREA)
    return corrcoef(template_edge, ce)


def derive_rows(bg):
    ys = bg["ys"]
    mg = bg["mg"]
    if bg["n"] == 8:
        centers = [0.5 * (ys[i] + ys[i + 1]) for i in range(7)]
        y_top, y_bot = ys[0], ys[-1]
    else:
        centers = [ys[i] - 0.5 * mg for i in range(7)]
        y_top, y_bot = ys[0] - mg, ys[-1] + 0.2 * mg
    return centers, y_top, y_bot, mg


def find_cols(bin_img, x0, x1, y_top, y_bot):
    h, w = bin_img.shape
    rx0 = int(max(0.50 * w, x0 + 0.50 * (x1 - x0)))
    rx1 = int(min(w, x1))
    ry0 = int(max(0, y_top))
    ry1 = int(min(h, y_bot))

    if rx1 - rx0 < 40 or ry1 - ry0 < 40:
        return [int(rx0 + f * (rx1 - rx0)) for f in (0.18, 0.44, 0.70, 0.90)]

    roi = bin_img[ry0:ry1, rx0:rx1]
    xsum = roi.sum(axis=0).astype(np.float32)
    xsum_s = cv2.GaussianBlur(xsum.reshape(1, -1), (1, 101), 0).ravel()

    peaks = []
    for i in range(1, len(xsum_s) - 1):
        if xsum_s[i] > xsum_s[i - 1] and xsum_s[i] > xsum_s[i + 1]:
            peaks.append((xsum_s[i], i))

    peaks.sort(reverse=True)
    chosen = []
    min_dist = max(12, int(0.05 * w))
    for val, idx in peaks:
        if all(abs(idx - c) >= min_dist for c in chosen):
            chosen.append(idx)
        if len(chosen) == 4:
            break

    if len(chosen) < 4:
        chosen = list(np.linspace(int(0.18 * (rx1 - rx0)), int(0.88 * (rx1 - rx0)), 4).astype(int))

    chosen = sorted(chosen)
    return [int(rx0 + c) for c in chosen]


def make_boxes(centers, cols, mg, shape):
    h, w = shape
    side = int(max(18, min(0.70 * mg, 0.08 * w)))
    half = side // 2

    boxes = []
    for y in centers:
        cy = int(round(y))
        row = []
        for cx in cols:
            cx = int(round(cx))
            x0 = max(0, cx - half)
            x1 = min(w, cx + half)
            y0 = max(0, cy - half)
            y1 = min(h, cy + half)
            row.append((x0, y0, x1, y1))
        boxes.append(row)
    return boxes


def masks_from_page(gray, boxes):
    masks = []
    for c in range(4):
        patches = []
        for r in range(7):
            x0, y0, x1, y1 = boxes[r][c]
            patch = cv2.resize(gray[y0:y1, x0:x1], (64, 64), interpolation=cv2.INTER_AREA)
            patches.append(patch)

        med = np.median(np.stack(patches, 0), axis=0).astype(np.uint8)
        b = binarize_adaptive(med)

        cnts, _ = cv2.findContours(b, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        mask = np.zeros_like(b)

        if cnts:
            cx0, cy0 = 32, 32
            best = None
            for ct in cnts:
                x, y, ww, hh = cv2.boundingRect(ct)
                cx = x + ww / 2
                cy = y + hh / 2
                d = (cx - cx0) ** 2 + (cy - cy0) ** 2
                if best is None or d < best[0]:
                    best = (d, ct)

            cv2.drawContours(mask, [best[1]], -1, 255, -1)
            mask = cv2.dilate(mask, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (9, 9)), 1)

        masks.append(mask)
    return masks


def parse_gad_page(bgr, qc_path=None):
    gray = cv2.cvtColor(bgr, cv2.COLOR_BGR2GRAY)
    bin_img = binarize_adaptive(gray)

    bg = best_group(detect_horiz_lines(bin_img), bin_img.shape)
    if bg is None:
        return None

    centers, y_top, y_bot, mg = derive_rows(bg)
    x0 = int(min(g[0] for g in bg["grp"]))
    x1 = int(max(g[1] for g in bg["grp"]))

    cols = find_cols(bin_img, x0, x1, y_top, y_bot)
    boxes = make_boxes(centers, cols, mg, gray.shape)
    masks = masks_from_page(gray, boxes)

    picks = []
    row_scores = []
    for r in range(7):
        sc = []
        for c in range(4):
            bx0, by0, bx1, by1 = boxes[r][c]
            patch = cv2.resize(gray[by0:by1, bx0:bx1], (64, 64), interpolation=cv2.INTER_AREA)
            b = binarize_adaptive(patch)
            resid = cv2.bitwise_and(b, cv2.bitwise_not(masks[c]))
            sc.append(int(np.sum(resid > 0)))
        row_scores.append(sc)
        picks.append(int(np.argmax(sc)))

    if qc_path is not None:
        dbg = bgr.copy()
        # draw detected lines (blue)
        for g in bg["grp"]:
            cv2.line(dbg, (x0, int(round(g[2]))), (x1, int(round(g[2]))), (255, 0, 0), 2)
        # draw boxes (green) + chosen (red)
        for r in range(7):
            for c in range(4):
                bx0, by0, bx1, by1 = boxes[r][c]
                cv2.rectangle(dbg, (bx0, by0), (bx1, by1), (0, 255, 0), 2)
            csel = picks[r]
            bx0, by0, bx1, by1 = boxes[r][csel]
            cv2.rectangle(dbg, (bx0, by0), (bx1, by1), (0, 0, 255), 3)
            cv2.putText(dbg, str(csel), (bx0, max(0, by0 - 6)), cv2.FONT_HERSHEY_SIMPLEX, 0.9, (0, 0, 255), 2)

        cv2.imwrite(str(qc_path), dbg)

    return {"picks": picks, "scores": row_scores, "cols": cols, "line_score": float(bg["score"])}


def process_pdf(pdf_path, template_edge, qc_dir):
    doc = fitz.open(pdf_path)
    best = None

    for pidx in range(doc.page_count):
        bgr0 = render_page_bgr(doc, pidx, zoom=0.95)
        h, w = bgr0.shape[:2]
        if w > 850:
            s = 850 / w
            bgr0 = cv2.resize(bgr0, (int(w * s), int(h * s)), interpolation=cv2.INTER_AREA)

        for rot in (0, 90, 180, 270):
            bgr = rotate_bgr(bgr0, rot)
            gray = cv2.cvtColor(bgr, cv2.COLOR_BGR2GRAY)
            bin_img = binarize_adaptive(gray)
            bg = best_group(detect_horiz_lines(bin_img), bin_img.shape)
            if bg is None:
                continue

            pk = col_peak_score(bin_img, bg)
            hm = header_match(gray, bg, template_edge)
            comb = pk * 6000 + bg["score"] + hm * 800

            if best is None or comb > best["comb"]:
                best = {"comb": comb, "pidx": pidx, "rot": rot, "pk": pk, "hm": hm}

    doc.close()
    if best is None:
        return None, {"pdf": str(pdf_path), "status": "no_table_found"}

    # high-res parse
    doc2 = fitz.open(pdf_path)
    bgr_hi0 = render_page_bgr(doc2, best["pidx"], zoom=1.8)
    doc2.close()
    bgr_hi = rotate_bgr(bgr_hi0, best["rot"])

    qc_name = f"{pdf_path.parent.parent.name}__{pdf_path.parent.name}__{pdf_path.stem}__p{best['pidx']+1}_rot{best['rot']}.png"
    qc_path = qc_dir / qc_name

    parsed = parse_gad_page(bgr_hi, qc_path=qc_path)
    if parsed is None:
        return None, {"pdf": str(pdf_path), "status": "parse_failed", "page": best["pidx"] + 1, "rot": best["rot"]}

    rec = {
        "study": pdf_path.parent.parent.name,
        "participant": pdf_path.parent.name,
        "source_pdf": pdf_path.name,
        "pdf": str(pdf_path),
        "gad7_page_in_pdf_1based": best["pidx"] + 1,
        "rotation": best["rot"],
        "selection_col_peak_score": float(best["pk"]),
        "selection_header_match": float(best["hm"]),
        "col_centers_json": json.dumps(parsed["cols"]),
        "row_scores_json": json.dumps(parsed["scores"]),
        "qc_image": qc_name,
    }
    for i in range(7):
        rec[f"item{i+1}"] = int(parsed["picks"][i])

    return rec, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zip", required=True, help="Path to questionnaireFiles.zip")
    ap.add_argument("--outdir", required=True, help="Output directory")
    args = ap.parse_args()

    zip_path = Path(args.zip)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    workdir = outdir / "_work"
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(workdir)

    root = workdir / "questionnaireFiles"
    pdfs = sorted(root.rglob("*.pdf"))
    if not pdfs:
        raise RuntimeError("No PDFs found under questionnaireFiles/")

    qc_dir = outdir / "gad7_qc_images"
    qc_dir.mkdir(parents=True, exist_ok=True)

    # Build a template edge from a likely-good file (first PDF)
    template_edge = make_template_edge(pdfs[0])

    records = []
    failures = []
    for pdf in pdfs:
        rec, fail = process_pdf(pdf, template_edge, qc_dir)
        if rec is not None:
            records.append(rec)
        if fail is not None:
            failures.append(fail)

    long_df = pd.DataFrame(records).sort_values(["study", "participant", "source_pdf"]).reset_index(drop=True)

    summary_rows = []
    for _, r in long_df.iterrows():
        row = {
            "study": r["study"],
            "participant": r["participant"],
            "source_pdf": r["source_pdf"],
            "gad7_page_in_pdf_1based": int(r["gad7_page_in_pdf_1based"]),
            "rotation": int(r["rotation"]),
        }
        for j, col in enumerate(ITEM_COLS, start=1):
            row[col] = int(r[f"item{j}"])
        summary_rows.append(row)
    summary_df = pd.DataFrame(summary_rows).sort_values(["study", "participant", "source_pdf"]).reset_index(drop=True)

    fail_df = pd.DataFrame(failures)

    summary_df.to_csv(outdir / "gad7_responses.csv", index=False)
    long_df.to_csv(outdir / "gad7_extractions_long.csv", index=False)
    fail_df.to_csv(outdir / "gad7_failures.csv", index=False)

    # zip QC images
    qc_zip = outdir / "gad7_qc_images.zip"
    if qc_zip.exists():
        qc_zip.unlink()
    with zipfile.ZipFile(qc_zip, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(qc_dir.rglob("*.png")):
            z.write(p, arcname=p.name)

    print("Wrote:", outdir / "gad7_responses.csv")
    print("Wrote:", outdir / "gad7_extractions_long.csv")
    print("Wrote:", outdir / "gad7_failures.csv")
    print("Wrote:", qc_zip)


if __name__ == "__main__":
    main()
