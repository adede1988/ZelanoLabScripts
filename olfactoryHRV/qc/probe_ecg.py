import glob, os, h5py, numpy as np

ROOT = os.environ.get("OHRV_DATA", r"R:\Neurology\Zelano_Lab\Lab_Common\Dupi")

def top_group(h):
    tops = [k for k in h.keys() if k not in ("#refs#", "#subsystem#")]
    for want in ("chanDat", "outDat", "out"):
        if want in tops:
            return want, h[want]
    return tops[0], h[tops[0]]

def read_labels(h, g):
    if "labels" not in g:
        return []
    out = []
    ds = g["labels"]
    try:
        refs = np.array(ds).ravel()
        for r in refs:
            try:
                arr = np.array(h[r]).ravel()
                out.append("".join(chr(int(c)) for c in arr if int(c) > 0))
            except Exception:
                pass
    except Exception:
        pass
    return out

pats = [("breathing", "*breathingPreProc.mat"),
        ("cue",       "*cueTaskPreProc.mat"),
        ("thresh",    "*PEA_threshold_preproc.mat"),
        ("O15",       "*O15preproc.mat")]

print(f"{'session':<26} {'task':<10} {'fin':<4} {'nECG':>5} {'RRint':>6} {'hBeats':>7} {'nSamp':>9}  ecg labels")
print("-" * 118)
rows = []
for task, pat in pats:
    for f in sorted(glob.glob(os.path.join(ROOT, "*", "preProc", pat))):
        base = os.path.basename(f)
        sess = base.split("_preproc")[0].split("_PEA")[0]
        for suf in ("_breathingPreProc.mat", "_cueTaskPreProc.mat", "_O15preproc.mat"):
            sess = sess.replace(suf, "")
        sess = sess.replace(".mat", "")
        if "separate" in base:
            continue
        try:
            with h5py.File(f, "r") as h:
                name, g = top_group(h)
                fin = "moreThan1" in g
                labs = read_labels(h, g)
                ecg = [l for l in labs if "ecg" in l.lower()]
                rr  = [l for l in labs if "rrint" in l.lower()]
                hb  = "heartBeats" in g
                nsamp = g["data"].shape[0] if "data" in g else 0
                print(f"{sess:<26} {task:<10} {str(fin):<4} {len(ecg):>5} {str(len(rr)>0):>6} {str(hb):>7} {nsamp:>9}  {','.join(ecg[:4])}")
                rows.append((sess, task, fin, len(ecg), len(rr) > 0, hb))
        except Exception as e:
            print(f"{sess:<26} {task:<10} ERROR {str(e)[:50]}")

print()
for task, _ in pats:
    sub = [r for r in rows if r[1] == task and r[2]]
    print(f"{task:<10} finished={len(sub):>3}  with ECG={sum(1 for r in sub if r[3] > 0):>3}  "
          f"with RRint={sum(1 for r in sub if r[4]):>3}  with heartBeats={sum(1 for r in sub if r[5]):>3}")
