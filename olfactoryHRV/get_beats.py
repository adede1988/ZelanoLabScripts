import glob, os, h5py, numpy as np

ROOT = os.environ.get("OHRV_DATA", r"R:\Neurology\Zelano_Lab\Lab_Common\Dupi")
OUT  = os.path.join(os.path.dirname(os.path.abspath(__file__)), "work")

files = sorted(glob.glob(os.path.join(ROOT, "*", "preProc", "*breathingPreProc.mat")))
files = [f for f in files if "separate" not in os.path.basename(f)]

print(f"{'session':<26} {'nBeats':>7} {'dur_min':>8} {'medHR':>7} {'pctPhys':>8}")
print("-" * 62)
for f in files:
    sess = os.path.basename(f).replace("_breathingPreProc.mat", "")
    try:
        with h5py.File(f, "r") as h:
            tops = [k for k in h.keys() if k not in ("#refs#", "#subsystem#")]
            name = next((t for t in ("chanDat", "outDat", "out") if t in tops), tops[0])
            g = h[name]
            if "heartBeats" not in g or "moreThan1" not in g:
                print(f"{sess:<26}   SKIP (unfinished or no heartBeats)")
                continue
            hb = np.array(g["heartBeats"]).ravel().astype(float)
            fs = float(np.array(g["fs"]).ravel()[0])
            nsamp = g["data"].shape[0]
    except Exception as e:
        print(f"{sess:<26}   ERROR {str(e)[:44]}")
        continue

    hb = np.sort(hb[np.isfinite(hb)])
    nn = np.diff(hb) / fs                      # NN intervals, seconds
    phys = (nn > 0.3) & (nn < 2.0)             # physiologically plausible
    np.savez(os.path.join(OUT, sess + "_beats.npz"), hb=hb, fs=fs, nsamp=nsamp)
    print(f"{sess:<26} {len(hb):>7} {nsamp/fs/60:>8.1f} "
          f"{60/np.median(nn[phys]) if phys.any() else float('nan'):>7.1f} "
          f"{100*phys.mean():>7.1f}%")
print("BEATS DONE")
