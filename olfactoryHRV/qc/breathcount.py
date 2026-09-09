import glob, os, h5py, numpy as np

pat = r"R:\Neurology\Zelano_Lab\Lab_Common\Dupi\*\preProc\*breathingPreProc.mat"
rows = []
for f in sorted(glob.glob(pat)):
    base = os.path.basename(f)
    if "separate" in base:
        continue
    rec = {"file": base, "nBreaths": None, "finished": False, "nSeg": None, "fs": None, "nSampMin": None}
    try:
        with h5py.File(f, "r") as h:
            top = [k for k in h.keys() if k not in ("#refs#", "#subsystem#")]
            g = h[top[0]] if len(top) == 1 else h[[t for t in top if t in ("chanDat", "outDat", "out")][0]]
            keys = set(g.keys())
            rec["finished"] = "moreThan1" in keys
            if "bmFeatures" in keys:
                bf = g["bmFeatures"]
                if "inhaleVolumes" in bf:
                    rec["nBreaths"] = int(np.prod(bf["inhaleVolumes"].shape))
                if "nBreathsSegmented" in bf:
                    rec["nSeg"] = int(np.array(bf["nBreathsSegmented"]).ravel()[0])
            if "bmObj" in keys:
                rec["nSampMin"] = list(g["bmObj"].shape)
            if "fs" in keys:
                rec["fs"] = float(np.array(g["fs"]).ravel()[0])
    except Exception as e:
        rec["file"] = base + "  ERROR: " + str(e)[:60]
    rows.append(rec)

print(f"{'file':<52} {'fin':<4} {'nBreath':>8} {'nSeg':>6} {'bmObj':>12}")
print("-" * 90)
for r in rows:
    print(f"{r['file']:<52} {str(r['finished']):<4} {str(r['nBreaths']):>8} {str(r['nSeg']):>6} {str(r['nSampMin']):>12}")
