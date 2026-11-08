"""Post-process saved T' matrices: every 1-norm leg, local slopes, Z sweep."""
import glob, math, re, sys
import numpy as np


def load(lmax=3, c=3, npts=24000):
    out = {}
    for f in sorted(glob.glob(f"debug/data/p60_Tp_l{lmax}_n*_c{c:g}_N{npts}.npz")):
        d = np.load(f)
        n = int(d["nmax"])
        out[n] = dict(Tp=d["Tp"], Rnu=d["Rnu"], lab=d["lab"], rmax=float(d["rmax"]))
    return out


def legs(rec, Z=2.0):
    Tp, Rnu = rec["Tp"], rec["Rnu"]
    K = len(Rnu)
    dTp = np.abs(np.diag(Tp))
    off = float(np.abs(Tp).sum() - dTp.sum())
    Mdiag = float(np.abs(Z * Rnu + np.diag(Tp)).sum())
    return dict(K=K,
                M_total=Mdiag + off,
                M_diag=Mdiag,
                M_off=off,
                Tp_full=float(np.abs(Tp).sum()),
                Tp_diag=float(dTp.sum()),
                T0=float(Z * Rnu.sum()))


ORDER = ["M_total", "M_diag", "M_off", "Tp_full", "Tp_diag", "T0"]


def show(data, Z=2.0, label="", kwin=None):
    ns = sorted(data)
    rows = [legs(data[n], Z) for n in ns]
    if kwin:
        keep = [i for i, r in enumerate(rows) if kwin[0] <= r["K"] <= kwin[1]]
        ns = [ns[i] for i in keep]; rows = [rows[i] for i in keep]
    K = np.array([r["K"] for r in rows], float)
    print(f"\n=== {label}  Z={Z}  K={int(K[0])}..{int(K[-1])} ===")
    print(f"{'nmax':>4} {'K':>5} " + " ".join(f"{k:>10}" for k in ORDER)
          + "  ||  local slopes: " + " ".join(f"{k:>8}" for k in ORDER))
    for i, r in enumerate(rows):
        vals = " ".join(f"{r[k]:>10.4f}" for k in ORDER)
        if i == 0:
            sl = " ".join(f"{'':>8}" for _ in ORDER)
        else:
            sl = " ".join(
                f"{(math.log(r[k])-math.log(rows[i-1][k]))/(math.log(K[i])-math.log(K[i-1])):>8.4f}"
                for k in ORDER)
        print(f"{ns[i]:>4} {int(K[i]):>5} {vals}  ||  {' '*17}{sl}")
    g = {k: float(np.polyfit(np.log(K), np.log([r[k] for r in rows]), 1)[0]) for k in ORDER}
    print("global fit: " + "  ".join(f"{k}={g[k]:.4f}" for k in ORDER))
    return rows, g


if __name__ == "__main__":
    lmax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    data = load(lmax=lmax)
    print("loaded nmax:", sorted(data), "boxes:", {n: data[n]["rmax"] for n in sorted(data)})
    show(data, 2.0, f"CONVERGED (c=3, graded N=24000), lmax={lmax}, full range")
    show(data, 2.0, f"CONVERGED, lmax={lmax}, PAPER WINDOW", kwin=(74, 164))
