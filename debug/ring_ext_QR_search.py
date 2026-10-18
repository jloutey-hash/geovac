"""QR step 3b: bounded numerical annihilator search, float64.
For f in {G (triangle), C3 (chain of three), C4 (open chain of four), S (ring of four)} and each
(order r, degree d in t^2): columns t^{2i} theta^j f, i <= d, j <= r.  Train on the 40 points
with t <= 1.6 (t in [0.25, 4], log-spaced), take the smallest right singular vector v, and test it
by EXTRAPOLATION on the points with t >= 2.4:
   res = max_k |sum_c v_c A_kc| / sum_c |v_c A_kc|.
A genuine operator extrapolates at the data's float64 level; a near-null vector of smooth columns
(the generic case once unknowns exceed what 16 digits can separate) does not.
Positive controls: G's known (theta-2) L_G (r = 4, d = 2, 15 unknowns), C3's W (r = 5, d = 3, 24)."""
import json
from pathlib import Path
import numpy as np

D = np.load(Path(__file__).parent / "data" / "ring_ext_QR_data.npz")
ts = D["ts"]
tr = np.where(ts <= 1.6)[0]; te = np.where(ts >= 2.4)[0]


def search(V, r, d):
    cols = [(i, j) for i in range(d + 1) for j in range(r + 1)]
    A = np.array([[ts[k]**(2 * i) * V[k, j] for (i, j) in cols] for k in range(len(ts))])
    A = A / np.sqrt(np.mean(A[tr]**2, axis=0))
    _, s, Vt = np.linalg.svd(A[tr], full_matrices=False)
    v = Vt[-1]
    res = np.max(np.abs(A[te] @ v) / (np.abs(A[te]) @ np.abs(v)))
    return float(s[-1] / s[0]), float(s[-1] / s[-2]), float(res)


def main():
    rec = {}
    for name in ["G", "C3", "C4", "S"]:
        rows = []
        for r in range(1, 11):
            for d in range(0, 8):
                n = (r + 1) * (d + 1)
                if n > min(44, len(tr) - 2):
                    continue
                rows.append([r, d, n, *search(D[name], r, d)])
        rec[name] = rows
        print(name, "best five (r, d, unknowns, smin/smax, smin/snext, extrapolation residual):", flush=True)
        for b in sorted(rows, key=lambda x: x[5])[:5]:
            print(f"    r={b[0]} d={b[1]} n={b[2]}  {b[3]:.1e} {b[4]:.1e} {b[5]:.1e}", flush=True)
    (Path(__file__).parent / "data" / "ring_ext_QR_search.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
