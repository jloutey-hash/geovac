"""Local-slope + global-fit table for any p60 ladder JSON."""
import json, math, sys
import numpy as np

LEGS = [("M_total", "p_total"), ("M_diag", "p_Mdiag"), ("M_off", "p_off"),
        ("Tp_full", "p_Tp_full"), ("Tp_diag", "p_Tpdiag"), ("T0", "p_T0")]


def table(rows, kmin=None, kmax=None, label=""):
    rows = sorted(rows, key=lambda r: r["K"])
    if kmin is not None:
        rows = [r for r in rows if kmin <= r["K"] <= kmax]
    K = np.array([r["K"] for r in rows], float)
    print(f"\n=== {label}  (K = {int(K[0])}..{int(K[-1])}, {len(K)} points) ===")
    hdr = f"{'K':>5} " + " ".join(f"{n:>10}" for _, n in LEGS) + "   | values"
    print(hdr)
    for i in range(len(rows)):
        cells = []
        for key, _ in LEGS:
            if i == 0:
                cells.append(f"{'':>10}")
            else:
                s = (math.log(rows[i][key]) - math.log(rows[i - 1][key])) / \
                    (math.log(K[i]) - math.log(K[i - 1]))
                cells.append(f"{s:>10.4f}")
        vals = " ".join(f"{rows[i][k]:9.4f}" for k, _ in LEGS)
        print(f"{int(K[i]):>5} " + " ".join(cells) + f"   | {vals}")
    gl = []
    for key, name in LEGS:
        y = np.array([r[key] for r in rows], float)
        gl.append(f"{name}={float(np.polyfit(np.log(K), np.log(y), 1)[0]):.4f}")
    print("global: " + "  ".join(gl))


if __name__ == "__main__":
    rows = json.load(open(sys.argv[1]))
    lab = sys.argv[1]
    table(rows, label=lab)
    if len(rows) > 4:
        table(rows, 74, 164, lab + " [paper window]")
