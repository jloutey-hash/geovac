"""nosplit step 5: independent cross-check of the full no-split through the order-4 factor.
U = ker P = Sol(R3' o (theta - 2)) with R3' = ann(S1(Sol R3)).  W splits  <=>  U splits
(V subset U complementing t^2)  <=>  exists S with (theta - 2) S = 1 mod R3'.
Same proved bound + exhaustive exact search as step 2, different operator pair (R3' has
apparent singularities at x = 1 and at the roots of 20x^2 - 7x + 5).
Optional argv[1] = dps: also re-measure alpha(C) = (P C)/h at t = 1.3 at that precision."""
from __future__ import annotations
import json, sys
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import nosplit_core as N
import nosplit_search as NS
import nosplit_order4 as O4

x = N.x


def main():
    ops = N.chain_ops()
    R3, R1, T2 = ops["R3"], ops["R1"], ops["T2"]
    S1 = O4.S1_exact()
    R3p = NS.annihilator(R3, S1)
    rec = {"R3prime": [str(c) for c in R3p]}
    E = N.compose(R1, T2)
    P = [sp.cancel(a - b) for a, b in zip(T2 + [0] * 3, N.compose(S1, E))]
    # P and R3' o (theta-2) have the same kernel: check P = c(x) * R3' o (theta-2)
    Wr = N.compose(R3p, T2)
    ratio = [sp.cancel(a / b) for a, b in zip(P, Wr)]
    print("P / (R3' o (theta-2)) coefficientwise:", set(ratio), flush=True)
    rec["P_is_left_multiple"] = len(set(ratio)) == 1
    rec["search"] = NS.run("remaining extension R3' o (theta-2)", R3p, T2, pads=(0, 2))
    if len(sys.argv) > 1:
        import mpmath as mp
        import density_loop_B2_inhom as B2
        mp.mp.dps = int(sys.argv[1])
        QP, den = O4.to_Q_theta(P)
        for ts in ["1.3"]:
            tv = mp.mpf(ts)
            dv = mp.mpf(sp.Rational(den.subs(x, sp.Rational(ts) ** 2)))
            a = B2.LC(QP, ts) / dv / (tv**3 / (1 - tv**2) ** 2)
            print(f"dps {mp.mp.dps} t {ts}: alpha(C) = {mp.nstr(a, 8)}", flush=True)
            rec[f"alpha_dps{mp.mp.dps}_t{ts}"] = mp.nstr(a, 8)
    (N.HERE / "data" / "nosplit_remaining.json").write_text(json.dumps(rec, indent=1, default=str))


if __name__ == "__main__":
    main()
