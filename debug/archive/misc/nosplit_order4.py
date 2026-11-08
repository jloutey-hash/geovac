"""nosplit step 3: the order-4 right factor of W, and whether the chain C lies on it.

Step 2 found the exact partial section S1 (R1 o S1 = 1 mod R3).  Then
  U = {y in Sol(W): (theta-2) y in S1(Sol R3)}  (dim 4, contains t^2)
is the unique 4-dim invariant subspace of Sol(W), and  P = (theta-2) - S1 o E  (E = R1 o (theta-2))
satisfies P y = alpha(y) h for y in Sol(W), h = t^3/(1-t^2)^2 the solution of R1, with U = ker P.
So C lies on the order-4 factor iff alpha(C) = 0.  alpha(C) = (P C)/h must be the same at every t
(consistency check).  Control: P applied with S1 perturbed must give a t-dependent ratio.
Also: R3' = ann(S1(Sol R3)) and the remaining extension R3' o (theta-2) -- searched with the step-2
bound machinery when its singular points are rational (cross-check of the full no-split)."""
from __future__ import annotations
import json, sys
from pathlib import Path
import mpmath as mp
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import nosplit_core as N
import nosplit_search as NS

x, t, u = N.x, N.t, sp.Symbol("u")
HERE = N.HERE


def S1_exact():
    S1 = [(12*x**3 + 229*x**2 - 31*x + 12)/(3*x*(x - 1)**2*(4*x + 11)),
          -2*(20*x**3 + 619*x**2 - 143*x + 80)/(15*x*(x - 1)**2*(4*x + 11)),
          -4*(x - 4)*(20*x**2 - 7*x + 5)/(15*x*(x - 1)**2*(4*x + 11))]
    return [sp.cancel(s) for s in S1]


def to_Q_theta(Pd):
    """delta-form (in x) -> polynomial Q-form sum_i t^{2i} Q_i(theta) after clearing a left
    denominator (returned)."""
    # delta^k = (theta/2)^k
    den = sp.Integer(1)
    for c in Pd:
        den = sp.lcm(den, sp.denom(sp.together(c)))
    Q = {}
    for k, c in enumerate(Pd):
        pc = sp.Poly(sp.expand(sp.cancel(c * den)), x)
        for (e,), v in zip(pc.monoms(), pc.coeffs()):
            Q[e] = Q.get(e, 0) + v * (u / 2) ** k
    return [sp.expand(Q.get(i, 0)) for i in range(max(Q) + 1)], den


def main():
    ops = N.chain_ops()
    R3, R1, T2 = ops["R3"], ops["R1"], ops["T2"]
    E = N.compose(R1, T2)
    S1 = S1_exact()
    assert NS.verify(R3, R1, S1)
    P = [sp.cancel(a - b) for a, b in zip(T2 + [0] * 3, N.compose(S1, E))]
    print("ord P", len(P) - 1, flush=True)
    # P o y for y in Sol(W): exact check that P maps Sol(W) into Sol(R1): R1 o P = 0 mod W
    _, remPW = N.rdiv(N.compose(R1, P), N.compose(R3, E))
    print("R1 o P mod W == 0:", all(sp.cancel(c) == 0 for c in remPW), flush=True)
    # P kills t^2 and S1-image: P o (theta-2)^{-1}... check P(t^2)=0 exactly
    Pt2 = sum(c * (sp.Integer(1) ** k) for k, c in enumerate(P))  # delta^k x = x
    print("P(x) = P(t^2):", sp.cancel(sum(c * x for c in P)), flush=True)
    import density_loop_B2_inhom as B2
    QP, den = to_Q_theta(P)
    S1bad = [S1[0] * sp.Rational(101, 100), S1[1], S1[2]]
    Pbad = [sp.cancel(a - b) for a, b in zip(T2 + [0] * 3, N.compose(S1bad, E))]
    QPb, denb = to_Q_theta(Pbad)
    rec = {"P_order": len(P) - 1, "R1P_mod_W_zero": all(sp.cancel(c) == 0 for c in remPW)}
    for dps in [int(a) for a in (sys.argv[1:] or ["20", "30"])]:
        mp.mp.dps = dps
        rows = {}
        for ts in ["0.5", "1.3", "2.5", "3.5"]:
            tv = mp.mpf(ts)
            dv = mp.mpf(sp.Rational(den.subs(x, sp.Rational(ts) ** 2)))
            dvb = mp.mpf(sp.Rational(denb.subs(x, sp.Rational(ts) ** 2)))
            h = tv**3 / (1 - tv**2) ** 2
            PC = B2.LC(QP, ts) / dv
            PCb = B2.LC(QPb, ts) / dvb
            C = B2.LC([sp.Integer(1)], ts)
            rows[ts] = {"C": mp.nstr(C, 15), "alpha": mp.nstr(PC / h, dps - 5), "alpha_control": mp.nstr(PCb / h, 12)}
            print(f"dps {dps} t {ts}: C {mp.nstr(C, 12)}  alpha = PC/h = {mp.nstr(PC / h, dps - 5)}"
                  f"   control (S1 coeff x1.01) {mp.nstr(PCb / h, 10)}", flush=True)
        rec[f"dps{dps}"] = rows
    # remaining extension: R3' = ann(S1(Sol R3)), W_U = R3' o (theta - 2)
    R3p = NS.annihilator(R3, S1)
    print("R3' =", R3p, flush=True)
    print("singular factors of R3':", N.singular_points(R3p), flush=True)
    rec["R3prime"] = [str(c) for c in R3p]
    try:
        res = NS.run("remaining extension: R3' o (theta-2)", R3p, T2, pads=(0, 2))
        rec["remaining"] = res
    except ValueError as ex:
        print("remaining-extension search skipped:", ex)
        rec["remaining"] = {"skipped": str(ex)}
    (HERE / "data" / "nosplit_order4.json").write_text(json.dumps(rec, indent=1, default=str))


if __name__ == "__main__":
    main()
