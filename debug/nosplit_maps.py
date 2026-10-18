"""nosplit step 6: explicit maps for the memo.
(a) W = Q1 o P exactly (P = (theta-2) - S1 o E of order 4, Q1 of order 1).
(b) B' = (S1 o A) mod L3 in the t-variable: (theta - 2) of any element of U is t^2 + B' f for
    f in Sol(L3); so every element of U is  c t^2 + t^2 int^t s^-3 (B' f)(s) ds."""
from __future__ import annotations
import json, sys
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import nosplit_core as N
import nosplit_order4 as O4
import density_loop_D10_k3_exact_check as D10

x, t = N.x, N.t


def main():
    ops = N.chain_ops()
    R3, R1, T2 = ops["R3"], ops["R1"], ops["T2"]
    E = N.compose(R1, T2)
    W = N.compose(R3, E)
    S1 = O4.S1_exact()
    P = [sp.cancel(a - b) for a, b in zip(T2 + [0] * 3, N.compose(S1, E))]
    Q1, rem = N.rdiv(W, P)
    print("W = Q1 o P, Q1 =", [sp.factor(c) for c in Q1], " remainder zero:", all(c == 0 for c in rem), flush=True)
    Pm = [sp.factor(sp.cancel(c / P[-1])) for c in P]
    print("P (monic, delta-form in x):", Pm, flush=True)
    # (b) t-variable, theta-form
    S1t = [sp.cancel(s.subs(x, t**2) / 2**k) for k, s in enumerate(S1)]
    L3t = [-t**2, -3 * t**2, -3 * t**2, 4 - t**2]
    Bp = D10.rdiv(D10.compose(S1t, D10.A), L3t)[1]
    Bp = [sp.factor(c) for c in Bp]
    print("B' = (S1 o A) mod L3 (theta-form, t):", Bp, flush=True)
    rec = {"Q1": [str(sp.factor(c)) for c in Q1], "W_eq_Q1oP": all(c == 0 for c in rem),
           "P_monic_delta_x": [str(c) for c in Pm], "Bprime_theta_t": [str(c) for c in Bp]}
    (N.HERE / "data" / "nosplit_maps.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
