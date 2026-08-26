"""Rung 3c -- BROADER guarded closure probe of the integrated collinear T2.

Extends rung3b (which tested only the central tau=i single-fibre ring + classical
polylogs, all negative) with the bases that search actually skipped:
  - CROSS-fibre CM Gamma-values (disc -4 AND disc -8 periods) -- because the
    integrated V spans the whole Legendre family, not one fibre;
  - the lemniscate constant and Gamma(1/8)Gamma(3/8) products;
  - weight-2 CROSS products K(m1)*K(m2), K*Catalan, Catalan^2.

Discipline (curve-fit audit, CLAUDE.md [[feedback_audit_numerical_claims]]):
  1. a CLOSURE requires rel[0] (V-coeff) != 0 AND small height;
  2. rel[0]==0 is a basis-internal identity, NOT a closure (the rung3b trap);
  3. every basis is ALSO run against a DECOY constant (random, same magnitude) --
     if the basis fits the decoy at the same height/precision, the basis is too
     rich to trust and any V-hit is meaningless.

V is pasted from the high-precision fast_evaluator run (dps>=40).
"""
from __future__ import annotations
import mpmath as mp

mp.mp.dps = 30
TOL = mp.mpf(10) ** -14   # matched to the ~17-digit value; higher-dps V => tighten

# --- integrated collinear T2 (fast_evaluator, ~17 stable digits; dps>=40 run pending) ---
V = mp.mpf('0.39535576590171392')

pi = mp.pi
G = mp.catalan
ln2 = mp.log(2)

# disc -4 fibre (tau=i, rho=1/2): K(1/2) = Gamma(1/4)^2/(4 sqrt pi)
K1 = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(pi))
# disc -8 fibre (tau=i sqrt2): period (1+sqrt2)^(1/2) Gamma(1/8)Gamma(3/8)/(2^(13/4) sqrt pi)
K2 = (1 + mp.sqrt(2)) ** mp.mpf('0.5') * mp.gamma(mp.mpf(1) / 8) * mp.gamma(mp.mpf(3) / 8) \
     / (2 ** mp.mpf('3.25') * mp.sqrt(pi))
# lemniscate constant varpi_L = Gamma(1/4)^2 / (2 sqrt(2 pi))
lem = mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * pi))


def guarded(name, target, basis, maxcoeff):
    rel = mp.pslq([target] + basis, tol=TOL, maxcoeff=maxcoeff, maxsteps=10 ** 6)
    if rel is None:
        return "(none)"
    if rel[0] == 0:
        return f"{rel}  (V-coeff 0 => basis-internal identity, NOT a closure)"
    if max(abs(c) for c in rel) <= 60:
        return f"{rel}  <<< SMALL, V-coeff nonzero => CANDIDATE -- AUDIT REQUIRED"
    return f"{rel}  (high-height => not a closure)"


BASES = {
    "cross-fibre CM K w<=1":
        [mp.mpf(1), pi, K1, K2],
    "cross-fibre CM K w<=2":
        [mp.mpf(1), pi, K1, K2, pi ** 2, K1 ** 2, K2 ** 2, K1 * K2, K1 * pi, K2 * pi],
    "K + Catalan w<=2":
        [mp.mpf(1), pi, K1, G, pi ** 2, K1 ** 2, G ** 2, K1 * G, K1 * pi, G * pi],
    "lemniscate + ln2 w<=2":
        [mp.mpf(1), pi, lem, ln2, pi ** 2, lem ** 2, ln2 ** 2, lem * ln2, lem * pi],
}

# a decoy of the same magnitude as V, algebraically unrelated to the bases
DECOY = mp.mpf(V) + mp.sqrt(mp.mpf(2)) / 1000 - mp.mpf('0.0007071067811865')  # ~ V + tiny irrational


def main():
    print(f"Rung 3c -- broader guarded closure probe.  dps={mp.mp.dps}")
    print(f"  V = {mp.nstr(V, mp.mp.dps - 3)}\n")
    for name, basis in BASES.items():
        hit = guarded(name, V, basis, maxcoeff=10 ** 4)
        decoy = guarded(name + " [DECOY]", DECOY, basis, maxcoeff=10 ** 4)
        print(f"  {name}")
        print(f"      V     : {hit}")
        print(f"      decoy : {decoy}")
    print("\nRead: only a V-line tagged CANDIDATE with a decoy line of (none) is worth")
    print("auditing.  Any CANDIDATE that the decoy ALSO fits is a basis-richness artifact.")


if __name__ == '__main__':
    main()
