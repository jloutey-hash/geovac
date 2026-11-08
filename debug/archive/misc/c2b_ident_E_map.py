"""EC3 stage E (docs/predictions/2026-10-22_c2b_identification.md): an exact map from the twisted
elliptic period operator to C2b.

Source: H_J = 2F1(1/12, 5/12; 1; 1728/J(x)) pulled back, then gauged by phi = x^(1/2) P^(-1/4) S^(1/2),
S in {1, (4x-1)(4x-25)}; phi'/phi is rational, so the source has polynomial coefficients b2, b1, b0.
Map: A = (N0 + N1 d/dx) / E, E = (4x - 1)^k (4x - 25)^k (4x + 5)  [source apparent point -5/4], deg N0 <= dN, deg N1 <= dN + 1.
Condition: C2b(A y) = 0 for y in Sol(source); with y'' = -(b1 y' + b0 y)/b2 this is U y + V y' = 0, i.e.
U = V = 0 as rational functions: an exact linear system over Q in the coefficients of N0, N1 (flint).

Usage: python debug/c2b_ident_E_map.py <solution json> <S: 1|q> <k> <dN>
<solution json> holds P, lam, widths (exact rationals), as written by stage D for a matching J."""
import json
import sys
from pathlib import Path

import flint
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
import c2b_ident_ops as O  # noqa: E402
import c2b_ident_C_jsolve as C  # noqa: E402
from c2b_ident_D_compare import pullback_op  # noqa: E402

x = sp.Symbol("x")
DATA = Path(__file__).parent / "data"


def fq(e):
    """sympy polynomial in x -> flint fmpq_poly"""
    cs = sp.Poly(sp.expand(e), x).all_coeffs()[::-1]
    return flint.fmpq_poly([flint.fmpq(int(sp.Rational(c).p), int(sp.Rational(c).q)) for c in cs])


def int_nullspace(rows, ncols):
    """exact nullspace of a rational matrix: rows scaled to integers, flint fmpz_mat.nullspace"""
    from math import lcm
    irows = []
    for row in rows:
        d = 1
        for v in row:
            d = lcm(d, int(flint.fmpq(v).q))
        irows.append([int(flint.fmpq(v) * d) for v in row])
    NS, nul = flint.fmpz_mat(irows).nullspace()
    return NS, nul


def gauge(ops, logder):
    """y = phi w, phi'/phi = logder: from a2 w'' + a1 w' + a0 w = 0 to the operator for y (polynomial form)"""
    a0, a1, a2 = ops
    # w = y/phi; w' = (y' - L y)/phi; w'' = (y'' - 2 L y' + (L^2 - L') y)/phi
    L = logder
    b2 = a2
    b1 = a1 - 2 * L * a2
    b0 = a0 - L * a1 + (L**2 - sp.diff(L, x)) * a2
    den = sp.lcm(*[sp.fraction(sp.cancel(v))[1] for v in (b0, b1, b2)])
    return [sp.expand(sp.cancel(v * den)) for v in (b0, b1, b2)]


def source_op(Pexpr, lam, widths, S):
    D = C.D_of([sp.Integer(0), sp.Rational(1, 4), sp.Rational(25, 4)], widths[:3])
    J = Pexpr**3 / (lam * D)
    ops = pullback_op(J)
    logder = sp.Rational(1, 2) / x - sp.Rational(1, 4) * sp.diff(Pexpr, x) / Pexpr
    if S == "q":
        logder += sp.Rational(1, 2) * (4 / (4 * x - 1) + 4 / (4 * x - 25))
    b = gauge(ops, sp.together(logder))
    g = sp.gcd_list([sp.expand(v) for v in b])
    return [sp.expand(sp.cancel(v / g)) for v in b]


def apply_c2b_to_image(a, b, n0, n1, E):
    """U, V (polynomials, common denominator cleared) with C2b((n0 y + n1 y')/E) = (U y + V y')/den.
    a: C2b coefficients [a0, a1, a2] (fmpq_poly), b: source [b0, b1, b2] (fmpq_poly)."""
    b0, b1, b2 = b
    # represent f y + g y' as pair (f, g) of rational functions with denominator tracked as power of b2 and E
    # z = (n0 y + n1 y')/E. Work with numerators over the common denominator E^m b2^k.
    # derivative of (f y + g y') / (E^m b2^k): use the rule on numerators.
    def deriv(f, g, m, k):
        # d/dx [ (f y + g y') / (E^m b2^k) ],  y'' = -(b1 y' + b0 y)/b2
        # numerator over E^(m+1) b2^(k+1):
        Ep = E.derivative()
        b2p = b2.derivative()
        F = (f.derivative() * E * b2 - f * (m * Ep * b2 + k * E * b2p)) - g * b0 * E
        G = (g.derivative() * E * b2 - g * (m * Ep * b2 + k * E * b2p)) + f * E * b2 - g * b1 * E
        return F, G, m + 1, k + 1
    f0, g0 = n0, n1
    f1, g1, m1, k1 = deriv(f0, g0, 1, 0)
    f2, g2, m2, k2 = deriv(f1, g1, m1, k1)
    a0, a1, a2 = a
    # common denominator E^3 b2^2: z = (f0, g0)/E, z' = (f1, g1)/(E^2 b2), z'' = (f2, g2)/(E^3 b2^2)
    U = a2 * f2 + a1 * f1 * E * b2 + a0 * f0 * E**2 * b2**2
    V = a2 * g2 + a1 * g1 * E * b2 + a0 * g0 * E**2 * b2**2
    return U, V


def main(solfile, S, k, dN):
    sol = json.load(open(solfile))
    Pexpr = sp.sympify(sol["P"])
    lam = sp.Rational(sol["lam"])
    widths = sol["widths"]
    src = source_op(Pexpr, lam, widths, S)
    b = [fq(v) for v in src]
    ac, _ = O.c2b()
    a = [flint.fmpq_poly([flint.fmpq(int(c)) for c in aj]) for aj in ac]
    E = fq(((4 * x - 1) * (4 * x - 25)) ** k * (4 * x + 5))
    cols = []
    labels = []
    X = flint.fmpq_poly([0, 1])
    zero = flint.fmpq_poly([])
    for j in range(dN + 1):
        U, V = apply_c2b_to_image(a, b, X**j, zero, E)
        cols.append((U, V)); labels.append(("N0", j))
    for j in range(dN + 2):
        U, V = apply_c2b_to_image(a, b, zero, X**j, E)
        cols.append((U, V)); labels.append(("N1", j))
    degU = max(max(c[0].degree(), c[1].degree()) for c in cols) + 1
    rows = []
    for r in range(degU):
        rows.append([c[0].coeffs()[r] if r < len(c[0].coeffs()) else 0 for c in cols])
        rows.append([c[1].coeffs()[r] if r < len(c[1].coeffs()) else 0 for c in cols])
    NS, nul = int_nullspace(rows, len(cols))
    print(f"S = {S}, k = {k}, dN = {dN}: {len(cols)} unknowns, {len(rows)} equations, nullity {nul}", flush=True)
    out = {"S": S, "k": k, "dN": dN, "nullity": nul, "source": [str(v) for v in src]}
    if nul >= 1:
        v = [NS[i, 0] for i in range(len(cols))]
        den = 1
        v = [int(c) for c in v]
        N0 = sum(c * x**j for (nm, j), c in zip(labels, v) if nm == "N0")
        N1 = sum(c * x**j for (nm, j), c in zip(labels, v) if nm == "N1")
        out.update({"N0": str(N0), "N1": str(N1), "E": str(sp.factor(((4 * x - 1) * (4 * x - 25)) ** k * (4 * x + 5)))})
        # exact re-check and the 1/1000 perturbation
        n0, n1 = fq(N0), fq(N1)
        U, V = apply_c2b_to_image(a, b, n0, n1, E)
        out["exact_identity"] = U.is_zero() and V.is_zero()
        lead = max(range(len(n0.coeffs())), key=lambda i: abs(n0.coeffs()[i]))
        cs = list(n0.coeffs()); cs[lead] = cs[lead] * flint.fmpq(1001, 1000)
        Up, Vp = apply_c2b_to_image(a, b, flint.fmpq_poly(cs), n1, E)
        out["perturbed_breaks"] = not (Up.is_zero() and Vp.is_zero())
        print("   exact identity:", out["exact_identity"], " 1/1000 perturbation breaks it:", out["perturbed_breaks"])
    json.dump(out, open(DATA / f"c2b_ident_E_map_{S}_{k}_{dN}.json", "w"), indent=1)
    return out


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
