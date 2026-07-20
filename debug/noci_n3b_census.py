"""N3b -- genuine-integral sparsity census (NOCI sandbox, branch sandbox/noci).

Roadmap step N3b of debug/noci_sandbox_notes.md: compute the GENUINE bare
(S, h, g) structure for a real GeoVac fragment pair -- hydrogenic orbitals
with l > 0 on two centers (A: Z=3 'Li-like' at origin, B: Z=1 'H' at R zhat)
-- and census its sparsity against (i) the corpus builder's surrogate tensor
and (ii) the symmetry-allowed dense ceiling.

Machinery:
  * S and one-electron kernels 1/r_A, 1/r_B: EXACT Mulliken/Ruedenberg
    auxiliary-integral route (reuses compute_topos3_exact_meet's radial and
    angular assembly; adds the kernel factor, which cancels against the
    volume element (xi^2 - eta^2) = (xi+eta)(xi-eta) so the integrand stays
    polynomial -- asserted).  Zero-decidability: entry = e^{-p} (U e^q +
    V e^{-q}) with U, V exact rationals; zero <=> U = V = 0 (Lindemann).
  * h census: hydrogenic eigen-trick.  Orbital b on center B has orbital
    charge Z_B, so (T + V_B)|b> = E_b|b> and
        h_ab = E_b S_ab - Z_A <a|1/r_A|b>          (ket trick)
             = E_a S_ab - Z_B <a|1/r_B|b>          (bra trick).
    The two assemblies must agree EXACTLY (rational U/V equality) -- this is
    asserted for every cross entry and certifies the whole machinery.
  * Same-center h off-blocks (<a|V_other|a'>): geovac.shibuya_wulfman
    (corpus module; exact structural selection rules, float values).
  * g: STRUCTURAL census only.  Exact axial selection rule
    m_p + m_r = m_q + m_s (rotation invariance about the molecular axis),
    plus multipole feasibility: all-four-same-center entries need a shared
    Gaunt-allowed L (single-center multipole); same-center distributions in
    cross classes need a per-side Gaunt-allowed L (bipolar expansion couples
    (L1, L2) freely); two-center distributions have no L constraint.
    No numerical two-center ERIs are computed (that engine is N4's build);
    zeros-by-symmetry are exact, 'allowed' entries are generically nonzero
    (no accidental-zero claim is made).

Validation: mpmath prolate-spheroidal quadrature cross-checks on kernel and
overlap entries; SW cross-check vs 1D radial quadrature; the exact bra/ket
identity above.

Run from debug/ cwd:  python noci_n3b_census.py
Output: debug/data/noci_n3b_census_results.json
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction

sys.path.insert(0, ".")
sys.path.insert(0, "..")

import mpmath as mp
import sympy as sp

from compute_topos3_exact_meet import (
    G_lm,
    eta,
    radial_norm_sq,
    radial_poly_coeffs,
    theta_norm_sq,
    xi,
)
from compute_topos3_two_center_meet import R_nl_modern, angular
from geovac.shibuya_wulfman import compute_cross_center_vne

mp.mp.dps = 20


# ------------------------------------------------------- exact two-center

def two_center_UV(Z1: Fraction, n1: int, l1: int,
                  Z2: Fraction, n2: int, l2: int,
                  m: int, R: Fraction, kernel: str | None = None):
    """Exact unnormalized two-center integral
        <phi_{n1 l1 m}(0, Z1) | K | phi_{n2 l2 m}(R zhat, Z2)>
    for K in {1 (overlap), 1/r_A ('inv_ra'), 1/r_B ('inv_rb')}.

    Returns (p, q, U, V):  value = e^{-p} (U e^{q} + V e^{-q})  [q != 0]
    or (p, 0, U, None):    value = e^{-p} U                      [q == 0].
    Generalizes compute_topos3_exact_meet.overlap_UV by the kernel factor,
    which cancels against the volume element -- cancellation asserted.
    """
    m = abs(m)
    Rs = sp.Rational(R.numerator, R.denominator)
    half = Rs / 2
    c1, a = radial_poly_coeffs(Z1, n1, l1)
    c2, b = radial_poly_coeffs(Z2, n2, l2)
    ra = half * (xi + eta)
    rb = half * (xi - eta)
    ct_a = (1 + xi * eta) / (xi + eta)
    ct_b = (xi * eta - 1) / (xi - eta)
    Rad1 = sum(c * ra ** k for k, c in c1.items())
    Rad2 = sum(c * rb ** k for k, c in c2.items())
    s2 = (xi ** 2 - 1) * (1 - eta ** 2)
    ang = (s2 ** m / ((xi + eta) * (xi - eta)) ** m
           * G_lm(l1, m, ct_a) * G_lm(l2, m, ct_b))
    vol = half ** 3 * (xi ** 2 - eta ** 2)
    if kernel == "inv_ra":
        ker = 2 / (Rs * (xi + eta))
    elif kernel == "inv_rb":
        ker = 2 / (Rs * (xi - eta))
    elif kernel is None:
        ker = sp.Integer(1)
    else:
        raise ValueError(kernel)
    expr = sp.cancel(sp.together(Rad1 * Rad2 * ang * vol * ker))
    num, den = sp.fraction(expr)
    assert den == 1, f"classical cancellation FAILED ({kernel}): den = {den}"
    poly = sp.Poly(sp.expand(num), xi, eta)
    p = (a + b) * Rs / 2
    q = (a - b) * Rs / 2

    max_i = max(mon[0] for mon in poly.monoms())
    max_j = max(mon[1] for mon in poly.monoms())
    a_coef = {0: 1 / p}
    for i in range(1, max_i + 1):
        a_coef[i] = (1 + i * a_coef[i - 1]) / p
    if q != 0:
        u = {0: 1 / q}
        v = {0: -1 / q}
        for j in range(1, max_j + 1):
            u[j] = sp.Rational((-1) ** j, 1) / q + j * u[j - 1] / q
            v[j] = sp.Rational(-1, 1) / q + j * v[j - 1] / q
        U = sp.Integer(0)
        V = sp.Integer(0)
        for (i, j), c in zip(poly.monoms(), poly.coeffs()):
            U += c * a_coef[i] * u[j]
            V += c * a_coef[i] * v[j]
        return sp.nsimplify(p), sp.nsimplify(q), sp.nsimplify(U), sp.nsimplify(V)
    U = sp.Integer(0)
    for (i, j), c in zip(poly.monoms(), poly.coeffs()):
        if j % 2 == 0:
            U += c * a_coef[i] * sp.Rational(2, j + 1)
    return sp.nsimplify(p), sp.Integer(0), sp.nsimplify(U), None


def uv_is_zero(U, V) -> bool:
    return (U == 0) if V is None else (U == 0 and V == 0)


def uv_float(p, q, U, V) -> float:
    if V is None:
        return float(sp.N(sp.exp(-p) * U, 30))
    return float(sp.N(sp.exp(-p) * (U * sp.exp(q) + V * sp.exp(-q)), 30))


def norm_const(Z: Fraction, n: int, l: int, m: int) -> sp.Expr:
    """Normalization of the unnormalized orbital used in two_center_UV
    (azimuthal 1/sqrt(2pi) factors cancel in m-diagonal entries)."""
    return sp.sqrt(radial_norm_sq(Z, n, l) * theta_norm_sq(l, abs(m)))


# ------------------------------------------------ mpmath quadrature check

def quad_two_center(Z1, n1, l1, Z2, n2, l2, m, R, kernel=None):
    """Normalized two-center integral by prolate-spheroidal quadrature."""
    R = mp.mpf(R)
    half = R / 2

    def integrand(x, e):
        r1 = half * (x + e)
        r2 = half * (x - e)
        if r1 == 0 or r2 == 0:
            return mp.mpf(0)
        ct1 = (1 + x * e) / (x + e)
        ct2 = (x * e - 1) / (x - e)
        val = (R_nl_modern(Z1, n1, l1, r1) * angular(l1, m, ct1)
               * R_nl_modern(Z2, n2, l2, r2) * angular(l2, m, ct2))
        w = half ** 3 * (x ** 2 - e ** 2)
        if kernel == "inv_ra":
            w *= 2 / (R * (x + e))
        elif kernel == "inv_rb":
            w *= 2 / (R * (x - e))
        return val * w

    s = float(1 / R)
    xi_pts = [1, 1 + 1 * s, 1 + 4 * s, 1 + 16 * s, mp.inf]
    return mp.quad(lambda x: mp.quad(lambda e: integrand(x, e), [-1, 0, 1]),
                   xi_pts)


# ------------------------------------------------------------- S/h census

BASIS_NMAX = 2
Z_A = Fraction(3)     # 'Li-like' fragment at the origin
Z_B = Fraction(1)     # H fragment at R zhat
R_AB = Fraction(3)    # bohr (LiH-scale)


def basis(n_max: int):
    return [(n, l, mm) for n in range(1, n_max + 1)
            for l in range(n) for mm in range(-l, l + 1)]


def E_hyd(Z: Fraction, n: int) -> sp.Rational:
    return -sp.Rational(Z.numerator ** 2, 2 * n ** 2 * Z.denominator ** 2)


def cross_block_census():
    """Exact census of the cross-center S and h blocks (A bra, B ket)."""
    orbs_A = basis(BASIS_NMAX)
    orbs_B = basis(BASIS_NMAX)
    S_entries = {}
    h_entries = {}
    identity_failures = []

    for (n1, l1, m1) in orbs_A:
        NA = norm_const(Z_A, n1, l1, m1)
        for (n2, l2, m2) in orbs_B:
            key = f"A({n1},{l1},{m1})|B({n2},{l2},{m2})"
            if m1 != m2:
                S_entries[key] = {"zero": True, "why": "m-rule", "value": 0.0}
                h_entries[key] = {"zero": True, "why": "m-rule", "value": 0.0}
                continue
            NB = norm_const(Z_B, n2, l2, m2)
            NN = NA * NB
            m = m1
            pS, qS, US, VS = two_center_UV(Z_A, n1, l1, Z_B, n2, l2, m, R_AB)
            pa, qa, Ua, Va = two_center_UV(Z_A, n1, l1, Z_B, n2, l2, m, R_AB,
                                           kernel="inv_ra")
            pb, qb, Ub, Vb = two_center_UV(Z_A, n1, l1, Z_B, n2, l2, m, R_AB,
                                           kernel="inv_rb")
            assert (pS, qS) == (pa, qa) == (pb, qb)

            S_zero = uv_is_zero(US, VS)
            S_val = 0.0 if S_zero else uv_float(pS, qS, US, VS) / float(sp.N(NN, 30))
            S_entries[key] = {"zero": bool(S_zero), "why": "exact" if S_zero
                              else None, "value": S_val}

            # h = T + V_A + V_B:
            #   ket trick: h = E_b S - Z_A I[1/r_A]
            #   bra trick: h = E_a S - Z_B I[1/r_B]
            Ea, Eb = E_hyd(Z_A, n1), E_hyd(Z_B, n2)
            ZAr = sp.Rational(Z_A.numerator, Z_A.denominator)
            ZBr = sp.Rational(Z_B.numerator, Z_B.denominator)
            U_ket = Eb * US - ZAr * Ua
            V_ket = None if VS is None else Eb * VS - ZAr * Va
            U_bra = Ea * US - ZBr * Ub
            V_bra = None if VS is None else Ea * VS - ZBr * Vb
            dU = sp.simplify(U_ket - U_bra)
            dV = sp.Integer(0) if VS is None else sp.simplify(V_ket - V_bra)
            if dU != 0 or dV != 0:
                identity_failures.append((key, str(dU), str(dV)))

            h_zero = uv_is_zero(U_ket, V_ket)
            h_val = 0.0 if h_zero else uv_float(pS, qS, U_ket, V_ket) / float(sp.N(NN, 30))
            h_entries[key] = {"zero": bool(h_zero), "why": "exact" if h_zero
                              else None, "value": h_val}

    return S_entries, h_entries, identity_failures


def same_center_h_blocks():
    """Same-center h blocks: E_n diagonal + SW <a|V_other|a'> (corpus module)."""
    states = basis(BASIS_NMAX)
    out = {}
    for label, Zorb, Znuc, parity in (("A", float(Z_A), float(Z_B), +1),
                                      ("B", float(Z_B), float(Z_A), -1)):
        v = compute_cross_center_vne(Zorb, states, Znuc, float(R_AB),
                                     L_max=4, nuc_parity=parity)
        blk = {}
        for i, (n1, l1, m1) in enumerate(states):
            for j, (n2, l2, m2) in enumerate(states):
                val = v[i, j]
                if i == j:
                    val += float(sp.N(E_hyd(Fraction(int(Zorb)), n1), 30))
                blk[f"({n1},{l1},{m1})|({n2},{l2},{m2})"] = val
        out[label] = blk
    return out


# --------------------------------------------------------------- g census

def side_L_range(l1: int, l2: int, M: int):
    lo = max(abs(l1 - l2), abs(M))
    hi = l1 + l2
    return [L for L in range(lo, hi + 1) if (l1 + l2 + L) % 2 == 0]


def census_g(n_max: int):
    """Structural (symmetry-rule) census of the genuine molecular ERI tensor
    vs the builder's block-diagonal surrogate, over all ordered (pq|rs)."""
    orbs = ([("A", n, l, mm) for (n, l, mm) in basis(n_max)]
            + [("B", n, l, mm) for (n, l, mm) in basis(n_max)])
    M = len(orbs)
    total = M ** 4
    mrule = 0
    genuine = 0
    builder = 0
    classes = {}

    def side_label(o1, o2):
        c = {o1[0], o2[0]}
        return "AB" if len(c) == 2 else (o1[0] + o1[0])

    for p in orbs:
        for q in orbs:
            M1 = q[3] - p[3]
            lab1 = side_label(p, q)
            s1_same = lab1 != "AB"
            L1 = side_L_range(p[2], q[2], M1) if s1_same else None
            for r in orbs:
                for s_ in orbs:
                    if p[3] + r[3] != q[3] + s_[3]:
                        continue
                    mrule += 1
                    lab2 = side_label(r, s_)
                    s2_same = lab2 != "AB"
                    M2 = s_[3] - r[3]  # = -M1 by the m-rule... (see note)
                    if s1_same and s2_same and p[0] == r[0]:
                        # all four on one center: shared-L single-center multipole
                        L2 = side_L_range(r[2], s_[2], M2)
                        ok = bool(set(L1) & set(L2))
                        if ok:
                            builder += 1
                    else:
                        ok1 = bool(L1) if s1_same else True
                        ok2 = bool(side_L_range(r[2], s_[2], M2)) if s2_same else True
                        ok = ok1 and ok2
                    if ok:
                        genuine += 1
                        cls = "|".join(sorted([lab1, lab2]))
                        classes[cls] = classes.get(cls, 0) + 1

    return {"M": M, "total_dense": total, "m_rule_allowed": mrule,
            "genuine_allowed": genuine, "builder_allowed": builder,
            "classes": classes}


# ------------------------------------------------------------ validations

def validations():
    out = {}

    # V1: p-p overlap, production config, exact vs quadrature
    p_, q_, U_, V_ = two_center_UV(Z_A, 2, 1, Z_B, 2, 1, 0, R_AB)
    NN = float(sp.N(norm_const(Z_A, 2, 1, 0) * norm_const(Z_B, 2, 1, 0), 30))
    ex = uv_float(p_, q_, U_, V_) / NN
    qd = float(quad_two_center(3, 2, 1, 1, 2, 1, 0, 3.0))
    out["V1_overlap_2p0A_2p0B"] = {"exact": ex, "quad": qd,
                                   "rel_err": abs(ex - qd) / abs(qd)}

    # V2: same pair, kernel 1/r_A
    p_, q_, U_, V_ = two_center_UV(Z_A, 2, 1, Z_B, 2, 1, 0, R_AB, kernel="inv_ra")
    ex = uv_float(p_, q_, U_, V_) / NN
    qd = float(quad_two_center(3, 2, 1, 1, 2, 1, 0, 3.0, kernel="inv_ra"))
    out["V2_inv_ra_2p0A_2p0B"] = {"exact": ex, "quad": qd,
                                  "rel_err": abs(ex - qd) / abs(qd)}

    # V3: 1s-1s, kernel 1/r_B
    p_, q_, U_, V_ = two_center_UV(Z_A, 1, 0, Z_B, 1, 0, 0, R_AB, kernel="inv_rb")
    NN2 = float(sp.N(norm_const(Z_A, 1, 0, 0) * norm_const(Z_B, 1, 0, 0), 30))
    ex = uv_float(p_, q_, U_, V_) / NN2
    qd = float(quad_two_center(3, 1, 0, 1, 1, 0, 0, 3.0, kernel="inv_rb"))
    out["V3_inv_rb_1sA_1sB"] = {"exact": ex, "quad": qd,
                                "rel_err": abs(ex - qd) / abs(qd)}

    # V4: SW same-center element vs 1D radial quadrature
    # <1s_A| -Z_B/r_B |1s_A> at Z_orb=3, Z_nuc=1: only L=0 survives (s-s).
    from geovac.shibuya_wulfman import compute_cross_center_vne_element
    sw = compute_cross_center_vne_element(3.0, 1, 0, 0, 1, 0, 0, 1.0, 3.0,
                                          L_max=4)
    f = lambda r: R_nl_modern(3, 1, 0, r) ** 2 * r ** 2 / max(r, mp.mpf(3))
    qd = -float(mp.quad(f, [0, 3, mp.inf]))
    out["V4_SW_1sA_VB_1sA"] = {"sw": sw, "quad": qd,
                               "rel_err": abs(sw - qd) / abs(qd)}
    return out


# ------------------------------------------------------------------ main

def main():
    results = {"config": {"Z_A": str(Z_A), "Z_B": str(Z_B), "R_AB": str(R_AB),
                          "basis_n_max": BASIS_NMAX,
                          "orbitals_per_center": len(basis(BASIS_NMAX))}}

    print("=== N3b genuine-integral sparsity census ===")
    print(f"A: Z={Z_A} at origin; B: Z={Z_B} at R={R_AB} zhat; "
          f"basis n_max={BASIS_NMAX} per center "
          f"({len(basis(BASIS_NMAX))} orbitals each)\n")

    print("--- validations ---")
    val = validations()
    results["validations"] = val
    for k, v in val.items():
        print(f"  {k}: rel_err = {v['rel_err']:.2e}")

    print("\n--- cross-center S and h blocks (exact) ---")
    S_cross, h_cross, id_fail = cross_block_census()
    results["identity_failures"] = id_fail
    print(f"  bra/ket eigen-trick exact identity: "
          f"{'PASS (all entries)' if not id_fail else f'FAIL {id_fail}'}")

    nS = sum(1 for v in S_cross.values() if not v["zero"])
    nh = sum(1 for v in h_cross.values() if not v["zero"])
    n_tot = len(S_cross)
    print(f"  S cross-block: {nS}/{n_tot} nonzero "
          f"(all {n_tot - nS} zeros are m-rule)")
    print(f"  h cross-block: {nh}/{n_tot} nonzero "
          f"(all {n_tot - nh} zeros are m-rule)")
    results["S_cross"] = S_cross
    results["h_cross"] = h_cross
    print("  cross-S magnitudes (nonzero entries):")
    for k, v in S_cross.items():
        if not v["zero"]:
            print(f"    {k}: S = {v['value']:+.6f}   h = {h_cross[k]['value']:+.6f}")

    print("\n--- same-center h blocks (SW corpus module) ---")
    sc = same_center_h_blocks()
    results["h_same_center"] = sc
    for lbl, blk in sc.items():
        nz = sum(1 for x in blk.values() if abs(x) > 0.0)
        print(f"  {lbl}-block: {nz}/{len(blk)} nonzero")

    # full-matrix counts (M=10)
    Msz = 2 * len(basis(BASIS_NMAX))
    S_nz = Msz + 2 * nS                      # identity blocks + 2 cross blocks
    h_same_nz = sum(sum(1 for x in blk.values() if abs(x) > 0.0)
                    for blk in sc.values())
    h_nz = h_same_nz + 2 * nh
    results["matrix_counts"] = {
        "M": Msz,
        "S_genuine_nonzero": S_nz, "S_builder_nonzero": Msz,
        "h_genuine_nonzero": h_nz, "h_builder_nonzero": h_same_nz,
        "dense": Msz * Msz}
    print(f"\n  full {Msz}x{Msz} matrices: "
          f"S genuine {S_nz} vs builder {Msz} (identity) vs dense {Msz * Msz}")
    print(f"                         "
          f"h genuine {h_nz} vs builder {h_same_nz} (no cross block, W1d) "
          f"vs dense {Msz * Msz}")

    print("\n--- g structural census (symmetry rules, ordered (pq|rs)) ---")
    results["g_census"] = {}
    for n_max in (2, 3):
        g = census_g(n_max)
        results["g_census"][f"n_max={n_max}"] = g
        M = g["M"]
        print(f"  n_max={n_max} (M={M}): dense {g['total_dense']:,} | "
              f"m-rule {g['m_rule_allowed']:,} "
              f"({g['m_rule_allowed'] / g['total_dense']:.1%}) | "
              f"genuine {g['genuine_allowed']:,} "
              f"({g['genuine_allowed'] / g['total_dense']:.1%}) | "
              f"builder {g['builder_allowed']:,} "
              f"({g['builder_allowed'] / g['total_dense']:.1%})")
        print(f"    classes: " + ", ".join(
            f"{k}: {v:,}" for k, v in sorted(g["classes"].items())))
        print(f"    genuine/builder inflation: "
              f"{g['genuine_allowed'] / g['builder_allowed']:.1f}x")

    with open("data/noci_n3b_census_results.json", "w") as fh:
        json.dump(results, fh, indent=1, default=str)
    print("\nwrote data/noci_n3b_census_results.json")


if __name__ == "__main__":
    main()
