"""Sprint Hodge-SL2 — Steps 2-3: polarization + Mumford-Tate group.

Step 1 established that GeoVac carries the complex structure J (J^2 = -1,
the Kramers / SU(2)=S^3 quaternionic structure). This driver builds the
polarized weight-1 Hodge structure on the fundamental doublet V = Q^2 and
computes its Mumford-Tate group, deciding the sprint's headline question:

  POSITIVE   MT = SL_2   (generic, non-CM): GeoVac's SL_2 IS the MT group.
  BORDERLINE MT = torus  (CM): GeoVac realizes a CM point; MT is a proper
             1-dim torus; identify the CM field.

Canonical modeling choice (stated honestly): V_fund = Q^2 is realized as
the j=1/2 fundamental doublet (m_j = +/-1/2) of the lowest Dirac shell,
and the complex structure is the Kramers J restricted to that doublet.
This is the natural GeoVac realization: the SL_2 comes from S^3 = SU(2),
whose fundamental IS the spin-1/2 doublet, and the Kramers J restricted
there is the SU(2) quaternionic structure.

Run: python debug/sprint_hodge_sl2_step23_mumford_tate.py
"""

from __future__ import annotations

import json
from typing import Dict, List, Tuple

import sympy as sp
from sympy import Integer, Matrix, Rational, eye, sqrt

from geovac.spectral_triple import FockSpectralTriple
from geovac.tannakian import sl2_standard_action


def extract_fundamental_J(n_max: int = 2) -> Tuple[Matrix, Dict]:
    """Restrict the Kramers J to a j=1/2 fundamental doublet of shell n=1.

    Returns (J_2x2, info). The doublet is the (n_fock=1, kappa, m_j=+/-1/2)
    pair for the first kappa carrying both m_j signs. J only connects
    m_j <-> -m_j within a fixed (n,kappa), so the 2x2 block is the full
    action on that subspace (off-block entries are zero by construction).
    """
    st = FockSpectralTriple(n_max=n_max, j_type="kramers")
    J = st.real_structure
    labels = st.labels

    # Group n_fock=1 states by kappa -> {two_m_j: index}
    by_kappa: Dict[int, Dict[int, int]] = {}
    for i, lab in enumerate(labels):
        if lab.n_fock == 1:
            by_kappa.setdefault(lab.kappa, {})[lab.two_m_j] = i

    # Pick the first kappa carrying both m_j = +/- 1/2 (two_m_j = +/-1)
    chosen_kappa = None
    for kappa, mp in sorted(by_kappa.items()):
        if 1 in mp and -1 in mp:
            chosen_kappa = kappa
            break
    if chosen_kappa is None:
        raise RuntimeError("No j=1/2 doublet found in shell n=1.")

    i_plus = by_kappa[chosen_kappa][1]    # m_j = +1/2
    i_minus = by_kappa[chosen_kappa][-1]  # m_j = -1/2
    idx = [i_plus, i_minus]

    # Verify the doublet is J-invariant (no leakage outside the block)
    leakage = Integer(0)
    for a in idx:
        for b in range(st.dim_H):
            if b not in idx and (J[b, a] != 0):
                leakage += sp.Abs(J[b, a])

    J2 = Matrix([[J[idx[r], idx[c]] for c in range(2)] for r in range(2)])

    info = {
        "n_max": n_max,
        "chosen_kappa": int(chosen_kappa),
        "idx_plus_minus": [int(i_plus), int(i_minus)],
        "J_block": [[str(J2[r, c]) for c in range(2)] for r in range(2)],
        "J_invariant_block": (leakage == 0),
        "J_block_squared_is_minus_I": (J2 * J2 == -eye(2)),
    }
    return J2, info


def mumford_tate(J2: Matrix) -> Dict:
    """Build the polarized weight-1 HS and compute its Mumford-Tate group.

    V = Q^2, complex structure J2 (J2^2 = -I), polarization Q = symplectic
    form (sign chosen for Riemann positivity). MT(V,J,Q) = smallest
    Q-algebraic subgroup of SL_2 whose R-points contain the Hodge circle
    {cos t * I + sin t * J2}. Dichotomy for rank-2 weight-1:
      J2 rational  ==>  MT = norm-1 torus of Q[J2]  (CM, dim 1).
      J2 irrational ==> generically MT = SL_2 (dim 3).
    GeoVac's J2 is rational, so we expect a CM torus; identify the field.
    """
    # --- minimal polynomial / CM field ---
    t = sp.symbols("t")
    charpoly = (J2 - t * eye(2)).det()        # t^2 - tr*t + det
    minpoly = sp.factor(charpoly)
    tr = J2.trace()
    det = J2.det()
    # J^2 = -I  =>  minpoly t^2 + 1  =>  Q[J] = Q(i)
    cm_field = "Q(i)" if (tr == 0 and det == 1) else f"Q[t]/({sp.factor(charpoly)})"

    j2_rational = all(
        J2[r, c].is_rational for r in range(2) for c in range(2)
    )

    # --- polarization: symplectic form, sign fixed for Riemann positivity ---
    eps = Matrix([[0, 1], [-1, 0]])
    # Riemann form S(x,y) = x^T Q J y must be symmetric positive-definite.
    for sign in (Integer(1), Integer(-1)):
        Q = sign * eps
        # J symplectic: J^T Q J = Q
        j_symplectic = (J2.T * Q * J2 == Q)
        Rmat = Q * J2  # x^T (Q J) x  -> symmetric part is Q(x, Jx)
        sym = (Rmat == Rmat.T)
        pos_def = sym and all(Rmat.eigenvals(multiple=True)[k] > 0 for k in range(2)) \
            if sym else False
        if j_symplectic and pos_def:
            polarization = {
                "Q": [[int(Q[r, c]) for c in range(2)] for r in range(2)],
                "J_symplectic": True,
                "riemann_form_QJ": [[int(Rmat[r, c]) for c in range(2)] for r in range(2)],
                "riemann_positive_definite": True,
                "sign": int(sign),
            }
            break
    else:
        polarization = {"riemann_positive_definite": False}

    # --- Mumford-Tate group ---
    # The Hodge circle h(t) = cos t I + sin t J generates, over Q, the torus
    # T_J = { a I + b J : a,b in field, det = 1 } = norm-1 torus of Q[J].
    # det(a I + b J) = a^2 + b^2 (since tr J = 0, det J = 1) = N_{Q(i)/Q}(a+bi).
    # T_J preserves the definite Riemann form (it equals SO(2) as a Q-group),
    # hence is a PROPER 1-dim subgroup of SL_2 (which preserves no definite
    # form). So MT = T_J  !=  SL_2 whenever J is rational: a CM torus.
    a, b = sp.symbols("a b")
    g = a * eye(2) + b * J2
    det_torus = sp.expand(g.det())  # should be a^2 + b^2 for J = [[0,-1],[1,0]]

    mt_is_sl2 = not j2_rational  # rational J  <=>  CM torus (proper) <=> not SL_2
    mt = {
        "J2_rational": j2_rational,
        "cm_field": cm_field,
        "minpoly_J": str(sp.factor(charpoly)),
        "torus_det_form": str(det_torus),
        "torus_dim": 1,
        "sl2_dim": 3,
        "MT_is_full_SL2": bool(mt_is_sl2),
        "MT_is_CM_torus": bool(not mt_is_sl2),
        "verdict": "POSITIVE (MT = SL_2)" if mt_is_sl2 else "BORDERLINE (MT = CM torus)",
    }
    return {"mumford_tate": mt, "polarization": polarization}


def torus_acts_in_sl2() -> Dict:
    """Sanity: a rational torus element is a genuine SL_2(Q) point acting
    standardly via GeoVac's sl2_standard_action. Use a Pythagorean rotation
    g = [[3/5,-4/5],[4/5,3/5]] in SO(2)(Q) ⊂ SL_2(Q)."""
    g = Matrix([[Rational(3, 5), Rational(-4, 5)], [Rational(4, 5), Rational(3, 5)]])
    in_sl2 = (g.det() == 1)
    acts_standard = (sl2_standard_action(g) == g)
    # on the CM torus?  g = a I + b J with J = [[0,-1],[1,0]], a=3/5,b=4/5,
    # a^2+b^2 = 1  -> yes, a norm-1 element of Q(i).
    on_torus = (Rational(3, 5) ** 2 + Rational(4, 5) ** 2 == 1)
    return {
        "pyth_rotation_in_SL2": bool(in_sl2),
        "acts_via_standard_rep": bool(acts_standard),
        "is_norm1_Qi_element": bool(on_torus),
    }


def main() -> None:
    results = {}

    J2, info = extract_fundamental_J(n_max=2)
    results["fundamental_J"] = info
    print("=== Step 2-3: fundamental complex structure ===", flush=True)
    print(f"  chosen kappa = {info['chosen_kappa']};  J block = {info['J_block']}", flush=True)
    print(f"  J-invariant block: {info['J_invariant_block']};  "
          f"J^2 = -I: {info['J_block_squared_is_minus_I']}", flush=True)

    mt = mumford_tate(J2)
    results.update(mt)
    pol = mt["polarization"]
    mtg = mt["mumford_tate"]
    print("\n=== polarization ===", flush=True)
    print(f"  Q = {pol.get('Q')};  J symplectic: {pol.get('J_symplectic')};  "
          f"Riemann pos-def: {pol.get('riemann_positive_definite')}", flush=True)
    print("\n=== Mumford-Tate group ===", flush=True)
    print(f"  J rational: {mtg['J2_rational']};  minpoly(J) = {mtg['minpoly_J']};  "
          f"CM field = {mtg['cm_field']}", flush=True)
    print(f"  torus det form = {mtg['torus_det_form']};  "
          f"dim(MT) = {mtg['torus_dim']} vs dim(SL_2) = {mtg['sl2_dim']}", flush=True)
    print(f"  >>> VERDICT: {mtg['verdict']}", flush=True)

    sanity = torus_acts_in_sl2()
    results["torus_in_sl2_sanity"] = sanity
    print("\n=== sanity: CM torus inside GeoVac SL_2 ===", flush=True)
    print(f"  Pythagorean rotation in SL_2(Q): {sanity['pyth_rotation_in_SL2']};  "
          f"acts via standard rep: {sanity['acts_via_standard_rep']};  "
          f"norm-1 Q(i) element: {sanity['is_norm1_Qi_element']}", flush=True)

    # Triangulation flag (reported, NOT asserted as proven-connected)
    results["triangulation_note"] = {
        "hodge_cm_field": mtg["cm_field"],
        "level4_field": "Q(i) = Q(mu_4), the field of G_4 = MT(Z[i,1/2])",
        "period_content_field": "Q(i): Catalan G mixed-Tate over Q(i) (EMN 2025)",
        "claim": "three roads (Hodge CM field, motivic-Galois level, period content) "
                 "single out Q(i); whether this is forced or coincidental is itself "
                 "an audit target, NOT asserted here as a proven identity.",
    }

    out = "debug/data/sprint_hodge_sl2_step23.json"
    with open(out, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nwrote {out}", flush=True)


if __name__ == "__main__":
    main()
