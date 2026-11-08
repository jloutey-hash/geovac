"""Driver for Sprint Q5'-TC-1f:\ SL_2 / Peter-Weyl substrate + full-panel injectivity.

Closes the explicit U^*_Levi = G_a^{3N} x SL_2 inclusion into Aut^otimes(omega)
on the abelian primitive substrate of the Tannakian closure arc.

Three blocks:
    1. SL_2 axioms on Peter-Weyl reps (invertibility, tensor, group hom).
    2. G_a x SL_2 commutativity on combined reps.
    3. Full-panel injectivity: at n_max=2, build all 15 = 3 * N(2) test reps
       (one per primitive generator) and verify Phi(t) is non-trivial on V_g
       iff t_g != 0.

See debug/sprint_q5p_tc1f_sl2_inclusion_memo.md.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

from sympy import Integer, Matrix, Rational

PKG_DIR = Path(__file__).resolve().parents[1]
if str(PKG_DIR) not in sys.path:
    sys.path.insert(0, str(PKG_DIR))

from geovac.pro_system import primitive_generators
from geovac.tannakian import (  # noqa: E402
    FinDimRep,
    PWRep,
    _pw_standard_rep,
    _pw_sym2_rep,
    _pw_trivial_rep,
    primitive_generator_rep,
    sl2_standard_action,
    verify_ga_sl2_commute,
    verify_injectivity_at_generator,
    verify_sl2_group_homomorphism,
    verify_sl2_invertibility,
    verify_sl2_tensor,
)

OUT = Path(__file__).resolve().parent / "data" / "q5p_tc1f_sl2_inclusion.json"


def _sl2_test_panel():
    """A panel of SL_2(Q) elements:\ identity, unipotent, torus, generic."""
    e = Matrix([[1, 0], [0, 1]])
    u = Matrix([[1, 1], [0, 1]])
    t = Matrix([[2, 0], [0, Rational(1, 2)]])
    # generic: det = 5 * 3 - 2 * 7 = 15 - 14 = 1
    s = Matrix([[5, 2], [7, 3]])
    # involution: det = 0 * 0 - 1 * (-1) = 1
    w = Matrix([[0, 1], [-1, 0]])
    return [
        ("identity", e),
        ("unipotent", u),
        ("torus", t),
        ("generic", s),
        ("weyl", w),
    ]


def _t_panel():
    """Panel of t-vectors over primitive generators of n_max=2."""
    gens = primitive_generators(2)
    return [
        {},  # zero
        {gens[0]: Rational(1)},
        {gens[5]: Rational(-2, 3)},  # slot k=1
        {gens[10]: Rational(7, 4)},  # slot k=2
        {gens[0]: Rational(1, 2), gens[5]: Rational(-1, 5), gens[10]: Rational(3)},
    ]


def _block_sl2_axioms():
    """Block 1:\ SL_2 axioms on PW reps."""
    reps = [
        _pw_trivial_rep(dim=1),
        _pw_standard_rep(),
        _pw_sym2_rep(),
    ]
    panel = _sl2_test_panel()
    out = {"invertibility": [], "tensor": [], "group_hom": []}
    for name, g in panel:
        for V in reps:
            v = verify_sl2_invertibility(g, V)
            v["g_label"] = name
            out["invertibility"].append(v)
        for V in reps:
            for W in reps:
                v = verify_sl2_tensor(g, V, W)
                v["g_label"] = name
                out["tensor"].append(v)
    # group hom on a few products
    for (n1, g1) in panel:
        for (n2, g2) in panel:
            for V in reps:
                v = verify_sl2_group_homomorphism(g1, g2, V)
                v["g1_label"] = n1
                v["g2_label"] = n2
                out["group_hom"].append(v)
    n_inv = len(out["invertibility"])
    n_inv_pass = sum(1 for v in out["invertibility"] if v["bit_exact"])
    n_tensor = len(out["tensor"])
    n_tensor_pass = sum(1 for v in out["tensor"] if v["bit_exact"])
    n_gh = len(out["group_hom"])
    n_gh_pass = sum(1 for v in out["group_hom"] if v["bit_exact"])
    return out, {
        "invertibility": (n_inv, n_inv_pass),
        "tensor": (n_tensor, n_tensor_pass),
        "group_hom": (n_gh, n_gh_pass),
    }


def _block_ga_sl2_commute():
    """Block 2:\ G_a x SL_2 commutativity on V (x) V_PW."""
    panel = _sl2_test_panel()
    t_panel = _t_panel()
    V_panel = [
        FinDimRep(
            n_max=2,
            dim=2,
            endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])},
            label="J2",
        ),
        FinDimRep(
            n_max=2,
            dim=3,
            endos={(1, 0, 0): Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])},
            label="J3",
        ),
    ]
    V_pw_panel = [_pw_standard_rep(), _pw_sym2_rep()]
    out = []
    for t in t_panel:
        for name, g in panel:
            for V in V_panel:
                for V_pw in V_pw_panel:
                    v = verify_ga_sl2_commute(t, g, V, V_pw)
                    v["g_label"] = name
                    out.append(v)
    n = len(out)
    n_pass = sum(1 for v in out if v["bit_exact"])
    return out, (n, n_pass)


def _block_full_panel_injectivity():
    """Block 3:\ full-panel injectivity at n_max=2 (15 generators)."""
    gens = primitive_generators(2)
    n_max = 2
    # Build all 15 test reps V_g.
    test_reps = [primitive_generator_rep(n_max, g) for g in gens]

    out = []
    # Case A: zero t -> all eta_{V_g}(t) = I.
    t = {}
    for g in gens:
        v = verify_injectivity_at_generator(t, g, n_max)
        v["case"] = "zero_t"
        out.append(v)
    # Case B: single-generator t -> only matching V_g is non-trivial.
    for g_active in gens:
        t = {g_active: Rational(1)}
        for g_test in gens:
            v = verify_injectivity_at_generator(t, g_test, n_max)
            v["case"] = "single_t"
            v["g_active"] = list(g_active)
            out.append(v)
    # Case C: a generic non-zero t on 3 generators.
    t = {
        gens[0]: Rational(1, 2),
        gens[5]: Rational(-1, 3),
        gens[10]: Rational(7, 4),
    }
    for g in gens:
        v = verify_injectivity_at_generator(t, g, n_max)
        v["case"] = "generic_t"
        out.append(v)
    n = len(out)
    n_pass = sum(1 for v in out if v["bit_exact"])
    n_gens = len(gens)
    return out, (n, n_pass, n_gens)


def main():
    print("Sprint Q5'-TC-1f: SL_2 / Peter-Weyl substrate + full-panel injectivity")
    print("=" * 72)

    print("\n[Block 1] SL_2 axioms on PW reps...")
    block1, totals1 = _block_sl2_axioms()
    print(f"  invertibility: {totals1['invertibility'][1]}/{totals1['invertibility'][0]} bit-exact")
    print(f"  tensor:        {totals1['tensor'][1]}/{totals1['tensor'][0]} bit-exact")
    print(f"  group hom:     {totals1['group_hom'][1]}/{totals1['group_hom'][0]} bit-exact")

    print("\n[Block 2] G_a x SL_2 commute on combined reps...")
    block2, totals2 = _block_ga_sl2_commute()
    print(f"  G_a x SL_2:    {totals2[1]}/{totals2[0]} bit-exact")

    print("\n[Block 3] Full-panel injectivity at n_max=2...")
    block3, totals3 = _block_full_panel_injectivity()
    print(f"  generators:    {totals3[2]} = 3 * N(2) = 3 * 5")
    print(f"  injectivity:   {totals3[1]}/{totals3[0]} bit-exact")

    n_total = (
        totals1["invertibility"][0]
        + totals1["tensor"][0]
        + totals1["group_hom"][0]
        + totals2[0]
        + totals3[0]
    )
    n_pass = (
        totals1["invertibility"][1]
        + totals1["tensor"][1]
        + totals1["group_hom"][1]
        + totals2[1]
        + totals3[1]
    )
    print(f"\nTOTAL: {n_pass}/{n_total} bit-exact residuals")
    assert n_pass == n_total, f"FAIL: {n_total - n_pass} non-bit-exact residuals"

    OUT.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "sprint": "Q5'-TC-1f",
        "n_total": n_total,
        "n_pass": n_pass,
        "block1_sl2_axioms": {
            "invertibility": {"n": totals1["invertibility"][0], "pass": totals1["invertibility"][1]},
            "tensor": {"n": totals1["tensor"][0], "pass": totals1["tensor"][1]},
            "group_hom": {"n": totals1["group_hom"][0], "pass": totals1["group_hom"][1]},
        },
        "block2_ga_sl2_commute": {"n": totals2[0], "pass": totals2[1]},
        "block3_injectivity_at_nmax2": {
            "n": totals3[0],
            "pass": totals3[1],
            "n_generators": totals3[2],
        },
    }
    with OUT.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=str)
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
