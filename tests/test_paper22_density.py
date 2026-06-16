"""
Paper 22 (Potential-Independent Angular Sparsity) — Theorem 3 density backing.

Paper 22 §III Theorem 3 reports TWO angular ERI densities over the single-shell
angular index set {(l, m) : l <= l_max}, by exact enumeration of the Gaunt
selection rules:

  * D(l_max)  — the PHYSICAL Coulomb density.  An ERI <ab|cd> survives iff
        (i)  global M_L conservation   m_a + m_b = m_c + m_d, AND
        (ii) there exists a multipole k satisfying, for BOTH pairs, the
             triangle inequality, the l+l'+k-even parity rule, the
             (l k l'; 0 0 0) reduced 3j, and the bottom-row 3j
             (l k l'; -m, m-m', m').
     Paper values: 14.8438 / 8.5200 / 6.0608 / 4.8274 / 3.9863 % at l_max=1..5.
     This is the headline density and the abstract/conclusion number.

  * D_pd(l_max) — the STRICTER pair-diagonal restriction (m_a = m_c AND
     m_b = m_d per pair), applicable only to axially-symmetric / m-decoupled
     bases.  Paper values: 7.8125 / 2.7587 / 1.4404 / 0.9024 / 0.6190 %.
     This is the density the composed pipeline realizes (Table I).

This test pins BOTH columns to the paper's stated values by independent exact
enumeration, EACH via its own correct convention:

  * D   via GLOBAL m-conservation + the full bottom-row Gaunt 3j (q = m - m').
  * D_pd via the pair-diagonal restriction (q = 0 on each pair).

Two enumerations are run and cross-checked:
  1. an exact sympy.wigner_3j zero/nonzero decision (no floats) — the
     authoritative count, run at l_max=1..3 fast and l_max=4..5 slow;
  2. the production float Gaunt coefficient geovac.casimir_ci._gaunt_ck — the
     same routine geovac.casimir_ci.two_electron_integral uses to build the
     physical ERIs — run at l_max=1..5 (fast), tying the claim to production
     code.

Potential-independence (Theorem 2) is additionally asserted at l_max<=2 by
reusing the production grid solver in geovac.nuclear.potential_sparsity: the
angular zero pattern is bit-identical across Coulomb / harmonic / Woods-Saxon /
square-well / Yukawa bases.

Provenance: SYMBOLIC PROOF (exact integer nonzero counts via sympy 3j) for the
density columns; PANEL-VERIFIED (five potentials, bit-identical) for
potential-independence.

NOTE (FLAG, do not "fix" here): the production enumerator
geovac.nuclear.potential_sparsity.angular_zero_count computes D_pd, not the
headline global-M_L D, despite its docstring describing global m-conservation.
This is documented as a regression assertion in
test_angular_zero_count_computes_pair_diagonal_not_global below and corresponds
to the open group4 carry-forward CF-1.  The composed pipeline may intend the
pair-diagonal density, so this is flagged, not changed.
"""
from __future__ import annotations

import numpy as np
import pytest
from sympy.physics.wigner import wigner_3j

from geovac.casimir_ci import _gaunt_ck


# ---------------------------------------------------------------------------
# Reference values from Paper 22 Theorem 3 (exact rational nonzero counts and
# their 4-decimal percentage renderings).
# ---------------------------------------------------------------------------

# l_max -> (N_orb^4 total, exact nonzero count, percentage string in the paper)
GLOBAL_D = {
    1: (256, 38, "14.8438"),
    2: (6561, 559, "8.5200"),
    3: (65536, 3972, "6.0608"),
    4: (390625, 18857, "4.8274"),
    5: (1679616, 66954, "3.9863"),
}

PAIR_DIAGONAL_D = {
    1: (256, 20, "7.8125"),
    2: (6561, 181, "2.7587"),
    3: (65536, 944, "1.4404"),
    4: (390625, 3525, "0.9024"),
    5: (1679616, 10396, "0.6190"),
}


# ---------------------------------------------------------------------------
# Single-shell angular index set and the two exact enumerations.
# ---------------------------------------------------------------------------

def _shell(l_max: int):
    """The single-shell angular orbital set {(l, m) : l <= l_max}."""
    return [(l, m) for l in range(l_max + 1) for m in range(-l, l + 1)]


def _global_nonzero(la, ma, lc, mc, lb, mb, ld, md, *, use_float: bool) -> bool:
    """True iff <ab|cd> survives the PHYSICAL Coulomb (global-M_L) Gaunt rule.

    Requires global m-conservation m_a + m_b = m_c + m_d, then a shared
    multipole k for which both pairs carry a nonzero Gaunt coefficient with
    the bottom-row magnetic projection q_pair = m - m' (q_ac + q_bd = 0 by
    m-conservation).
    """
    if ma + mb != mc + md:
        return False
    k_min = max(abs(la - lc), abs(lb - ld))
    k_max = min(la + lc, lb + ld)
    for k in range(k_min, k_max + 1):
        if (la + lc + k) % 2 or (lb + ld + k) % 2:
            continue
        if use_float:
            if abs(_gaunt_ck(la, ma, lc, mc, k)) < 1e-15:
                continue
            if abs(_gaunt_ck(lb, mb, ld, md, k)) < 1e-15:
                continue
        else:
            if wigner_3j(la, k, lc, 0, 0, 0) == 0:
                continue
            if wigner_3j(lb, k, ld, 0, 0, 0) == 0:
                continue
            if wigner_3j(la, k, lc, -ma, ma - mc, mc) == 0:
                continue
            if wigner_3j(lb, k, ld, -mb, mb - md, md) == 0:
                continue
        return True
    return False


def _pair_diagonal_nonzero(la, ma, lc, mc, lb, mb, ld, md, *, use_float: bool) -> bool:
    """True iff <ab|cd> survives the STRICTER pair-diagonal rule (m_a=m_c, m_b=m_d)."""
    if ma != mc or mb != md:
        return False
    k_min = max(abs(la - lc), abs(lb - ld))
    k_max = min(la + lc, lb + ld)
    for k in range(k_min, k_max + 1):
        if (la + lc + k) % 2 or (lb + ld + k) % 2:
            continue
        if use_float:
            if abs(_gaunt_ck(la, ma, lc, mc, k)) < 1e-15:
                continue
            if abs(_gaunt_ck(lb, mb, ld, md, k)) < 1e-15:
                continue
        else:
            if wigner_3j(la, k, lc, 0, 0, 0) == 0:
                continue
            if wigner_3j(lb, k, ld, 0, 0, 0) == 0:
                continue
            if wigner_3j(la, k, lc, -ma, 0, mc) == 0:
                continue
            if wigner_3j(lb, k, ld, -mb, 0, md) == 0:
                continue
        return True
    return False


def _count(l_max: int, predicate, *, use_float: bool):
    """Count surviving four-index combinations over the single shell."""
    orbs = _shell(l_max)
    total = len(orbs) ** 4
    nz = 0
    for (la, ma) in orbs:
        for (lb, mb) in orbs:
            for (lc, mc) in orbs:
                for (ld, md) in orbs:
                    if predicate(la, ma, lc, mc, lb, mb, ld, md, use_float=use_float):
                        nz += 1
    return nz, total


def _pct_string(nz: int, total: int) -> str:
    return f"{100.0 * nz / total:.4f}"


# ---------------------------------------------------------------------------
# Headline density D (global M_L) — exact sympy counts.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("l_max", [1, 2, 3])
def test_global_D_exact_sympy(l_max):
    """D(l_max): exact integer nonzero count + percentage, global-M_L rule, sympy 3j."""
    total_ref, nz_ref, pct_ref = GLOBAL_D[l_max]
    nz, total = _count(l_max, _global_nonzero, use_float=False)
    assert total == total_ref
    assert nz == nz_ref, f"global D nonzero count l_max={l_max}: got {nz}, want {nz_ref}"
    assert _pct_string(nz, total) == pct_ref


@pytest.mark.slow
@pytest.mark.parametrize("l_max", [4, 5])
def test_global_D_exact_sympy_high_lmax(l_max):
    """D(l_max=4,5): exact integer nonzero count + percentage (slow, sympy 3j)."""
    total_ref, nz_ref, pct_ref = GLOBAL_D[l_max]
    nz, total = _count(l_max, _global_nonzero, use_float=False)
    assert total == total_ref
    assert nz == nz_ref
    assert _pct_string(nz, total) == pct_ref


@pytest.mark.parametrize("l_max", [1, 2, 3, 4, 5])
def test_global_D_production_gaunt(l_max):
    """D(l_max=1..5) via the production casimir_ci._gaunt_ck float coefficient.

    Ties the headline density to the exact routine that builds the physical
    two_electron_integral, across the full l_max range stated in the paper.
    """
    total_ref, nz_ref, pct_ref = GLOBAL_D[l_max]
    nz, total = _count(l_max, _global_nonzero, use_float=True)
    assert total == total_ref
    assert nz == nz_ref
    assert _pct_string(nz, total) == pct_ref


# ---------------------------------------------------------------------------
# Pair-diagonal density D_pd — exact sympy counts.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("l_max", [1, 2, 3])
def test_pair_diagonal_D_exact_sympy(l_max):
    """D_pd(l_max): exact integer nonzero count + percentage, pair-diagonal rule."""
    total_ref, nz_ref, pct_ref = PAIR_DIAGONAL_D[l_max]
    nz, total = _count(l_max, _pair_diagonal_nonzero, use_float=False)
    assert total == total_ref
    assert nz == nz_ref, f"D_pd nonzero count l_max={l_max}: got {nz}, want {nz_ref}"
    assert _pct_string(nz, total) == pct_ref


@pytest.mark.slow
@pytest.mark.parametrize("l_max", [4, 5])
def test_pair_diagonal_D_exact_sympy_high_lmax(l_max):
    """D_pd(l_max=4,5): exact integer nonzero count + percentage (slow, sympy 3j)."""
    total_ref, nz_ref, pct_ref = PAIR_DIAGONAL_D[l_max]
    nz, total = _count(l_max, _pair_diagonal_nonzero, use_float=False)
    assert total == total_ref
    assert nz == nz_ref
    assert _pct_string(nz, total) == pct_ref


@pytest.mark.parametrize("l_max", [1, 2, 3, 4, 5])
def test_pair_diagonal_D_production_gaunt(l_max):
    """D_pd(l_max=1..5) via the production casimir_ci._gaunt_ck float coefficient."""
    total_ref, nz_ref, pct_ref = PAIR_DIAGONAL_D[l_max]
    nz, total = _count(l_max, _pair_diagonal_nonzero, use_float=True)
    assert total == total_ref
    assert nz == nz_ref
    assert _pct_string(nz, total) == pct_ref


# ---------------------------------------------------------------------------
# The two densities are genuinely distinct (D ~ 4x D_pd) — guards against the
# two conventions silently collapsing.
# ---------------------------------------------------------------------------

def test_global_strictly_denser_than_pair_diagonal():
    """D > D_pd at every l_max >= 1, and the pair-diagonal set is a subset of D."""
    for l_max in (1, 2, 3):
        _, nz_g, _ = GLOBAL_D[l_max]
        _, nz_pd, _ = PAIR_DIAGONAL_D[l_max]
        assert nz_g > nz_pd
        # Pair-diagonal survivors must also survive the global rule (subset).
        orbs = _shell(l_max)
        for (la, ma) in orbs:
            for (lb, mb) in orbs:
                for (lc, mc) in orbs:
                    for (ld, md) in orbs:
                        if _pair_diagonal_nonzero(la, ma, lc, mc, lb, mb, ld, md,
                                                  use_float=False):
                            assert _global_nonzero(la, ma, lc, mc, lb, mb, ld, md,
                                                   use_float=False)


# ---------------------------------------------------------------------------
# Potential-independence (Theorem 2) at l_max <= 2, reusing the production
# grid solver.  The angular zero PATTERN is bit-identical across five
# potentials; radial integrals only reweight surviving entries.
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_potential_independence_lmax2():
    """Theorem 2: angular zero pattern bit-identical across five potentials.

    Reuses geovac.nuclear.potential_sparsity (production grid solver) at
    n_max=1, l_max=2 across Coulomb / harmonic / Woods-Saxon / square-well /
    Yukawa.  Asserts the nonzero ERI MASK is bit-identical across all five
    (the paper's "every element zero for Coulomb is zero for all others"
    claim) and that the realized density is uniform across potentials.

    The comparison is taken over the COMMON set of (n_r, l) orbitals that all
    five potentials bind (the shallow default Yukawa well does not support the
    (n_r=1, l=0) state, so it carries 17 orbitals vs 18 for the others); the
    bit-identical-mask statement is the actual content of Theorem 2 and holds
    at matched orbitals, which is precisely what restricting to the common
    set enforces.
    """
    from geovac.nuclear.potential_sparsity import (
        DEFAULT_POTENTIALS,
        compute_eri_tensor,
        radial_wavefunctions_for_potential,
    )

    n_max, l_max = 1, 2
    wfs_all = {
        name: radial_wavefunctions_for_potential(name, params, n_max, l_max)
        for name, params in DEFAULT_POTENTIALS.items()
    }
    common_keys = set.intersection(*(set(w.keys()) for w in wfs_all.values()))
    assert len(common_keys) >= 5, common_keys

    masks = {}
    densities = {}
    for name, wfs in wfs_all.items():
        wfs_common = {k: v for k, v in wfs.items() if k in common_keys}
        eri, stats = compute_eri_tensor(wfs_common, threshold=1e-10)
        masks[name] = (np.abs(eri) > 1e-10)
        densities[name] = stats["eri_density_pct"]

    ref_mask = masks["coulomb"]
    for name, mask in masks.items():
        assert mask.shape == ref_mask.shape, name
        assert np.array_equal(mask, ref_mask), (
            f"angular zero pattern for {name} differs from coulomb"
        )

    # Radial zeros contribute nothing: the realized density is identical
    # across all five potentials (bit-identical to the Coulomb value).
    ref_density = densities["coulomb"]
    for name, dens in densities.items():
        assert abs(dens - ref_density) < 1e-9, (
            f"{name}: density {dens}% != coulomb {ref_density}%"
        )


# ---------------------------------------------------------------------------
# FLAG (regression assertion, NOT a fix): the production enumerator
# angular_zero_count computes D_pd, while its docstring + the Paper 22 footnote
# describe the global-M_L D.  This is the open group4 carry-forward CF-1.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("l_max", [1, 2, 3])
def test_angular_zero_count_computes_pair_diagonal_not_global(l_max):
    """FLAG/CF-1: angular_zero_count returns D_pd despite a global-m docstring.

    geovac.nuclear.potential_sparsity.angular_zero_count enforces
    m_a+m_b = m_c+m_d textually, but its (a,c)/(b,d) k-set intersection drops
    the q-pairing that a true global-M_L sum requires, so the count collapses
    to the pair-diagonal D_pd.  This test documents the CURRENT behavior (it
    equals D_pd, NOT the headline global D) so any future change to that
    routine is caught.  Do not "fix" angular_zero_count here: the composed
    pipeline may intend the pair-diagonal density (CF-1 decision pending).
    """
    from geovac.nuclear.potential_sparsity import angular_zero_count

    zeros, total = angular_zero_count(n_max=0, l_max=l_max)
    nz = total - zeros
    total_pd, nz_pd, _ = PAIR_DIAGONAL_D[l_max]
    total_g, nz_g, _ = GLOBAL_D[l_max]

    assert total == total_pd == total_g
    # It matches the PAIR-DIAGONAL count ...
    assert nz == nz_pd, (
        f"angular_zero_count l_max={l_max}: got {nz}, "
        f"expected pair-diagonal {nz_pd}"
    )
    # ... and is strictly below the headline global-M_L count it documents.
    assert nz < nz_g
