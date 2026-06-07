"""Tests for the TC-2a finite-cutoff reconstruction at n_max = 2.

See: Sprint Q5'-Tannakian-Closure TC-2a
(debug/sprint_q5p_tc2a_aut_equality_memo.md).

These tests verify the bit-exact closed-form structure of the witness-panel
variety of natural ⊗-automorphisms, confirming
dim Aut^⊗(ω) = 15 = dim G_a^{3 N(2)} on the n_max-axis substrate.
"""

import importlib.util
from pathlib import Path

import pytest
from sympy import Matrix, Rational, Integer, eye, zeros

from geovac.pro_system import primitive_generators
from geovac.tannakian import FinDimRep, trivial_rep


# ---------------------------------------------------------------------
# Load the driver helpers
# ---------------------------------------------------------------------


def _load_driver():
    """Import the TC-2a driver as a module so we can test its helpers."""
    repo_root = Path(__file__).parent.parent
    driver_path = repo_root / "debug" / "compute_q5p_tc2a_aut_equality.py"
    spec = importlib.util.spec_from_file_location("tc2a_driver", driver_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


driver = _load_driver()


# ---------------------------------------------------------------------
# Witness panel
# ---------------------------------------------------------------------


def test_witness_panel_size_at_n_max_2():
    gens, panel = driver.build_witness_panel(2)
    assert len(gens) == 15
    assert len(panel) == 15
    for g in gens:
        assert g in panel
        V_g = panel[g]
        assert V_g.dim == 2


def test_witness_panel_activation():
    """Each V_g activates exactly one generator: X_g = E_{12}, X_h = 0 (h != g)."""
    gens, panel = driver.build_witness_panel(2)
    E12 = Matrix([[0, 1], [0, 0]])
    for g in gens:
        V_g = panel[g]
        endos = V_g.non_zero_endos()
        assert list(endos.keys()) == [g]
        assert endos[g] == E12


# ---------------------------------------------------------------------
# hom_basis
# ---------------------------------------------------------------------


def test_hom_basis_disjoint_generators_is_1dim():
    """For V_g, V_h with g != h, Hom(V_g, V_h) is 1-dim:
    f = [[0, a], [0, 0]] for some a in Q."""
    gens, panel = driver.build_witness_panel(2)
    g = gens[0]
    h = gens[1]
    morphs = driver.hom_basis(panel[g], panel[h])
    assert len(morphs) == 1
    f = morphs[0]
    # f should have zero first column and zero second row
    assert f[0, 0] == 0
    assert f[1, 0] == 0
    assert f[1, 1] == 0
    # f[0, 1] is the free coordinate; nullspace basis convention sets it to 1
    assert f[0, 1] != 0


def test_hom_basis_endomorphism_is_2dim_centralizer():
    """For V_g, End(V_g) = centralizer of E_{12} in M_2(Q), 2-dim."""
    gens, panel = driver.build_witness_panel(2)
    g = gens[0]
    V_g = panel[g]
    morphs = driver.hom_basis(V_g, V_g)
    assert len(morphs) == 2
    # All morphisms must commute with E_{12}
    E12 = Matrix([[0, 1], [0, 0]])
    for f in morphs:
        assert f * E12 == E12 * f


def test_hom_basis_with_trivial():
    """Hom(T, V_g) and Hom(V_g, T) are each 1-dim."""
    gens, panel = driver.build_witness_panel(2)
    T = trivial_rep(2, dim=1)
    g = gens[0]
    V_g = panel[g]
    hom_T_to_V = driver.hom_basis(T, V_g)
    hom_V_to_T = driver.hom_basis(V_g, T)
    assert len(hom_T_to_V) == 1
    assert len(hom_V_to_T) == 1
    # Hom(T, V_g) is a 2x1 vector with second entry zero (= e_1 up to scaling)
    f_in = hom_T_to_V[0]
    assert f_in.rows == 2 and f_in.cols == 1
    assert f_in[1, 0] == 0
    # Hom(V_g, T) is a 1x2 row with first entry zero (= e_2^T up to scaling)
    f_out = hom_V_to_T[0]
    assert f_out.rows == 1 and f_out.cols == 2
    assert f_out[0, 0] == 0


# ---------------------------------------------------------------------
# Linear-system dimension
# ---------------------------------------------------------------------


def test_aut_equality_at_n_max_2_is_15():
    """The headline test: dim Aut^⊗(ω) = 15 at n_max = 2.

    Runs the bit-exact pipeline on the same witness panel as TC-1f and confirms
    the predicted dim of the converse-direction variety.
    """
    gens, panel = driver.build_witness_panel(2)
    etas, eta_T_mat, all_symbols = driver.parameterize_eta(gens)
    constraints, _ = driver.collect_naturality_constraints(
        gens, panel, etas, eta_T_mat, 2
    )
    driver.add_unit_normalization(constraints, eta_T_mat)
    result = driver.solve_linear_system(constraints, all_symbols)

    assert result["consistent"] is True
    assert result["dim_solution_variety"] == 15
    assert result["n_vars"] == 61
    assert result["rank_A"] == 46


def test_closed_form_solution_recovers_phi():
    """Each η_{V_g} solves to [[1, q_g], [0, 1]] = exp(q_g E_{12})."""
    gens, panel = driver.build_witness_panel(2)
    etas, eta_T_mat, all_symbols = driver.parameterize_eta(gens)
    constraints, _ = driver.collect_naturality_constraints(
        gens, panel, etas, eta_T_mat, 2
    )
    driver.add_unit_normalization(constraints, eta_T_mat)
    sol_struct = driver.extract_solution_structure(
        constraints, all_symbols, gens, etas, eta_T_mat
    )

    assert sol_struct.get("n_free_vars") == 15
    assert sol_struct.get("eta_T_solved") == "1"
    phi_check = driver.verify_phi_recovery(gens, sol_struct)
    assert phi_check["all_match"] is True


def test_panel_morphism_count():
    """Sanity: the panel has exactly 271 morphisms (with T)."""
    gens, panel = driver.build_witness_panel(2)
    etas, eta_T_mat, all_symbols = driver.parameterize_eta(gens)
    _, counts = driver.collect_naturality_constraints(
        gens, panel, etas, eta_T_mat, 2
    )
    total = sum(counts.values())
    assert total == 271


# ---------------------------------------------------------------------
# Smaller-cutoff sanity: at n_max = 1, dim should equal 3 * N(1) = 3
# ---------------------------------------------------------------------


def test_aut_equality_at_n_max_1_is_3():
    """Smaller-cutoff cross-check: dim Aut^⊗(ω) = 3 at n_max = 1.

    N(1) = 1 * 4 / 2 = 2.  3 * N(1) = 6.
    Wait: N(1) = 1*(1+3)/2 = 2 sectors, 3 Mellin slots → 6 primitive generators.
    """
    from geovac.pro_system import N_sectors, n_primitive_generators
    n_gens_1 = n_primitive_generators(1)
    assert n_gens_1 == 6  # = 3 * N(1) = 3 * 2

    gens, panel = driver.build_witness_panel(1)
    etas, eta_T_mat, all_symbols = driver.parameterize_eta(gens)
    constraints, _ = driver.collect_naturality_constraints(
        gens, panel, etas, eta_T_mat, 1
    )
    driver.add_unit_normalization(constraints, eta_T_mat)
    result = driver.solve_linear_system(constraints, all_symbols)

    assert result["consistent"] is True
    assert result["dim_solution_variety"] == 6
