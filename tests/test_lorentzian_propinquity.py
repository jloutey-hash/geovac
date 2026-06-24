"""Tests for the K^+-restricted weak-form Lorentzian propinquity assembly.

Sprint L3b-2 Sub-Sprint D (2026-05-18): verification suite for the
K^+-restricted weak-form assembly on the compact-temporal Lorentzian
truncated Krein spectral triple. (The "propinquity convergence theorem"
framing is RETRACTED -- Paper 45 Remark rem:history; the K^+ compression
annihilates the spatial Dirac. The c3_joint_panel_sup assertions below
test the WITHDRAWN sqrt-envelope arithmetic, not the L3 constant; the
actual L3 comparison constant is C_3 = 1, Paper 38 L3, with the
operator-system refutation in test_paper46_c3_operator_system.py.)

Per CLAUDE.md §13.4a, every equation in the proof memo has a
corresponding unit test.  Coverage:

  - LorentzianTunnelingPair construction at panel cells
  - K^+ restriction correctness (commutativity with J)
  - reach_B, height_B, height_P bounds
  - Main theorem bound at panel cells
  - Riemannian-limit recovery at N_t = 1 (load-bearing falsifier)
  - Monotone decrease in n_max at fixed N_t
  - Asymptotic rate consistency
  - Five-Lemma roadmap status
"""

from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from geovac.joint_berezin_compact_temporal import (
    joint_axisymmetric_positive,
    joint_constant_function,
    joint_panel,
    joint_separable_single_mode,
)
from geovac.lorentzian_propinquity_compact_temporal import (
    C_LIPSCHITZ_JOINT_ASYMPTOTIC,
    LorentzianFiveLemmaStatus,
    LorentzianPropinquityBound,
    LorentzianTunnelingPair,
    asymptotic_rate_ratio,
    c3_joint_panel_sup,
    compute_lorentzian_propinquity_bound,
    lorentzian_gh_convergence_table,
    lorentzian_theorem_statement,
    verify_monotone_decrease_in_n_max,
    verify_riemannian_limit_at_N_t_1,
)


# ---------------------------------------------------------------------------
# Constants and helpers
# ---------------------------------------------------------------------------


class TestConstants:
    """Tests for module-level constants."""

    def test_C_LIPSCHITZ_JOINT_ASYMPTOTIC_is_one(self):
        """The joint Lichnerowicz constant asymptote is 1 (Paper 38 §L3 inherited)."""
        assert C_LIPSCHITZ_JOINT_ASYMPTOTIC == 1.0

    # NOTE: c3_joint_panel_sup computes the WITHDRAWN sqrt-envelope formula
    # (NOT the L3 comparison constant, which is C_3 = 1; Paper 38 L3). These
    # tests pin only the withdrawn-formula arithmetic for documentation; the
    # operator-system refutation of the sqrt-story is in
    # test_paper46_c3_operator_system.py.
    def test_c3_joint_panel_sup_n_max_2_withdrawn_formula(self):
        """WITHDRAWN formula arithmetic: sup_{N<=2} (N-1)/sqrt(N^2-1) = 1/sqrt(3)."""
        val = c3_joint_panel_sup(2)
        expected = 1.0 / np.sqrt(3.0)
        assert abs(val - expected) < 1e-12

    def test_c3_joint_panel_sup_n_max_3_withdrawn_formula(self):
        """WITHDRAWN formula arithmetic: sup = 2/sqrt(8) = 1/sqrt(2)."""
        val = c3_joint_panel_sup(3)
        expected = 2.0 / np.sqrt(8.0)
        assert abs(val - expected) < 1e-12

    def test_c3_joint_panel_sup_monotone_withdrawn_formula(self):
        """WITHDRAWN formula arithmetic: monotone in n_max. (NOT the L3 constant,
        which is the flat C_3 = 1; this sqrt-sup is op-norm-false as a constant.)"""
        vals = [c3_joint_panel_sup(n) for n in [2, 3, 4, 5, 6, 10, 100]]
        for a, b in zip(vals[:-1], vals[1:]):
            assert b > a, f"non-monotone: {a} >= {b}"
        assert vals[-1] < 1.0
        assert vals[-1] > 0.99


# ---------------------------------------------------------------------------
# LorentzianTunnelingPair construction
# ---------------------------------------------------------------------------


class TestTunnelingPairConstruction:
    """Tests for LorentzianTunnelingPair.build."""

    def test_build_n_max_2_N_t_3(self):
        """Build at the (2, 3) panel cell; verify metadata."""
        pair = LorentzianTunnelingPair.build(n_max=2, N_t=3)
        assert pair.n_max == 2
        assert pair.N_t == 3
        assert pair.T == pytest.approx(2 * np.pi)
        assert pair.cb_norm_joint == sp.Rational(2, 3)
        # c_lipschitz_joint is the L3 comparison constant C_3 = 1 (Paper 38 L3,
        # gradient seminorm). The withdrawn sqrt-envelope formula is no longer
        # wired here (op-norm-false; see c3_joint_panel_sup docstring).
        assert pair.c_lipschitz_joint == pytest.approx(1.0)
        assert pair.gamma_joint_su2 > 0
        assert pair.gamma_joint_u1 > 0
        assert pair.gamma_joint_L1 > pair.gamma_joint_su2

    def test_build_n_max_3_N_t_5(self):
        """Build at (3, 5); cb_norm = 1/2."""
        pair = LorentzianTunnelingPair.build(n_max=3, N_t=5)
        assert pair.cb_norm_joint == sp.Rational(1, 2)

    def test_build_invalid_n_max(self):
        with pytest.raises(ValueError):
            LorentzianTunnelingPair.build(n_max=0, N_t=3)

    def test_build_invalid_N_t(self):
        with pytest.raises(ValueError):
            LorentzianTunnelingPair.build(n_max=2, N_t=0)

    def test_build_invalid_T(self):
        with pytest.raises(ValueError):
            LorentzianTunnelingPair.build(n_max=2, N_t=3, T=-1.0)


# ---------------------------------------------------------------------------
# K^+ restriction correctness
# ---------------------------------------------------------------------------


class TestKPlusRestriction:
    """Tests for the K^+ commutativity property (Sub-sprint C §6)."""

    def test_K_plus_compatibility_constant(self):
        """B^joint(constant) commutes with J bit-exact."""
        pair = LorentzianTunnelingPair.build(n_max=2, N_t=3)
        f = joint_constant_function()
        preserves, residual = pair.verify_K_plus_compatibility(f)
        assert preserves
        assert residual < 1e-10

    def test_K_plus_compatibility_single_harmonic(self):
        """B^joint(Y_{2,0,0} ⊗ e^{iqt}) commutes with J bit-exact."""
        pair = LorentzianTunnelingPair.build(n_max=2, N_t=3)
        f = joint_separable_single_mode(2, 0, 0, 0)
        preserves, residual = pair.verify_K_plus_compatibility(f)
        assert preserves
        assert residual < 1e-10

    def test_K_plus_compatibility_panel(self):
        """K^+ preservation at every entry of the joint panel (Sub-sprint C §6)."""
        pair = LorentzianTunnelingPair.build(n_max=2, N_t=3)
        for f in joint_panel(2, 3):
            preserves, residual = pair.verify_K_plus_compatibility(f)
            assert preserves, f"K^+ failed for {f.name}, residual={residual}"

    def test_krein_state_space_J_eigendecomp(self):
        """The K^+ state space has J = +I on K^+ subspace, J = -I on K^-."""
        pair = LorentzianTunnelingPair.build(n_max=2, N_t=3)
        ok, details = pair.krein_state_space.verify_J_eigendecomp()
        assert ok
        assert details["K_plus_dim"] + details["K_minus_dim"] == details["dim_K"]


# ---------------------------------------------------------------------------
# Constituent bounds: reach_B, height_B, height_P
# ---------------------------------------------------------------------------


class TestConstituentBounds:
    """Tests for the four propinquity constituents (memo §3)."""

    def test_height_P_zero(self):
        """P^joint is an orthogonal projection: height_P = 0 exactly."""
        pair = LorentzianTunnelingPair.build(n_max=2, N_t=3)
        assert pair.height_P() == 0.0

    def test_reach_B_finite_and_nonneg(self):
        """reach_B(f) is finite and non-negative for every f."""
        pair = LorentzianTunnelingPair.build(n_max=2, N_t=3)
        for f in joint_panel(2, 3):
            r = pair.reach_B(f)
            assert r >= 0.0
            assert np.isfinite(r)

    def test_height_B_finite_and_nonneg(self):
        """height_B(f) is finite and non-negative for every f."""
        pair = LorentzianTunnelingPair.build(n_max=2, N_t=3)
        for f in joint_panel(2, 3):
            h = pair.height_B(f)
            assert h >= 0.0
            assert np.isfinite(h)


# ---------------------------------------------------------------------------
# Main theorem panel
# ---------------------------------------------------------------------------


class TestMainTheoremPanel:
    """Tests for compute_lorentzian_propinquity_bound at panel cells."""

    def test_bound_at_2_3(self):
        """Lambda^L bound at (n_max, N_t) = (2, 3) is finite."""
        b = compute_lorentzian_propinquity_bound(n_max=2, N_t=3)
        assert isinstance(b, LorentzianPropinquityBound)
        assert b.n_max == 2
        assert b.N_t == 3
        assert b.propinquity_bound > 0
        assert np.isfinite(b.propinquity_bound)

    def test_bound_at_3_5(self):
        """Lambda^L bound at (3, 5) is finite."""
        b = compute_lorentzian_propinquity_bound(n_max=3, N_t=5)
        assert b.propinquity_bound > 0

    def test_bound_qualitative_rate_only(self):
        """Sub-sprint D bound is qualitative-rate (Track C constant not yet)."""
        b = compute_lorentzian_propinquity_bound(n_max=2, N_t=3)
        assert b.qualitative_rate_only

    def test_bound_cb_norm_matches(self):
        """The Lambda^L bound carries the correct cb-norm = 2/(n_max+1)."""
        b = compute_lorentzian_propinquity_bound(n_max=3, N_t=5)
        assert b.cb_norm_joint == pytest.approx(2.0 / 4.0)  # 1/2

    def test_bound_to_dict_serializable(self):
        """The bound serializes to a JSON-compatible dict."""
        b = compute_lorentzian_propinquity_bound(n_max=2, N_t=3)
        d = b.to_dict()
        assert "propinquity_bound" in d
        assert "n_max" in d
        assert "N_t" in d
        # Check JSON-serializable
        import json
        json_str = json.dumps(d, default=str)
        assert len(json_str) > 0


# ---------------------------------------------------------------------------
# Riemannian-limit recovery (load-bearing falsifier)
# ---------------------------------------------------------------------------


class TestRiemannianLimitRecovery:
    """Load-bearing falsifier: at N_t = 1 the bound reduces to Paper 38 bit-exactly."""

    def test_riemannian_limit_n_max_2(self):
        """At n_max=2, N_t=1: gamma matches Paper 38 bit-exact."""
        match, details = verify_riemannian_limit_at_N_t_1(n_max=2)
        assert match
        assert details["gamma_residual"] == 0.0  # bit-exact
        assert details["cb_match"]

    def test_riemannian_limit_n_max_3(self):
        """At n_max=3, N_t=1: gamma matches Paper 38 bit-exact."""
        match, details = verify_riemannian_limit_at_N_t_1(n_max=3)
        assert match
        assert details["gamma_residual"] == 0.0
        assert details["cb_match"]

    def test_riemannian_limit_n_max_4(self):
        """At n_max=4, N_t=1: gamma matches Paper 38 bit-exact."""
        match, details = verify_riemannian_limit_at_N_t_1(n_max=4)
        assert match
        assert details["gamma_residual"] == 0.0
        assert details["cb_match"]

    def test_riemannian_limit_pair_method(self):
        """The pair-level reduces_to_paper38_at_N_t_1 also returns True at N_t=1."""
        pair = LorentzianTunnelingPair.build(n_max=2, N_t=1)
        match, details = pair.reduces_to_paper38_at_N_t_1()
        assert match
        assert details["gamma_residual_F"] == 0.0


# ---------------------------------------------------------------------------
# Convergence properties
# ---------------------------------------------------------------------------


class TestConvergence:
    """Tests for monotone decrease and asymptotic rate consistency."""

    def test_monotone_decrease_at_N_t_1(self):
        """Lambda^L decreases monotonically with n_max at fixed N_t=1."""
        bounds = lorentzian_gh_convergence_table([(2, 1), (3, 1), (4, 1)])
        is_monotone, violations = verify_monotone_decrease_in_n_max(bounds)
        assert is_monotone, f"violations: {violations}"

    def test_asymptotic_rate_ratio_decreases(self):
        """Lambda^L at higher n_max is smaller than at lower n_max."""
        # Use (2, 1) and (4, 1) to keep memory budget reasonable
        bounds = lorentzian_gh_convergence_table([(2, 1), (4, 1)])
        ratio = asymptotic_rate_ratio(bounds, (2, 1), (4, 1))
        assert ratio < 1.0  # the bound decreases

    def test_asymptotic_rate_ratio_invalid_cell(self):
        """asymptotic_rate_ratio raises KeyError for unknown cells."""
        bounds = lorentzian_gh_convergence_table([(2, 3)])
        with pytest.raises(KeyError):
            asymptotic_rate_ratio(bounds, (2, 3), (99, 99))


# ---------------------------------------------------------------------------
# Five-lemma roadmap status
# ---------------------------------------------------------------------------


class TestFiveLemmaStatus:
    """Tests for the L3b-2 roadmap status."""

    def test_all_done(self):
        """All five lemmas L1', L2, L3, L4, L5 are DONE."""
        status = LorentzianFiveLemmaStatus()
        assert status.all_done()

    def test_status_to_dict(self):
        """to_dict returns a dict with all lemma keys."""
        status = LorentzianFiveLemmaStatus()
        d = status.to_dict()
        assert "L1_prime" in d
        assert "L2" in d
        assert "L3" in d
        assert "L4" in d
        assert "L5" in d
        assert d["all_done"]

    def test_status_l5_2026_05_18(self):
        """L5 is dated 2026-05-18 (this sprint)."""
        status = LorentzianFiveLemmaStatus()
        assert "2026-05-18" in status.L5
        assert "Sub-sprint D" in status.L5


# ---------------------------------------------------------------------------
# Theorem statement
# ---------------------------------------------------------------------------


class TestTheoremStatement:
    """Tests for the CURRENT Paper 45 headline statement (the K^+ annihilation
    theorem). Guards against the WITHDRAWN weak-form convergence theorem
    re-entering the returned string (the v4.43.5 code-side-zombie pattern)."""

    def test_statement_is_the_annihilation_theorem(self):
        """The current headline is the K^+ annihilation (degeneracy) theorem."""
        stmt = lorentzian_theorem_statement()
        assert "annihilat" in stmt.lower()
        assert "vanishes identically" in stmt or "P_+ (D_GV (x) I) P_+ = 0" in stmt
        assert "DEGENERATE" in stmt or "degenerate" in stmt.lower()
        # the surviving convergence is van Suijlekom STATE-SPACE GH, not propinquity
        assert "state-space" in stmt.lower()

    def test_withdrawn_convergence_claim_is_flagged_not_live(self):
        """If the retracted weak-form convergence form is mentioned, it must
        carry the WITHDRAWN/retracted flag -- never presented as live."""
        stmt = lorentzian_theorem_statement()
        if "C_3^joint" in stmt or "-> 0" in stmt:
            assert "WITHDRAWN" in stmt or "retracted" in stmt.lower()

    def test_statement_names_honest_scope(self):
        """The statement explicitly names the strong-form open gap."""
        stmt = lorentzian_theorem_statement()
        assert "Honest scope" in stmt or "strong-form" in stmt


# ---------------------------------------------------------------------------
# Slow tests (full panel at larger n_max)
# ---------------------------------------------------------------------------


@pytest.mark.slow
class TestFullPanel:
    """Slow tests on the full Sub-sprint D panel."""

    def test_full_panel_at_2_3_3_5_4_7(self):
        """Compute Lambda^L at three panel cells."""
        bounds = lorentzian_gh_convergence_table([(2, 3), (3, 5), (4, 7)])
        for (n_max, N_t), b in bounds.items():
            assert b.propinquity_bound > 0
            assert np.isfinite(b.propinquity_bound)

    def test_full_panel_monotone_su2(self):
        """gamma_su2 decreases monotonically with n_max (Paper 38 § L2)."""
        bounds = lorentzian_gh_convergence_table([(2, 3), (3, 5), (4, 7)])
        gammas = [bounds[c].gamma_joint_su2 for c in [(2, 3), (3, 5), (4, 7)]]
        assert gammas[0] > gammas[1] > gammas[2]
