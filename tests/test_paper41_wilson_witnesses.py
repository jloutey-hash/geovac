# -*- coding: utf-8 -*-
"""Pinning tests for Paper 41's seven-witness Rule B Wilson U(1) verdict.

Backing drivers live in tests/wilson_rule_b_support/ (2026-07-03 durability
migration from debug/archive/{qed_arc,misc,chemistry_qc_arc}/ -- the wh7_support
pattern, CLAUDE.md v4.49.1).  Archived production JSON lives in
tests/wilson_rule_b_support/data/.

Pin strategy (per witness):
  W1 spectral dimension   -- RECOMPUTED (heat-kernel d_s + Weyl d_W, <1 s).
  W2 Migdal-Kadanoff RG   -- RECOMPUTED (k_eff, no fixed point, monotone flow).
  W3 Gaussian Wilson loop -- RECOMPUTED at n_max=2 (full enumeration,
                             deterministic); JSON band-pin at n_max=3,4,5.
  W4 strong-coupling      -- RECOMPUTED (sigma(beta) = -ln(I1/I0) positivity +
                             both asymptotic limits) + JSON cross-check.
  W5 monopole density     -- structural facts RECOMPUTED (k_delta=3, 60
                             elementary sites, P=44); production MC fit values
                             (c=9.40/8.95, R^2=0.981) pinned from the archived
                             seeded JSON; reduced-stat seeded MC smoke verifies
                             the exponential-suppression signal + determinism.
  W6 full-MC Wilson loops -- structural setup + archived seeded-JSON pins
                             (sigma_ens > 0 monotone; WEAK-PASS flags).  The
                             production MC (n_sample=2000 x 11 beta) is too
                             heavy to recompute in a test; the archived JSON
                             carries fixed seeds 42/43.
  W7 Polyakov loop        -- homological class RECOMPUTED (H^1 defect
                             E - V + 1 - rank(d_1): 0 at N_t=2, 1 at N_t>=3,
                             product-graph counts); confinement verdict pinned
                             from the archived seeded JSON; reduced-stat seeded
                             MC smoke verifies <P> ~ 0 + determinism.

Paper anchors: papers/group5_qed_gauge/paper_41_rule_b_wilson_u1.tex
(sec:spectral_dim, sec:mk_rg, sec:gaussian_wilson, sec:strong_coupling,
sec:monopole_density, sec:full_mc_wilson, sec:polyakov_loop,
sec:seven_witness).
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse.csgraph as csgraph

SUPPORT = Path(__file__).resolve().parent / "wilson_rule_b_support"
DATA = SUPPORT / "data"
sys.path.insert(0, str(SUPPORT))

from geovac.ihara_zeta_dirac import build_dirac_s3_graph  # noqa: E402

import xcwg_u1_wilson_rule_b_pilot as pilot  # noqa: E402
import xcwg_wilson_loop_scaling as wls  # noqa: E402
import xcwg_wls_v3 as wls_v3  # noqa: E402
import xcwg_strong_coupling_wilson as scw  # noqa: E402
import xcwg_nlo_character_expansion as nlo  # noqa: E402
import xcwg_monopole_density as mono  # noqa: E402
import xcwg_mk_blockspin_rule_b as mk  # noqa: E402
import xcwg_polyakov_loop as poly  # noqa: E402
import xcwg_rule_b_spectral_dim as sdim  # noqa: E402


def _load(name: str) -> dict:
    with open(DATA / name) as f:
        return json.load(f)


# =====================================================================
# Rule B graph structural facts (foundation for all seven witnesses)
# =====================================================================

# (V, E) for the Dirac Rule B graph, Paper 41 Tables (foundation + W1)
RULE_B_COUNTS = {2: (10, 20), 3: (28, 106), 4: (60, 312), 5: (110, 692)}


@pytest.mark.parametrize("n_max", [2, 3, 4, 5])
def test_rule_b_graph_counts_connected_bipartite(n_max):
    """V, E, beta_1 = E - V + 1, connectedness, and bipartiteness (Delta-l = +-1
    parity flip => even cycles only => the L=4 plaquette complex is consistent)."""
    A, labels, deg, _ = build_dirac_s3_graph(n_max, "B")
    V_exp, E_exp = RULE_B_COUNTS[n_max]
    V = int(A.shape[0])
    E = int(A.sum()) // 2
    assert V == V_exp
    assert E == E_exp
    n_comp, _ = csgraph.connected_components((A > 0).astype(int), directed=False)
    assert n_comp == 1
    # bipartite: BFS 2-coloring never conflicts
    color = -np.ones(V, dtype=int)
    color[0] = 0
    stack = [0]
    while stack:
        u = stack.pop()
        for v in np.nonzero(A[u])[0]:
            v = int(v)
            if color[v] == -1:
                color[v] = 1 - color[u]
                stack.append(v)
            else:
                assert color[v] != color[u], "Rule B graph must be bipartite"


@pytest.mark.parametrize("n_max,ker_L1", [(2, 11), (3, 79)])
def test_hodge_identity(n_max, ker_L1):
    """Paper 41 sec 'Hodge identity': nonzero spec(L0) = nonzero spec(L1),
    dim ker L0 = 1 (connected), dim ker L1 = beta_1.  Pilot recompute."""
    A, _, _, _ = build_dirac_s3_graph(n_max, "B")
    B_inc, _ = pilot.signed_incidence(A)
    h = pilot.hodge_identity_check(B_inc @ B_inc.T, B_inc.T @ B_inc)
    assert h["hodge_identity_passed"]
    assert h["kernel_dim_L0"] == 1
    assert h["kernel_dim_L1"] == ker_L1
    assert h["shared_max_diff"] < 1e-9


@pytest.mark.parametrize("n_max,n_plaq", [(2, 44), (3, 994)])
def test_l4_plaquette_counts(n_max, n_plaq):
    """Primitive non-backtracking 4-cycle (plaquette) counts feeding the
    Wilson 2-complex: 44 at n_max=2, 994 at n_max=3 (Paper 41 / archived JSON)."""
    A, _, _, _ = build_dirac_s3_graph(n_max, "B")
    adj = wls.adjacency_list(A)
    plaqs = list(wls_v3.enumerate_primitive_closed_walks_streaming(adj, 4))
    assert len(plaqs) == n_plaq


# =====================================================================
# Witness 1: spectral dimension (RECOMPUTED)
# =====================================================================

T_GRID = np.logspace(-2.5, 2.5, 200)  # production grid of the driver

# Paper 41 seven-witness table: d_s 1.86 -> 2.27 -> 2.54, d_W 2.48 -> 3.46 -> 3.72
RULE_B_DIMS = {3: (1.86, 2.48), 4: (2.27, 3.46), 5: (2.54, 3.72)}
SCALAR_DS = {3: 1.68, 4: 1.76, 5: 1.79}


@pytest.mark.parametrize("n_max", [3, 4, 5])
def test_witness1_spectral_dimension_rule_b(n_max):
    ds_exp, dw_exp = RULE_B_DIMS[n_max]
    L, A, V = sdim.build_rule_b_laplacian(n_max)
    c = sdim.count_components(A)
    ev = np.clip(np.linalg.eigvalsh(L), 0.0, None)
    _, P = sdim.heat_kernel_return_probability(ev, V, T_GRID)
    d_s = sdim.extract_spectral_dimension(T_GRID, P, c, V)["d_s_window"]
    d_W = sdim.weyl_dimension(ev, V, c)["d_W"]
    assert abs(d_s - ds_exp) < 0.02, f"d_s({n_max}) = {d_s}"
    assert abs(d_W - dw_exp) < 0.02, f"d_W({n_max}) = {d_W}"


def test_witness1_ds_monotone_toward_3_scalar_plateaus():
    """The witness content: Rule B d_s climbs monotonically toward 3 while the
    scalar Fock graph plateaus near 2."""
    ds_b, ds_s = [], []
    for n_max in (3, 4, 5):
        for build, acc in ((sdim.build_rule_b_laplacian, ds_b),
                           (sdim.build_scalar_laplacian, ds_s)):
            L, A, V = build(n_max)
            c = sdim.count_components(A)
            ev = np.clip(np.linalg.eigvalsh(L), 0.0, None)
            _, P = sdim.heat_kernel_return_probability(ev, V, T_GRID)
            acc.append(sdim.extract_spectral_dimension(T_GRID, P, c, V)["d_s_window"])
    assert ds_b[0] < ds_b[1] < ds_b[2] < 3.0
    for n_max, ds in zip((3, 4, 5), ds_s):
        assert abs(ds - SCALAR_DS[n_max]) < 0.02
    assert ds_s[2] < 2.0  # scalar plateaus below 2
    assert ds_b[2] > ds_s[2] + 0.5  # Rule B decisively above scalar


# =====================================================================
# Witness 2: Migdal-Kadanoff RG (RECOMPUTED)
# =====================================================================

def test_witness2_mk_k_eff_and_flow():
    """k_eff (mean plaquettes per edge) = 8.8 at n_max=2, 37.509... at n_max=3
    (archived JSON k_eff_canonical); MK flow monotone to beta=0; no fixed point
    in (0, 100); strong-coupling limit matched to 0.25% at beta=0.03 and
    weak-coupling limit to 2% at beta=200 (8.6% at beta=50 -- the O(1/beta)
    approach; matches the archived JSON beta_eff_table bit-for-bit)."""
    for n_max, k_exp in ((2, 8.8), (3, 994 * 4 / 106)):
        A, _, _, _ = build_dirac_s3_graph(n_max, "B")
        cycles = mk.primitive_4_plaquettes(A)
        _, stats = mk.plaquette_edge_incidence(cycles, A.shape[0])
        assert abs(stats["k_eff_mean"] - k_exp) < 1e-9
    k_eff = 8.8
    # monotone flow to zero from several starting couplings
    for beta_0 in (0.5, 1.0, 2.0, 5.0, 10.0):
        flow = mk.iterate_mk(beta_0, k_eff, n_steps=12)
        assert all(flow[i + 1] < flow[i] or flow[i] == 0.0
                   for i in range(len(flow) - 1))
        assert flow[-1] < 1e-10
    # no non-trivial fixed point (g(beta) = beta_eff - beta has no sign change)
    fp = mk.find_fixed_points(k_eff, beta_lo=1e-3, beta_hi=100.0, n_grid=301)
    assert fp["n_sign_changes"] == 0
    assert fp["fixed_points"] == []
    # strong-coupling limit: <= 0.25% at beta = 0.03 (paper sec MK asymptotics)
    assert abs(mk.beta_eff_mk(0.03, k_eff)
               / mk.strong_coupling_limit(0.03, k_eff) - 1.0) < 0.0025
    # weak-coupling limit: O(1/beta) approach -- 8.6% at beta=50 (the archived
    # JSON value), crossing <= 2% at beta = 200
    err_50 = abs(mk.beta_eff_mk(50.0, k_eff)
                 / mk.weak_coupling_limit(50.0, k_eff) - 1.0)
    assert abs(err_50 - 0.0861) < 0.001
    assert abs(mk.beta_eff_mk(200.0, k_eff)
               / mk.weak_coupling_limit(200.0, k_eff) - 1.0) < 0.02


# =====================================================================
# Witness 3: Gaussian Wilson-loop scaling (RECOMPUTED at n_max=2)
# =====================================================================

def test_witness3_gaussian_alpha_nmax2_recompute():
    """<S(L)> ~ L^alpha with alpha = 1.177 at n_max=2 (full enumeration --
    deterministic: 44/144/1412 walks at L=4/6/8).  Paper table row 1:
    0.25 (JSON) / 0.3843 / 0.5687, alpha = 1.18 +- 0.08, in [0.89, 1.18]."""
    A, _, _, _ = build_dirac_s3_graph(2, "B")
    V = int(A.shape[0])
    E = int(A.sum()) // 2
    _, edges, edge_idx = wls.signed_incidence(A)
    adj = wls.adjacency_list(A)
    plaqs = list(wls_v3.enumerate_primitive_closed_walks_streaming(adj, 4))
    _, K = wls.build_K_from_L4_plaquettes(plaqs, edge_idx, E)
    K_pinv = np.linalg.pinv(K)
    assert np.linalg.matrix_rank(K) == 11  # = beta_1 at n_max=2
    means = {}
    counts_exp = {4: 44, 6: 144, 8: 1412}
    means_exp = {4: 0.2500, 6: 0.38426, 8: 0.56870}
    for L in (4, 6, 8):
        r = wls_v3.measure_S(adj, edges, edge_idx, K_pinv, L, V, cap=10_000_000)
        assert not r["capped"]
        assert r["count"] == counts_exp[L]
        assert abs(r["mean"] - means_exp[L]) < 5e-4
        means[L] = r["mean"]
    fit = wls.fit_scaling_exponent([4, 6, 8], [means[L] for L in (4, 6, 8)])
    assert abs(fit["alpha"] - 1.1775) < 0.01
    assert 0.89 <= fit["alpha"] <= 1.18
    assert fit["r_squared"] > 0.99


def test_witness3_alpha_band_all_nmax_from_archive():
    """Archived consolidated JSON: alpha in [0.89, 1.18] at every
    n_max = 2, 3, 4, 5 (paper: 1.18 / 1.03 / 0.89 / 0.99)."""
    d = _load("xcwg_wilson_loop_scaling.json")
    alpha_exp = {"n_max_2": 1.18, "n_max_3": 1.03, "n_max_4": 0.89, "n_max_5": 0.99}
    for key, a_exp in alpha_exp.items():
        alpha = d[key]["scaling_fit_full"]["alpha"]
        assert 0.89 - 0.005 <= alpha <= 1.18 + 0.005, f"{key}: alpha={alpha}"
        assert abs(alpha - a_exp) < 0.01, f"{key}: alpha={alpha} vs paper {a_exp}"


# =====================================================================
# Witness 4: strong-coupling area law (RECOMPUTED)
# =====================================================================

def test_witness4_strong_coupling_sigma():
    """sigma(beta) = -ln(I1(beta)/I0(beta)) > 0 at every beta > 0, monotone
    decreasing, with both asymptotic limits: sigma -> ln(2/beta) (small beta)
    and sigma -> 1/(2 beta) (large beta, 5e-4 rel. at beta = 1000)."""
    betas = np.geomspace(0.005, 1000.0, 60)
    sig = np.array([scw.sigma_beta(b) for b in betas])
    assert (sig > 0).all()
    assert (np.diff(sig) < 0).all()
    # strong-coupling limit
    assert abs(scw.sigma_beta(0.01) / np.log(2 / 0.01) - 1.0) < 1e-5
    # weak-coupling limit (paper: 5e-4 at the JSON's beta = 1000 check)
    assert abs(scw.sigma_beta(1000.0) / (1 / 2000.0) - 1.0) < 6e-4
    # values pinned in the archived JSON verdict
    d = _load("xcwg_strong_coupling_wilson.json")["verdict_witness_2"]
    assert d["sigma_all_positive"] and d["sigma_monotone_decreasing_with_beta"]
    for key, beta in (("sigma_at_small_beta_0p3", 0.3),
                      ("sigma_at_moderate_beta_3p0", 3.0),
                      ("sigma_at_large_beta_30", 30.0)):
        assert abs(scw.sigma_beta(beta) - d[key]) < 1e-12
    # LO area law: <W> = (I1/I0)^A depends only on area A
    for beta in (0.3, 1.0, 3.0):
        for A_ in (1, 2, 5):
            assert abs(scw.wilson_strong_coupling(beta, A_)
                       - np.exp(-A_ * scw.sigma_beta(beta))) < 1e-12


# =====================================================================
# Witness 5: monopole density (structural RECOMPUTE + archived MC pins)
# =====================================================================

def test_witness5_monopole_structural_recompute():
    """k_delta = 3 (smallest closed 2-cycle = triangle-prism) and 60 elementary
    monopole sites at n_max=2; Dirac-string charges are exact integers and
    vanish at theta = 0 (DeGrand-Toussaint on Rule B)."""
    g = nlo.build_d1_and_plaquettes(2, plaq_cap=200)
    assert (g["V"], g["E"], g["P"], g["beta_1_graph"]) == (10, 20, 44, 11)
    sm = nlo.find_smallest_closed_2cycle(g["d_1"], g["plaq_adj"],
                                         max_size=8, time_budget_sec=120)
    assert sm is not None and sm["size"] == 3
    sites = nlo.enumerate_closed_2cycles_up_to_size(
        g["d_1"], g["plaq_adj"], size_target=3,
        max_count=100_000, time_budget_sec=120)
    assert len(sites) == 60
    # Dirac-string charge: integer-valued, zero at theta = 0
    assert np.all(mono.dirac_string_charge(g["d_1"], np.zeros(g["E"])) == 0)
    rng = np.random.default_rng(11)
    theta = rng.uniform(-np.pi, np.pi, size=g["E"])
    m_S = mono.monopole_charges(g["d_1"], theta, sites)
    assert m_S.dtype.kind == "i"  # Stokes: enclosed charge is an exact integer


def test_witness5_monopole_dilute_gas_fit_archive():
    """Archived seeded MC (seeds 42/43): Polyakov dilute-gas fit
    rho_M = A exp(-c beta) with c = 9.40 (n_max=2) / 8.95 (n_max=3),
    R^2 = 0.981 both; fifth witness PASS."""
    d = _load("xcwg_monopole_density.json")
    f2 = d["n_max_2"]["polyakov_dilute_gas_fit"]
    f3 = d["n_max_3"]["polyakov_dilute_gas_fit"]
    assert abs(f2["c"] - 9.40) < 0.01
    assert abs(f3["c"] - 8.95) < 0.01
    assert f2["r_squared"] > 0.98 and f3["r_squared"] > 0.98
    assert d["n_max_2"]["rng_seed"] == 42 and d["n_max_3"]["rng_seed"] == 43
    assert d["n_max_2"]["k_delta"] == 3 and d["n_max_3"]["k_delta"] == 3
    assert d["n_max_2"]["n_elementary_monopole_sites"] == 60
    assert d["n_max_3"]["n_elementary_monopole_sites"] == 80
    assert d["composite_verdict"]["fifth_witness_pass"] is True


def test_witness5_monopole_mc_smoke_seeded():
    """Reduced-stat seeded MC reproduces the qualitative witness: rho_M large
    in the condensed regime (beta = 0.3), exponentially suppressed at
    beta = 2.  Fixed seed => deterministic; verified bit-identical on rerun."""
    g = nlo.build_d1_and_plaquettes(2, plaq_cap=200)
    sites = nlo.enumerate_closed_2cycles_up_to_size(
        g["d_1"], g["plaq_adj"], size_target=3,
        max_count=100_000, time_budget_sec=120)
    edge_data = mono.plaquettes_containing_edge(g["d_1"])

    def run(beta):
        rng = np.random.default_rng(777)
        return mono.mc_run_one_beta(g["d_1"], edge_data, sites, beta,
                                    n_therm=400, n_sample=300,
                                    sample_interval=5, rng=rng, verbose=False)

    r_lo, r_hi = run(0.3), run(2.0)
    assert r_lo["rho_M_mean"] > 0.05          # condensed phase
    assert r_hi["rho_M_mean"] < 0.005         # dilute-gas suppression
    assert r_lo["rho_M_mean"] > 10 * max(r_hi["rho_M_mean"], 1e-9)
    # determinism (stability across runs at fixed seed)
    assert run(0.3)["rho_M_mean"] == r_lo["rho_M_mean"]


# =====================================================================
# Witness 6: full-MC Wilson loops (structural setup + archived pins)
# =====================================================================

def test_witness6_full_mc_archive_pins():
    """Archived seeded MC (seed 42 at n_max=2 / 43 at n_max=3): ensemble-fit
    sigma_ens(beta) > 0 and monotone decreasing (0.377 -> 0.011 -> 0.003 across
    beta = 0.1 -> 3 -> 10 at n_max=2); verdict = WEAK PASS, perimeter-dominated
    finite-volume regime (the paper's nuanced sixth witness)."""
    d = _load("xcwg_full_mc_wilson_loops.json")
    n2 = d["n_max_2"]
    assert (n2["V"], n2["E"], n2["P"], n2["beta_1_graph"]) == (10, 20, 44, 11)
    assert n2["rng_seed"] == 42 and d["n_max_3"]["rng_seed"] == 43
    rows = {r["beta"]: r["sigma_MC"] for r in n2["sigma_results"]}
    assert abs(rows[0.1] - 0.377) < 0.001
    assert abs(rows[3.0] - 0.011) < 0.001
    assert abs(rows[10.0] - 0.003) < 0.001
    sig = [r["sigma_MC"] for r in n2["sigma_results"]]
    assert all(s > 0 for s in sig)
    assert all(sig[i + 1] < sig[i] for i in range(len(sig) - 1))
    v = d["verdict"]
    assert v["ens_fit_sigma_all_positive"] and v["ens_fit_monotone_in_beta"]
    assert v["perimeter_dominated_finite_volume_regime"] is True
    assert v["sixth_witness_clean_pass"] is False      # honest nuance preserved
    assert v["sixth_witness_weak_pass_perimeter_dominated"] is True
    assert v["sixth_witness_overall"] is True


# =====================================================================
# Witness 7: Polyakov loop (homological RECOMPUTE + archived MC pins)
# =====================================================================

def test_witness7_polyakov_homological_class():
    """H^1 cocycle defect  E - V + 1 - rank(d_1)  on G_B x C_{N_t} at n_max=2:
    0 at N_t=2 (kinematic degeneracy, P == 1), exactly 1 at N_t >= 3 (the
    forced Polyakov class).  Product-graph counts pinned vs the archived JSON."""
    counts_exp = {2: (20, 50, 108), 4: (40, 120, 256), 8: (80, 240, 512)}
    defect_exp = {2: 0, 3: 1, 4: 1, 8: 1}
    for N_t, dfct in defect_exp.items():
        g = poly.build_product_graph(2, N_t)
        rank = int(np.linalg.matrix_rank(g["d_1"].astype(np.float64)))
        assert g["E"] - g["V"] + 1 - rank == dfct, f"N_t={N_t}"
        if N_t in counts_exp:
            assert (g["V"], g["E"], g["P"]) == counts_exp[N_t], f"N_t={N_t}"
    # cross-check the archived JSON graph summaries
    d = _load("xcwg_polyakov_loop.json")
    for gs in d["graph_summaries"]:
        assert (gs["V"], gs["E"], gs["P"]) == counts_exp[gs["N_t"]]


def test_witness7_polyakov_confinement_archive():
    """Archived seeded MC: |<P>| = 0 within noise at all 6 non-degenerate
    (N_t, beta) cells; N_t=2 flagged DEGENERATE (P == 1 identically);
    seventh witness PASS."""
    d = _load("xcwg_polyakov_loop.json")
    assert d["seventh_witness_pass"] is True
    n_confined = n_degenerate = 0
    for row in d["verdict_summary"]:
        if row["degenerate_N_t_2"]:
            assert row["N_t"] == 2
            assert row["abs_P_mean"] == 1.0  # P == 1 kinematic identity
            n_degenerate += 1
        else:
            assert row["verdict_label"] == "CONFINED"
            n_confined += 1
    assert n_confined == 6 and n_degenerate == 3


def test_witness7_polyakov_mc_smoke_seeded():
    """Reduced-stat seeded MC at (N_t=4, beta=1.0): Re<P>, Im<P> within 3 sigma
    of zero (persistent confinement).  Fixed seed => deterministic; verified
    bit-identical on rerun."""
    g = poly.build_product_graph(2, 4)

    def run():
        return poly.mc_run_one_point(g, beta=1.0, n_therm=300, n_sample=200,
                                     sample_interval=5, rng_seed=12345,
                                     verbose=False)

    r = run()
    assert abs(r["re_P_mean"]) <= 3 * r["re_P_se"]
    assert abs(r["im_P_mean"]) <= 3 * r["im_P_se"]
    assert r["abs_P_mean"] < 0.2  # far from the deconfined |<P>| ~ 1
    r2 = run()
    assert r2["re_P_mean"] == r["re_P_mean"]
    assert r2["abs_P_mean"] == r["abs_P_mean"]


# =====================================================================
# Archived spectral-dim JSON is faithful to the recompute (W1 archive)
# =====================================================================

def test_witness1_archive_matches_recompute():
    d = _load("xcwg_rule_b_spectral_dim.json")
    for n_max in (3, 4, 5):
        L, A, V = sdim.build_rule_b_laplacian(n_max)
        c = sdim.count_components(A)
        ev = np.clip(np.linalg.eigvalsh(L), 0.0, None)
        _, P = sdim.heat_kernel_return_probability(ev, V, T_GRID)
        d_s = sdim.extract_spectral_dimension(T_GRID, P, c, V)["d_s_window"]
        d_s_arch = d["rule_b"][str(n_max)]["spectral_dimension"]["d_s_window"]
        assert abs(d_s - d_s_arch) < 1e-9
