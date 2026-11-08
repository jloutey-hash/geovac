# -*- coding: utf-8 -*-
"""Pinning tests for Paper 41's seven-witness Rule B Wilson U(1) verdict.

Backing drivers live in tests/wilson_rule_b_support/ (2026-07-03 durability
migration from debug/archive/{qed_arc,misc,chemistry_qc_arc}/ -- the wh7_support
pattern, CLAUDE.md v4.49.1).  Archived production JSON lives in
tests/wilson_rule_b_support/data/.

Pin strategy (per witness):
  W1 spectral dimension   -- RECOMPUTED (heat-kernel d_s, <1 s).
  W2 Migdal-Kadanoff RG   -- RECOMPUTED (k_eff, no fixed point, monotone flow).
  W3 Gaussian Wilson loop -- RECOMPUTED at n_max=2 (full enumeration,
                             deterministic); JSON band-pin at n_max=3,4,5.
  W4 strong-coupling      -- RECOMPUTED (sigma(beta) = -ln(I1/I0) positivity +
                             the small-beta limit) + JSON cross-check.
  W5 monopole density     -- structural facts RECOMPUTED (k_delta=3, 60
                             elementary sites, P=44); production MC fit values
                             (n_max=2: c=9.40, R^2=0.981; n_max=3 full
                             complex: c=22.9, R^2=0.935) pinned from the
                             archived seeded JSON; reduced-stat seeded MC smoke verifies
                             the exponential-suppression signal + determinism.
  W6 full-MC Wilson loops -- structural setup + archived seeded-JSON pins
                             (sigma_ens > 0 monotone; WEAK-PASS flags).  The
                             production MC (n_sample=2000 x 11 beta) is too
                             heavy to recompute in a test; the archived JSON
                             carries fixed seeds 42/43.
  W7 Polyakov loop        -- homological class RECOMPUTED (H^1 defect
                             E - V + 1 - rank(d_1): 0 at N_t=2, 1 at N_t>=3,
                             product-graph counts); centre symmetry (Delta S = 0,
                             P -> e^{i phi} P) RECOMPUTED, so <P> = 0 is forced;
                             the order parameter <|P_lat|> pinned from the
                             archived seeded JSON (the 'confined' reading is
                             withdrawn, 2026-11-07).
  W1 note (2026-11-07): the sliding-window d_s passes 3 at n_max = 8; the
                             fixed-scale gap over the scalar graph is the
                             witness content, not a trend toward 3.

Paper anchors: papers/group5_qed_gauge/paper_30_su2_wilson.tex sec:ruleb
("Seven diagnostics"), which absorbed the archived
papers/archive/paper_41_rule_b_wilson_u1.tex.  Asserts backing Paper 41
claims that Paper 30 did not keep (the Weyl dimension d_W, the MK asymptotic
accuracies, the large-beta limit of sigma, the dual-graph counts) were
removed 2026-11-08.  Archived Paper 41 section labels:
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
RULE_B_DS = {3: 1.86, 4: 2.27, 5: 2.54}
SCALAR_DS = {3: 1.68, 4: 1.76, 5: 1.79}


@pytest.mark.parametrize("n_max", [3, 4, 5])
def test_witness1_spectral_dimension_rule_b(n_max):
    ds_exp = RULE_B_DS[n_max]
    L, A, V = sdim.build_rule_b_laplacian(n_max)
    c = sdim.count_components(A)
    ev = np.clip(np.linalg.eigvalsh(L), 0.0, None)
    _, P = sdim.heat_kernel_return_probability(ev, V, T_GRID)
    d_s = sdim.extract_spectral_dimension(T_GRID, P, c, V)["d_s_window"]
    assert abs(d_s - ds_exp) < 0.02, f"d_s({n_max}) = {d_s}"


def _return_probability(build, n_max):
    L, A, V = build(n_max)
    c = sdim.count_components(A)
    ev = np.clip(np.linalg.eigvalsh(L), 0.0, None)
    _, P = sdim.heat_kernel_return_probability(ev, V, T_GRID)
    return P, c, V


def test_witness1_paper_estimator_passes_3():
    """Rebuilt 2026-11-07 (group5 guard pass; the former test asserted a
    monotone climb TOWARD 3, which Paper 41 withdraws as an estimator
    artifact).  With the paper's own sliding-window estimator, Rule B's d_s
    keeps rising past 3: 2.74, 2.89, 3.04, 3.15 at n_max = 6..9.  Rejects
    'd_s saturates at 3'."""
    printed = {6: 2.74, 7: 2.89, 8: 3.04, 9: 3.15}
    for n_max, val in printed.items():
        P, c, V = _return_probability(sdim.build_rule_b_laplacian, n_max)
        d_s = sdim.extract_spectral_dimension(T_GRID, P, c, V)["d_s_window"]
        assert abs(d_s - val) < 0.01, (n_max, d_s)
    assert d_s > 3.1


def _ds_fixed_window(P, c, V, lo=0.1, hi=0.5):
    """d_s fitted over a FIXED return-probability band lo < P - c/V < hi
    (the same diffusion scale at every graph size)."""
    dP = P - c / V
    m = (dP > lo) & (dP < hi)
    assert m.sum() >= 10
    return -2.0 * np.polyfit(np.log(T_GRID[m]), np.log(dP[m]), 1)[0]


def test_witness1_fixed_scale_gap_no_trend():
    """Paper 41 W1 reading: at fixed diffusion scale Rule B is nearly flat
    (~2.15 -> 2.30 for n_max = 3..8) and sits about 0.45 above the scalar
    Fock graph (~1.7 -> 1.8).  Rejects a trend toward 3 (Rule B would rise
    by far more than 0.2) and a vanishing or reversed gap.  The band is the
    one Paper 41 states (sec:spectral_dim, 'Reading')."""
    ds_b, ds_s = [], []
    for n_max in range(3, 9):
        ds_b.append(_ds_fixed_window(*_return_probability(sdim.build_rule_b_laplacian, n_max)))
        ds_s.append(_ds_fixed_window(*_return_probability(sdim.build_scalar_laplacian, n_max)))
    ds_b, ds_s = np.array(ds_b), np.array(ds_s)
    assert np.all((ds_b > 2.10) & (ds_b < 2.35)), ds_b
    assert ds_b.max() - ds_b.min() < 0.2, ds_b
    assert np.all((ds_s > 1.65) & (ds_s < 1.85)), ds_s
    gap = float(np.mean(ds_b - ds_s))
    assert 0.35 < gap < 0.55, gap
    # scalar paper-estimator values (Paper 41 W1 table)
    for n_max in (3, 4, 5):
        P, c, V = _return_probability(sdim.build_scalar_laplacian, n_max)
        d_s = sdim.extract_spectral_dimension(T_GRID, P, c, V)["d_s_window"]
        assert abs(d_s - SCALAR_DS[n_max]) < 0.02


# =====================================================================
# Witness 2: Migdal-Kadanoff RG (RECOMPUTED)
# =====================================================================

def test_witness2_mk_k_eff_and_flow():
    """Paper 30 sec:ruleb (2): k = <q_e> (mean plaquettes per edge) = 8.8 at
    n_max=2, 37.51 at n_max=3 (tab. of Rule-B counts); MK flow monotone to
    beta=0; no fixed point in (0, 100)."""
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
    decreasing, with the strong-coupling limit sigma -> ln(2/beta) (Paper 30
    sec:ruleb (4))."""
    betas = np.geomspace(0.005, 1000.0, 60)
    sig = np.array([scw.sigma_beta(b) for b in betas])
    assert (sig > 0).all()
    assert (np.diff(sig) < 0).all()
    # strong-coupling limit
    assert abs(scw.sigma_beta(0.01) / np.log(2 / 0.01) - 1.0) < 1e-5
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
    """Paper 30 sec:ruleb (5), archived seeded MC.  n_max = 2 (seed 42, all
    60 sites): rho_M = A exp(-c beta) with c = 9.40, R^2 = 0.981.  n_max = 3
    on the full complex (seed 43, 994 plaquettes, all 3342 sites;
    debug/p30_full_complex_rerun.py): c = 22.9, R^2 = 0.935 on
    0.05 <= beta <= 0.5, zero from beta = 0.7.  The retired n_max = 3 fit
    (c = 8.95, R^2 = 0.981) came from 200 plaquettes and 80 sites."""
    d = _load("xcwg_monopole_density.json")
    f2 = d["n_max_2"]["polyakov_dilute_gas_fit"]
    assert abs(f2["c"] - 9.40) < 0.01
    assert f2["r_squared"] > 0.98
    assert d["n_max_2"]["rng_seed"] == 42
    assert d["n_max_2"]["k_delta"] == 3
    assert d["n_max_2"]["n_elementary_monopole_sites"] == 60
    d3 = _load("p30_full_complex_mono3.json")
    assert (d3["P"], d3["n_elementary_monopole_sites"], d3["k_delta"],
            d3["rng_seed"]) == (994, 3342, 3, 43)
    f3 = d3["polyakov_dilute_gas_fit"]
    assert abs(f3["c"] - 22.86) < 0.01
    assert abs(f3["r_squared"] - 0.935) < 0.001
    assert (f3["beta_min_fit"], f3["beta_max_fit"]) == (0.05, 0.5)
    rho = {r["beta"]: r["rho_M_mean"] for r in d3["scan_results"]}
    assert rho[0.5] > 0.0 and all(rho[b] == 0.0 for b in (0.7, 1.0, 2.0, 5.0, 30.0))
    assert rho[0.1] / rho[0.2] > 50.0   # the two-decade drop between 0.1 and 0.2
    assert abs(f3["c"] - 8.95) > 5.0   # the retired truncated-complex fit


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


def test_witness7_centre_symmetry_forces_P_zero():
    """Rebuilt 2026-11-07 (group5 guard pass; the former test pinned the
    driver's 'CONFINED' labels, a reading Paper 41 withdraws).  Shifting
    every temporal link on one time slice by phi leaves the Wilson action
    unchanged and multiplies every Polyakov loop by e^{i phi}; hence
    <P> = 0 at every beta in finite volume and carries no phase information.
    Rejects '<P> = 0 => confined'."""
    rng = np.random.default_rng(3)
    for N_t in (3, 4, 8):
        g = poly.build_product_graph(2, N_t)
        d_1, path = g["d_1"], g["polyakov_path"]
        theta = rng.uniform(-np.pi, np.pi, d_1.shape[1])
        phi = 0.731
        shifted = theta.copy()
        for v in range(len(path)):
            e_idx, sign = path[v][0]          # the temporal edge on slice 0
            shifted[e_idx] += sign * phi
        S0 = np.sum(np.cos(poly.plaquette_angles(d_1, theta)))
        S1 = np.sum(np.cos(poly.plaquette_angles(d_1, shifted)))
        assert abs(S1 - S0) < 1e-10, (N_t, S1 - S0)
        P0 = poly.measure_polyakov(theta, path)
        P1 = poly.measure_polyakov(shifted, path)
        assert np.allclose(P1, np.exp(1j * phi) * P0)


def test_witness7_order_parameter_rises_with_beta():
    """Paper 41's stored order parameter <|P_lat|> (archived seeded MC):
    N_t = 4: 0.285, 0.498, 0.833; N_t = 8: 0.285, 0.341, 0.691 at
    beta = 0.3, 1, 3.  The beta = 0.3 value is the independent-random-phase
    value on V_B = 10 sites (~sqrt(pi/40)), recomputed here.  Rejects
    'persistent confinement' (<|P_lat|> pinned at the random-phase value)."""
    d = _load("xcwg_polyakov_loop.json")
    printed = {(4, 0.3): 0.285, (4, 1.0): 0.498, (4, 3.0): 0.833,
               (8, 0.3): 0.285, (8, 1.0): 0.341, (8, 3.0): 0.691}
    got = {(r["N_t"], r["beta"]): r["abs_lat_avg_mean"] for r in d["all_results"]
           if r["N_t"] in (4, 8)}
    for key, val in printed.items():
        assert abs(got[key] - val) < 5e-4, (key, got[key])
    V_B = len(poly.build_product_graph(2, 4)["polyakov_path"])
    assert V_B == 10
    rng = np.random.default_rng(5)
    rand = np.abs(np.exp(1j * rng.uniform(0, 2 * np.pi, (200000, V_B))).mean(axis=1)).mean()
    assert abs(rand - got[(4, 0.3)]) < 0.012, rand
    for N_t in (4, 8):
        assert got[(N_t, 0.3)] < got[(N_t, 1.0)] < got[(N_t, 3.0)]
        assert got[(N_t, 3.0)] > rand + 0.3


def test_witness7_polyakov_mc_smoke_seeded():
    """Reduced-stat seeded MC at (N_t=4, beta=1.0): Re<P>, Im<P> within 3 sigma
    of zero, as centre symmetry requires (not a confinement signal).  Fixed
    seed => deterministic; verified bit-identical on rerun."""
    g = poly.build_product_graph(2, 4)

    def run():
        return poly.mc_run_one_point(g, beta=1.0, n_therm=300, n_sample=200,
                                     sample_interval=5, rng_seed=12345,
                                     verbose=False)

    r = run()
    assert abs(r["re_P_mean"]) <= 3 * r["re_P_se"]
    assert abs(r["im_P_mean"]) <= 3 * r["im_P_se"]
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


# =====================================================================
# Cert-run additions (2026-07-04): numeric re-derivations + v5 pins
# =====================================================================

def test_witness6_perimeter_dominance_from_numeric_fields():
    """W6 re-derived from the archived JSON's NUMERIC fields (not the
    driver's own verdict flags): at every beta the joint-fit area
    coefficient is tiny while the perimeter coefficient is large and
    positive -- the honest weak-pass 'perimeter-dominated' reading.
    Added 2026-07-04 (cert run: the verdict-label asserts alone could
    not catch a silently flipped fit)."""
    d = _load("xcwg_full_mc_wilson_loops.json")
    rows = d["n_max_2"]["sigma_results"]
    assert len(rows) >= 5
    for r in rows:
        s, mu = r["sigma_MC_combined"], r["mu_MC_combined"]
        assert abs(s) < 0.05          # area coefficient tiny at every beta
        assert mu > 0.0               # perimeter coefficient positive
        assert mu > 5 * abs(s)        # perimeter dominance, per row
    # strong-coupling rows carry the large perimeter signal (mu decays
    # with beta as <W> -> 1; the flat signal lives at small beta)
    strong = {r["beta"]: r["mu_MC_combined"] for r in rows if r["beta"] <= 0.3}
    assert all(mu > 0.1 for mu in strong.values()) and len(strong) == 3


def test_witness7_confinement_from_numeric_fields():
    """W7 table from the archived JSON's numeric fields: every
    non-degenerate (N_t, beta) cell has Re<P> and Im<P> zero within
    2 sigma (paper: all 6 cells within 2 sigma), as centre symmetry
    forces; this carries no confinement information.  Added 2026-07-04."""
    d = _load("xcwg_polyakov_loop.json")
    checked = 0
    for row in d["verdict_summary"]:
        if row["degenerate_N_t_2"]:
            continue
        assert abs(row["re_P_mean"]) <= 2.0 * row["re_P_se"], row
        assert abs(row["im_P_mean"]) <= 2.0 * row["im_P_se"], row
        checked += 1
    assert checked == 6


def test_nmax4_full_complex_sigma_comb_pins():
    """Paper 30 sec:ruleb (6): at n_max = 4 on all 5047 plaquettes the joint
    fit gives sigma = -0.0043(15), -0.0006(7), -0.0005(1) against
    mu = 0.022, 0.0058, 0.0020 at beta = 0.3, 1, 3 (archived seeded run,
    debug/p30_full_complex_rerun.py, seed 44).  The retired -0.080 at
    beta = 0.3 came from a run on 1500 of the 5047 plaquettes."""
    d = _load("p30_full_complex_wilson4.json")
    assert (d["V"], d["E"], d["beta_1_graph"], d["P"]) == (60, 312, 253, 5047)
    rows = {r["beta"]: r for r in d["sigma_results"]}
    for beta, sig, se, mu in ((0.3, -0.0043, 0.0015, 0.0220),
                              (1.0, -0.0006, 0.0007, 0.0058),
                              (3.0, -0.0005, 0.0001, 0.0020)):
        r = rows[beta]
        assert abs(r["sigma_MC_combined"] - sig) < 5e-5, (beta, r)
        assert abs(r["sigma_MC_combined_se"] - se) < 5e-5, (beta, r)
        assert abs(r["mu_MC_combined"] - mu) < 5e-5, (beta, r)
    # the retired truncated-complex value is far outside the full-complex fit
    assert abs(rows[0.3]["sigma_MC_combined"] - (-0.080)) > 0.05


def test_v5_polyakov_rate_resolution_pins():
    """Paper 30 sec:ruleb (6): at n_max = 2 the exponential rate of the
    area-only slope, c_sigma_ens = 1.196, matches that of the perimeter
    coefficient, c_mu_comb = 1.242 (4% agreement)."""
    d = _load("xcwg_polyakov_rate_refinement.json")
    m = d["n_max_2_what_measured"]
    assert abs(m["c_sigma_ens"] - 1.196) < 0.001
    assert abs(m["c_mu_comb"] - 1.242) < 0.001
    assert abs(m["c_sigma_ens"] / m["c_mu_comb"] - 1.0) < 0.05
