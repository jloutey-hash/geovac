# -*- coding: utf-8 -*-
"""B3 Phase 3, Sprint 3: cone-graded admissibility, cutoff stability of the
interval-penalty structure, Bures positivity adversarial (2026-06-10).

Executes the three named Sprint-3 items from the Sprint-2 memo
(debug/sprint_wh7_b3_phase3_sprint2_memo.md):

  T1  Cone-graded admissibility: fold the Phase-2 multiplier classes onto the
      window wedge (m' -> -m' reflection, plain-swap convention matching
      HemisphericWedge), build the unfolded K_W = V^dag |K| V, and run
      per-class penalty ensembles on BOTH substrates (full window + folded
      wedge). Decide: does the causal grading q_F = (2m'^2 - b(b+1))/(b(b+1))
      or the K-weight transfer 2|m'| organize the linear-penalty magnitudes?
      (Spearman rank correlation of ensemble-median |D+| with each predictor.)

  T2  Cutoff stability (convergence groundwork, geovac wedge n_max = 2..5):
      fixed band-limited reference perturbation + kicks (n_fock <= 2 labels,
      identical at every cutoff). Structural point: the wedge KMS state has NO
      trace-norm limit (Z grows, the state flattens -- the finite-cutoff
      shadow of type III), so naive trace-distance objects decay ~ 1/Z; but
      D_max penalties of band-limited kicks converge because the untouched
      bulk cancels in the ratio operator sigma^{-1/2} rho sigma^{-1/2}.
      Band-relative convergence is the Phase-3 convergence mode. Also pins:
      evenness dichotomy at EVERY cutoff (cutoff-uniform exact symmetry), and
      the PERIOD-ATTRIBUTION CORRECTION: U_pi = -1 already on the FULL
      Dirac space (all two_m_j odd => weight differences even), so the
      operator/orbit period is pi upstairs AND on the wedge -- the halving is
      the spinor spin-statistics grading (Phase-1 sigma_pi(F) = (-1)^{2b} F
      on a pure half-integer sector), NOT a wedge-folding effect. The
      Sprint-2 "folding halves the period" attribution is corrected here.

  T3  Bures positivity adversarial: the Sprint-2 panel observation (Bures
      excess positive on every ensemble cell for low-transfer classes) is
      attacked with random flow-commuting kicks ([G, K] = 0 by construction:
      random Hermitian blocks within K-weight eigenspaces), 2400 cells over
      25 seeds x 4 perturbation strengths x 3 flow times x 8 generators with
      sampled eps. D_max counted on the same panel for contrast.

Companion: debug/wh7_b3_phase3_sprint2.py (Sprint 2, untouched).
Frozen falsifier: tests/test_wh7_b3_phase3_sprint3.py.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(REPO))
import wh7_b1_joint_product_gh as b1  # noqa: E402
import wh7_b3_phase3_state_intervals as s1  # noqa: E402
import wh7_b3_phase3_sprint2 as p32  # noqa: E402
from geovac.modular_hamiltonian import for_bisognano_wichmann  # noqa: E402

BETA = 1.0
NW = b1.NW
E0 = 1e-4                                   # derivative rung

# Phase-2 exact causal ratios and K-weight transfers per class
QF = {"(0.5,0.5) spacelike": -1.0 / 3, "(1.0,0.0) spacelike": -1.0,
      "(2.0,0.0) spacelike": -1.0, "(2.0,1.0) spacelike": -2.0 / 3,
      "(1.0,1.0) null": 0.0, "(1.5,1.5) timelike": 1.0 / 5,
      "(2.0,2.0) timelike": 1.0 / 3}
TRANSFER = {"(0.5,0.5) spacelike": 1.0, "(1.0,0.0) spacelike": 0.0,
            "(2.0,0.0) spacelike": 0.0, "(2.0,1.0) spacelike": 2.0,
            "(1.0,1.0) null": 2.0, "(1.5,1.5) timelike": 3.0,
            "(2.0,2.0) timelike": 4.0}
BVAL = {n: float(n.split(",")[0][1:]) for n in QF}


# ---------- T1: window wedge (m' -> -m' folding) -------------------------------
def window_reflection():
    """Plain-swap reflection R: (j, a, c) -> (j, 2j-a, c), i.e. m' -> -m'
    (HemisphericWedge convention: no phase factors)."""
    idx = {lab: i for i, lab in enumerate(b1.LABELS)}
    R = np.zeros((NW, NW))
    for i, (j, a, c) in enumerate(b1.LABELS):
        dim = int(round(2 * j + 1))
        R[idx[(j, dim - 1 - a, c)], i] = 1.0
    return R


def wedge_isometry():
    """Isometry V onto the +1 eigenspace of R (sym combinations + m'=0 states);
    returns V (NW x dW) and the unfolded weights w = 2|m'| per wedge vector."""
    cols, weights, done = [], [], set()
    idx = {lab: i for i, lab in enumerate(b1.LABELS)}
    for i, (j, a, c) in enumerate(b1.LABELS):
        if i in done:
            continue
        dim = int(round(2 * j + 1))
        mp = j - a
        if mp < 0:
            continue
        v = np.zeros(NW)
        if mp == 0:
            v[i] = 1.0
        else:
            ii = idx[(j, dim - 1 - a, c)]
            v[i] = v[ii] = 1.0 / np.sqrt(2)
            done.add(ii)
        cols.append(v)
        weights.append(2.0 * mp)
    return np.array(cols).T, np.array(weights)


def class_ensemble(cfg_factory, gens, K, n_seeds=12):
    """Per-class small-eps linear coefficients + evenness + production excess
    across an ensemble of reference states. cfg_factory(seed) -> cfg dict."""
    stats = {name: {"Dp": [], "even": [], "E02": []} for name in gens}
    for seed in range(2000, 2000 + n_seeds):
        cfg = cfg_factory(seed)
        base12 = s1.d_max(cfg["s1"], cfg["mid"])
        base23 = s1.d_max(cfg["mid"], cfg["s3"])
        base = base12 + base23 - s1.d_max(cfg["s1"], cfg["s3"])
        for name, G in gens.items():
            def legs(e):
                s2_ = p32.kicked(cfg, G, e)
                return s1.d_max(cfg["s1"], s2_), s1.d_max(s2_, cfg["s3"])
            lp, lm = legs(E0), legs(-E0)
            Ep = (lp[0] - base12) + (lp[1] - base23)
            stats[name]["Dp"].append(Ep / E0)
            stats[name]["even"].append(
                max(abs(lp[0] - lm[1]), abs(lm[0] - lp[1])))
            stats[name]["E02"].append(
                p32.deficit(s1.d_max, cfg, p32.kicked(cfg, G, 0.2)) - base)
    out = {}
    for name, d in stats.items():
        out[name] = {
            "median_abs_Dp": float(np.median(np.abs(d["Dp"]))),
            "median_even": float(np.median(d["even"])),
            "frac_E02_positive": float(np.mean([e > 1e-12 for e in d["E02"]])),
            "K_commutator": None,  # filled by caller
        }
    return out


def organization_test(med_by_class):
    """Spearman rank correlation of median |D+| against the three candidate
    organizers, on the five transfer>0 classes (the linear regime)."""
    names = [n for n in QF if TRANSFER[n] > 0]
    y = [med_by_class[n]["median_abs_Dp"] for n in names]
    out = {}
    for tag, pred in (("transfer", TRANSFER), ("qF", QF), ("b", BVAL)):
        r, p = spearmanr([pred[n] for n in names], y)
        out[tag] = {"spearman_r": float(r), "p_value": float(p)}
    return out


def run_t1(n_seeds=12):
    out = {}
    R = window_reflection()
    out["R_involution_residual"] = float(np.linalg.norm(R @ R - np.eye(NW)))
    V, w = wedge_isometry()
    dW = V.shape[1]
    out["dim_wedge"] = int(dW)
    out["weights"] = sorted(float(x) for x in w)
    KW = np.diag(w.astype(complex))
    rho_W = np.diag((np.exp(-BETA * w) / np.sum(np.exp(-BETA * w)))
                    .astype(complex))

    # folded class generators: fold norms, weight-transfer witness
    gens_full = p32.class_gens()
    gens_W, fold_info = {}, {}
    for name, G in gens_full.items():
        GW = V.T @ G @ V
        GW = (GW + GW.conj().T) / 2
        nrm = float(np.linalg.norm(GW, 2))
        fold_info[name] = {"fold_norm": nrm}
        if nrm < 1e-12:
            fold_info[name]["folds_to_zero"] = True
            continue
        GW = GW / nrm
        gens_W[name] = GW
        fold_info[name]["K_commutator_wedge"] = float(
            np.linalg.norm(KW @ GW - GW @ KW, 2))
    out["fold_info"] = fold_info

    # ensembles on both substrates
    def cfg_wedge(seed):
        rng = np.random.default_rng(seed)
        H = rng.normal(size=(dW, dW)) + 1j * rng.normal(size=(dW, dW))
        H = (H + H.conj().T) / 2
        om0 = s1.conj(expm(0.3j * H), rho_W)
        return {"K": KW, "s1": om0,
                "mid": s1.conj(np.diag(np.exp(-0.5j * w)), om0),
                "s3": s1.conj(np.diag(np.exp(-1.0j * w)), om0)}

    def cfg_full(seed):
        return p32.make_config(seed, theta=0.3)

    ens_W = class_ensemble(cfg_wedge, gens_W, KW, n_seeds)
    for name in ens_W:
        ens_W[name]["K_commutator"] = fold_info[name]["K_commutator_wedge"]
    ens_F = class_ensemble(cfg_full, gens_full, s1.boost_K(), n_seeds)
    K_full = s1.boost_K()
    for name, G in gens_full.items():
        ens_F[name]["K_commutator"] = float(
            np.linalg.norm(K_full @ G - G @ K_full, 2))
    out["ensemble_wedge"] = ens_W
    out["ensemble_full"] = ens_F
    out["organization_full"] = organization_test(ens_F)
    if all(n in ens_W for n in QF if TRANSFER[n] > 0):
        out["organization_wedge"] = organization_test(ens_W)
    return out


# ---------- T2: geovac wedge cutoff sweep --------------------------------------
def run_t2(n_list=(2, 3, 4, 5), t_total=1.0):
    out = {"per_n": {}}
    fixed_band_keys = None
    Hb_fixed = None
    kick_pairs = None
    for n_max in n_list:
        mh = for_bisognano_wichmann(n_max)
        wi = mh.wedge_basis_indices()
        labs = [mh.basis[i] for i in wi]
        w = np.array([float(l.two_m_j) for l in labs])
        dW = len(wi)

        # upstairs period correction: all two_m_j odd => U_pi = -1 on H_full
        diag_full = np.real(np.diag(mh.K_geometric))
        u_pi_full = float(np.max(np.abs(np.exp(-1j * np.pi * diag_full) + 1)))

        # canonical band (n_fock <= 2), identical labels at every cutoff
        band = [k for k in range(dW) if labs[k].n_fock <= 2]
        band.sort(key=lambda k: (labs[k].n_fock, labs[k].l,
                                 labs[k].two_m_j, labs[k].chirality))
        keys = tuple((labs[k].n_fock, labs[k].l, labs[k].two_m_j,
                      labs[k].chirality) for k in band)
        if fixed_band_keys is None:
            fixed_band_keys = keys
            nb = len(band)
            rng = np.random.default_rng(7)
            Hb = rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb))
            Hb_fixed = (Hb + Hb.conj().T) / 2
            Hb_fixed /= np.linalg.norm(Hb_fixed, 2)
            # fixed kick pairs in canonical band order: first same-weight pair
            # (Dw = 0) and first Dw = 2 pair
            kick_pairs = {}
            for dw in (0, 2):
                for i in range(nb):
                    for j in range(i + 1, nb):
                        if abs(keys[i][2] - keys[j][2]) == dw:
                            kick_pairs[dw] = (i, j)
                            break
                    if dw in kick_pairs:
                        break
        assert keys == fixed_band_keys, f"band labels differ at n_max={n_max}"
        nb = len(band)

        H = np.zeros((dW, dW), dtype=complex)
        H[np.ix_(band, band)] = Hb_fixed
        gens = {}
        for dw, (i, j) in kick_pairs.items():
            G = np.zeros((dW, dW), dtype=complex)
            G[band[i], band[j]] = 1.0
            G[band[j], band[i]] = 1.0
            gens[dw] = G / np.linalg.norm(G, 2)

        Z = float(np.sum(np.exp(-BETA * w)))
        rho = np.diag((np.exp(-BETA * w) / Z).astype(complex))
        om0 = s1.conj(expm(0.3j * H), rho)
        cfg = {"K": np.diag(w.astype(complex)), "s1": om0,
               "mid": s1.conj(np.diag(np.exp(-0.5j * t_total * w)), om0),
               "s3": s1.conj(np.diag(np.exp(-1j * t_total * w)), om0)}
        base12 = s1.d_max(cfg["s1"], cfg["mid"])
        base23 = s1.d_max(cfg["mid"], cfg["s3"])
        base = base12 + base23 - s1.d_max(cfg["s1"], cfg["s3"])

        def legs(G, e):
            s2_ = p32.kicked(cfg, G, e)
            return s1.d_max(cfg["s1"], s2_), s1.d_max(s2_, cfg["s3"])

        lp0, lm0 = legs(gens[0], E0), legs(gens[0], -E0)
        lp2 = legs(gens[2], E0)
        row = {
            "dim_full": len(mh.basis), "dim_W": dW, "n_band": nb, "Z": Z,
            "U_pi_plus_identity_full": u_pi_full,
            "baseline_deficit": float(base),
            "E02_dw0": float(p32.deficit(s1.d_max, cfg,
                                         p32.kicked(cfg, gens[0], 0.2)) - base),
            "E02_dw2": float(p32.deficit(s1.d_max, cfg,
                                         p32.kicked(cfg, gens[2], 0.2)) - base),
            "Dp_dw2": float(((lp2[0] - base12) + (lp2[1] - base23)) / E0),
            "evenness_dw0": float(max(abs(lp0[0] - lm0[1]),
                                      abs(lm0[0] - lp0[1]))),
        }
        # naive trace-norm orbit scale (decays ~ 1/Z) vs D_max-matching
        # interval recovery (band-stable)
        ts = np.linspace(0.05, np.pi - 0.05, 20)
        prof = [s1.trace_dist(s1.conj(np.diag(np.exp(-1j * t * w)), om0), om0)
                for t in ts]
        row["orbit_scale_trace"] = float(np.max(prof))

        dt_true = 0.7
        om_a = s1.conj(np.diag(np.exp(-0.4j * w)), om0)
        om_b = s1.conj(np.diag(np.exp(-1j * (0.4 + dt_true) * w)), om0)
        taus = np.linspace(0.0, np.pi, 314, endpoint=False)
        dvals = [s1.d_max(s1.conj(np.diag(np.exp(-1j * t * w)), om_a), om_b)
                 for t in taus]
        i0 = int(np.argmin(dvals))
        im, ip = (i0 - 1) % 314, (i0 + 1) % 314
        den = dvals[im] - 2 * dvals[i0] + dvals[ip]
        sh = 0.5 * (dvals[im] - dvals[ip]) / den if abs(den) > 1e-300 else 0.0
        row["interval_recovery_err_dmax"] = float(
            abs((taus[i0] + sh * (np.pi / 314)) % np.pi - dt_true))
        out["per_n"][str(n_max)] = row
    # stabilization against the largest cutoff
    ref = out["per_n"][str(n_list[-1])]
    for n_max in n_list[:-1]:
        r = out["per_n"][str(n_max)]
        r["delta_E02_dw0_vs_ref"] = float(abs(r["E02_dw0"] - ref["E02_dw0"]))
        r["delta_E02_dw2_vs_ref"] = float(abs(r["E02_dw2"] - ref["E02_dw2"]))
        r["delta_Dp_dw2_vs_ref"] = float(abs(r["Dp_dw2"] - ref["Dp_dw2"]))
    return out


# ---------- T3: Bures positivity adversarial ------------------------------------
def commuting_gen(rng, K):
    """Random Hermitian commuting with K: blocks within K-weight eigenspaces."""
    w = np.real(np.diag(K))
    G = np.zeros((NW, NW), dtype=complex)
    for val in np.unique(np.round(w, 9)):
        idx = np.where(np.abs(w - val) < 1e-9)[0]
        m = len(idx)
        B = rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m))
        G[np.ix_(idx, idx)] = (B + B.conj().T) / 2
    return G / np.linalg.norm(G, 2)


def run_t3():
    rng = np.random.default_rng(20260615)
    n_cells, neg_b, neg_d = 0, 0, 0
    worst_b, worst_d = np.inf, np.inf
    worst_params = None
    max_comm = 0.0
    class_cells, class_neg = 0, 0
    m0_gens = {n: g for n, g in p32.class_gens().items()
               if TRANSFER[n] == 0.0}
    for seed in range(25):
        for theta in (0.1, 0.3, 0.6, 1.0):
            for t_tot in (0.3, 1.0, 2.5):
                cfg = p32.make_config(3000 + seed, theta=theta,
                                      t_total=t_tot)
                base_b = p32.deficit(p32.bures_angle, cfg, cfg["mid"])
                base_d = p32.deficit(s1.d_max, cfg, cfg["mid"])
                for _ in range(8):
                    G = commuting_gen(rng, cfg["K"])
                    max_comm = max(max_comm, float(np.linalg.norm(
                        cfg["K"] @ G - G @ cfg["K"], 2)))
                    eps = float(rng.choice([0.05, 0.2, 0.5]))
                    s2 = p32.kicked(cfg, G, eps)
                    Eb = p32.deficit(p32.bures_angle, cfg, s2) - base_b
                    Ed = p32.deficit(s1.d_max, cfg, s2) - base_d
                    n_cells += 1
                    if Eb < worst_b:
                        worst_b, worst_params = Eb, (seed, theta, t_tot, eps)
                    worst_d = min(worst_d, Ed)
                    if Eb < -1e-10:
                        neg_b += 1
                    if Ed < -1e-10:
                        neg_d += 1
                # continuity cells: the two m'=0 class generators
                for name, G in m0_gens.items():
                    Eb = p32.deficit(p32.bures_angle, cfg,
                                     p32.kicked(cfg, G, 0.2)) - base_b
                    class_cells += 1
                    if Eb < -1e-10:
                        class_neg += 1
                    worst_b = min(worst_b, Eb)
    return {"n_cells": n_cells, "bures_negative": neg_b,
            "dmax_negative": neg_d, "worst_bures_excess": float(worst_b),
            "worst_dmax_excess": float(worst_d),
            "worst_params": worst_params,
            "max_K_commutator": max_comm,
            "class_gen_cells": class_cells, "class_gen_negative": class_neg}


def run():
    return {"T1_cone_admissibility": run_t1(),
            "T2_cutoff_stability": run_t2(),
            "T3_bures_adversarial": run_t3()}


if __name__ == "__main__":
    res = run()
    (ROOT / "data" / "wh7_b3_phase3_sprint3.json").write_text(
        json.dumps(res, indent=1, default=str), encoding="utf-8")

    t1 = res["T1_cone_admissibility"]
    print(f"T1 window wedge: dim {t1['dim_wedge']}, R involution "
          f"{t1['R_involution_residual']:.1e}")
    print("   class                     fold_nrm  [K,G]_W   med|D+|_W "
          f" med|D+|_full  evenW")
    for name in QF:
        fi = t1["fold_info"][name]
        if name in t1["ensemble_wedge"]:
            ew, ef = t1["ensemble_wedge"][name], t1["ensemble_full"][name]
            print(f"   {name:24s}: {fi['fold_norm']:.3f}   "
                  f"{ew['K_commutator']:.3f}    {ew['median_abs_Dp']:.4f}   "
                  f"{ef['median_abs_Dp']:.4f}      {ew['median_even']:.1e}")
        else:
            print(f"   {name:24s}: {fi['fold_norm']:.3f}   FOLDS TO ZERO")
    for tag in ("organization_full", "organization_wedge"):
        if tag in t1:
            o = t1[tag]
            print(f"   {tag}: " + "  ".join(
                f"{k}: r={v['spearman_r']:+.2f} (p={v['p_value']:.3f})"
                for k, v in o.items()))

    t2 = res["T2_cutoff_stability"]
    print("T2 cutoff sweep (geovac wedge):")
    for n, r in t2["per_n"].items():
        extra = ""
        if "delta_E02_dw0_vs_ref" in r:
            extra = (f" | dE0={r['delta_E02_dw0_vs_ref']:.2e} "
                     f"dE2={r['delta_E02_dw2_vs_ref']:.2e} "
                     f"dDp={r['delta_Dp_dw2_vs_ref']:.2e}")
        print(f"   n_max={n}: dim {r['dim_full']}->{r['dim_W']} Z={r['Z']:.2f} "
              f"| U_pi_full {r['U_pi_plus_identity_full']:.1e} | "
              f"E02(dw0)={r['E02_dw0']:+.4f} E02(dw2)={r['E02_dw2']:+.4f} "
              f"D+(dw2)={r['Dp_dw2']:+.4f} even(dw0)={r['evenness_dw0']:.1e} | "
              f"orbit_scale={r['orbit_scale_trace']:.4f} "
              f"rec_err={r['interval_recovery_err_dmax']:.1e}{extra}")

    t3 = res["T3_bures_adversarial"]
    print(f"T3 Bures adversarial: {t3['bures_negative']}/{t3['n_cells']} "
          f"negative (worst {t3['worst_bures_excess']:+.2e}); D_max contrast "
          f"{t3['dmax_negative']}/{t3['n_cells']} "
          f"(worst {t3['worst_dmax_excess']:+.2e}); max [K,G] "
          f"{t3['max_K_commutator']:.1e}; class-gen cells "
          f"{t3['class_gen_negative']}/{t3['class_gen_cells']} negative")
