# -*- coding: utf-8 -*-
"""B3 Phase 3, Sprint 2: ensemble sign statistics, cost-functional comparison,
bimodal-p mechanism, and the wedge substrate (2026-06-10).

Executes the four named Sprint-2 items from the Sprint-1 memo
(debug/sprint_wh7_b3_phase3_sprint1_memo.md):

  T1  Ensemble statistics over reference states: the excess-sign question as a
      distributional statement (seeds x perturbation strengths x kick classes),
      for ALL cost functionals in the comparison set.
  T2  Cost-functional comparison: D_max vs Umegaki vs Bures angle vs trace
      distance vs Jeffreys (symmetrized Umegaki) -- chain-inequality status on
      the Sprint-1 96-cell panel design + which functional (if any) produces a
      state-independent positive excess.
  T3  Mechanism of the bimodal p: one-sided derivatives of the excess at
      eps -> 0+/-, per-leg first-order terms (numerical + analytic rho-slot
      perturbation formula), top-eigenvalue gaps, local log-log slope profile.
      Decides kink (|eps|) vs quadratic and locates the non-smoothness.
  T4  Wedge substrate stand-up: the Sprint-1 anchors (KMS / flow-invariance /
      period / orbit injectivity / D_max chain) on the PROPER Phase-3 substrate
      (HemisphericWedge of geovac/modular_hamiltonian.py, unfolded K_W with odd
      positive integer spectrum), kick grading by K-weight transfer
      Dw in {0, 2, 4} (the wedge analog of the Sprint-1 |m'| split), and the
      OPERATIONAL interval functional ell = flow parameter recovered from state
      pairs by trace-distance matching (well-definedness + additivity).

Substrate notes documented here, not load-bearing: beta = 1 on the wedge (the
Sprint-1 convention); on the wedge the weight differences are even integers, so
the operator-level orbit period is pi (U_pi = -I), HALF the full-space 2*pi --
the positive-weight folding halves the thermal-time period.

Companion: debug/wh7_b3_phase3_state_intervals.py (Sprint 1, untouched).
Frozen falsifier: tests/test_wh7_b3_phase3_sprint2.py.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent.parent  # project root (moved debug/ -> tests/wh7_support/, 2026-06-26)
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(REPO))
import wh7_b1_joint_product_gh as b1  # noqa: E402
import wh7_b3_phase3_state_intervals as s1  # noqa: E402
from geovac.modular_hamiltonian import for_bisognano_wichmann  # noqa: E402

BETA = 1.0
NW = b1.NW
EPS_PROD = [0.05, 0.1, 0.2, 0.4]            # Sprint-1 production ladder
EPS_SMALL = [1e-4, 1e-3, 1e-2]              # mechanism ladder (T3)


# ---------- cost functionals --------------------------------------------------
def _psd_power(A, p):
    ev, V = np.linalg.eigh(A)
    ev = np.maximum(ev, 0.0)
    return V @ np.diag(ev ** p) @ V.conj().T


def bures_angle(rho, sigma):
    """A(rho, sigma) = arccos sqrt(F); F = (tr sqrt(rho^1/2 sigma rho^1/2))^2.
    A true metric (triangle inequality is a theorem)."""
    rh = _psd_power(rho, 0.5)
    M = rh @ sigma @ rh
    ev = np.maximum(np.linalg.eigvalsh(M), 0.0)
    f = float(np.sum(np.sqrt(ev)))
    return float(np.arccos(np.clip(f, 0.0, 1.0)))


def jeffreys(rho, sigma):
    """Symmetrized Umegaki divergence J = D(rho||sigma) + D(sigma||rho)."""
    return s1.umegaki(rho, sigma) + s1.umegaki(sigma, rho)


COSTS = {
    "dmax": s1.d_max,
    "umegaki": s1.umegaki,
    "bures": bures_angle,
    "trace": s1.trace_dist,
    "jeffreys": jeffreys,
}
METRICS = ("bures", "trace")               # true metrics: chain is a theorem


# ---------- shared configuration machinery ------------------------------------
def make_config(seed, theta=0.3, t_total=1.0):
    """Reference orbit triple (s1, mid, s3) for a seeded perturbed KMS state."""
    K = s1.boost_K()
    rho = s1.kms_state(K)
    rng = np.random.default_rng(seed)
    H = rng.normal(size=(NW, NW)) + 1j * rng.normal(size=(NW, NW))
    H = (H + H.conj().T) / 2
    om0 = s1.conj(expm(1j * theta * H), rho)
    return {
        "K": K,
        "s1": om0,
        "mid": s1.conj(s1.flow_U(K, t_total / 2), om0),
        "s3": s1.conj(s1.flow_U(K, t_total), om0),
    }


def deficit(cost, cfg, s2):
    return (cost(cfg["s1"], s2) + cost(s2, cfg["s3"])
            - cost(cfg["s1"], cfg["s3"]))


def kicked(cfg, G, eps):
    return s1.conj(expm(1j * eps * G), cfg["mid"])


def class_gens():
    return {name: s1.class_generator(b, r, c)
            for name, (b, r, c) in s1.CLASSES.items()}


# ---------- T1 + T2: ensemble sign statistics, all functionals ----------------
def run_ensemble(gens, n_seeds=24, thetas=(0.1, 0.3, 0.6), eps=0.2):
    cells = []
    for seed in range(1000, 1000 + n_seeds):
        for theta in thetas:
            cfg = make_config(seed, theta=theta)
            base = {k: deficit(c, cfg, cfg["mid"]) for k, c in COSTS.items()}
            for name, G in gens.items():
                s2 = kicked(cfg, G, eps)
                row = {"seed": seed, "theta": theta, "class": name}
                for k, c in COSTS.items():
                    row[f"excess_{k}"] = deficit(c, cfg, s2) - base[k]
                cells.append(row)
    summary = {}
    for k in COSTS:
        per_class = {}
        for name in gens:
            vals = [c[f"excess_{k}"] for c in cells if c["class"] == name]
            per_class[name] = {
                "frac_positive": float(np.mean([v > 1e-12 for v in vals])),
                "min": float(np.min(vals)), "median": float(np.median(vals)),
            }
        allv = [c[f"excess_{k}"] for c in cells]
        per_theta = {}
        for th in thetas:
            tv = [c[f"excess_{k}"] for c in cells if c["theta"] == th]
            per_theta[str(th)] = float(np.mean([v > 1e-12 for v in tv]))
        summary[k] = {
            "frac_positive_overall": float(np.mean([v > 1e-12 for v in allv])),
            "frac_positive_per_theta": per_theta,
            "per_class": per_class,
        }
    return summary, cells


# ---------- T2: chain-inequality panel (Sprint-1 96-cell design) ---------------
def run_chain_panel(gens):
    cfg = make_config(20260613, theta=0.3)   # the Sprint-1 primary seed
    rng = np.random.default_rng(20260613)
    s2_list = []
    baselines = {k: deficit(c, cfg, cfg["mid"]) for k, c in COSTS.items()}
    for name, G in gens.items():
        Gp = G - (float(np.real(np.vdot(1j * (G @ cfg["mid"] - cfg["mid"] @ G),
                                        1j * (cfg["K"] @ cfg["mid"]
                                              - cfg["mid"] @ cfg["K"])))) /
                  max(np.linalg.norm(1j * (cfg["K"] @ cfg["mid"]
                                           - cfg["mid"] @ cfg["K"])) ** 2,
                      1e-300)) * cfg["K"]
        Gp = Gp / np.linalg.norm(Gp, 2)
        for gen in (G, Gp):
            for eps in EPS_PROD:
                s2_list.append(kicked(cfg, gen, eps))
    # 7 classes x 2 generators x 4 eps = 56; pad with 40 random kicks -> 96
    for _ in range(40):
        Hr = rng.normal(size=(NW, NW)) + 1j * rng.normal(size=(NW, NW))
        Hr = (Hr + Hr.conj().T) / 2
        Hr /= np.linalg.norm(Hr, 2)
        s2_list.append(kicked(cfg, Hr, rng.uniform(0.05, 0.5)))
    out = {}
    for k, c in COSTS.items():
        ds = [deficit(c, cfg, s2) for s2 in s2_list]
        viol = sum(1 for d in ds if d < -1e-10)
        out[k] = {"n_cells": len(s2_list), "chain_violations": viol,
                  "baseline_deficit": baselines[k],
                  "min_deficit": float(np.min(ds))}
    return out


# ---------- T3: mechanism of the bimodal p -------------------------------------
def leg_costs(cfg, s2):
    """The two kicked legs separately (D_max), plus the relevant M operators."""
    c12 = s1.d_max(cfg["s1"], s2)            # s2 in the sigma slot
    c23 = s1.d_max(s2, cfg["s3"])            # s2 in the rho slot
    return c12, c23


def analytic_rho_slot_slope(cfg, G):
    """d/d eps log lmax(s3^{-1/2} s2(eps) s3^{-1/2}) at eps = 0 for simple lmax:
    = <v| s3^{-1/2} (i [G, mid]) s3^{-1/2} |v> / lmax."""
    ev3, V3 = np.linalg.eigh(cfg["s3"])
    ev3 = np.maximum(ev3, 1e-300)
    sih = V3 @ np.diag(ev3 ** -0.5) @ V3.conj().T
    M0 = sih @ cfg["mid"] @ sih
    M0 = (M0 + M0.conj().T) / 2
    ew, W = np.linalg.eigh(M0)
    lmax, v = ew[-1], W[:, -1]
    gap = float(ew[-1] - ew[-2])
    dmid = 1j * (G @ cfg["mid"] - cfg["mid"] @ G)
    dM = sih @ dmid @ sih
    slope = float(np.real(v.conj() @ dM @ v)) / lmax
    return slope, gap


def run_mechanism(gens):
    cfg = make_config(20260613, theta=0.3)
    K = cfg["K"]
    base12, base23 = leg_costs(cfg, cfg["mid"])
    out = {"base_leg_equality": float(abs(base12 - base23))}
    eps_all = EPS_SMALL + EPS_PROD
    table = {}
    for name, G in gens.items():
        entry = {}
        # weight-transfer witness: [K, G] = 0 iff the kick conserves K-weight
        entry["K_commutator_norm"] = float(np.linalg.norm(K @ G - G @ K, 2))
        # signed leg costs on the two-sided ladder
        legs = {}
        for eps in eps_all:
            for sgn in (+1.0, -1.0):
                c12, c23 = leg_costs(cfg, kicked(cfg, G, sgn * eps))
                legs[(sgn, eps)] = (c12, c23)
        # EVENNESS IDENTITY: for [G, K] = 0, conjugation by e^{-i eps G}
        # commutes with the flow, giving c12(eps) = c23(-eps) exactly for ANY
        # unitarily invariant cost; hence E(eps) = E(-eps) (even function,
        # quadratic leading order). Residual measures the broken symmetry.
        even_res = max(
            max(abs(legs[(+1.0, e)][0] - legs[(-1.0, e)][1]),
                abs(legs[(-1.0, e)][0] - legs[(+1.0, e)][1])) for e in eps_all)
        entry["leg_evenness_residual"] = float(even_res)
        # signed excess ladder (both signs)
        entry["E_rows"] = [
            {"eps": float(sgn * e),
             "E": float((legs[(sgn, e)][0] - base12)
                        + (legs[(sgn, e)][1] - base23))}
            for e in eps_all for sgn in (+1.0, -1.0)]
        # one-sided derivative estimates from the smallest rung
        e0 = EPS_SMALL[0]
        Ep = (legs[(+1.0, e0)][0] - base12) + (legs[(+1.0, e0)][1] - base23)
        Em = (legs[(-1.0, e0)][0] - base12) + (legs[(-1.0, e0)][1] - base23)
        Dp, Dm = Ep / e0, Em / (-e0)
        entry["D_plus"] = float(Dp)
        entry["D_minus"] = float(Dm)
        scale = max(abs(Dp), abs(Dm), 1e-300)
        if scale < 1e-2:
            kind = "quadratic"
        elif abs(Dp + Dm) / scale < 0.2:
            kind = "kink_abs"               # E ~ a |eps|
        elif abs(Dp - Dm) / scale < 0.2:
            kind = "smooth_linear"          # E ~ a eps (sign-carrying)
        else:
            kind = "mixed"
        entry["kind"] = kind
        entry["E_over_eps2_smallest"] = float(Ep / e0 ** 2)
        # analytic rho-slot first-order term + top gap: per-leg slopes do NOT
        # individually vanish for the quadratic classes -- they cancel between
        # the legs (the evenness identity), which is the selection mechanism
        slope23, gap23 = analytic_rho_slot_slope(cfg, G)
        entry["analytic_rho_slot_slope"] = float(slope23)
        entry["rho_slot_top_gap"] = float(gap23)
        entry["numerical_rho_slot_slope"] = float(
            (legs[(+1.0, e0)][1] - base23) / e0)
        table[name] = entry
    out["classes"] = table
    return out


# ---------- T4: wedge substrate -------------------------------------------------
def wedge_setup(n_max=3):
    mh = for_bisognano_wichmann(n_max)
    KW = np.real(np.diag(mh.restrict_K_alpha_to_wedge()))   # odd positive ints
    dW = len(KW)
    w = np.exp(-BETA * KW)
    rho = np.diag((w / np.sum(w)).astype(complex))
    return {"mh": mh, "KW": KW, "dim_full": len(mh.basis), "dim_W": dW,
            "rho": rho}


def wedge_flow_U(KW, t):
    return np.diag(np.exp(-1j * t * KW))


def wedge_kick_gens(KW, rng, n_per=3):
    """Hermitized matrix-unit generators graded by K-weight transfer Dw."""
    dW = len(KW)
    gens = {}
    for dw in (0, 2, 4):
        pairs = [(i, j) for i in range(dW) for j in range(i + 1, dW)
                 if abs(KW[i] - KW[j]) == dw]
        rng.shuffle(pairs)
        reps = []
        for (i, j) in pairs[:n_per]:
            G = np.zeros((dW, dW), dtype=complex)
            ph = np.exp(1j * rng.uniform(0, 2 * np.pi))
            G[i, j] = ph
            G[j, i] = np.conj(ph)
            reps.append(G / np.linalg.norm(G, 2))
        gens[dw] = reps
    return gens


def run_wedge(n_max=3, t_total=1.0):
    W = wedge_setup(n_max)
    KWd = np.diag(W["KW"]).astype(complex)
    dW = W["dim_W"]
    rng = np.random.default_rng(20260614)
    out = {"dim_full": W["dim_full"], "dim_W": dW,
           "KW_spectrum_min": float(np.min(W["KW"])),
           "KW_spectrum_max": float(np.max(W["KW"]))}

    # anchors: KMS at beta = 1, flow invariance, period (predicted pi on wedge)
    kms_res = 0.0
    wts = np.exp(-BETA * W["KW"])
    for _ in range(6):
        A = rng.normal(size=(dW, dW)) + 1j * rng.normal(size=(dW, dW))
        B = rng.normal(size=(dW, dW)) + 1j * rng.normal(size=(dW, dW))
        Bi = (wts[:, None] * B) * (1.0 / wts)[None, :]
        lhs = np.trace(W["rho"] @ A @ Bi)
        rhs = np.trace(W["rho"] @ B @ A)
        kms_res = max(kms_res, abs(lhs - rhs) / max(1.0, abs(rhs)))
    out["kms_residual"] = float(kms_res)
    out["flow_invariance"] = float(np.linalg.norm(
        s1.conj(wedge_flow_U(W["KW"], 1.234), W["rho"]) - W["rho"]))

    H0 = rng.normal(size=(dW, dW)) + 1j * rng.normal(size=(dW, dW))
    H0 = (H0 + H0.conj().T) / 2
    om0 = s1.conj(expm(0.3j * H0), W["rho"])
    out["orbit_period_pi_residual"] = float(s1.trace_dist(
        s1.conj(wedge_flow_U(W["KW"], np.pi), om0), om0))
    out["U_pi_is_minus_identity"] = float(np.max(np.abs(
        wedge_flow_U(W["KW"], np.pi) + np.eye(dW))))
    ts = np.linspace(0.05, np.pi - 0.05, 40)
    prof = [s1.trace_dist(s1.conj(wedge_flow_U(W["KW"], t), om0), om0)
            for t in ts]
    out["min_orbit_distance"] = float(np.min(prof))

    # D_max chain + excess scaling by weight-transfer grade
    cfg = {"K": KWd, "s1": om0,
           "mid": s1.conj(wedge_flow_U(W["KW"], t_total / 2), om0),
           "s3": s1.conj(wedge_flow_U(W["KW"], t_total), om0)}
    base = deficit(s1.d_max, cfg, cfg["mid"])
    out["baseline_deficit_dmax"] = float(base)
    gens = wedge_kick_gens(W["KW"], rng)
    base12 = s1.d_max(cfg["s1"], cfg["mid"])
    base23 = s1.d_max(cfg["mid"], cfg["s3"])
    table = {}
    viol, n_cells = 0, 0
    e0 = EPS_SMALL[0]
    for dw, reps in gens.items():
        rows = []
        for G in reps:
            ys = []
            for eps in EPS_PROD:
                d = deficit(s1.d_max, cfg, kicked(cfg, G, eps))
                n_cells += 1
                if d < -1e-10:
                    viol += 1
                ys.append(float(d - base))
            # small-eps one-sided derivatives + evenness residual (mechanism)
            def _legs(e):
                s2 = kicked(cfg, G, e)
                return s1.d_max(cfg["s1"], s2), s1.d_max(s2, cfg["s3"])
            lp, lm = _legs(e0), _legs(-e0)
            Ep = (lp[0] - base12) + (lp[1] - base23)
            Em = (lm[0] - base12) + (lm[1] - base23)
            even_res = max(abs(lp[0] - lm[1]), abs(lm[0] - lp[1]))
            row = {"K_commutator_norm": float(np.linalg.norm(
                       KWd @ G - G @ KWd, 2)),
                   "leg_evenness_residual": float(even_res),
                   "D_plus": float(Ep / e0), "D_minus": float(Em / (-e0)),
                   "excess_prod": ys}
            if all(y > 0 for y in ys):
                p, _ = np.polyfit(np.log(EPS_PROD), np.log(ys), 1)
                row["fit_power"] = float(p)
            rows.append(row)
        table[str(dw)] = rows
    # random kicks to mirror the panel design
    for _ in range(40):
        Hr = rng.normal(size=(dW, dW)) + 1j * rng.normal(size=(dW, dW))
        Hr = (Hr + Hr.conj().T) / 2
        Hr /= np.linalg.norm(Hr, 2)
        d = deficit(s1.d_max, cfg, kicked(cfg, Hr, rng.uniform(0.05, 0.5)))
        n_cells += 1
        if d < -1e-10:
            viol += 1
    out["chain"] = {"n_cells": n_cells, "dmax_chain_violations": viol}
    out["transfer_table"] = table

    # operational interval functional: recover Delta-tau from state pairs by
    # trace-distance matching; verify well-definedness + additivity (mod pi)
    def ell_hat(om_a, om_b, n_grid=600):
        taus = np.linspace(0.0, np.pi, n_grid, endpoint=False)
        d = [s1.trace_dist(s1.conj(wedge_flow_U(W["KW"], t), om_a), om_b)
             for t in taus]
        i0 = int(np.argmin(d))
        # parabolic refinement on the grid minimum
        im, ip = (i0 - 1) % n_grid, (i0 + 1) % n_grid
        a, b, c = d[im], d[i0], d[ip]
        denom = (a - 2 * b + c)
        shift = 0.5 * (a - c) / denom if abs(denom) > 1e-300 else 0.0
        return float((taus[i0] + shift * (np.pi / n_grid)) % np.pi)

    rec_errs, add_errs = [], []
    for _ in range(8):
        ta, dt1, dt2 = rng.uniform(0, np.pi, 3) * np.array([1.0, 0.3, 0.3])
        om_a = s1.conj(wedge_flow_U(W["KW"], ta), om0)
        om_b = s1.conj(wedge_flow_U(W["KW"], ta + dt1), om0)
        om_c = s1.conj(wedge_flow_U(W["KW"], ta + dt1 + dt2), om0)
        l_ab, l_bc, l_ac = ell_hat(om_a, om_b), ell_hat(om_b, om_c), \
            ell_hat(om_a, om_c)
        rec_errs.append(abs(l_ab - dt1))
        add_errs.append(abs((l_ab + l_bc - l_ac) % np.pi))
    out["interval_recovery_max_err"] = float(np.max(rec_errs))
    out["interval_additivity_max_err"] = float(np.max(
        [min(e, np.pi - e) for e in add_errs]))
    return out


def run():
    gens = class_gens()
    out = {}
    summary, _cells = run_ensemble(gens)
    out["T1_ensemble"] = summary
    out["T2_chain_panel"] = run_chain_panel(gens)
    out["T3_mechanism"] = run_mechanism(gens)
    out["T4_wedge"] = run_wedge()
    return out


if __name__ == "__main__":
    res = run()
    (ROOT / "data" / "wh7_b3_phase3_sprint2.json").write_text(
        json.dumps(res, indent=1, default=float), encoding="utf-8")

    print("T1 ensemble (frac positive excess overall | per class):")
    for k, sm in res["T1_ensemble"].items():
        pc = " ".join(f"{n.split()[0]}={d['frac_positive']:.2f}"
                      for n, d in sm["per_class"].items())
        print(f"   {k:9s}: overall {sm['frac_positive_overall']:.3f} | {pc}")
        print(f"             per theta: {sm['frac_positive_per_theta']}")
    print("T2 chain panel (violations / cells | baseline deficit, min deficit):")
    for k, d in res["T2_chain_panel"].items():
        print(f"   {k:9s}: {d['chain_violations']}/{d['n_cells']} | "
              f"base {d['baseline_deficit']:+.4f}, min {d['min_deficit']:+.4f}")
    print(f"T3 mechanism (base12 = base23 residual "
          f"{res['T3_mechanism']['base_leg_equality']:.2e}):")
    for name, e in res["T3_mechanism"]["classes"].items():
        print(f"   {name:24s}: kind={e['kind']:13s} D+={e['D_plus']:+.4f} "
              f"D-={e['D_minus']:+.4f} | [K,G]={e['K_commutator_norm']:.2e} "
              f"evenness={e['leg_evenness_residual']:.2e} | rho-slot "
              f"an={e['analytic_rho_slot_slope']:+.4f} "
              f"num={e['numerical_rho_slot_slope']:+.4f}")
    w = res["T4_wedge"]
    print(f"T4 wedge: dim {w['dim_full']} -> {w['dim_W']}; "
          f"KMS {w['kms_residual']:.2e}; flow-inv {w['flow_invariance']:.2e}; "
          f"period-pi {w['orbit_period_pi_residual']:.2e} "
          f"(U_pi+I {w['U_pi_is_minus_identity']:.2e})")
    print(f"   chain {w['chain']['dmax_chain_violations']}/{w['chain']['n_cells']}; "
          f"baseline deficit {w['baseline_deficit_dmax']:.3f}; "
          f"min orbit dist {w['min_orbit_distance']:.3f}")
    for dw, rows in w["transfer_table"].items():
        descr = ", ".join(
            (f"p={r['fit_power']:.2f}" if "fit_power" in r else "signed")
            + f"(D+={r['D_plus']:+.3f},ev={r['leg_evenness_residual']:.1e})"
            for r in rows)
        print(f"   Dw={dw}: {descr}")
    print(f"   interval recovery max err {w['interval_recovery_max_err']:.2e}; "
          f"additivity max err {w['interval_additivity_max_err']:.2e}")
