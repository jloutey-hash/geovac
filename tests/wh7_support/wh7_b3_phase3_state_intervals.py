# -*- coding: utf-8 -*-
"""B3 Phase 3, Sprint 1: state-level thermal-time cost structure on modular orbits
(2026-06-10, v3.115.0).

Design fixed by Phase 2's frozen negative: the reverse-triangle (twin-paradox)
content cannot come from seminorm legs (timelike sector positive-definite); it must
enter at the STATE level along modular orbits, with costs from entropy production.

Substrate: j <= 1 Peter-Weyl window (dim 14), boost K = 2 J_z (Paper 42 two_m_j
convention), wedge KMS state rho = e^{-beta K}/Z (beta = 1 in these units; the
KMS property itself is verified as T0, conventions documented in the memo).

Cost functional: Datta max-divergence D_max(rho||sigma) = log ||sigma^{-1/2} rho
sigma^{-1/2}||_op -- per the standing Paper 49 lesson (memory: Umegaki chain fails
~5% on cocycle triples; D_max chain D_max(1||3) <= D_max(1||2) + D_max(2||3) is a
theorem). Umegaki computed as comparison column.

Checks:
  T0  KMS anchor: <A sigma_{i beta}(B)>_rho = <B A>_rho bit-check; flow-invariance
      of rho; orbit period 2 pi.
  T1  orbit well-definedness: trace-distance profile t -> d(omega_0, omega_t) > 0
      on (0, 2 pi) for a generic perturbed state (intervals well-defined mod 2 pi).
  T2  D_max chain inequality on orbit triples with off-orbit kicks: deficit
      Delta = c12 + c23 - c13 >= 0 ALWAYS (theorem anchor); Umegaki violation
      count > 0 expected (reproduces the Paper 49 Datta fix on this substrate).
  T3  twin deficit: kick-attributable excess Delta(eps) - Delta(0) >= 0, scaling
      ~ eps^p with p approx 2.
  T4  cone-class organization: the quadratic coefficient C_class of the deficit,
      tabulated by the causal class (b, |m'|) of the kick generator -- does the
      Phase-2 grading (timelike / null / spacelike) organize the state-level costs?

Frozen falsifier: tests/test_wh7_b3_phase3.py.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import wh7_b1_joint_product_gh as b1  # noqa: E402
import wh7_b3_boost_seminorm_probe as b3  # noqa: E402

RNG = np.random.default_rng(20260613)
BETA = 1.0
NW = b1.NW


def boost_K():
    return 2.0 * b1.A_GENS[2]                       # K = 2 J_z, diagonal in PW basis


def kms_state(K, beta=BETA):
    w = np.exp(-beta * np.real(np.diag(K)))
    return np.diag(w / np.sum(w)).astype(complex)


def flow_U(K, t):
    return np.diag(np.exp(-1j * t * np.real(np.diag(K))))


def conj(U, rho):
    return U @ rho @ U.conj().T


def d_max(rho, sigma):
    """Datta max-divergence: log lambda_max(sigma^{-1/2} rho sigma^{-1/2})."""
    ev, V = np.linalg.eigh(sigma)
    ev = np.maximum(ev, 1e-300)
    s_inv_half = V @ np.diag(ev ** -0.5) @ V.conj().T
    M = s_inv_half @ rho @ s_inv_half
    return float(np.log(np.max(np.linalg.eigvalsh((M + M.conj().T) / 2))))


def umegaki(rho, sigma):
    def mlog(A):
        ev, V = np.linalg.eigh(A)
        ev = np.maximum(ev, 1e-300)
        return V @ np.diag(np.log(ev)) @ V.conj().T
    return float(np.real(np.trace(rho @ (mlog(rho) - mlog(sigma)))))


def trace_dist(a, b):
    ev = np.linalg.eigvalsh(a - b)
    return 0.5 * float(np.sum(np.abs(ev)))


def class_generator(b, row, col):
    """Hermitized compression of D^b_{m'm}: G = C + C^dagger, normalized."""
    D = b3.dcache(b)
    f = np.sqrt(2 * b + 1) * D[:, row, col]
    C = b1.compress_spatial(f)
    G = C + C.conj().T
    return G / np.linalg.norm(G, 2)


# representative kick generators by causal class (b, |m'|); col chosen mid-band
CLASSES = {
    "(0.5,0.5) spacelike": (0.5, 0, 0),
    "(1.0,0.0) spacelike": (1.0, 1, 1),
    "(2.0,0.0) spacelike": (2.0, 2, 2),
    "(2.0,1.0) spacelike": (2.0, 1, 2),
    "(1.0,1.0) null":      (1.0, 0, 1),    # the op-norm null ray C^1_{1,0}
    "(1.5,1.5) timelike":  (1.5, 0, 1),
    "(2.0,2.0) timelike":  (2.0, 0, 2),
}


def run():
    out = {}
    K = boost_K()
    rho = kms_state(K)

    # T0 — KMS anchor + flow invariance + period
    kms_res = 0.0
    for _ in range(6):
        A = RNG.normal(size=(NW, NW)) + 1j * RNG.normal(size=(NW, NW))
        B = RNG.normal(size=(NW, NW)) + 1j * RNG.normal(size=(NW, NW))
        # sigma_{i*beta}(B) = e^{-beta K} B e^{+beta K} for rho = e^{-beta K}/Z
        w = np.exp(-BETA * np.real(np.diag(K)))
        Bi = (w[:, None] * B) * (1.0 / w)[None, :]
        lhs = np.trace(rho @ A @ Bi)
        rhs = np.trace(rho @ B @ A)
        kms_res = max(kms_res, abs(lhs - rhs) / max(1.0, abs(rhs)))
    inv_res = float(np.linalg.norm(conj(flow_U(K, 1.234), rho) - rho))
    per_res = 0.0
    H0 = RNG.normal(size=(NW, NW)) + 1j * RNG.normal(size=(NW, NW))
    H0 = (H0 + H0.conj().T) / 2
    omega0 = conj(expm(0.3j * H0), rho)
    per_res = trace_dist(conj(flow_U(K, 2 * np.pi), omega0), omega0)
    out["T0"] = {"kms_residual": kms_res, "flow_invariance": inv_res,
                 "orbit_period_2pi_residual": per_res}

    # T1 — orbit injectivity profile
    ts = np.linspace(0.1, 2 * np.pi - 0.1, 40)
    prof = [trace_dist(conj(flow_U(K, t), omega0), omega0) for t in ts]
    out["T1"] = {"min_orbit_distance": float(np.min(prof)),
                 "max_orbit_distance": float(np.max(prof))}

    # T2/T3/T4 — twin deficits with class kicks
    t_total = 1.0
    s1 = omega0
    s3 = conj(flow_U(K, t_total), omega0)
    mid = conj(flow_U(K, t_total / 2), omega0)

    def deficits(s2):
        dmax = (d_max(s1, s2) + d_max(s2, s3) - d_max(s1, s3))
        du = (umegaki(s1, s2) + umegaki(s2, s3) - umegaki(s1, s3))
        return dmax, du

    base_max, base_u = deficits(mid)               # eps = 0 baseline (on-orbit midpoint)
    eps_list = [0.05, 0.1, 0.2, 0.4]
    # orbit tangent at the midpoint state (HS geometry): T_op = i [K, mid]
    T_op = 1j * (K @ mid - mid @ K)
    T_nrm2 = float(np.real(np.vdot(T_op, T_op)))

    def state_velocity(G):
        return 1j * (G @ mid - mid @ G)

    table = {}
    chain_viol_dmax, chain_viol_umegaki, n_cells = 0, 0, 0
    for name, (b, r, c) in CLASSES.items():
        G = class_generator(b, r, c)
        dω = state_velocity(G)
        lam = float(np.real(np.vdot(T_op, dω))) / T_nrm2
        overlap = abs(np.vdot(T_op, dω)) / max(
            1e-300, np.linalg.norm(dω) * np.sqrt(T_nrm2))
        Gp = G - lam * K                       # flow-tangent component removed
        Gp = Gp / np.linalg.norm(Gp, 2)
        entry = {"tangent_overlap": float(overlap)}
        for tag, gen in (("raw", G), ("perp", Gp)):
            rows = []
            for eps in eps_list:
                s2 = conj(expm(1j * eps * gen), mid)
                dmax, du = deficits(s2)
                n_cells += 1
                if dmax < -1e-10:
                    chain_viol_dmax += 1
                if du < -1e-10:
                    chain_viol_umegaki += 1
                rows.append({"eps": eps, "excess_dmax": dmax - base_max,
                             "excess_umegaki": du - base_u})
            xs = np.log(eps_list)
            ys = np.log(np.maximum([r_["excess_dmax"] for r_ in rows], 1e-300))
            p, logC = np.polyfit(xs, ys, 1)
            entry[tag] = {"rows": rows, "fit_power": float(p),
                          "fit_coeff": float(np.exp(logC))}
        table[name] = entry
    # random off-class kicks for the Umegaki-violation hunt (broader panel)
    for _ in range(40):
        Hr = RNG.normal(size=(NW, NW)) + 1j * RNG.normal(size=(NW, NW))
        Hr = (Hr + Hr.conj().T) / 2
        Hr /= np.linalg.norm(Hr, 2)
        s2 = conj(expm(1j * RNG.uniform(0.05, 0.5) * Hr), mid)
        dmax, du = deficits(s2)
        n_cells += 1
        if dmax < -1e-10:
            chain_viol_dmax += 1
        if du < -1e-10:
            chain_viol_umegaki += 1
    out["T2"] = {"n_cells": n_cells, "dmax_chain_violations": chain_viol_dmax,
                 "umegaki_chain_violations": chain_viol_umegaki,
                 "baseline_deficit_dmax": base_max, "baseline_deficit_umegaki": base_u}
    out["T3_T4_class_table"] = table

    # T5 — robustness of the p-split across reference states (3 fresh seeds,
    # two representative classes from each scaling group)
    rob = {}
    for seed in (101, 202, 303):
        rng = np.random.default_rng(seed)
        Hs = rng.normal(size=(NW, NW)) + 1j * rng.normal(size=(NW, NW))
        Hs = (Hs + Hs.conj().T) / 2
        om0 = conj(expm(0.3j * Hs), rho)
        s1s, s3s = om0, conj(flow_U(K, t_total), om0)
        mids = conj(flow_U(K, t_total / 2), om0)
        b0 = (d_max(s1s, mids) + d_max(mids, s3s) - d_max(s1s, s3s))
        for name in ("(2.0,0.0) spacelike", "(2.0,2.0) timelike"):
            b_, r_, c_ = CLASSES[name]
            G = class_generator(b_, r_, c_)
            ys = []
            for eps in eps_list:
                s2 = conj(expm(1j * eps * G), mids)
                d = (d_max(s1s, s2) + d_max(s2, s3s) - d_max(s1s, s3s)) - b0
                ys.append(max(d, 1e-300))
            p, _ = np.polyfit(np.log(eps_list), np.log(ys), 1)
            rob.setdefault(name, []).append(float(p))
    out["T5_robustness"] = rob
    return out


if __name__ == "__main__":
    res = run()
    (ROOT / "data" / "wh7_b3_phase3_state_intervals.json").write_text(
        json.dumps(res, indent=1, default=float), encoding="utf-8")
    t0 = res["T0"]
    print(f"T0 KMS residual {t0['kms_residual']:.2e}; flow-inv {t0['flow_invariance']:.2e}; "
          f"period {t0['orbit_period_2pi_residual']:.2e}")
    t1 = res["T1"]
    print(f"T1 orbit distance in ({t1['min_orbit_distance']:.3f}, "
          f"{t1['max_orbit_distance']:.3f}) on (0, 2pi)  -> injective orbit")
    t2 = res["T2"]
    print(f"T2 chain: D_max violations {t2['dmax_chain_violations']}/{t2['n_cells']}, "
          f"Umegaki violations {t2['umegaki_chain_violations']}/{t2['n_cells']}; "
          f"baseline deficit dmax {t2['baseline_deficit_dmax']:.3e}, "
          f"umegaki {t2['baseline_deficit_umegaki']:.3e}")
    print("T3/T4 kick-class table (excess D_max deficit ~ C * eps^p):")
    for name, d in res["T3_T4_class_table"].items():
        print(f"   {name:24s}: overlap = {d['tangent_overlap']:.3f} | "
              f"raw p = {d['raw']['fit_power']:+.2f} C = {d['raw']['fit_coeff']:.3e} | "
              f"perp p = {d['perp']['fit_power']:+.2f} C = {d['perp']['fit_coeff']:.3e}")
