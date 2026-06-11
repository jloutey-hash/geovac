# -*- coding: utf-8 -*-
"""Shared substrate library for the B3 Phase-3 band-exhaustion program
(2026-06-10). Built in main session; probe agents import from here and must
NOT modify this file.

Substrate: the FOLDED Peter-Weyl wedge at growing window j <= j_max.
 - Window labels (j, m', m), j = 0, 1/2, ..., j_max.
 - Multiplier compressions C^b_{mu',mu} built from exact Clebsch-Gordan
   coefficients (float-evaluated, construction validated bit-level against
   the Sprint-3b exact-arithmetic driver and the b1 quadrature machinery).
 - Plain-swap wedge reflection R: m' -> -m'; folding isometry V onto the +1
   eigenspace; unfolded wedge boost K_W = diag(2|m'|) >= 0 (positive
   generator -- the reason the wedge is the right substrate for growing
   windows: the unfolded K = 2 J_z is unbounded below and the KMS state
   would concentrate on the window edge).
 - Wedge KMS state rho = e^{-beta K_W}/Z (beta = 1 default).
 - Band exhaustion = growing j_max with FIXED multiplier-class data
   (reference perturbation and kicks are window compressions of fixed
   functions on SU(2) ~ S^3), i.e. the Paper-38 compression family at the
   state level.

Conventions: window label order is canonical (j ascending; m' descending
j..-j; m descending). Wedge label order: (j, m', m) with m' >= 0, same sort.
All costs/states are numpy complex matrices in the wedge basis.
"""
from functools import lru_cache
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import expm
from sympy.physics.quantum.cg import CG as _CG

ROOT = Path(__file__).resolve().parent
HALF = sp.Rational(1, 2)

# the seven Phase-2 causal classes as (b, mu', mu); QF/TRANSFER for reference
CLASSES = {
    "(0.5,0.5) spacelike": (HALF, HALF, HALF),
    "(1.0,0.0) spacelike": (sp.Integer(1), sp.Integer(0), sp.Integer(0)),
    "(2.0,0.0) spacelike": (sp.Integer(2), sp.Integer(0), sp.Integer(0)),
    "(2.0,1.0) spacelike": (sp.Integer(2), sp.Integer(1), sp.Integer(0)),
    "(1.0,1.0) null":      (sp.Integer(1), sp.Integer(1), sp.Integer(0)),
    "(1.5,1.5) timelike":  (sp.Rational(3, 2), sp.Rational(3, 2), HALF),
    "(2.0,2.0) timelike":  (sp.Integer(2), sp.Integer(2), sp.Integer(0)),
}

# standard window ladder for the program (agents may add 7/2 if runtime allows)
JMAX_LADDER = [sp.Integer(1), sp.Rational(3, 2), sp.Integer(2),
               sp.Rational(5, 2), sp.Integer(3)]


def jlist(jmax):
    out, j = [], sp.Integer(0)
    while j <= jmax:
        out.append(j)
        j += HALF
    return out


@lru_cache(maxsize=None)
def labels(jmax):
    labs = []
    for j in jlist(jmax):
        ms = [j - k for k in range(int(2 * j) + 1)]
        for mp in ms:
            for m in ms:
                labs.append((j, mp, m))
    return tuple(labs)


@lru_cache(maxsize=None)
def _cgf(a, al, b, be, c, ga):
    return float(_CG(a, al, b, be, c, ga).doit())


@lru_cache(maxsize=None)
def build_C(b, mu_p, mu, jmax):
    """Dense float matrix of C^b_{mu',mu} on the FULL window j <= jmax.
    Entry = sqrt((2j2+1)(2b+1)/(2j1+1)) <b mu'; j2 m2'|j1 m1'> <b mu; j2 m2|j1 m1>."""
    labs = labels(jmax)
    idx = {l: i for i, l in enumerate(labs)}
    N = len(labs)
    A = np.zeros((N, N))
    js = jlist(jmax)
    for (j2, m2p, m2) in labs:
        m1p, m1 = m2p + mu_p, m2 + mu
        for j1 in js:
            if abs(m1p) > j1 or abs(m1) > j1:
                continue
            if not (abs(b - j2) <= j1 <= b + j2):
                continue
            c1 = _cgf(b, mu_p, j2, m2p, j1, m1p)
            if c1 == 0.0:
                continue
            c2 = _cgf(b, mu, j2, m2, j1, m1)
            if c2 == 0.0:
                continue
            Nf = float(sp.sqrt((2 * j2 + 1) * (2 * b + 1) / (2 * j1 + 1)))
            A[idx[(j1, m1p, m1)], idx[(j2, m2p, m2)]] = Nf * c1 * c2
    return A


@lru_cache(maxsize=None)
def hermitized_C(b, mu_p, mu, jmax):
    """G = C + C^T on the full window (real symmetric), UNNORMALIZED."""
    A = build_C(b, mu_p, mu, jmax)
    return A + A.T


@lru_cache(maxsize=None)
def wedge(jmax):
    """Folding data: (wedge_labels, V, w) with V the (N x dW) isometry onto
    the R: m' -> -m' +1-eigenspace and w the unfolded weights 2|m'| >= 0."""
    labs = labels(jmax)
    idx = {l: i for i, l in enumerate(labs)}
    wlabs, cols, w = [], [], []
    for l in labs:
        j, mp, m = l
        if mp < 0:
            continue
        v = np.zeros(len(labs))
        if mp == 0:
            v[idx[l]] = 1.0
        else:
            v[idx[l]] = v[idx[(j, -mp, m)]] = 1.0 / np.sqrt(2)
        wlabs.append(l)
        cols.append(v)
        w.append(2.0 * float(mp))
    return tuple(wlabs), np.array(cols).T, np.array(w)


def fold(M, jmax):
    """V^dag M V: wedge-block compression of a full-window operator."""
    _, V, _ = wedge(jmax)
    return V.T @ M @ V


@lru_cache(maxsize=None)
def class_gen_folded(name, jmax, normalized=False):
    """Folded Hermitized class generator (raw by default; the raw compression
    is the honest fixed-continuum-object choice -- record norms separately)."""
    b, mp, mu = CLASSES[name]
    GW = fold(hermitized_C(b, mp, mu, jmax), jmax)
    if normalized:
        nrm = np.linalg.norm(GW, 2)
        return GW / nrm if nrm > 1e-14 else GW
    return GW


def kms_state(jmax, beta=1.0):
    _, _, w = wedge(jmax)
    z = np.exp(-beta * w)
    return np.diag((z / np.sum(z)).astype(complex)), float(np.sum(z))


def flow_U(jmax, t):
    _, _, w = wedge(jmax)
    return np.diag(np.exp(-1j * t * w))


def conj(U, rho):
    return U @ rho @ U.conj().T


def d_max(rho, sigma):
    """Datta max-divergence log lmax(sigma^{-1/2} rho sigma^{-1/2})."""
    ev, V = np.linalg.eigh(sigma)
    ev = np.maximum(ev, 1e-300)
    sih = V @ np.diag(ev ** -0.5) @ V.conj().T
    M = sih @ rho @ sih
    return float(np.log(np.max(np.linalg.eigvalsh((M + M.conj().T) / 2))))


def trace_dist(a, b):
    ev = np.linalg.eigvalsh(a - b)
    return 0.5 * float(np.sum(np.abs(ev)))


def make_config(jmax, H_ref, theta=0.3, t_total=1.0, beta=1.0):
    """Orbit triple from the FIXED reference perturbation H_ref (wedge-basis
    Hermitian, typically a folded multiplier): om0 = e^{i theta H} rho e^{-...}."""
    rho, Z = kms_state(jmax, beta)
    om0 = conj(expm(1j * theta * H_ref), rho)
    return {"jmax": jmax, "Z": Z,
            "s1": om0,
            "mid": conj(flow_U(jmax, t_total / 2), om0),
            "s3": conj(flow_U(jmax, t_total), om0)}


def deficit(cfg, s2):
    return (d_max(cfg["s1"], s2) + d_max(s2, cfg["s3"])
            - d_max(cfg["s1"], cfg["s3"]))


def kicked(cfg, G, eps):
    return conj(expm(1j * eps * G), cfg["mid"])


def reference_H(jmax, which="null"):
    """The two standard FIXED reference perturbations (folded multipliers):
    'null'  = folded Hermitized C^1_{1,0} (breaks flow-invariance at every
              window: nonzero K_W commutator);
    'mixed' = folded Hermitized C^{1/2}_{1/2,1/2}."""
    if which == "null":
        return class_gen_folded("(1.0,1.0) null", jmax)
    if which == "mixed":
        return class_gen_folded("(0.5,0.5) spacelike", jmax)
    raise ValueError(which)


def selftest():
    """Validate against the Sprint-3b EXACT results at j_max = 1 and the
    window-edge theorem at 3/2. Returns dict of residuals/booleans."""
    out = {}
    jm = sp.Integer(1)
    exact_ratios = {"(1.0,0.0) spacelike": 6 / 19, "(2.0,0.0) spacelike": 5 / 6,
                    "(2.0,1.0) spacelike": 0.0, "(1.0,1.0) null": 11 / 19,
                    "(2.0,2.0) timelike": 0.5,
                    "(0.5,0.5) spacelike": 3 / 8, "(1.5,1.5) timelike": 0.25}
    dev = 0.0
    for name, target in exact_ratios.items():
        b, mp, mu = CLASSES[name]
        G = hermitized_C(b, mp, mu, jm)
        GW = fold(G, jm)
        r = np.sum(GW ** 2) / np.sum(G ** 2)
        dev = max(dev, abs(r - target))
    out["fold_ratio_max_dev_vs_3b_exact"] = float(dev)
    # commutation table: mu'=0 commute at 1 and 3/2; (2,2) only at 1
    def comm(name, jmx):
        _, _, w = wedge(jmx)
        GW = class_gen_folded(name, jmx)
        K = np.diag(w)
        return float(np.linalg.norm(K @ GW - GW @ K, 2))
    out["mu0_commute_1"] = max(comm("(1.0,0.0) spacelike", jm),
                               comm("(2.0,0.0) spacelike", jm))
    jm32 = sp.Rational(3, 2)
    out["mu0_commute_32"] = max(comm("(1.0,0.0) spacelike", jm32),
                                comm("(2.0,0.0) spacelike", jm32))
    out["c22_commute_1"] = comm("(2.0,2.0) timelike", jm)
    out["c22_commute_32"] = comm("(2.0,2.0) timelike", jm32)   # nonzero
    out["c21_norm_1"] = float(np.linalg.norm(
        class_gen_folded("(2.0,1.0) spacelike", jm), 2))        # zero
    out["c21_norm_32"] = float(np.linalg.norm(
        class_gen_folded("(2.0,1.0) spacelike", jm32), 2))      # nonzero
    out["dims"] = {str(j): (len(labels(j)), len(wedge(j)[0]))
                   for j in JMAX_LADDER}
    return out


if __name__ == "__main__":
    st = selftest()
    print(f"fold ratios vs 3b exact : {st['fold_ratio_max_dev_vs_3b_exact']:.2e}")
    print(f"mu'=0 [K,G] at 1 / 3/2  : {st['mu0_commute_1']:.2e} / "
          f"{st['mu0_commute_32']:.2e}")
    print(f"(2,2) [K,G] at 1 / 3/2  : {st['c22_commute_1']:.2e} / "
          f"{st['c22_commute_32']:.3f}  (zero then NONZERO)")
    print(f"(2,1) norm at 1 / 3/2   : {st['c21_norm_1']:.2e} / "
          f"{st['c21_norm_32']:.3f}  (zero then NONZERO)")
    print("dims (full, wedge):", st["dims"])
