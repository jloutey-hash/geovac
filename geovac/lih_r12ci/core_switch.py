r"""Parameter-free Li-core switch w(r) for the LiH electron-nucleus (e-n) columns.

Backs (Paper 19, LiH explicit-r12 marriage): the core switch used by en_factors.py (analytic) and by
the VMC linear-method comparator (debug/en_factors_stage0_lm.py).  PI RULING (2026-10-05): the switch
carries NO free radius; it is DEFINED by the reference itself:

    w(r) = |phi_c(r)|^2 / |phi_c(R_Li)|^2,

phi_c = the core natural orbital of the truncated M=16 pi reference = the largest-occupation natural
orbital that is (i) sigma (m = 0) and (ii) Li-centred (|phi(R_Li)| > |phi(R_H)|).  The NOs come
ordered by occupation from pi_m16.build_ref_from_prims, so phi_c is NO 0 in practice (occupation
~1.999); the selection is computed, not chosen.  w = 1 at the Li nucleus by construction, it is the
core density shape (no radius), axially symmetric (phi_c is sigma), so it carries no circulation.

Closed-form derivatives from the NO (phi_c real at phi_azi = 0, m = 0):
    grad w = 2 phi_c grad phi_c / phi0^2,     lap w = 2 (phi_c lap phi_c + |grad phi_c|^2) / phi0^2.
lap w carries the 1/r_Li cusp singularity of lap phi_c (integrable, like V_ne).

Light module: numpy + geovac.lih_vmc only (no grid / kernel imports).
"""
from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np

import geovac.lih_vmc as V


def nucleus_values(prims: Sequence[dict], T_no: np.ndarray, R: float) -> Tuple[np.ndarray, np.ndarray]:
    """phi_p at the Li (z=-R/2) and H (z=+R/2) nuclei for every NO column p of T_no (real parts)."""
    pts = np.array([[0.0, 0.0, -R / 2.0], [0.0, 0.0, R / 2.0]])
    Phi = V.eval_values_batch(prims, pts, R)                 # (Nprim, 2)
    no = np.einsum('ap,ae->pe', T_no, Phi)                   # (Mk, 2)
    return no[:, 0].real, no[:, 1].real


def select_core_no(occ: np.ndarray, m_no: np.ndarray, prims: Sequence[dict], T_no: np.ndarray,
                   R: float) -> int:
    """Index of the largest-occupation sigma, Li-centred natural orbital (the core NO)."""
    vLi, vH = nucleus_values(prims, T_no, R)
    order = np.argsort(-np.asarray(occ), kind='stable')
    for p in order:
        if int(m_no[p]) == 0 and abs(vLi[p]) > abs(vH[p]):
            return int(p)
    raise RuntimeError("no sigma Li-centred natural orbital found")


class CoreSwitch:
    """w(r) = phi_c(r)^2 / phi_c(R_Li)^2 from the primitive expansion tcore (Nprim,) of the core NO."""

    def __init__(self, prims: Sequence[dict], tcore: np.ndarray, R: float):
        self.prims = list(prims)
        self.tcore = np.asarray(tcore, dtype=float)
        self.R = float(R)
        vLi, _ = nucleus_values(self.prims, self.tcore[:, None], self.R)
        self.phi0 = float(vLi[0])
        self.inv0 = 1.0 / self.phi0 ** 2

    @classmethod
    def from_reference(cls, prims, T_no, occ, m_no, R) -> "CoreSwitch":
        p = select_core_no(occ, m_no, prims, T_no, R)
        sw = cls(prims, T_no[:, p], R)
        sw.no_index = p
        return sw

    def value(self, r: np.ndarray) -> np.ndarray:
        Phi = V.eval_values_batch(self.prims, r, self.R)
        ph = np.tensordot(self.tcore, Phi, axes=(0, 0)).real
        return ph * ph * self.inv0

    def vgl(self, r: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """w, grad w (r.shape), lap w at points r (..., 3)."""
        ph = np.zeros(r.shape[:-1]); gph = np.zeros(r.shape); lph = np.zeros(r.shape[:-1])
        for c, p in zip(self.tcore, self.prims):
            if c == 0.0:
                continue
            v, g, l = V.orbital_vgl(p, r, self.R)
            ph += c * v.real; gph += c * g.real; lph += c * l.real
        w = ph * ph * self.inv0
        gw = 2.0 * ph[..., None] * gph * self.inv0
        lw = 2.0 * (ph * lph + np.sum(gph * gph, axis=-1)) * self.inv0
        return w, gw, lw


class CoreFractionSwitch:
    """VARIANT (comparison only): w(r) = n_c |phi_c|^2 / sum_p n_p |phi_p|^2 -- the core NO's share of
    the full-CI reference density (n_p = NO occupations).  Parameter-free; w -> ~1 inside the core,
    -> 0 in the valence.  Quotient-rule derivatives from the orbital value/grad/lap."""

    def __init__(self, prims: Sequence[dict], T_no: np.ndarray, occ: np.ndarray, core: int, R: float):
        self.prims = list(prims); self.T = np.asarray(T_no, float); self.occ = np.asarray(occ, float)
        self.core = int(core); self.R = float(R)

    def _orbs(self, r):
        M = len(self.prims)
        Phi = np.empty((M,) + r.shape[:-1], complex); G = np.empty((M,) + r.shape, complex)
        L = np.empty((M,) + r.shape[:-1], complex)
        for a, p in enumerate(self.prims):
            Phi[a], G[a], L[a] = V.orbital_vgl(p, r, self.R)
        ph = np.tensordot(self.T, Phi, axes=(0, 0)); gph = np.tensordot(self.T, G, axes=(0, 0))
        lph = np.tensordot(self.T, L, axes=(0, 0))
        return ph, gph, lph

    def vgl(self, r: np.ndarray):
        ph, gph, lph = self._orbs(r)
        dens = np.abs(ph) ** 2                                            # (Mk, ...)
        gd = 2.0 * np.real(np.conj(ph)[..., None] * gph)
        ld = 2.0 * np.real(np.conj(ph) * lph) + 2.0 * np.sum(np.abs(gph) ** 2, axis=-1)
        n = self.occ.reshape((-1,) + (1,) * (ph.ndim - 1))
        rho = np.sum(n * dens, 0); grho = np.sum(n[..., None] * gd, 0); lrho = np.sum(n * ld, 0)
        a = self.occ[self.core] * dens[self.core]; ga = self.occ[self.core] * gd[self.core]
        la = self.occ[self.core] * ld[self.core]
        w = a / rho
        gw = (ga - w[..., None] * grho) / rho[..., None]
        lw = (la - 2.0 * np.sum(gw * grho, -1) - w * lrho) / rho
        return w, gw, lw

    def value(self, r: np.ndarray) -> np.ndarray:
        return self.vgl(r)[0]


def grid_fraction_switch(nval: np.ndarray, ngrad: np.ndarray, nlap: np.ndarray, occ: np.ndarray,
                         core: int):
    """grid version of CoreFractionSwitch (w, lap w) from the pi_m16 NO grid rows (phi_azi = 0)."""
    ph = np.asarray(nval); gph = np.asarray(ngrad); lph = np.asarray(nlap)
    dens = np.abs(ph) ** 2
    gd = 2.0 * np.real(np.conj(ph)[..., None] * gph)
    ld = 2.0 * np.real(np.conj(ph) * lph) + 2.0 * np.sum(np.abs(gph) ** 2, axis=-1)
    n = np.asarray(occ, float)[:, None]
    rho = np.sum(n * dens, 0); grho = np.sum(n[..., None] * gd, 0); lrho = np.sum(n * ld, 0)
    a = occ[core] * dens[core]; ga = occ[core] * gd[core]; la = occ[core] * ld[core]
    w = a / rho
    gw = (ga - w[:, None] * grho) / rho[:, None]
    lw = (la - 2.0 * np.sum(gw * grho, -1) - w * lrho) / rho
    return w, lw


def grid_switch(nval_c: np.ndarray, ngrad_c: np.ndarray, nlap_c: np.ndarray, phi0: float):
    """w and lap w on the kernels grid from the core NO's grid rows (pi_m16 info nval/ngrad/nlap)."""
    ph = np.asarray(nval_c).real
    g2 = np.sum(np.abs(ngrad_c) ** 2, axis=-1)
    w = ph * ph / phi0 ** 2
    lw = 2.0 * (ph * np.asarray(nlap_c).real + g2) / phi0 ** 2
    return w, lw
