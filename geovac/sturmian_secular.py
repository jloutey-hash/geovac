"""Isoenergetic generalized-Sturmian secular equation for He, extended to l>0 (s,p,d,f).

Tracked port of the ``debug/sturmian_he_lmax.py`` driver (Paper 60, group2). This
module backs three load-bearing Paper 60 claims, regression-protected by
``tests/test_sturmian_secular.py``:

  (i)   single-config 1s^2  ->  E = -2.847 Ha  (textbook variational He);
  (ii)  the entrywise 1-norm ``||M||_1`` grows SUBLINEARLY with the config count K,
        approaching ``~K^0.84`` for the full s+p+d+f basis (Paper 60 ``eq:sublinear``);
  (iii) restoring the L2 overlap metric S (generalized eigenproblem ``M B = p S B``)
        is ill-conditioned -- ``cond(S)`` climbs from ~4 into the thousands (the
        paper's "4 -> 3673");
  (iv)  ``[eq:secular]`` the interelectron matrix T' is a matrix of PURE NUMBERS,
        independent of the nuclear charge Z.

Physics
-------
Secular equation [Avery BK6 6.35]::

    [ diag(Z * R_nu) + T' - p_kappa * I ] B = 0 ,     E = -p_kappa^2 / 2 .

``R_nu = sqrt(sum_j 1/n_j^2)``. The weighted charge is built at the reference value
``p_kappa = 1``, so each configuration's hydrogenic radial orbitals sit at charge
``Q_nu = 1 / R_nu`` -- NOT at Z. Consequently T' is a matrix of pure numbers,
independent of both ``p_kappa`` and Z [BK6 boxed remark]::

    T'_{ij} = -<Psi_i | 1/r12 | Psi_j>   (mixed-scale, non-orthogonal).

The metric-free *standard* eigenproblem is the paper's atomic construction (no
overlap metric S in the primary method); the L2 overlap metric is reintroduced only
diagnostically to exhibit its ill-conditioning (claim iii).

Configurations are Goscinskian two-electron singlets (S=0) of hydrogenic orbitals
``(n_a, l)(n_b, l)`` coupled to total L=0 (equal l on the two electrons is required
for L=0), spatially symmetric.

Interelectron matrix elements use the Slater-Condon Legendre multipole expansion of
``1/r12``: angular factors are products of Gaunt integrals (Wigner 3j), radial
factors are Slater integrals ``R^k``.

Grid note: the radial grid (``R_MAX``, ``N_GRID``) is calibrated so the
cumulative-trapezoid (O(dr^2)) Slater potentials are converged -- ``(5/8) Q`` and
``E(1s^2) = -2.84766`` to 5-6 digits. Do not change it without re-validating the
gate numbers.
"""
from __future__ import annotations

import math
import warnings
from itertools import combinations_with_replacement
from math import factorial as fac
from typing import Dict, List, Tuple, Union

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.linalg import eigh
from scipy.special import genlaguerre

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------------------
# Radial grid.  Small-r resolution matters (high-n high-charge configs contract).
# Slater potentials use cumulative-TRAPEZOID (O(dr^2)); validated exact on (5/8)Q
# (grid-converged already at N~12000: (5/8)/sqrt2 to 6 digits, E(1s^2)=-2.84766).
# --------------------------------------------------------------------------------------
R_MAX: float = 60.0
N_GRID: int = 18000
r: np.ndarray = np.linspace(1e-7, R_MAX, N_GRID)
dr: float = r[1] - r[0]
r2: np.ndarray = r * r

# Default nuclear charge (helium).  Isoenergetic secular construction; Z enters ONLY
# through the diagonal diag(Z * R_nu), never through T' (see build_Tprime).
Z_HE: float = 2.0

# Module caches / id counter (reset by build_configs before each solve).
_RK_CACHE: Dict[tuple, float] = {}
_GAUNT_CACHE: Dict[tuple, float] = {}
_NEXT_RID: int = 0


def reset_caches() -> None:
    """Clear the radial Slater-integral cache and reset the radial-id counter.

    The Gaunt cache is pure angular data (no configuration dependence) and is left
    persistent across solves. Called by :func:`build_configs` before each build so
    that radial ids are deterministic and the ``R^k`` cache does not leak between
    unrelated configuration sets.
    """
    global _NEXT_RID
    _RK_CACHE.clear()
    _NEXT_RID = 0


def _ctrap_fwd(y: np.ndarray) -> np.ndarray:
    """Forward cumulative integral ``int_{r0}^{r_i} y dr`` (same length as ``y``, 0 at i=0)."""
    return np.concatenate(([0.0], cumulative_trapezoid(y, dx=dr)))


def _ctrap_rev(y: np.ndarray) -> np.ndarray:
    """Reverse cumulative integral ``int_{r_i}^{rmax} y dr`` (same length as ``y``)."""
    return _ctrap_fwd(y[::-1])[::-1]


# --------------------------------------------------------------------------------------
# Angular machinery: self-contained Wigner-3j (Racah) + real-Y Gaunt integral.
# --------------------------------------------------------------------------------------
def wigner3j(j1: int, j2: int, j3: int, m1: int, m2: int, m3: int) -> float:
    """Wigner 3j symbol via the Racah single-sum formula (integer arguments)."""
    if m1 + m2 + m3 != 0:
        return 0.0
    if not (abs(j1 - j2) <= j3 <= j1 + j2):
        return 0.0
    if any(abs(m) > j for m, j in ((m1, j1), (m2, j2), (m3, j3))):
        return 0.0
    if (j1 + j2 + j3) < 0:
        return 0.0
    delta = math.sqrt(
        fac(j1 + j2 - j3) * fac(j1 - j2 + j3) * fac(-j1 + j2 + j3)
        / fac(j1 + j2 + j3 + 1)
    )
    pref = math.sqrt(
        fac(j1 + m1) * fac(j1 - m1) * fac(j2 + m2) * fac(j2 - m2)
        * fac(j3 + m3) * fac(j3 - m3)
    )
    tmin = max(0, j2 - j3 - m1, j1 - j3 + m2)
    tmax = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    s = 0.0
    for t in range(tmin, tmax + 1):
        denom = (fac(t) * fac(j1 + j2 - j3 - t) * fac(j1 - m1 - t)
                 * fac(j2 + m2 - t) * fac(j3 - j2 + m1 + t) * fac(j3 - j1 - m2 + t))
        s += (-1) ** t / denom
    return (-1) ** (j1 - j2 - m3) * delta * pref * s


def gaunt(l1: int, l2: int, l3: int, m1: int, m2: int, m3: int) -> float:
    """Gaunt integral ``int Y_{l1 m1} Y_{l2 m2} Y_{l3 m3} dOmega`` (no conjugation)."""
    key = (l1, l2, l3, m1, m2, m3)
    v = _GAUNT_CACHE.get(key)
    if v is not None:
        return v
    w0 = wigner3j(l1, l2, l3, 0, 0, 0)
    if w0 == 0.0:
        _GAUNT_CACHE[key] = 0.0
        return 0.0
    val = (math.sqrt((2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) / (4 * math.pi))
           * w0 * wigner3j(l1, l2, l3, m1, m2, m3))
    _GAUNT_CACHE[key] = val
    return val


def cg_L0(l: int, m: int) -> float:
    """Clebsch-Gordan ``<l m; l -m | 0 0> = (-1)^{l-m} / sqrt(2l+1)``."""
    return (-1) ** (l - m) / math.sqrt(2 * l + 1)


# --------------------------------------------------------------------------------------
# Radial orbitals and Slater integrals.
# --------------------------------------------------------------------------------------
def hyd_radial(n: int, l: int, Q: float) -> np.ndarray:
    """Hydrogenic radial ``R_{nl}`` at charge ``Q``, L2-normalized on the grid.

    Normalized so that ``int R^2 r^2 dr = 1`` on the module grid.
    """
    a = Q / n
    f = (2 * a * r) ** l * np.exp(-a * r) * genlaguerre(n - l - 1, 2 * l + 1)(2 * a * r)
    nrm = np.sqrt(np.trapezoid(f * f * r2, r))
    return f / nrm


def radial_overlap(Pa: np.ndarray, Pb: np.ndarray) -> float:
    """Radial overlap ``int Pa Pb r^2 dr`` on the module grid."""
    return float(np.trapezoid(Pa * Pb * r2, r))


def slater_Rk(ida: int, idc: int, idb: int, idd: int,
              Pa: np.ndarray, Pc: np.ndarray, Pb: np.ndarray, Pd: np.ndarray,
              k: int) -> float:
    """Radial Slater integral ``R^k(ac ; bd)`` with the multipole potential of the
    electron-2 density built by cumulative trapezoid.

    ``R^k = int int Pa(r1)Pc(r1) (r_<^k / r_>^{k+1}) Pb(r2)Pd(r2) r1^2 r2^2 dr1 dr2``.
    Cached on a canonical key that folds the ``(a,c)``/``(b,d)`` and electron-swap
    symmetries so radial functions shared across m-values are deduplicated.
    """
    p1 = (ida, idc) if ida <= idc else (idc, ida)
    p2 = (idb, idd) if idb <= idd else (idd, idb)
    key = (p1, p2, k) if p1 <= p2 else (p2, p1, k)
    v = _RK_CACHE.get(key)
    if v is not None:
        return v
    g = Pb * Pd * r2                                   # electron-2 density * r2
    inner = _ctrap_fwd(g * r ** k)                     # int_0^{r1} dens2 r2^{k+2} dr2
    outer = _ctrap_rev(g * r ** (-(k + 1)))            # int_{r1}^inf dens2 r2^{1-k} dr2
    Uk = inner * r ** (-(k + 1)) + outer * r ** k
    val = float(np.trapezoid(Pa * Pc * Uk * r2, r))
    _RK_CACHE[key] = val
    return val


# --------------------------------------------------------------------------------------
# Orbital representation and coupled two-electron configuration wavefunctions.
# An "orbital" is a dict {rid, l, m, P (radial array)}.  Radial id encodes (config,
# which-orb) so the Slater-integral cache dedups across m (radial depends only on
# n, l, Q, not m).  A configuration term list is [(coeff, orb_e1, orb_e2), ...]
# representing Psi^unnorm(1, 2).
# --------------------------------------------------------------------------------------
class Config:
    """A Goscinskian two-electron ``(n_a, l)(n_b, l)`` singlet coupled to total L=0.

    The hydrogenic radial orbitals are built at the reference weighted charge
    ``Q = pk_ref / R_nu`` (NOT at Z), so the configuration -- and hence T' -- carries
    no nuclear-charge dependence. Z enters the secular matrix only through the
    diagonal ``Z * R_nu`` assembled in :func:`build_M`.
    """

    def __init__(self, l: int, na: int, nb: int, pk_ref: float = 1.0) -> None:
        global _NEXT_RID
        self.l, self.na, self.nb = l, na, nb
        self.Rnu: float = math.sqrt(1.0 / na ** 2 + 1.0 / nb ** 2)
        self.Q: float = pk_ref / self.Rnu
        # two radial orbitals (distinct radial ids); if na==nb they are the same fn
        self.Pa: np.ndarray = hyd_radial(na, l, self.Q)
        self.rid_a: int = _NEXT_RID
        _NEXT_RID += 1
        if nb == na:
            self.Pb: np.ndarray = self.Pa
            self.rid_b: int = self.rid_a
        else:
            self.Pb = hyd_radial(nb, l, self.Q)
            self.rid_b = _NEXT_RID
            _NEXT_RID += 1
        self.terms: List[Tuple[float, dict, dict]] = self._build_terms()
        self.norm: float = 1.0 / math.sqrt(self._self_overlap())

    def _orb(self, which: str, m: int) -> dict:
        if which == 'a':
            return dict(rid=self.rid_a, l=self.l, m=m, P=self.Pa)
        return dict(rid=self.rid_b, l=self.l, m=m, P=self.Pb)

    def _build_terms(self) -> List[Tuple[float, dict, dict]]:
        """``Psi^unnorm(1, 2)`` coupled to L=0, spatially symmetric.

        Same orbital (na==nb): ``Phi_aa = sum_m cg * phi_{a m}(1) phi_{a -m}(2)``
        (already symmetric). Distinct: symmetric combination ``Phi_ab + Phi_ba``.
        """
        l = self.l
        terms: List[Tuple[float, dict, dict]] = []
        for m in range(-l, l + 1):
            c = cg_L0(l, m)
            terms.append((c, self._orb('a', m), self._orb('b', -m)))   # Phi_ab
        if self.na != self.nb:
            for m in range(-l, l + 1):
                c = cg_L0(l, m)
                terms.append((c, self._orb('b', m), self._orb('a', -m)))  # Phi_ba
        return terms

    def _self_overlap(self) -> float:
        return overlap_terms(self.terms, self.terms)


# --------------------------------------------------------------------------------------
# Two-electron pair Coulomb primitive <phi_a(1) phi_b(2)| 1/r12 |phi_c(1) phi_d(2)>.
#   = sum_k (4 pi / (2k+1)) R^k(ac, bd) sum_q (-1)^{m_a+q+m_b}
#         gaunt(l_a, k, l_c, -m_a, -q, m_c) gaunt(l_b, k, l_d, -m_b, q, m_d)
# --------------------------------------------------------------------------------------
def pair_coulomb(oa: dict, ob: dict, oc: dict, od: dict) -> float:
    """Pair Coulomb primitive ``<phi_a(1) phi_b(2) | 1/r12 | phi_c(1) phi_d(2)>``."""
    la, ma = oa['l'], oa['m']
    lb, mb = ob['l'], ob['m']
    lc, mc = oc['l'], oc['m']
    ld, md = od['l'], od['m']
    if (ma + mb) != (mc + md):
        return 0.0
    kmax = min(la + lc, lb + ld)
    kmin = max(abs(la - lc), abs(lb - ld))
    total = 0.0
    for k in range(kmin, kmax + 1):
        # parity from (l_a k l_c;000): needs la+k+lc even; gaunt handles the zeros
        # NOTE (2026-08-29 wrong-sign-q audit): this q = mc - ma is CORRECT
        # here -- it is NEGATED at the call site (gaunt(..., -ma, -q, mc), so
        # the m-arguments sum to -ma - (mc-ma) + mc = 0).  Do NOT "fix" it in
        # a mechanical sweep; flipping it would INTRODUCE the bug this audit
        # removed elsewhere.  See debug/sprint_eri_evaluator_defects_memo.md.
        q = mc - ma
        g1 = gaunt(la, k, lc, -ma, -q, mc)
        if g1 == 0.0:
            continue
        g2 = gaunt(lb, k, ld, -mb, q, md)
        if g2 == 0.0:
            continue
        ang = (-1) ** (ma + q + mb) * g1 * g2 * (4 * math.pi / (2 * k + 1))
        Rk = slater_Rk(oa['rid'], oc['rid'], ob['rid'], od['rid'],
                       oa['P'], oc['P'], ob['P'], od['P'], k)
        total += ang * Rk
    return total


def overlap_terms(termsA: List[Tuple[float, dict, dict]],
                  termsB: List[Tuple[float, dict, dict]]) -> float:
    """``<Psi_A^unnorm | Psi_B^unnorm> = sum W_A W_B <u_A|u_B><v_A|v_B>``."""
    tot = 0.0
    for (wa, ua, va) in termsA:
        for (wb, ub, vb) in termsB:
            if ua['l'] != ub['l'] or ua['m'] != ub['m']:
                continue
            if va['l'] != vb['l'] or va['m'] != vb['m']:
                continue
            su = 1.0 if ua['rid'] == ub['rid'] else radial_overlap(ua['P'], ub['P'])
            sv = 1.0 if va['rid'] == vb['rid'] else radial_overlap(va['P'], vb['P'])
            tot += wa * wb * su * sv
    return tot


def repulsion_terms(termsA: List[Tuple[float, dict, dict]],
                    termsB: List[Tuple[float, dict, dict]]) -> float:
    """``<Psi_A^unnorm | 1/r12 | Psi_B^unnorm>``."""
    tot = 0.0
    for (wa, ua, va) in termsA:
        for (wb, ub, vb) in termsB:
            g = pair_coulomb(ua, va, ub, vb)
            if g != 0.0:
                tot += wa * wb * g
    return tot


# --------------------------------------------------------------------------------------
# Secular-matrix assembly.  M = diag(Z * R_nu) + T',  T'_{ij} = -<Psi_i|1/r12|Psi_j>.
# T' is separated out so its Z-independence (Paper 60 eq:secular) is directly testable.
# --------------------------------------------------------------------------------------
def build_Tprime(configs: List[Config]) -> np.ndarray:
    """Interelectron matrix ``T'_{ij} = -<Psi_i | 1/r12 | Psi_j>`` -- a matrix of
    PURE NUMBERS, independent of the nuclear charge Z.

    Z-independence is structural: the configurations' radial orbitals are built at
    the reference weighted charge ``Q_nu = 1 / R_nu`` (see :class:`Config`), so no Z
    enters here. This is Paper 60 ``eq:secular``.
    """
    K = len(configs)
    Tp = np.zeros((K, K))
    for i in range(K):
        ci = configs[i]
        for j in range(i, K):
            cj = configs[j]
            g = ci.norm * cj.norm * repulsion_terms(ci.terms, cj.terms)
            Tp[i, j] = Tp[j, i] = -g
    return Tp


def build_M(configs: List[Config], Z: float = Z_HE) -> np.ndarray:
    """Assemble the isoenergetic secular matrix ``M = diag(Z * R_nu) + T'``.

    Z enters ONLY the diagonal; the interelectron block ``T'`` (from
    :func:`build_Tprime`) is Z-independent.
    """
    M = build_Tprime(configs)
    for i, ci in enumerate(configs):
        M[i, i] += Z * ci.Rnu
    return M


def build_S(configs: List[Config]) -> np.ndarray:
    """L2 overlap metric between the (L2-normalized) mixed-scale configurations.

    Reintroduced only for the diagnostic generalized eigenproblem (claim iii); the
    primary metric-free method does not use S.
    """
    K = len(configs)
    S = np.zeros((K, K))
    for i in range(K):
        ci = configs[i]
        for j in range(i, K):
            cj = configs[j]
            s = ci.norm * cj.norm * overlap_terms(ci.terms, cj.terms)
            S[i, j] = S[j, i] = s
    return S


def gen_configs(lmax: int,
                nmax_per_l: Union[int, Dict[int, int]]) -> List[Tuple[int, int, int]]:
    """Enumerate ``(l, n_a, n_b)`` configuration tuples up to ``lmax``.

    ``nmax_per_l`` is either an int (same principal-qn ceiling for every l) or a dict
    ``l -> nmax``. For each l, ``n`` ranges over ``l+1 .. nmax`` and pairs
    ``(n_a, n_b)`` are taken with replacement (n_a <= n_b).
    """
    if isinstance(nmax_per_l, int):
        nmax_per_l = {l: nmax_per_l for l in range(lmax + 1)}
    configs: List[Tuple[int, int, int]] = []
    for l in range(lmax + 1):
        nmx = nmax_per_l.get(l, 0)
        ns = list(range(l + 1, nmx + 1))
        for na, nb in combinations_with_replacement(ns, 2):
            configs.append((l, na, nb))
    return configs


def build_configs(config_tuples: List[Tuple[int, int, int]]) -> List[Config]:
    """Reset the module caches and construct :class:`Config` objects from tuples."""
    reset_caches()
    return [Config(l, na, nb) for (l, na, nb) in config_tuples]


def solve(config_tuples: List[Tuple[int, int, int]],
          Z: float = Z_HE) -> Tuple[float, float, int, np.ndarray]:
    """Solve the metric-free (standard) isoenergetic secular equation.

    Returns ``(E, one_norm, K, M)`` where ``E = -p_kappa^2 / 2`` for the largest
    root ``p_kappa`` (deepest binding), ``one_norm = ||M||_1`` (entrywise sum of
    ``|M_ij|``), ``K`` the config count, and ``M`` the secular matrix.
    """
    cfgs = build_configs(config_tuples)
    M = build_M(cfgs, Z)
    p = np.sort(eigh(M, eigvals_only=True))[-1]      # largest root = deepest binding
    E = -p ** 2 / 2
    onenorm = float(np.abs(M).sum())
    return E, onenorm, len(cfgs), M


def solve_with_metric(config_tuples: List[Tuple[int, int, int]],
                      Z: float = Z_HE) -> Tuple[float, float, float, int]:
    """Diagnostic: reintroduce the L2 overlap metric -> generalized eigenproblem
    ``M B = p S B``.

    Returns ``(E_metricfree, E_with_S, cond_S, K)``. The L2 framing (Paper 60 Sec.2)
    is ill-conditioned: ``cond(S)`` grows into the thousands as the basis grows,
    whereas the metric-free standard eigenproblem stays well-behaved.
    """
    cfgs = build_configs(config_tuples)
    M = build_M(cfgs, Z)
    S = build_S(cfgs)
    p_std = np.sort(eigh(M, eigvals_only=True))[-1]
    p_gen = np.sort(eigh(M, S, eigvals_only=True))[-1]
    return -p_std ** 2 / 2, -p_gen ** 2 / 2, float(np.linalg.cond(S)), len(cfgs)
