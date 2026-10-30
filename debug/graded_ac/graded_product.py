"""Graded almost-commutative product on the true Camporesi-Higuchi S^3 window.

PROTOTYPE INSTRUMENT (debug/, 2026-10-29).  Builds the pieces a later
pre-registered test needs to pose the Chamseddine-Connes-Marcolli relation
Y_2(S) = 4 g^2 at finite n_max (debug/lit_scan/ccm_spectral_action_yukawa_memo.md
Sec. 5).  It does NOT form the scored quantity: no function here is ever
called by this prototype with the S^3 outer factor to produce K_A, K_Phi or
c(n_max).  See graded_product_memo.md.

Three product Dirac operators on H_GV (x) H_F are provided:

  * ``dirac_dd``    D = D_GV (x) gamma_F + 1 (x) D_F   (Dabrowski-Dossena
                    eq. (18), the product their Sec. 4.2 prescribes when the
                    first factor is odd).  D^2 = D_GV^2 (x) 1 + 1 (x) D_F^2
                    because {gamma_F, D_F} = 0.  RECOMMENDED.
  * ``dirac_swap``  D = D_GV (x) 1 + gamma_s (x) D_F with gamma_s the
                    left/right (inversion) swap between the +lambda and
                    -lambda Dirac eigenspaces of the true window.  gamma_s^2 = 1
                    and {gamma_s, D_GV} = 0, so D^2 also factorizes, but
                    gamma_s does NOT commute with the operator system (the true
                    window is odd, Paper 32 prop:odd(c)), so one-forms do not
                    split into gauge + Higgs.  Built to name that failure.
  * ``dirac_sign``  D = D_GV (x) 1 + sign(D_GV) (x) D_F  (the production
                    geovac/almost_commutative.py product).  Negative control:
                    cross term 2|D_GV| (x) D_F.

Outer factor: the TRUE window (complete Dirac eigenspaces |lambda| <= n_max+1/2)
of the round-S^3 spinor Dirac, realized on L^2(SU(2)) (x) C^2 with
D = 2 sigma.J_R + 3/2 (port of tests/test_p38_true_ch_window.py, Route A).
Here the window eigenbasis is built EXACTLY from Clebsch-Gordan coupling of
(J_R, spin) -- no eigh -- so D_GV is diagonal with rational entries and the
real structure J_GV comes out as an exact signed permutation.

Exact arithmetic: ``SparseQ`` holds Gaussian-rational sparse matrices
(Fractions).  D_GV, D_F, gamma_F, gamma_s, J_GV, J_F and every product Dirac
are exact; Tr(D^{2k}) of the unfluctuated products is exact.  The compressed
multipliers carry Clebsch-Gordan square roots, so fluctuated quantities on
S^3 are float64 (the toy outer factor below is exact throughout).

Second-order (epsilon^2) machinery:
  * ``eps2_heat_coefficient(D, F, t)``  (1/2) d^2/d eps^2 Tr exp(-t (D+eps F)^2)
    at eps = 0, by the exact second-order Duhamel formula in the eigenbasis
    of D^2 (no fit, no finite difference).
  * ``kinetic_duhamel(D, C, t)``  (t^2/2) int_0^1 Tr(C e^{-s t D^2} C^dag
    e^{-(1-s) t D^2}) ds  -- the Higgs kinetic (kin-kin) functional.
  * ``eps2_trace_power(D, F, n)``  eps^2 coefficient of Tr (D + eps F)^n
    (exact when D, F are SparseQ).

Toy outer factor (control (d)): the flat 2-torus with antiperiodic spin
structure, window |k| <= R (complete Dirac eigenspaces), plane-wave
multipliers -- see ``TorusWindow`` and ``torus_continuum_*``.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from typing import Dict, Iterable, List, Optional, Tuple, Union

import numpy as np
from sympy import Rational
from sympy.physics.wigner import clebsch_gordan

Num = Union[int, Fraction]


# =============================================================================
# Exact Gaussian-rational sparse matrices
# =============================================================================


class QI:
    """Gaussian rational re + i*im with Fraction parts."""

    __slots__ = ("re", "im")

    def __init__(self, re: Num = 0, im: Num = 0) -> None:
        self.re = Fraction(re)
        self.im = Fraction(im)

    def __add__(self, o: "QI") -> "QI":
        return QI(self.re + o.re, self.im + o.im)

    def __sub__(self, o: "QI") -> "QI":
        return QI(self.re - o.re, self.im - o.im)

    def __mul__(self, o: "QI") -> "QI":
        return QI(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)

    def __neg__(self) -> "QI":
        return QI(-self.re, -self.im)

    def conj(self) -> "QI":
        return QI(self.re, -self.im)

    def is_zero(self) -> bool:
        return self.re == 0 and self.im == 0

    def __eq__(self, o: object) -> bool:
        return isinstance(o, QI) and self.re == o.re and self.im == o.im

    def __repr__(self) -> str:
        return f"QI({self.re}, {self.im})"

    def __complex__(self) -> complex:
        return complex(float(self.re), float(self.im))


class SparseQ:
    """Exact sparse square matrix over Q(i): dict (row, col) -> QI."""

    def __init__(self, n: int, data: Optional[Dict[Tuple[int, int], QI]] = None) -> None:
        self.n = n
        self.d: Dict[Tuple[int, int], QI] = {}
        if data:
            for k, v in data.items():
                if not v.is_zero():
                    self.d[k] = v

    # -- construction -------------------------------------------------------
    @staticmethod
    def identity(n: int) -> "SparseQ":
        return SparseQ(n, {(i, i): QI(1) for i in range(n)})

    @staticmethod
    def diag(vals: Iterable[Num]) -> "SparseQ":
        vals = list(vals)
        return SparseQ(len(vals), {(i, i): QI(v) for i, v in enumerate(vals)})

    @staticmethod
    def from_gaussian_integer_array(arr: np.ndarray, tol: float = 1e-12) -> "SparseQ":
        """Exact reconstruction of an array whose entries are (to tol) Gaussian
        integers; raises if any entry is not."""
        n = arr.shape[0]
        out: Dict[Tuple[int, int], QI] = {}
        re_r = np.rint(arr.real)
        im_r = np.rint(arr.imag)
        dev = max(np.abs(arr.real - re_r).max(), np.abs(arr.imag - im_r).max())
        if dev > tol:
            raise ValueError(f"not a Gaussian-integer matrix (deviation {dev:.3e})")
        for i, j in zip(*np.nonzero((re_r != 0) | (im_r != 0))):
            out[(int(i), int(j))] = QI(int(re_r[i, j]), int(im_r[i, j]))
        return SparseQ(n, out)

    # -- algebra ------------------------------------------------------------
    def __add__(self, o: "SparseQ") -> "SparseQ":
        r = dict(self.d)
        for k, v in o.d.items():
            r[k] = r[k] + v if k in r else v
        return SparseQ(self.n, r)

    def __sub__(self, o: "SparseQ") -> "SparseQ":
        return self + o.scale(QI(-1))

    def scale(self, c: QI) -> "SparseQ":
        return SparseQ(self.n, {k: v * c for k, v in self.d.items()})

    def __matmul__(self, o: "SparseQ") -> "SparseQ":
        rows: Dict[int, List[Tuple[int, QI]]] = {}
        for (i, j), v in o.d.items():
            rows.setdefault(i, []).append((j, v))
        r: Dict[Tuple[int, int], QI] = {}
        for (i, k), v in self.d.items():
            for j, w in rows.get(k, ()):
                key = (i, j)
                p = v * w
                r[key] = r[key] + p if key in r else p
        return SparseQ(self.n, r)

    def kron(self, o: "SparseQ") -> "SparseQ":
        r: Dict[Tuple[int, int], QI] = {}
        for (i, j), v in self.d.items():
            for (k, l), w in o.d.items():
                r[(i * o.n + k, j * o.n + l)] = v * w
        return SparseQ(self.n * o.n, r)

    def conj(self) -> "SparseQ":
        return SparseQ(self.n, {k: v.conj() for k, v in self.d.items()})

    def dag(self) -> "SparseQ":
        return SparseQ(self.n, {(j, i): v.conj() for (i, j), v in self.d.items()})

    def trace(self) -> QI:
        t = QI(0)
        for (i, j), v in self.d.items():
            if i == j:
                t = t + v
        return t

    def is_zero(self) -> bool:
        return len(self.d) == 0

    def nnz(self) -> int:
        return len(self.d)

    def max_abs(self) -> float:
        return max((abs(complex(v)) for v in self.d.values()), default=0.0)

    def to_numpy(self) -> np.ndarray:
        a = np.zeros((self.n, self.n), complex)
        for (i, j), v in self.d.items():
            a[i, j] = complex(v)
        return a


def sq_power(A: SparseQ, k: int) -> SparseQ:
    R = SparseQ.identity(A.n)
    for _ in range(k):
        R = R @ A
    return R


# =============================================================================
# True Camporesi-Higuchi window (exact coupled eigenbasis)
# =============================================================================

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)
ISY = np.array([[0, 1], [-1, 0]], complex)
HALF = Fraction(1, 2)


@lru_cache(maxsize=None)
def _cg(j1: Fraction, m1: Fraction, j2: Fraction, m2: Fraction, J: Fraction, M: Fraction) -> float:
    if m1 + m2 != M or not (abs(j1 - j2) <= J <= j1 + j2) or abs(M) > J:
        return 0.0
    R = lambda x: Rational(x.numerator, x.denominator)  # noqa: E731
    return float(clebsch_gordan(R(j1), R(j2), R(J), R(m1), R(m2), R(M)))


def _ms(j: Fraction) -> List[Fraction]:
    return [j - k for k in range(int(2 * j) + 1)]


def _spin_mats(j: Fraction) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    ms = _ms(j)
    d = len(ms)
    jz = np.diag([float(m) for m in ms]).astype(complex)
    jp = np.zeros((d, d), complex)
    for a, m in enumerate(ms):
        if a > 0:
            jp[a - 1, a] = np.sqrt(float(j * (j + 1) - m * (m + 1)))
    return (jp + jp.conj().T) / 2, (jp - jp.conj().T) / (2j), jz


@dataclass
class TrueWindow:
    """True CH window at cutoff n_max, in the exact coupled eigenbasis.

    Window labels: (sign, n, m, M) with sign = +1: j = n/2, J = n/2 + 1/2,
    eigenvalue +(n + 3/2); sign = -1: j = (n+1)/2, J = n/2, eigenvalue
    -(n + 3/2); m in j (left label), M in J (right label coupled to spin).
    """

    n_max: int
    amb_labels: List[Tuple[Fraction, Fraction, Fraction, int]]
    amb_idx: Dict[Tuple[Fraction, Fraction, Fraction, int], int]
    D_amb: np.ndarray
    win_labels: List[Tuple[int, int, Fraction, Fraction]]
    win_idx: Dict[Tuple[int, int, Fraction, Fraction], int]
    V: np.ndarray                 # ambient x window, real (CG entries)
    eig: List[Fraction]           # exact window eigenvalues

    @property
    def dim(self) -> int:
        return len(self.win_labels)

    # -- exact operators ----------------------------------------------------
    def D_exact(self) -> SparseQ:
        return SparseQ.diag(self.eig)

    def D_float(self) -> np.ndarray:
        return np.diag([float(x) for x in self.eig]).astype(complex)

    def sign_exact(self) -> SparseQ:
        return SparseQ.diag([1 if x > 0 else -1 for x in self.eig])

    def J_unitary_float(self) -> np.ndarray:
        """J_GV = U_J K on the window: U_J = V^T U V (V real)."""
        U = np.zeros((len(self.amb_labels),) * 2, complex)
        for (j, m, mp, s) in self.amb_labels:
            for s2 in (0, 1):
                if ISY[s2, s] != 0:
                    U[self.amb_idx[(j, -m, -mp, s2)], self.amb_idx[(j, m, mp, s)]] = \
                        (-1) ** int(m - mp) * ISY[s2, s]
        return self.V.T @ U @ self.V

    def J_exact(self) -> SparseQ:
        return SparseQ.from_gaussian_integer_array(self.J_unitary_float())

    def swap_grading_exact(self, phase_rule: str = "none") -> SparseQ:
        """gamma_s: |+, n, m, M> <-> |-, n, m' = M, M' = m> (left/right swap).

        phase_rule 'none': all +1; 'mM': (-1)^(m - M) on the + -> - leg and its
        transpose on the way back (keeps gamma_s a real symmetric involution).
        """
        d: Dict[Tuple[int, int], QI] = {}
        for (sg, n, m, M), i in self.win_idx.items():
            if sg != +1:
                continue
            k = self.win_idx[(-1, n, M, m)]
            ph = 1
            if phase_rule == "mM":
                ph = (-1) ** int(m - M + HALF)   # m - M is a half-integer here
            d[(k, i)] = QI(ph)
            d[(i, k)] = QI(ph)
        return SparseQ(self.dim, d)

    # -- multipliers -------------------------------------------------------
    def multiplier_ambient(self, J: Fraction, M: Fraction, Mp: Fraction) -> np.ndarray:
        """Multiplication by sqrt(2J+1)-normalized D^J_{M M'} on the ambient
        space (Clebsch-Gordan on (m, m'), identity on spin)."""
        n = len(self.amb_labels)
        A = np.zeros((n, n), complex)
        for (j1, m1, m1p, s), i1 in self.amb_idx.items():
            m2, m2p = m1 + M, m1p + Mp
            for j2 in {jj for (jj, _, _, _) in self.amb_labels}:
                key = (j2, m2, m2p, s)
                if key not in self.amb_idx:
                    continue
                c = _cg(j1, m1, J, M, j2, m2) * _cg(j1, m1p, J, Mp, j2, m2p)
                if c != 0.0:
                    A[self.amb_idx[key], i1] = np.sqrt(float((2 * j1 + 1) / (2 * j2 + 1))) * c
        return A

    def multiplier(self, J: Fraction, M: Fraction, Mp: Fraction) -> np.ndarray:
        """Compression P f P of f = D^J_{M M'} onto the window (an element of
        the operator system O^CH_{n_max})."""
        A = self.multiplier_ambient(J, M, Mp)
        return self.V.T @ A @ self.V

    def operator_system_basis(self) -> List[np.ndarray]:
        out = []
        for N in range(1, 2 * self.n_max + 1):
            J = Fraction(N - 1, 2)
            for M, Mp in itertools.product(_ms(J), _ms(J)):
                out.append(self.multiplier(J, M, Mp))
        return out


@lru_cache(maxsize=None)
def true_window(n_max: int) -> TrueWindow:
    js = [Fraction(k, 2) for k in range(0, n_max + 1)]
    amb: List[Tuple[Fraction, Fraction, Fraction, int]] = []
    for j in js:
        for m, mp, s in itertools.product(_ms(j), _ms(j), (0, 1)):
            amb.append((j, m, mp, s))
    aidx = {lab: i for i, lab in enumerate(amb)}
    D = np.zeros((len(amb), len(amb)), complex)
    for j in js:
        d = int(2 * j) + 1
        jx, jy, jz = _spin_mats(j)
        blk_r = 2 * (np.kron(jx, SX) + np.kron(jy, SY) + np.kron(jz, SZ)) + 1.5 * np.eye(2 * d)
        blk = np.kron(np.eye(d), blk_r)
        ids = [aidx[(j, m, mp, s)] for m in _ms(j) for mp in _ms(j) for s in (0, 1)]
        D[np.ix_(ids, ids)] = blk
    wl: List[Tuple[int, int, Fraction, Fraction]] = []
    eig: List[Fraction] = []
    cols: List[np.ndarray] = []
    spin = {0: HALF, 1: -HALF}
    for sg in (+1, -1):
        for n in range(n_max):
            if sg == +1:
                j, Jt = Fraction(n, 2), Fraction(n, 2) + HALF
            else:
                j, Jt = Fraction(n + 1, 2), Fraction(n, 2)
            for m in _ms(j):
                for M in _ms(Jt):
                    v = np.zeros(len(amb))
                    for mp in _ms(j):
                        for s in (0, 1):
                            c = _cg(j, mp, HALF, spin[s], Jt, M)
                            if c != 0.0:
                                v[aidx[(j, m, mp, s)]] = c
                    wl.append((sg, n, m, M))
                    eig.append(sg * (Fraction(n) + Fraction(3, 2)))
                    cols.append(v)
    V = np.array(cols).T
    widx = {lab: i for i, lab in enumerate(wl)}
    return TrueWindow(n_max, amb, aidx, D, wl, widx, V, eig)


# =============================================================================
# Finite factor (CCM doubled, KO-6), exact
# =============================================================================


@dataclass(frozen=True)
class FiniteExact:
    """Exact version of geovac.almost_commutative.ElectroweakFiniteTriple.

    H_F = matter (nu_L, e_L, nu_R, e_R) (+) antimatter (same), dim 8.
    D_F = diag(M, conj M), M = [[0, Y^dag], [Y, 0]], Y = diag(y_nu, y_e).
    gamma_F = diag(1,1,-1,-1, -1,-1,1,1);  J_F = U_F K, U_F = sigma_x (x) 1_4.
    """

    y_nu: Fraction = Fraction(0)
    y_e: Fraction = Fraction(0)

    def D(self) -> SparseQ:
        d: Dict[Tuple[int, int], QI] = {}
        for off, ys in ((0, False), (4, True)):
            for a, y in enumerate((self.y_nu, self.y_e)):
                q = QI(y)
                qc = q.conj()
                Y = qc if ys else q          # antimatter block = conj(M)
                Ydag = q if ys else qc
                d[(off + 2 + a, off + a)] = Y
                d[(off + a, off + 2 + a)] = Ydag
        return SparseQ(8, d)

    @staticmethod
    def gamma() -> SparseQ:
        return SparseQ.diag([1, 1, -1, -1, -1, -1, 1, 1])

    @staticmethod
    def U_J() -> SparseQ:
        d = {(i, i + 4): QI(1) for i in range(4)}
        d.update({(i + 4, i): QI(1) for i in range(4)})
        return SparseQ(8, d)

    @staticmethod
    def action(lam: QI, q: Tuple[QI, QI, QI, QI]) -> SparseQ:
        """pi_F(lambda, q) on matter only; q = q1 + q2 i s1 + q3 i s2 + q4 i s3."""
        q1, q2, q3, q4 = q
        i = QI(0, 1)
        # i*sigma_1 = [[0,i],[i,0]]; i*sigma_2 = [[0,1],[-1,0]]; i*sigma_3 = [[i,0],[0,-i]]
        m = {
            (0, 0): q1 + q4 * i, (1, 1): q1 - q4 * i,
            (0, 1): q2 * i + q3, (1, 0): q2 * i - q3,
            (2, 2): lam, (3, 3): lam.conj(),
        }
        return SparseQ(8, m)


# =============================================================================
# Product Dirac operators (exact and float)
# =============================================================================


def dirac_dd(Dgv: SparseQ, DF: SparseQ, gF: SparseQ) -> SparseQ:
    """D_GV (x) gamma_F + 1 (x) D_F."""
    return Dgv.kron(gF) + SparseQ.identity(Dgv.n).kron(DF)


def dirac_swap(Dgv: SparseQ, gs: SparseQ, DF: SparseQ) -> SparseQ:
    """D_GV (x) 1 + gamma_s (x) D_F."""
    return Dgv.kron(SparseQ.identity(DF.n)) + gs.kron(DF)


def dirac_sign(Dgv: SparseQ, sgn: SparseQ, DF: SparseQ) -> SparseQ:
    """D_GV (x) 1 + sign(D_GV) (x) D_F (ungraded, production module)."""
    return Dgv.kron(SparseQ.identity(DF.n)) + sgn.kron(DF)


def cross_term(D: SparseQ, Dgv: SparseQ, DF: SparseQ) -> SparseQ:
    """D^2 - D_GV^2 (x) 1 - 1 (x) D_F^2 (exact)."""
    return (D @ D) - (Dgv @ Dgv).kron(SparseQ.identity(DF.n)) \
        - SparseQ.identity(Dgv.n).kron(DF @ DF)


def exact_trace_D2k(D: SparseQ, k: int) -> QI:
    return sq_power(D @ D, k).trace()


def factorized_trace_D2k(Dgv: SparseQ, DF: SparseQ, k: int) -> QI:
    """sum_i C(k,i) Tr(D_GV^{2i}) Tr(D_F^{2(k-i)}) -- what Tr(D^{2k}) must be
    when D^2 = D_GV^2 (x) 1 + 1 (x) D_F^2."""
    tot = QI(0)
    A, B = Dgv @ Dgv, DF @ DF
    for i in range(k + 1):
        tot = tot + QI(math.comb(k, i)) * sq_power(A, i).trace() * sq_power(B, k - i).trace()
    return tot


# -- real structures ---------------------------------------------------------


def antiunitary_relations(U: SparseQ, D: SparseQ, gamma: Optional[SparseQ] = None) -> Dict[str, object]:
    """For J = U K: exact signs of J^2 = eps, J D = eps' D J, J gamma = eps'' gamma J.

    Returns 'eps', "eps'", "eps''" as +1 / -1 / None (None = neither sign holds),
    each decided by exact zero tests."""
    n = U.n
    I = SparseQ.identity(n)
    J2 = U @ U.conj()
    out: Dict[str, object] = {}
    out["eps"] = +1 if (J2 - I).is_zero() else (-1 if (J2 + I).is_zero() else None)
    UD, DU = U @ D.conj(), D @ U
    out["eps'"] = +1 if (UD - DU).is_zero() else (-1 if (UD + DU).is_zero() else None)
    if gamma is not None:
        Ug, gU = U @ gamma.conj(), gamma @ U
        out["eps''"] = +1 if (Ug - gU).is_zero() else (-1 if (Ug + gU).is_zero() else None)
    return out


# =============================================================================
# One-forms (float)
# =============================================================================


def one_form(D: np.ndarray, a: np.ndarray, alpha: np.ndarray, b: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """(a (x) alpha) [D, b (x) beta]."""
    A = np.kron(a, alpha)
    B = np.kron(b, beta)
    return A @ (D @ B - B @ D)


def gauge_one_form(D: np.ndarray, a: np.ndarray, b: np.ndarray, T: np.ndarray) -> np.ndarray:
    """(a (x) T) [D, b (x) 1]: for the DD product = a[D_GV, b] (x) T gamma_F."""
    return one_form(D, a, T, b, np.eye(T.shape[0]))


def higgs_one_form(D: np.ndarray, a: np.ndarray, b: np.ndarray, alpha: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """(a (x) alpha)[D, b (x) beta] - (a (x) alpha beta)[D, b (x) 1].

    For the DD product with beta even this is (a b) (x) alpha [D_F, beta]: the
    SAME GV multiplier pair (a, b) that builds the gauge one-form."""
    I = np.eye(alpha.shape[0])
    return one_form(D, a, alpha, b, beta) - one_form(D, a, alpha @ beta, b, I)


def selfadjoint_fluctuation(omega: np.ndarray, U: Optional[np.ndarray] = None,
                            eps_prime: int = +1) -> np.ndarray:
    """omega + omega^dag, plus eps' J (.) J^{-1} if U (J = U K) is given."""
    F = omega + omega.conj().T
    if U is not None:
        F = F + eps_prime * (U @ np.conj(F) @ U.conj().T)
    return F


# =============================================================================
# Second-order spectral-action machinery
# =============================================================================


def _phi(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """int_0^1 exp(-s x - (1-s) y) ds, stable on the diagonal."""
    dlt = x - y
    small = np.abs(dlt) < 1e-12
    safe = np.where(small, 1.0, dlt)
    val = np.exp(-y) * (-np.expm1(-safe)) / safe
    return np.where(small, np.exp(-0.5 * (x + y)), val)


def eps2_heat_coefficient(D: np.ndarray, F: np.ndarray, t: float) -> float:
    """(1/2) d^2/d eps^2 Tr exp(-t (D + eps F)^2) at eps = 0 (exact Duhamel).

    (D+eps F)^2 = H + eps X1 + eps^2 X2, H = D^2, X1 = {D,F}, X2 = F^2;
    coefficient = -t Tr(X2 e^{-tH}) + (t^2/2) int_0^1 Tr(X1 e^{-stH} X1 e^{-(1-s)tH}) ds.
    """
    H = D @ D
    X1 = D @ F + F @ D
    X2 = F @ F
    mu, Q = np.linalg.eigh(H)
    E = np.exp(-t * mu)
    X2q = Q.conj().T @ X2 @ Q
    X1q = Q.conj().T @ X1 @ Q
    first = -t * float(np.real(np.sum(np.diag(X2q) * E)))
    W = _phi(t * mu[:, None], t * mu[None, :])
    second = 0.5 * t * t * float(np.real(np.sum(np.abs(X1q) ** 2 * W)))
    return first + second


def kinetic_duhamel(D: np.ndarray, C: np.ndarray, t: float) -> float:
    """(t^2/2) int_0^1 Tr(C e^{-stD^2} C^dag e^{-(1-s)tD^2}) ds >= 0."""
    mu, Q = np.linalg.eigh(D @ D)
    Cq = Q.conj().T @ C @ Q
    W = _phi(t * mu[:, None], t * mu[None, :])
    return 0.5 * t * t * float(np.real(np.sum(np.abs(Cq) ** 2 * W)))


def eps2_trace_power(D, F, n: int):
    """eps^2 coefficient of Tr (D + eps F)^n (works for numpy arrays or SparseQ)."""
    if isinstance(D, SparseQ):
        M0, M1, M2 = SparseQ.identity(D.n), SparseQ(D.n), SparseQ(D.n)
        for _ in range(n):
            M0, M1, M2 = M0 @ D, (M1 @ D) + (M0 @ F), (M2 @ D) + (M1 @ F)
        return M2.trace()
    M0 = np.eye(D.shape[0], dtype=complex)
    M1 = np.zeros_like(M0)
    M2 = np.zeros_like(M0)
    for _ in range(n):
        M0, M1, M2 = M0 @ D, M1 @ D + M0 @ F, M2 @ D + M1 @ F
    return np.trace(M2)


# =============================================================================
# Toy outer factor: flat 2-torus, antiperiodic spin structure
# =============================================================================


@dataclass
class TorusWindow:
    """Plane-wave window k in (Z+1/2)^2, |k|^2 <= R2, spinor C^2.

    Coordinates are stored doubled (K = 2k, odd integers) so everything is
    exact.  D = sigma.k on each k, D^2 = |k|^2.  The window is a union of
    complete Dirac eigenspaces (both signs +-|k| share a k)."""

    R2: Fraction
    ks: List[Tuple[int, int]]          # doubled coordinates
    idx: Dict[Tuple[int, int, int], int]

    @property
    def dim(self) -> int:
        return 2 * len(self.ks)

    def D_exact(self) -> SparseQ:
        d: Dict[Tuple[int, int], QI] = {}
        for (K1, K2) in self.ks:
            k1, k2 = Fraction(K1, 2), Fraction(K2, 2)
            i0, i1 = self.idx[(K1, K2, 0)], self.idx[(K1, K2, 1)]
            # sigma.k = [[0, k1 - i k2], [k1 + i k2, 0]]
            d[(i0, i1)] = QI(k1, -k2)
            d[(i1, i0)] = QI(k1, k2)
        return SparseQ(self.dim, d)

    def shift_exact(self, P1: int, P2: int) -> SparseQ:
        """Compression of e^{i p.x}, p = (P1, P2) integer momentum."""
        d: Dict[Tuple[int, int], QI] = {}
        for (K1, K2) in self.ks:
            tgt = (K1 + 2 * P1, K2 + 2 * P2)
            for s in (0, 1):
                if (tgt[0], tgt[1], s) in self.idx:
                    d[(self.idx[(tgt[0], tgt[1], s)], self.idx[(K1, K2, s)])] = QI(1)
        return SparseQ(self.dim, d)


def torus_window(R2: Num) -> TorusWindow:
    R2 = Fraction(R2)
    lim = int(math.isqrt(int(4 * R2))) + 2
    ks = [(K1, K2) for K1 in range(-lim, lim + 1, 1) for K2 in range(-lim, lim + 1, 1)
          if K1 % 2 and K2 % 2 and Fraction(K1 * K1 + K2 * K2, 4) <= R2]
    idx: Dict[Tuple[int, int, int], int] = {}
    for a, (K1, K2) in enumerate(ks):
        idx[(K1, K2, 0)] = 2 * a
        idx[(K1, K2, 1)] = 2 * a + 1
    return TorusWindow(R2, ks, idx)


def torus_continuum_gauge(p: Tuple[int, int], q: Tuple[int, int], t: float) -> float:
    """Continuum (infinite-torus-lattice, Poisson-exact up to e^{-pi^2/t})
    eps^2 heat coefficient of D + eps A on the flat torus [0,2 pi]^2, with
    A = a[D,b] + h.c., a = e^{i q.x}, b = e^{i p.x}  =>  A = (sigma.p) e^{iP.x} + h.c.,
    P = p + q.  Independent route: symbolic Gaussian k-integration (sympy),
    then the Duhamel s-integral by mpmath quadrature.  No matrices."""
    return _torus_continuum("gauge", p, q, t)


def torus_continuum_higgs_kin(p: Tuple[int, int], q: Tuple[int, int], t: float) -> float:
    """Continuum kinetic functional (t^2/2) int Tr(C e C^dag e) for C = [D, e^{iP.x}]."""
    return _torus_continuum("higgs", p, q, t)


@lru_cache(maxsize=None)
def _torus_continuum_symbolic(kind: str, p: Tuple[int, int], q: Tuple[int, int]):
    import sympy as sp
    u1, u2, s, t = sp.symbols("u1 u2 s t", real=True)
    I = sp.I
    sx = sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Matrix([[0, -I], [I, 0]])
    one = sp.eye(2)
    sdot = lambda v: v[0] * sx + v[1] * sy  # noqa: E731
    P = (p[0] + q[0], p[1] + q[1])
    terms = []          # (kind, polynomial-in-u, delta) contributions to the s-integrand
    first = sp.Integer(0)
    if kind == "gauge":
        ap = sdot(p)                      # Fourier component at +P and at -P (both sigma.p)
        comps = {P: ap, (-P[0], -P[1]): ap}
        # Tr(A^2 e^{-tD^2}) = (pi/t) sum_delta tr(a_{-delta} a_delta)
        first = -t * (sp.pi / t) * sum((comps[(-d[0], -d[1])] * comps[d]).trace() for d in comps)
        for d, ad in comps.items():
            amd = comps[(-d[0], -d[1])]
            # k = u - s*delta
            k = (u1 - s * d[0], u2 - s * d[1])
            kd = (k[0] + d[0], k[1] + d[1])
            b_fwd = sdot(kd) * ad + ad * sdot(k)          # <k+d|B|k>
            b_back = sdot(k) * amd + amd * sdot(kd)        # <k|B|k+d>
            terms.append((sp.expand((b_back * b_fwd).trace()), d))
    else:
        cP = sdot(P)                      # <k+P|[D, e^{iPx}]|k> = sigma.P
        terms.append((sp.expand((cP * cP.H).trace()), P))
    # Gaussian moments int u^n e^{-t u^2} du
    def mom(nn: int):
        if nn % 2:
            return sp.Integer(0)
        return sp.gamma(sp.Rational(nn + 1, 2)) / t ** sp.Rational(nn + 1, 2)
    integrand_s = sp.Integer(0)
    for poly, d in terms:
        pl = sp.Poly(poly, u1, u2)
        val = sp.Integer(0)
        for (e1, e2), c in pl.terms():
            val += c * mom(e1) * mom(e2)
        dd = d[0] ** 2 + d[1] ** 2
        integrand_s += sp.simplify(val) * sp.exp(-t * s * (1 - s) * dd)
    second = t ** 2 / 2 * integrand_s
    # for "higgs" there is a single (P) term; the trace over k of C e C^dag e
    return sp.lambdify((s, t), second, "mpmath"), sp.lambdify((t,), first, "mpmath")


def _torus_continuum(kind: str, p: Tuple[int, int], q: Tuple[int, int], t: float) -> float:
    import mpmath as mp
    mp.mp.dps = 30
    f2, f1 = _torus_continuum_symbolic(kind, tuple(p), tuple(q))
    val = mp.quad(lambda s: f2(s, t), [0, 1])
    if kind == "gauge":
        val += f1(t)
    return float(mp.re(val))
