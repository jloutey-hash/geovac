"""N-electron s-only variational R12-CI engine (step 4).

Basis:  all N-electron determinants over ns Loewdin-orthonormal s orbitals,
        PLUS one correlated function  |G> = F |Phi_0>,  F = sum_{i<j} f(r_ij).

Every matrix element is a sum over permutations P in S_N of a spin-matched radial
contraction, dispatched on leg shape by debug/r12ci_ne_core.py (all shapes validated
against direct 3-electron quadrature).  For N = 3 the product F H F stays at most
3-body, so this is EXACT -- no RI, no CABS.  N >= 4 would need disjoint-pair (4-leg)
shapes, and the evaluator raises rather than silently dropping them.

Term inventory (verified in debug/r12ci_3e_vertex_rules.py):
  <Phi|F H F|Phi> = <Phi|F^2 V|Phi> + (1/2) sum_i INT |grad_i(F Phi)|^2
  and grad_i(F Phi) = F grad_i Phi + Phi grad_i F, giving
    F^2 (grad Phi)^2            -- legs from F^2
    2 F Phi (grad_i F . grad_i Phi) -- one projected leg + legs from F
    Phi^2 (grad_i F)^2          -- same-pair f'^2, and shared-vertex proj x proj
"""
from __future__ import annotations

import itertools
import os
import sys

import numpy as np
from scipy.linalg import eigh

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r12ci_ne_core import (  # noqa: E402
    c_none, c_one, c_shared, c_triangle, coul_moments, l0_kernel, l0_kernel_proj,
    moments,
)

try:
    from scipy.special import eval_genlaguerre
except ImportError:  # pragma: no cover
    raise


# ---------------------------------------------------------------------------
# radial basis (identical convention to debug/ctf12_r12ci_he.py)
# ---------------------------------------------------------------------------
def R_and_dR(n: int, r: np.ndarray, k: float):
    N = 2.0 * k ** 1.5 / n
    x = 2 * k * r
    L1 = eval_genlaguerre(n - 1, 1, x)
    R = N * np.exp(-k * r) * L1
    L2 = eval_genlaguerre(n - 2, 2, x) if n >= 2 else np.zeros_like(r)
    dR = N * np.exp(-k * r) * (-k * L1 - 2 * k * L2)
    return R, dR


def make_grid(k: float, Ng: int, r_max: float | None = None):
    if r_max is None:
        r_max = 42.0 / k
    t = np.linspace(0.0, 1.0, Ng)
    r = r_max * t ** 2
    r[0] = 1e-9
    wr = np.zeros(Ng)
    wr[1:-1] = (r[2:] - r[:-2]) / 2.0
    wr[0] = (r[1] - r[0]) / 2.0
    wr[-1] = (r[-1] - r[-2]) / 2.0
    return r, wr


# ---------------------------------------------------------------------------
# leg / term algebra
# ---------------------------------------------------------------------------
class Leg:
    """A pair factor on electrons (a, b).  kind: 'scalar' or 'proj' (projected on a)."""

    __slots__ = ("a", "b", "tag", "kind")

    def __init__(self, a, b, tag, kind="scalar"):
        self.a, self.b, self.tag, self.kind = a, b, tag, kind

    def key(self):
        return (min(self.a, self.b), max(self.a, self.b))


class Term:
    """coef * [per-electron density modifiers] * [legs]."""

    __slots__ = ("coef", "mods", "legs")

    def __init__(self, coef, mods, legs):
        self.coef = coef          # float
        self.mods = mods          # dict electron -> 'dd' | 'dr' | 'invr'
        self.legs = legs          # list[Leg]


class KernelBank:
    """Builds and caches L=0 / projected kernels and Legendre moments of PRODUCTS."""

    def __init__(self, r, gammas, Lmax=20, nx=300, exact_coul=True):
        """`gammas` is a LIST -- one entry per correlated basis function.

        Several geminals with different ranges let the CI mix them linearly, so the
        wavefunction can build a different effective correlation length for different
        orbital-pair types (1s1s vs 1s2s vs 2s2s) without gamma being tied to a pair --
        which it cannot be, since electrons are indistinguishable.  Matrix elements need
        different gamma on bra and ket; every product stays elementary because
        exp(-ga r) exp(-gb r) = exp(-(ga+gb) r).
        """
        if np.isscalar(gammas):
            gammas = [float(gammas)]
        self.r, self.gammas, self.Lmax, self.nx = r, list(gammas), Lmax, nx
        self.exact_coul = exact_coul
        self._k0, self._kp, self._mom = {}, {}, {}
        self.fn = {"coul": lambda x: 1.0 / x}
        for a, g in enumerate(self.gammas):
            self.fn[f"f{a}"] = (lambda gg: (lambda x: np.exp(-gg * x)))(g)
            self.fn[f"fp{a}"] = (lambda gg: (lambda x: -gg * np.exp(-gg * x)))(g)

    def _prod(self, tags):
        fns = [self.fn[t] for t in tags]

        def p(x):
            out = fns[0](x)
            for fn in fns[1:]:
                out = out * fn(x)
            return out
        return p

    def k0(self, tags):
        key = tuple(sorted(tags))
        if key not in self._k0:
            if key == ("coul",) and self.exact_coul:
                # EXACT: (1/2) INT_-1^1 dx / r_12 = 1 / max(r_i, r_j).
                # Gauss-Legendre on 1/r_12 converges only as ~1/nx because of the
                # |r_i - r_j| singularity at x = 1 (2.7e-3 rel at nx=160), and the
                # triangle path uses the analytic moments -- so quadrature here makes
                # the same kernel disagree with itself.  Caught by GATE 1.
                self._k0[key] = coul_moments(self.r, 0)[0]
            else:
                self._k0[key] = l0_kernel(self._prod(key), self.r, self.nx)
        return self._k0[key]

    def kproj(self, tags):
        key = tuple(sorted(tags))
        if key not in self._kp:
            self._kp[key] = l0_kernel_proj(self._prod(key), self.r, self.nx)
        return self._kp[key]

    def mom(self, tags):
        key = tuple(sorted(tags))
        if key not in self._mom:
            if key == ("coul",):
                self._mom[key] = coul_moments(self.r, self.Lmax)
            else:
                self._mom[key] = moments(self._prod(key), self.r, self.Lmax, self.nx)
        return self._mom[key]


# ---------------------------------------------------------------------------
# operator term generators.  N = number of electrons.
# ---------------------------------------------------------------------------
def pairs(N):
    return list(itertools.combinations(range(N), 2))


def terms_identity(N):
    return [Term(1.0, {}, [])]


def terms_F(N, a):
    """F_a = sum_{i<j} f_{gamma_a}(r_ij)."""
    return [Term(1.0, {}, [Leg(i, j, f"f{a}")]) for (i, j) in pairs(N)]


def terms_FF(N, a, b):
    """F_a F_b (bra geminal index a, ket geminal index b)."""
    out = []
    for (i, j) in pairs(N):
        for (k, l) in pairs(N):
            out.append(Term(1.0, {}, [Leg(i, j, f"f{a}"), Leg(k, l, f"f{b}")]))
    return out


def terms_V(N):
    out = [Term(1.0, {m: "invr"}, []) for m in range(N)]          # -Z/r_m (Z applied later)
    out += [Term(1.0, {}, [Leg(m, n, "coul")]) for (m, n) in pairs(N)]
    return out


def _mul(t1, t2):
    mods = dict(t1.mods)
    for e, m in t2.mods.items():
        if e in mods:
            raise ValueError("two density modifiers on one electron")
        mods[e] = m
    return Term(t1.coef * t2.coef, mods, t1.legs + t2.legs)


def terms_T_grad(N, bra_idx, ket_idx):
    """(1/2) sum_i INT grad_i(F_bra Phi) . grad_i(F_ket Phi).

    `bra_idx` / `ket_idx` are geminal indices, or None for a plain determinant side.
    """
    bra_F = bra_idx is not None
    ket_F = ket_idx is not None
    out = []
    for i in range(N):
        # (a) both gradients on the orbital part
        base = Term(0.5, {i: "dd"}, [])
        pre = []
        if bra_F and ket_F:
            pre = terms_FF(N, bra_idx, ket_idx)
        elif bra_F:
            pre = terms_F(N, bra_idx)
        elif ket_F:
            pre = terms_F(N, ket_idx)
        else:
            pre = [Term(1.0, {}, [])]
        for p in pre:
            out.append(_mul(base, p))
        # (b) cross terms.  With Bra_grad = F^a grad Phi_A + a Phi_A grad F and
        # likewise for the ket, the two cross pieces are
        #   beta : a F^b Phi_A (grad F . grad Phi_B)   -> derivative on the KET radial
        #   gamma: b F^a Phi_B (grad Phi_A . grad F)   -> derivative on the BRA radial
        # They are NOT equal once a permutation makes pa != pb.
        one = [Term(1.0, {}, [])]
        cross_specs = []
        if bra_F:   # beta: grad of the BRA's F, orbital gradient on the KET radial
            cross_specs.append(("rd", bra_idx,
                                terms_F(N, ket_idx) if ket_F else one))
        if ket_F:   # gamma: grad of the KET's F, orbital gradient on the BRA radial
            cross_specs.append(("dr", ket_idx,
                                terms_F(N, bra_idx) if bra_F else one))
        for mod, gidx, mults in cross_specs:
            for j in range(N):
                if j == i:
                    continue
                cross = Term(0.5, {i: mod}, [Leg(i, j, f"fp{gidx}", kind="proj")])
                for p in mults:
                    out.append(_mul(cross, p))
        # (c) both gradients on F (needs F on both sides)
        if bra_F and ket_F:
            for j in range(N):
                if j == i:
                    continue
                for kk in range(N):
                    if kk == i:
                        continue
                    if j == kk:
                        # same pair: rhat_ij . rhat_ij = 1, so this is SCALAR
                        out.append(Term(0.5, {}, [Leg(i, j, f"fp{bra_idx}"),
                                                  Leg(i, j, f"fp{ket_idx}")]))
                    else:
                        out.append(Term(0.5, {}, [Leg(i, j, f"fp{bra_idx}", kind="proj"),
                                                  Leg(i, kk, f"fp{ket_idx}", kind="proj")]))
    return out


def operator_terms(name, N, a=None, b=None):
    """`a`, `b` are the bra / ket geminal indices (None = plain determinant side)."""
    if name == "S":
        return terms_identity(N)
    if name == "S_FG":                       # <Phi| F_b |Phi_0>
        return terms_F(N, b)
    if name == "S_GG":                       # <Phi_0| F_a F_b |Phi_0>
        return terms_FF(N, a, b)
    if name == "H":
        return terms_V(N) + terms_T_grad(N, None, None)
    if name == "H_FG":                       # <Phi| H F_b |Phi_0>
        return [_mul(v, f) for v in terms_V(N) for f in terms_F(N, b)] \
               + terms_T_grad(N, None, b)
    if name == "H_GG":                       # <Phi_0| F_a H F_b |Phi_0>
        return [_mul(v, f) for v in terms_V(N) for f in terms_FF(N, a, b)] \
               + terms_T_grad(N, a, b)
    raise ValueError(name)


# ---------------------------------------------------------------------------
# evaluator
# ---------------------------------------------------------------------------
def eval_term(term, dens_plain, dens_dd, dens_dr, dens_rd, dens_invr, bank, Z):
    dens = list(dens_plain)
    coef = term.coef
    for e, m in term.mods.items():
        if m == "dd":
            dens[e] = dens_dd[e]
        elif m == "dr":
            dens[e] = dens_dr[e]
        elif m == "rd":
            dens[e] = dens_rd[e]
        elif m == "invr":
            dens[e] = dens_invr[e]
            coef *= -Z
    # group legs by pair
    groups = {}
    for lg in term.legs:
        groups.setdefault(lg.key(), []).append(lg)
    shapes = []
    for key, lgs in groups.items():
        tags = []
        kind = "scalar"
        vertex = None
        for lg in lgs:
            tags.append(lg.tag)
            if lg.kind == "proj":
                kind = "proj"
                vertex = lg.a
        shapes.append((key, tuple(sorted(tags)), kind, vertex))

    n = len(shapes)
    if n == 0:
        return coef * c_none(dens)
    if n == 1:
        key, tags, kind, vertex = shapes[0]
        a, b = key
        if kind == "proj":
            a, b = vertex, (key[0] if key[1] == vertex else key[1])
            return coef * c_one(bank.kproj(tags), a, b, dens)
        return coef * c_one(bank.k0(tags), a, b, dens)
    if n == 2:
        (k1, t1, kd1, v1), (k2, t2, kd2, v2) = shapes
        shared = set(k1) & set(k2)
        if len(shared) != 1:
            raise NotImplementedError("disjoint 2-leg shape (needs N>=4)")
        a = shared.pop()
        b = (set(k1) - {a}).pop()
        c = (set(k2) - {a}).pop()

        def oriented(tags, kind, vertex, other):
            """Kernel indexed [r_a, r_other].  kproj(t)[i,j] projects onto i, so a
            leg projected onto the NON-shared vertex needs the transpose:
            kproj.T[a,other] = <fun(r) (r_other - r_a x)/r>, i.e. projection onto
            `other`.  Both orientations occur: grad_i F puts the projection on
            electron i, which need not be the vertex the two legs share."""
            if kind != "proj":
                return bank.k0(tags)
            return bank.kproj(tags) if vertex == a else bank.kproj(tags).T

        K1 = oriented(t1, kd1, v1, b)
        K2 = oriented(t2, kd2, v2, c)
        return coef * c_shared(K1, K2, a, b, c, dens)
    if n == 3:
        for _, _, kd, _ in shapes:
            if kd == "proj":
                raise NotImplementedError("projected leg inside a triangle")
        verts = set()
        for key, _, _, _ in shapes:
            verts |= set(key)
        if len(verts) != 3:
            raise NotImplementedError("3-leg non-triangle shape")
        a, b, c = sorted(verts)
        mA = mB = mC = None
        for key, tags, _, _ in shapes:
            if key == (a, b):
                mA = bank.mom(tags)
            elif key == (a, c):
                mB = bank.mom(tags)
            elif key == (b, c):
                mC = bank.mom(tags)
        if mA is None or mB is None or mC is None:
            raise NotImplementedError("3-leg shape not a closed triangle")
        return coef * c_triangle(mA, mB, mC, a, b, c, dens)
    raise NotImplementedError(f"{n}-leg shape (needs N>=4)")


# ---------------------------------------------------------------------------
# determinants
# ---------------------------------------------------------------------------
def spatial(so):
    return so // 2


def spin(so):
    return so % 2


def make_dets(ns, N, ms2):
    """All determinants over 2*ns spin-orbitals with 2*Ms = ms2."""
    out = []
    for combo in itertools.combinations(range(2 * ns), N):
        s = sum(1 if spin(x) == 0 else -1 for x in combo)
        if s == ms2:
            out.append(combo)
    return out


def perm_sign(p):
    p = list(p)
    sgn, seen = 1, [False] * len(p)
    for i in range(len(p)):
        if seen[i]:
            continue
        j, ln = i, 0
        while not seen[j]:
            seen[j] = True
            j = p[j]
            ln += 1
        if ln % 2 == 0:
            sgn = -sgn
    return sgn


# ---------------------------------------------------------------------------
# the engine
# ---------------------------------------------------------------------------
def build(ns, N, Z, k, gamma, ms2, Ng=220, r_max=None, Lmax=20, nx=300,
          with_geminal=True, ref_det=None, lowdin=True, exact_coul=True, lams=None):
    """`lams`: optional per-function exponents (the multi-lambda / multi-zeta basis).

    Default (None) is the framework's SHARED scale -- every function at the same k,
    which is what Papers 8-9's structural theorem is about.  Passing a list gives each
    radial function its own lambda, so a tight core function and a diffuse valence
    function can coexist.  The angular labels are untouched, so this does not touch the
    prime directive; it is a radial-amplitude choice.
    """
    if lams is None:
        lams = [k] * ns
    lams = [float(x) for x in lams]
    r, wr = make_grid(min(lams), Ng, r_max)
    W = r * r * wr
    Wr = r * wr                                   # for -Z/r
    raw = {p: R_and_dR(p + 1, r, lams[p]) for p in range(ns)}
    S1 = np.array([[np.sum(raw[p][0] * raw[q][0] * W) for q in range(ns)]
                   for p in range(ns)])
    if lowdin:
        # Loewdin: free classically, same span.  Makes the determinant block exactly I.
        ev, U = eigh(S1)
        X = U @ np.diag(ev ** -0.5) @ U.T
    else:
        # raw (non-orthonormal) Sturmians.  melem computes S explicitly and solve()
        # handles the generalized problem, so this is legitimate -- and it is what
        # makes orbital 0 EXACTLY the 1s Sturmian, so |G> = f * R_1(r1) R_1(r2)
        # matches the reference engine's geminal function for the cross-check.
        X = np.eye(ns)
    R = np.array([sum(X[p, q] * raw[q][0] for q in range(ns)) for p in range(ns)])
    dR = np.array([sum(X[p, q] * raw[q][1] for q in range(ns)) for p in range(ns)])

    dets = make_dets(ns, N, ms2)
    if ref_det is None:
        ref_det = dets[0]
    didx = {d: i for i, d in enumerate(dets)}
    nd = len(dets)
    bank = KernelBank(r, gamma, Lmax=Lmax, nx=nx, exact_coul=exact_coul)
    n_gem = len(bank.gammas) if with_geminal else 0
    perms = list(itertools.permutations(range(N)))
    signs = [perm_sign(p) for p in perms]

    cache = {}

    def terms_for(opname, a, b):
        key = (opname, a, b)
        if key not in cache:
            cache[key] = operator_terms(opname, N, a, b)
        return cache[key]

    def melem(A, B, opname, a=None, b=None):
        terms = terms_for(opname, a, b)
        tot = 0.0
        for P, sg in zip(perms, signs):
            ket = [B[P[kk]] for kk in range(N)]
            if any(spin(A[kk]) != spin(ket[kk]) for kk in range(N)):
                continue
            pa = [spatial(A[kk]) for kk in range(N)]
            pb = [spatial(ket[kk]) for kk in range(N)]
            dp = [R[pa[kk]] * R[pb[kk]] * W for kk in range(N)]
            dd = [dR[pa[kk]] * dR[pb[kk]] * W for kk in range(N)]
            dr = [dR[pa[kk]] * R[pb[kk]] * W for kk in range(N)]
            rd = [R[pa[kk]] * dR[pb[kk]] * W for kk in range(N)]
            di = [R[pa[kk]] * R[pb[kk]] * Wr for kk in range(N)]
            sub = 0.0
            for t in terms:
                sub += eval_term(t, dp, dd, dr, rd, di, bank, Z)
            tot += sg * sub
        return tot

    nb = nd + n_gem
    Smat = np.zeros((nb, nb))
    Hmat = np.zeros((nb, nb))
    for i, A in enumerate(dets):
        for j, B in enumerate(dets):
            if j < i:
                continue
            Smat[i, j] = Smat[j, i] = melem(A, B, "S")
            Hmat[i, j] = Hmat[j, i] = melem(A, B, "H")
    for a in range(n_gem):
        ca = nd + a
        for i, A in enumerate(dets):
            Smat[i, ca] = Smat[ca, i] = melem(A, ref_det, "S_FG", None, a)
            Hmat[i, ca] = Hmat[ca, i] = melem(A, ref_det, "H_FG", None, a)
        for b in range(a, n_gem):
            cb = nd + b
            Smat[ca, cb] = Smat[cb, ca] = melem(ref_det, ref_det, "S_GG", a, b)
            Hmat[ca, cb] = Hmat[cb, ca] = melem(ref_det, ref_det, "H_GG", a, b)
    return Smat, Hmat, dets


def solve(Smat, Hmat, nd=None, res_tol=1e-10, thr=1e-10):
    """Solve the generalized problem, projecting the geminal against the determinant
    space first.

    The determinant block is exactly the identity (Loewdin-orthonormal orbitals), so
    the geminal's residual norm outside that span is the Schur complement
    S_GG - sum_I S_IG^2.  As gamma -> 0 the geminal collapses onto |Phi_0>, that
    residual -> 0, and solving the raw generalized problem picks up a spurious
    variational collapse from the near-null direction.  Projecting first (and
    dropping the column when the residual is numerically zero) makes the gamma -> 0
    control measure physics rather than conditioning.  Same congruence transform as
    the Q12 strong-orthogonality projection; here it is numerical hygiene.
    """
    n = Smat.shape[0]
    if nd is not None and nd < n:
        P, G = slice(0, nd), slice(nd, n)
        C = np.linalg.solve(Smat[P, P], Smat[P, G])
        T = np.eye(n)
        T[P, G] = -C
        Smat = T.T @ Smat @ T
        Hmat = T.T @ Hmat @ T
        resid = np.diag(Smat)[nd:] / np.maximum(np.diag(Smat)[nd:].max(), 1e-300)
        keepcols = list(range(nd)) + [nd + i for i, v in enumerate(np.diag(Smat)[nd:])
                                      if v > res_tol]
        Smat = Smat[np.ix_(keepcols, keepcols)]
        Hmat = Hmat[np.ix_(keepcols, keepcols)]
    ev, U = np.linalg.eigh(Smat)
    keep = ev > thr
    Xo = U[:, keep] / np.sqrt(ev[keep])
    return float(np.linalg.eigvalsh(Xo.T @ Hmat @ Xo)[0]), int(keep.sum())
