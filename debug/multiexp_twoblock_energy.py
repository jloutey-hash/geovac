r"""Two-block radial exponent for Paper 12's H2 prolate CI: conditioning, then energy.

THE TWO QUESTIONS, in the order they gate each other.

  PHASE 1 (cheap, gating).  Does an ORTHOGONAL re-basing absorb the two-block
  conditioning the way it absorbs the single-block case?  The monomial two-block
  overlap comes out at cond ~9.6e16 at (3,3) mu<=1 against one-block's 2.046e15
  (debug/multiexp_overlap_poc.py) -- a 47x handicap on a basis that is severely
  ill-conditioned either way, which is why geovac/prolate_recondition.py exists.

  PHASE 2 (only on GO).  At FIXED FUNCTION COUNT, does a two-block exponent beat
  the single-alpha optimum?  alpha_opt drifts upward with basis size (~1.20-1.25
  at (4,4) mu<=1, 1.40 at (5,5) mu<=1), which is the signature of one exponent
  straining to serve functions of different effective range.

THE STRUCTURAL PREDICTION UNDER TEST -- and it is a prediction, not an assumption.
A transform row is a polynomial against ONE exponential: mixing blocks means
mixing e^{-a1 xi} with e^{-a2 xi}, which leaves the orthogonal family entirely.
So the transform must be BLOCK-DIAGONAL in the degree, and block-diagonal
re-basing can only fix the WITHIN-block Hankel conditioning -- it has no row that
can address a near-degeneracy BETWEEN blocks.  If the measured 47x lives across
the blocks, re-basing cannot reach it.

WHY THE TABLE HAS FIVE COLUMNS AND NOT TWO.  A bare single-vs-two-block pair
cannot say WHERE a penalty comes from, and there are two candidate mechanisms,
not one.  Block-diagonalising the transform also RESTRICTS it: block 2 spans the
monomials xi^{j_split+1}..xi^{j_max} only, so its rows cannot be the plain
Laguerre L_n (which need the low-degree tail sitting in block 1's exponent).
That restriction costs conditioning by itself, at a1 == a2, with no second
exponent anywhere.  So:

    A  single alpha, production transform (full Laguerre family)   baseline
    B  single alpha, BLOCK-DIAGONAL transform (a1 == a2)           isolates the
                                                                  restriction
    C  two blocks (a1 != a2), block-diagonal transform            adds the
                                                                  cross-block part
    D  two blocks, monomial (no re-basing)                        the 47x handicap
    E  single alpha, monomial (no re-basing)                      its baseline

cond(B)/cond(A) is then the transform-restriction cost and cond(C)/cond(B) the
cross-block cost, measured separately rather than inferred from one ratio.

TWO SCHEMES for block 2's rows, because the choice is not forced and the cheaper
thing is to measure both:
  * "shift" -- xi^{d0} L_n(2 a2 (xi-1)), n = 0..(j_max-d0): the natural analogue,
    a triangular invertible transform of block 2's span.
  * "trunc" -- L_o(2 a2 (xi-1)) with its columns below d0 deleted: also
    triangular and invertible, but it is a DIFFERENT set of functions (the
    deleted tail is not re-added elsewhere), so it spans block 2 too.

FALSIFIERS (each must pass before any number here is quotable):
  FA  at a1 == a2, the block-diagonal pipeline must reproduce production's
      ENERGY -- identical span, different basis rows, so the generalized
      spectrum is invariant.  This is the end-to-end one; the v5.14.0 record is
      four matrix-level checks passing while a fixture solved the wrong
      Hamiltonian, so only the energy proves the wiring.
  FB  T1 is invertible (triangular, non-zero diagonal) -- printed, not assumed.
  FC  the one-body pencil's lowest eigenvalues agree between transforms at
      a1 == a2 (a congruence check that does not need V_ee).

THREE TRAPS, each already paid for once in this sprint:
  1. `pr.vee_mp` does NOT substitute a default at l_neumann <= 0; it only does
     min(l_neumann, q_max + 2max(s)), so 0 stays 0, while `recondition_energy`
     substitutes 2*l_max + 4*mu_max + 10.  Passing 0 compares different physics
     (it produced a 0.155-relative false failure).  Every call here passes a
     POSITIVE l_neumann explicitly.
  2. Three production sites read the exponent as `basis[0].alpha` -- a single
     scalar off the FIRST function (neumann_vee_general_m.py:482,
     prolate_general_m.py:259 and :370).  A split basis reaching them is
     silently evaluated at block 1's exponent with NO error raised.  So nothing
     here calls one_body_mp, vee_matrix or the general-m V_ee on a split basis:
     the one-body half goes through multiexp_overlap_poc and V_ee through
     multiexp_vee_xtable.
  3. The two-block X-table costs 8.0-11.4x the single-rate one (measured), so
     (4,4) and up is minutes, not seconds.

Run: python debug/multiexp_twoblock_energy.py {fals|cond|energy|sweep}
"""
from __future__ import annotations

import sys
import time
from typing import Callable, Dict, List, Sequence, Tuple

import mpmath as mp
import numpy as np

from geovac import prolate_recondition as pr

sys.path.insert(0, __file__.rsplit("multiexp_twoblock_energy.py", 1)[0])
import multiexp_overlap_poc as poc          # noqa: E402
import multiexp_vee_xtable as xt            # noqa: E402

R = pr.R_DEFAULT


# ---------------------------------------------------------------------------
# the block-diagonal analytic transform
# ---------------------------------------------------------------------------
def block_transform_1d(j_max: int, j_split: int, a1: float, a2: float, mu: int,
                       family: str = "laguerre_legendre",
                       scheme: str = "shift") -> np.ndarray:
    """1D radial transform T1[o, j] = coeff of monomial xi^j in orth. function o.

    BLOCK-DIAGONAL by construction: a row for a degree in block b has non-zero
    columns only inside block b, because a row is a polynomial times ONE
    exponential and block b's monomials are the only ones carrying e^{-a_b xi}.
    That is the structural constraint the whole measurement turns on.

    j_split >= j_max (or < 0) degenerates to ONE block at a1, in which case the
    "shift" scheme reproduces `pr._transforms_per_mu`'s radial family exactly --
    which is what makes column A of the table share this code path.
    """
    w = j_max + 1
    if family == "laguerre_legendre":
        def rc(n: int, a: float) -> List[mp.mpf]:
            return pr.laguerre_coeffs(n, a, w)
    elif family == "gegenbauer":
        def rc(n: int, a: float) -> List[mp.mpf]:
            return pr.assoc_laguerre_coeffs(n, mu, a, w)
    else:
        raise ValueError(f"unknown family {family!r}")

    if j_split < 0 or j_split >= j_max:
        blocks = [(0, j_max, a1)]
    else:
        blocks = [(0, j_split, a1), (j_split + 1, j_max, a2)]

    T1 = np.empty((w, w), object)
    for i in range(w):
        for j in range(w):
            T1[i, j] = mp.mpf(0)
    for (d0, d1, a) in blocks:
        for n, o in enumerate(range(d0, d1 + 1)):
            if scheme == "shift":
                row = rc(n, a)                      # degree n, columns 0..n
                for j in range(len(row) - d0):
                    if row[j] != 0:
                        T1[o, j + d0] = row[j]      # times xi^{d0}
            elif scheme == "trunc":
                row = rc(o, a)                      # degree o, columns 0..o
                for j in range(d0, min(len(row), w)):
                    T1[o, j] = row[j]
            elif scheme == "assoc":
                # xi^{d0} L_n^{(beta)}(2a(xi-1)) with beta = 2 d0 (+ mu for the
                # mu-adapted family).  THE POINT: block 2's reduced polynomial
                # part lives against the measure xi^{2 d0} (xi^2-1)^mu e^{-2a xi},
                # not e^{-2a xi}, so the plain Laguerre family is the wrong
                # orthogonal set there -- which is exactly what "shift" uses and
                # what FC measured as a 5600x conditioning cost at a1 == a2.
                # At d0 = 0 (block 1, or a degenerate single block) beta = mu and
                # this reduces to the production family.
                beta = 2 * d0 + (mu if family == "gegenbauer" else 0)
                row = pr.assoc_laguerre_coeffs(n, beta, a, w)
                for j in range(len(row) - d0):
                    if row[j] != 0:
                        T1[o, j + d0] = row[j]
            else:
                raise ValueError(f"unknown scheme {scheme!r}")
    return T1


def transforms_two_block(j_max: int, l_max: int, mu_max: int, j_split: int,
                         a1: float, a2: float,
                         family: str = "laguerre_legendre",
                         scheme: str = "shift"
                         ) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """Per-mu (T_radial, T_angular), mirroring `pr._transforms_per_mu`'s shapes.

    The two-electron radial transform is the Kronecker-style product of two
    copies of the SAME 1D block-diagonal matrix, which is what keeps the
    exponent consistent: orthogonal function (a_, c_) draws electron-1 monomials
    only from a_'s block, hence only at alpha_of(a_).  The angular side carries
    no exponential and is untouched.
    """
    Tr_list: List[np.ndarray] = []
    Ta_list: List[np.ndarray] = []
    rlist = [(a, c) for a in range(j_max + 1) for c in range(j_max + 1)]
    alist = [(l, m) for l in range(l_max + 1) for m in range(l_max + 1)
             if (l + m) % 2 == 0]
    for mu in range(mu_max + 1):
        T1 = block_transform_1d(j_max, j_split, a1, a2, mu, family, scheme)
        if family == "laguerre_legendre":
            Te = [pr.legendre_coeffs(nn, l_max + 1) for nn in range(l_max + 1)]
        else:
            Te = [pr.gegenbauer_coeffs(nn, mp.mpf(mu) + mp.mpf('0.5'), l_max + 1)
                  for nn in range(l_max + 1)]
        Tem = np.array(Te, object)
        Tr = np.empty((len(rlist), len(rlist)), object)
        for o, (a_, c_) in enumerate(rlist):
            for i, (j, k) in enumerate(rlist):
                Tr[o, i] = T1[a_, j] * T1[c_, k]
        Ta = np.empty((len(alist), len(alist)), object)
        for o, (b_, d_) in enumerate(alist):
            for i, (l, m) in enumerate(alist):
                Ta[o, i] = Tem[b_][l] * Tem[d_][m]
        Tr_list.append(Tr)
        Ta_list.append(Ta)
    return Tr_list, Ta_list


# ---------------------------------------------------------------------------
# conditioning of the normalized overlap -- the decision-gate metric
# ---------------------------------------------------------------------------
def _cond_norm(S: np.ndarray) -> Tuple[float, float, float, bool]:
    """(cond, min eig, max eig, SPD) of the UNIT-NORMALIZED overlap.

    Mirrors `pr._normalized_solve`'s `cond_norm` definition exactly -- wmax over
    the smallest eigenvalue above 1e-14*wmax of the normalized matrix -- so the
    numbers here are comparable to the recorded 4.90e10 at the (5,5)+delta
    headline.
    """
    Sf = np.asarray(S, dtype=float)
    d = np.sqrt(np.diag(Sf))
    Sh = (Sf / d[:, None]) / d[None, :]
    Sh = 0.5 * (Sh + Sh.T)
    w = np.linalg.eigvalsh(Sh)
    wmax = w[-1]
    cond = float(wmax / w[w > 1e-14 * wmax].min())
    return cond, float(w[0]), float(wmax), bool(w[0] > 0)


def _cond_raw(S: np.ndarray) -> float:
    """Plain cond of the un-normalized overlap (the poc's 9.6e16 / 2.046e15)."""
    Sf = np.asarray(S, dtype=float)
    Sf = 0.5 * (Sf + Sf.T)
    w = np.linalg.eigvalsh(Sf)
    return float(w[-1] / w[0]) if w[0] > 0 else float("inf")


def _n_mom(j_max: int, l_max: int, mu_max: int) -> int:
    """Same moment-table sizing as `build_one_body_direct`, deliberately generous.

    The poc measured cond EXACTLY flat in n_mom (9.610e16 at 48/96/192 terms), so
    a short table is not a candidate explanation for anything measured here.
    """
    return 6 * max(j_max, l_max) + 6 * (mu_max + 2) + 20


def _mats(j_max: int, l_max: int, mu_max: int, alpha_of: Callable[[int], float],
          a_ref: float) -> Tuple[List[pr.ProductFn], np.ndarray, int, int, int]:
    idx = pr._product_index(j_max, l_max, mu_max)
    fns = [pr.ProductFn(j, l, k, m, mu, a_ref) for (j, l, k, m, mu) in idx]
    S = poc.overlap_multi(fns, alpha_of, _n_mom(j_max, l_max, mu_max))
    Nr = (j_max + 1) ** 2
    Na = len([1 for l in range(l_max + 1) for m in range(l_max + 1)
              if (l + m) % 2 == 0])
    return fns, S, mu_max + 1, Nr, Na


def cond_row(j_max: int, l_max: int, mu_max: int, j_split: int, a1: float,
             a2: float, family: str = "laguerre_legendre",
             scheme: str = "shift", label: str = "") -> Dict[str, float]:
    """One row of the conditioning table: monomial and re-based cond, same basis."""
    t0 = time.time()
    with mp.workdps(pr.DEFAULT_DPS):
        aof = poc.alpha_of_split(j_split, a1, a2) if 0 <= j_split < j_max \
            else (lambda j: a1)
        fns, S, Nmu, Nr, Na = _mats(j_max, l_max, mu_max, aof, a1)
        S_mono = pr._to_f64(S)
        Tr, Ta = transforms_two_block(j_max, l_max, mu_max, j_split, a1, a2,
                                      family, scheme)
        S_o = pr._to_f64(pr._factored_cob(S, Nmu, Nr, Na, Tr, Ta))
    c_mono_n, e_mono, _mx, spd_m = _cond_norm(S_mono)
    c_mono_r = _cond_raw(S_mono)
    c_o_n, e_o, _mx2, spd_o = _cond_norm(S_o)
    out = dict(N=len(fns), cond_mono_norm=c_mono_n, cond_mono_raw=c_mono_r,
               cond_reb_norm=c_o_n, mineig_reb=e_o, spd_mono=spd_m, spd_reb=spd_o,
               secs=time.time() - t0)
    print(f"  {label:<34s} N={out['N']:<5d} cond(mono,norm)={c_mono_n:.3e}  "
          f"cond(mono,raw)={c_mono_r:.3e}  cond(rebased,norm)={c_o_n:.3e}  "
          f"SPD={spd_m}/{spd_o}  [{out['secs']:.0f}s]", flush=True)
    return out


# ---------------------------------------------------------------------------
# energy
# ---------------------------------------------------------------------------
def two_block_energy(j_max: int, l_max: int, mu_max: int, j_split: int,
                     a1: float, a2: float, family: str = "laguerre_legendre",
                     scheme: str = "shift", l_neumann: int = 0,
                     verbose: bool = True, label: str = "") -> pr.ReconditionResult:
    """Full two-block H2 energy: S + H1 (per-pair dispatch) + V_ee (per-rate-pair).

    Deliberately NOT routed through `recondition_energy`: its `direct` engine
    calls `build_one_body_direct` and `vee_mp`, both single-alpha, and trap 2
    above is that a split basis reaching a `basis[0].alpha` reader is evaluated
    at block 1's exponent in silence.  The Hamiltonian is assembled here from the
    two validated multi-rate halves only.
    """
    if l_neumann <= 0:
        l_neumann = 2 * l_max + 4 * mu_max + 10        # trap 1: never pass 0 on
    t0 = time.time()
    single = not (0 <= j_split < j_max)
    with mp.workdps(pr.DEFAULT_DPS):
        aof = (lambda j: a1) if single else poc.alpha_of_split(j_split, a1, a2)
        idx = pr._product_index(j_max, l_max, mu_max)
        fns = [pr.ProductFn(j, l, k, m, mu, a1) for (j, l, k, m, mu) in idx]
        nm = _n_mom(j_max, l_max, mu_max)
        S = poc.overlap_multi(fns, aof, nm)
        H1 = poc.h1_multi(fns, aof, nm)
        V = xt.vee_multi(fns, a1, R, l_neumann, alpha_of=aof, verbose=False)
        n = len(fns)
        bad = [(i, j) for i in range(n) for j in range(n) if V[i, j] is None]
        if bad:
            raise RuntimeError(f"V_ee gather left {len(bad)} entries unfilled, "
                               f"first {bad[:3]}")
        Sf = 1.0 / mp.mpf(R)
        H = np.empty((n, n), object)
        for i in range(n):
            for j in range(n):
                H[i, j] = H1[i, j] + V[i, j] + Sf * S[i, j]
        Nr = (j_max + 1) ** 2
        Na = len([1 for l in range(l_max + 1) for m in range(l_max + 1)
                  if (l + m) % 2 == 0])
        Tr, Ta = transforms_two_block(j_max, l_max, mu_max, j_split, a1, a2,
                                      family, scheme)
        S_o = pr._to_f64(pr._factored_cob(S, mu_max + 1, Nr, Na, Tr, Ta))
        H_o = pr._to_f64(pr._factored_cob(H, mu_max + 1, Nr, Na, Tr, Ta))
    E, cond, nk, sweep = pr._normalized_solve(S_o, H_o)
    de = 100.0 * (-1.0 - E) / pr.DE_EXACT
    err = (pr.E_EXACT - E) * 1000.0
    res = pr.ReconditionResult(E, de, err, n, nk, cond,
                               E > pr.E_EXACT - 5e-6, family,
                               (j_max, l_max, mu_max), sweep, "twoblock")
    if verbose:
        print(f"  {label:<34s} N={n:<5d} E={E:.7f}  D_e%={de:.3f}  "
              f"err={abs(err):.3f}mHa  cond={cond:.2e}  keep={nk}  "
              f"{'var' if res.variational else 'NON-VARIATIONAL'}  "
              f"[{time.time() - t0:.0f}s]", flush=True)
    return res


# ---------------------------------------------------------------------------
# falsifiers
# ---------------------------------------------------------------------------
def fb_invertible(j_max: int = 5, j_split: int = 2, a1: float = 1.60,
                  a2: float = 1.00) -> bool:
    """FB: the block-diagonal T1 is invertible, and block-diagonal in fact."""
    print("\n=== FB: T1 block-diagonal and invertible ===")
    ok = True
    with mp.workdps(pr.DEFAULT_DPS):
        for scheme in ("shift", "trunc"):
            T1 = block_transform_1d(j_max, j_split, a1, a2, 0,
                                    "laguerre_legendre", scheme)
            Tf = np.array([[float(T1[i, j]) for j in range(j_max + 1)]
                           for i in range(j_max + 1)])
            leak = 0
            for o in range(j_max + 1):
                b_o = 0 if o <= j_split else 1
                for j in range(j_max + 1):
                    b_j = 0 if j <= j_split else 1
                    if b_o != b_j and Tf[o, j] != 0.0:
                        leak += 1
            det = float(np.linalg.det(Tf))
            cond = float(np.linalg.cond(Tf))
            good = leak == 0 and det != 0.0
            ok = ok and good
            print(f"  {scheme:<6s} cross-block non-zeros={leak}  det={det:.4e}  "
                  f"cond(T1)={cond:.3e}  {'ok' if good else 'FAIL'}")
    print(f"  {'PASS' if ok else 'FAIL'}")
    return ok


def _ov_quad(p: int, q: int, mu: int, c: float) -> mp.mpf:
    """One electron's overlap factor by DIRECT quadrature, from the definition.

    int int (xi^2 - eta^2) xi^p eta^q (xi^2-1)^mu (1-eta^2)^mu e^{-c xi}, i.e.
    `pr._ov` with no monomial moment table and no `_mom_eta` rational sum -- the
    two pieces the per-pair dispatch actually routes.
    """
    cm = mp.mpf(c)
    xa = mp.quad(lambda x: x ** (p + 2) * (x * x - 1) ** mu * mp.e ** (-cm * x),
                 [1, mp.inf])
    xb = mp.quad(lambda x: x ** p * (x * x - 1) ** mu * mp.e ** (-cm * x),
                 [1, mp.inf])
    ya = mp.quad(lambda y: y ** q * (1 - y * y) ** mu, [-1, 1])
    yb = mp.quad(lambda y: y ** (q + 2) * (1 - y * y) ** mu, [-1, 1])
    return xa * ya - xb * yb


def fd_two_rate_overlap(j_max: int = 3, l_max: int = 3, mu_max: int = 1,
                        j_split: int = 1, a1: float = 1.60, a2: float = 1.00,
                        ncheck: int = 10) -> bool:
    """FD: MIXED-RATE overlap entries vs direct quadrature (the one open gap).

    Every existing check on the two-block overlap is either DEGENERATE (a1 == a2,
    where the mixed rate a1+a2 never arises) or structural (SPD, flat in n_mom).
    The conditioning result rests on entries whose electron-1 factor decays at
    a1 + a2 -- a rate no single-exponent code path ever builds -- so those need a
    route outside `_mono_moments`/`_mom_eta` before any cond number is quotable.
    Mirrors F3's logic in multiexp_vee_xtable.py for the V_ee side.

    Deliberately picks pairs that are CROSS-BLOCK on electron 1, since a dispatch
    that silently used 2*a1 (the `basis[0].alpha` failure mode, trap 2) would
    agree on every same-block entry and differ only here.
    """
    print(f"\n=== FD: mixed-rate ({a1}+{a2}) overlap entries vs direct "
          f"quadrature ===")
    mixed = mp.mpf(a1) + mp.mpf(a2)
    with mp.workdps(30):
        aof = poc.alpha_of_split(j_split, a1, a2)
        idx = pr._product_index(j_max, l_max, mu_max)
        fns = [pr.ProductFn(j, l, k, m, mu, a1) for (j, l, k, m, mu) in idx]
        S = poc.overlap_multi(fns, aof, _n_mom(j_max, l_max, mu_max))
        h6 = (mp.mpf(R) / 2) ** 6
        picks: List[Tuple[int, int]] = []
        for i in range(len(fns)):
            for jj in range(i, len(fns)):
                bi, bj = fns[i], fns[jj]
                if bi.mu != bj.mu:
                    continue
                if abs((aof(bi.j) + aof(bj.j)) - float(mixed)) > 1e-12:
                    continue
                picks.append((i, jj))
        step = max(1, len(picks) // ncheck)
        picks = picks[::step][:ncheck]
        ok = True
        for (i, jj) in picks:
            bi, bj = fns[i], fns[jj]
            mu = bi.mu
            c1 = aof(bi.j) + aof(bj.j)
            c2 = aof(bi.k) + aof(bj.k)
            ref = (h6 * pr._phi_cc(mu)
                   * _ov_quad(bi.j + bj.j, bi.l + bj.l, mu, c1)
                   * _ov_quad(bi.k + bj.k, bi.m + bj.m, mu, c2))
            got = S[i, jj]
            den = abs(ref) if ref != 0 else mp.mpf(1)
            rel = abs(got - ref) / den
            good = rel < mp.mpf('1e-18')
            ok = ok and good
            print(f"  ({i:>3d},{jj:>3d}) mu={mu} rates=({c1:.2f},{c2:.2f})  "
                  f"table={mp.nstr(got, 10)}  quad={mp.nstr(ref, 10)}  "
                  f"rel={mp.nstr(rel, 3)}  {'ok' if good else 'MISMATCH'}")
        print(f"  {len(picks)} mixed-rate entries checked of "
              f"{len(mixed and picks) and 'many'}; {'PASS' if ok else 'FAIL'}")
    return ok


def fc_one_body_pencil(j_max: int = 3, l_max: int = 3, mu_max: int = 1,
                       alpha: float = 1.40, scheme: str = "shift",
                       j_split: int = 1) -> bool:
    """FC: at a1 == a2 the one-body pencil spectrum is transform-invariant.

    Same span, different rows, so the generalized eigenvalues of (H1, S) must
    agree.  Needs no V_ee, which makes it the cheap congruence check.
    """
    print("\n=== FC: one-body pencil invariant under the block transform "
          f"(a1 == a2 == {alpha}) ===")
    with mp.workdps(pr.DEFAULT_DPS):
        idx = pr._product_index(j_max, l_max, mu_max)
        fns = [pr.ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
        nm = _n_mom(j_max, l_max, mu_max)
        S = poc.overlap_multi(fns, lambda j: alpha, nm)
        H1 = poc.h1_multi(fns, lambda j: alpha, nm)
        Nr = (j_max + 1) ** 2
        Na = len([1 for l in range(l_max + 1) for m in range(l_max + 1)
                  if (l + m) % 2 == 0])
        outs = {}
        for name, (Tr, Ta) in (
                ("production", pr._transforms_per_mu("laguerre_legendre", j_max,
                                                     l_max, mu_max, alpha)),
                (f"block[{scheme}]", transforms_two_block(
                    j_max, l_max, mu_max, j_split, alpha, alpha,
                    "laguerre_legendre", scheme))):
            S_o = pr._to_f64(pr._factored_cob(S, mu_max + 1, Nr, Na, Tr, Ta))
            H_o = pr._to_f64(pr._factored_cob(H1, mu_max + 1, Nr, Na, Tr, Ta))
            E, cond, nk, _sw = pr._normalized_solve(S_o, H_o)
            outs[name] = (E, cond, nk)
    (E1, c1, k1), (E2, c2, k2) = list(outs.values())
    d = abs(E1 - E2)
    ok = d < 1e-8
    for name, (E, c, k) in outs.items():
        print(f"  {name:<16s} E_1body={E:.10f}  cond={c:.2e}  keep={k}")
    print(f"  |dE| = {d:.3e}   {'PASS' if ok else 'FAIL'}")
    return ok


def fa_degenerate_energy(j_max: int = 3, l_max: int = 3, mu_max: int = 1,
                         alpha: float = 1.40, scheme: str = "shift",
                         j_split: int = 1) -> bool:
    """FA (end-to-end): at a1 == a2 the full energy must equal production's.

    Identical span, identical operator; only the basis rows differ.  The v5.14.0
    record is why this is not redundant with FB/FC -- four matrix-level checks
    passed there while a fixture solved the wrong Hamiltonian.
    """
    print("\n=== FA: full energy invariant at a1 == a2 (span identical) ===")
    ref = pr.recondition_energy(j_max, l_max, mu_max, alpha,
                                "laguerre_legendre", engine="direct")
    got = two_block_energy(j_max, l_max, mu_max, j_split, alpha, alpha,
                           scheme=scheme, label=f"block[{scheme}] a1=a2")
    print(f"  production      E={ref.energy:.10f}  D_e%={ref.de_pct:.4f}  "
          f"cond={ref.cond_norm:.2e}  keep={ref.n_kept}")
    d = abs(ref.energy - got.energy)
    ok = d < 1e-7
    print(f"  |dE| = {d:.3e} Ha   {'PASS' if ok else 'FAIL'}")
    return ok


# ---------------------------------------------------------------------------
# phase 1 / phase 2 drivers
# ---------------------------------------------------------------------------
_TRUNCS = ((3, 3, 1), (4, 4, 1))


def phase1(truncs: Sequence[Tuple[int, int, int]] = _TRUNCS,
           a_single: float = 1.40, a1: float = 1.60, a2: float = 1.00,
           scheme: str = "shift") -> None:
    print("=== PHASE 1: does block-diagonal re-basing absorb the two-block "
          "penalty? ===")
    print("cond(norm) is `_normalized_solve`'s metric (headline (5,5)+delta = "
          "4.90e10).\n")
    for (j, l, mu) in truncs:
        print(f"({j},{l}) mu<={mu}:")
        js = max(1, j // 2)
        cond_row(j, l, mu, -1, a_single, a_single, scheme=scheme,
                 label=f"A single a={a_single} prod-transform")
        cond_row(j, l, mu, js, a_single, a_single, scheme=scheme,
                 label=f"B single a={a_single} block-transform")
        cond_row(j, l, mu, js, a1, a2, scheme=scheme,
                 label=f"C two-block {a1}/{a2} j<={js}")
        print()


def phase2(truncs: Sequence[Tuple[int, int, int]] = _TRUNCS,
           scheme: str = "shift") -> None:
    print("=== PHASE 2: two-block vs the single-alpha optimum, FIXED N ===\n")
    for (j, l, mu) in truncs:
        print(f"({j},{l}) mu<={mu}:")
        for a in (1.20, 1.40):
            two_block_energy(j, l, mu, -1, a, a, scheme=scheme,
                             label=f"single a={a}")
        js = max(1, j // 2)
        for (a1, a2) in ((1.60, 1.00), (1.80, 1.10), (1.20, 1.60)):
            two_block_energy(j, l, mu, js, a1, a2, scheme=scheme,
                             label=f"two-block {a1}/{a2} j<={js}")
        print()


def schemes(truncs: Sequence[Tuple[int, int, int]] = ((3, 3, 1), (4, 4, 1)),
            a_single: float = 1.40, a1: float = 1.60, a2: float = 1.00) -> None:
    """Decompose the penalty: restriction (B/A) vs cross-block (C/B), per scheme.

    Run because FC measured the restriction cost FIRST and it was large (2.24e5
    -> 1.26e9 at a1 == a2, scheme "shift"), i.e. most of a naive two-block
    conditioning number may be nothing to do with having two exponents.  If a
    better block-2 family removes it, the cross-block part is what is left.
    """
    print("=== scheme comparison: what does BLOCK-DIAGONALISING cost, and what "
          "does the SECOND EXPONENT cost? ===")
    print("A = production family (no restriction).  B = block transform at "
          "a1 == a2 (restriction only).\nC = block transform, two exponents "
          "(restriction + cross-block).  cond is cond(norm).\n")
    for (j, l, mu) in truncs:
        js = max(1, j // 2)
        print(f"({j},{l}) mu<={mu}  j_split={js}")
        cond_row(j, l, mu, -1, a_single, a_single, scheme="shift",
                 label=f"A single a={a_single} prod-family")
        for sc in ("shift", "trunc", "assoc"):
            cond_row(j, l, mu, js, a_single, a_single, scheme=sc,
                     label=f"B single a={a_single} block[{sc}]")
        for sc in ("shift", "trunc", "assoc"):
            cond_row(j, l, mu, js, a1, a2, scheme=sc,
                     label=f"C two-block {a1}/{a2} block[{sc}]")
        print()


def split_sweep(j_max: int = 4, l_max: int = 4, mu_max: int = 1,
                a1: float = 1.60, a2: float = 1.00,
                scheme: str = "assoc") -> None:
    """cond(norm) across every admissible split point, at fixed exponents."""
    print(f"=== split-point sweep, ({j_max},{l_max}) mu<={mu_max}, "
          f"{a1}/{a2}, block[{scheme}] ===")
    cond_row(j_max, l_max, mu_max, -1, a1, a1, scheme=scheme,
             label=f"single a={a1} prod-family")
    for js in range(0, j_max):
        cond_row(j_max, l_max, mu_max, js, a1, a2, scheme=scheme,
                 label=f"j_split={js}")


if __name__ == "__main__":
    mp.mp.dps = pr.DEFAULT_DPS
    which = sys.argv[1] if len(sys.argv) > 1 else "fals"
    if which == "fals":
        fb_invertible()
        fc_one_body_pencil()
        fa_degenerate_energy()
    elif which == "fd":
        fd_two_rate_overlap()
    elif which == "fals2":
        for sc in ("trunc", "assoc"):
            fc_one_body_pencil(scheme=sc)
            fa_degenerate_energy(scheme=sc)
    elif which == "cond":
        phase1()
    elif which == "schemes":
        schemes()
    elif which == "split":
        split_sweep()
    elif which == "energy":
        phase2(scheme=sys.argv[2] if len(sys.argv) > 2 else "assoc")
    elif which == "one":
        # one J L MU JSPLIT A1 A2 [SCHEME] -- one energy point, so the phase-2
        # grid can be spread across processes (the only cheap parallelism here:
        # mpmath is single-threaded and the X-table dominates).
        _j, _l, _mu, _js = (int(x) for x in sys.argv[2:6])
        _a1, _a2 = float(sys.argv[6]), float(sys.argv[7])
        _sc = sys.argv[8] if len(sys.argv) > 8 else "assoc"
        two_block_energy(_j, _l, _mu, _js, _a1, _a2, scheme=_sc,
                         label=f"({_j},{_l})mu<={_mu} js={_js} {_a1}/{_a2} [{_sc}]")
    else:
        print(f"unknown phase {which!r}; expected fals, fals2, cond, schemes, "
              f"split or energy")
