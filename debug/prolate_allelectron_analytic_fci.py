r"""Route C / C3: the SIGMA-ONLY all-electron LiH FCI from validated analytic pieces.

Wires together, WITHOUT re-deriving any of them (all imported):
  * C1  debug/prolate_atomcentered_core.py   -- tight-core R-accuracy gate (M_xi/M_eta)
  * C2  debug/prolate_mixed_eri.py            -- analytic mixed core-valence sigma ERIs
        (eri_sigma; orbital constructors sto_orbital / valence_prolate_orbital)
  * FCI debug/prolate_allelectron_fci.py:385  fci_energy(h1, eri, M, nelec)

THE ONE NEW PIECE built here is the UNIFIED one-body engine ``one_body_sigma`` over
the SAME ``Orbital`` objects the ERIs use (Route-A: the core's eta-exponential is the
Legendre truncation, so core and valence share one representation
    phi = norm * xi^p * P_eta(eta) * e^{-alpha xi}
and one-body + ERI + FCI all see the identical basis).  Every sigma one-electron
matrix element is closed form in the C1 xi-moments A_k(c) = int_1^inf xi^k e^{-c xi}
and the elementary eta-moments int_{-1}^1 eta^k P(eta):

    S    = 2pi (R/2)^3 [ A_{P+2}<E>  - A_P <eta^2 E> ]
    <1/rA> = 2pi (R/2)^2 [ A_{P+1}<E> - A_P<eta E> ]     (rA=(R/2)(xi+eta))
    <1/rB> = 2pi (R/2)^2 [ A_{P+1}<E> + A_P<eta E> ]     (rB=(R/2)(xi-eta))
    T    = 1/2 (4/R^2)(R/2)^3 2pi [ (int (xi^2-1) f_i' f_j') <E>
                                    + A_P <(1-eta^2) g_i' g_j'> ]
with P = p_i+p_j, c = alpha_i+alpha_j, E = P_eta,i * P_eta,j (the xi^2-eta^2
Jacobian and the m=0 azimuthal 2pi are exact).  The heteronuclear V_ne combines the
two centres: -Z_A<1/rA> - Z_B<1/rB>.  Reduces to prolate_recondition's single-
electron _ov/_kin/_vne bit-for-bit at monomial-eta / single-alpha (validated below).

Run from root:
  python debug/prolate_allelectron_analytic_fci.py validate   # one-body + controls
  python debug/prolate_allelectron_analytic_fci.py h2         # C-A wiring control
  python debug/prolate_allelectron_analytic_fci.py lih        # the R_eq scan (default)
"""
from __future__ import annotations

import os
import sys
import time
from typing import List, Sequence, Tuple

import numpy as np
import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))                 # debug/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # root

from geovac import neumann_vee_general_m as ngm      # noqa: E402
from geovac import prolate_recondition as pr          # noqa: E402
import prolate_mixed_eri as pmx                        # noqa: E402
from prolate_mixed_eri import (                        # noqa: E402
    Orbital, sto_orbital, valence_prolate_orbital, eri_sigma, sto_eta_poly,
    _i_sph, ZC_LI,
)
from prolate_allelectron_fci import fci_energy, build_mo_integrals  # noqa: E402

mp.mp.dps = 60   # match C2 (core Neumann needs the forward Q_l headroom)

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


# ==========================================================================
# H-centered (centre B) 1s STO as a Route-A Orbital
#   e^{-zeta r_B} = e^{-alpha xi} e^{+alpha eta},  alpha = zeta R/2,
#   e^{+alpha eta} = sum_l (2l+1) i_l(alpha) P_l(eta)   (NO (-1)^l vs sto_eta_poly)
# ==========================================================================
def sto_eta_poly_plus(alpha, L: int = 24) -> List[mp.mpf]:
    alpha = mp.mpf(alpha)
    coeffs: List[mp.mpf] = [mp.mpf(0)]
    for l in range(L):
        a_l = (2 * l + 1) * _i_sph(l, alpha)          # e^{+a eta}: no (-1)^l
        coeffs = pmx._pa(coeffs, [a_l * c for c in ngm._leg_coeffs(l)])
    return coeffs


def sto_orbital_B(zeta, R, L: int = 24) -> Orbital:
    """Unit-normalized centre-B (H) 1s STO sqrt(zeta^3/pi) e^{-zeta r_B}."""
    zeta = mp.mpf(zeta)
    R = mp.mpf(R)
    alpha = zeta * R / 2
    norm = mp.sqrt(zeta ** 3 / mp.pi)
    return Orbital(0, sto_eta_poly_plus(alpha, L), alpha, norm, zeta=None, is_core=False)


# ==========================================================================
# The one new one-body piece: unified sigma S, h1 = T + V_ne over Orbital objects
# ==========================================================================
def _pder(poly: Sequence[mp.mpf]) -> List[mp.mpf]:
    if len(poly) <= 1:
        return [mp.mpf(0)]
    return [mp.mpf(k) * poly[k] for k in range(1, len(poly))]


def one_body_sigma(orbs: Sequence[Orbital], R, Z_A, Z_B
                   ) -> Tuple[np.ndarray, np.ndarray]:
    """Float S[M,M], h1[M,M] = T + V_ne over the mixed {core, valence} sigma set."""
    with mp.workdps(60):
        R = mp.mpf(R)
        hR = R / 2
        two_pi = 2 * mp.pi
        pref_S = two_pi * hR ** 3
        pref_V = two_pi * hR ** 2
        pref_T = mp.mpf('0.5') * (4 / R ** 2) * hR ** 3 * two_pi
        xi2m1 = pr._xi2m1(1)
        meta2 = pr._meta2(1)
        ZA, ZB = mp.mpf(Z_A), mp.mpf(Z_B)
        n = len(orbs)
        S = np.empty((n, n), object)
        H = np.empty((n, n), object)
        for i in range(n):
            oi = orbs[i]
            gi = list(oi.eta_poly)
            dgi = _pder(gi)
            for jj in range(i, n):
                oj = orbs[jj]
                gj = list(oj.eta_poly)
                c = oi.alpha + oj.alpha
                P = oi.xi_power + oj.xi_power
                A = ngm._mono_moments(c, P + 6)
                E = pr._pm(gi, gj)
                Ns = oi.norm * oj.norm
                momE = pr._mom_eta(E)
                momE1 = pr._mom_eta(pr._shift(E, 1))
                momE2 = pr._mom_eta(pr._shift(E, 2))
                # overlap
                Sij = Ns * pref_S * (A[P + 2] * momE - A[P] * momE2)
                # V_ne  (heteronuclear: rA=(R/2)(xi+eta), rB=(R/2)(xi-eta))
                term_xi = A[P + 1] * momE
                term_eta = A[P] * momE1
                vA = term_xi - term_eta
                vB = term_xi + term_eta
                Vne = Ns * pref_V * (-ZA * vA - ZB * vB)
                # kinetic (gradient form; xi^2-eta^2 Jacobian cancels)
                mxi = pr._mx_poly(oi.xi_power, oi.alpha)
                mxj = pr._mx_poly(oj.xi_power, oj.alpha)
                Kxi = pr._mom_xi(pr._pm(pr._pm(mxi, mxj), xi2m1), A)
                dgj = _pder(gj)
                Ke = pr._mom_eta(pr._pm(pr._pm(dgi, dgj), meta2))
                grad = Kxi * momE + A[P] * Ke
                Tij = Ns * pref_T * grad
                S[i, jj] = S[jj, i] = Sij
                H[i, jj] = H[jj, i] = Tij + Vne
        Sf = np.array([[float(S[i, j]) for j in range(n)] for i in range(n)])
        Hf = np.array([[float(H[i, j]) for j in range(n)] for i in range(n)])
        return Sf, Hf


# ==========================================================================
# ERI tensor via C2's eri_sigma, 8-fold-symmetry cached
# ==========================================================================
def _canon(p, q, r, s):
    a = (p, q) if p <= q else (q, p)
    b = (r, s) if r <= s else (s, r)
    return (a, b) if a <= b else (b, a)


def build_eri_tensor(orbs, R, verbose=False):
    M = len(orbs)
    eri = np.zeros((M, M, M, M))
    cache = {}
    t0 = time.time()
    for p in range(M):
        for q in range(M):
            for r in range(M):
                for s in range(M):
                    key = _canon(p, q, r, s)
                    v = cache.get(key)
                    if v is None:
                        v = float(eri_sigma(orbs[p], orbs[q], orbs[r], orbs[s], R))
                        cache[key] = v
                    eri[p, q, r, s] = v
    if verbose:
        print(f"    ERI: {len(cache)} unique of {M**4} in {time.time()-t0:.0f}s", flush=True)
    return eri


# ==========================================================================
# Assemble: one-body + ERI + Loewdin -> fci_energy
# ==========================================================================
def assemble_energy(orbs, R, Z_A, Z_B, nelec, Vnn, cond_tol=1e-10, verbose=False):
    S, h1 = one_body_sigma(orbs, R, Z_A, Z_B)
    eri = build_eri_tensor(orbs, R, verbose=verbose)
    w, U = np.linalg.eigh(S)
    wmax = w[-1]
    keep = w > cond_tol * wmax
    X = U[:, keep] / np.sqrt(w[keep])
    Mk = int(keep.sum())
    h1_o = X.T @ h1 @ X
    eri_o = np.einsum('ap,bq,cr,ds,abcd->pqrs', X, X, X, X, eri, optimize=True)
    E_elec, ndet = fci_energy(h1_o, eri_o, Mk, nelec)
    E_tot = E_elec + Vnn
    cond = wmax / max(w[keep].min(), 1e-300)
    if verbose:
        print(f"    M={len(orbs)} kept={Mk} ndet={ndet} cond(S)={cond:.1e} "
              f"E_elec={E_elec:.5f} E_tot={E_tot:.5f}", flush=True)
    return E_tot, E_elec, Mk, ndet, cond


# ==========================================================================
# Orbital-set builders
# ==========================================================================
def lih_orbitals(R, extended=False):
    """sigma-only LiH set: Li 1s core + valence spanning Li-2s / H-1s / bond."""
    orbs = [sto_orbital(ZC_LI, R, is_core=True)]        # Li 1s core, zeta=2.6875
    orbs.append(sto_orbital(mp.mpf('0.65'), R))          # Li 2s outer lobe (diffuse Li)
    orbs.append(sto_orbital_B(mp.mpf('1.0'), R))         # H 1s
    orbs.append(valence_prolate_orbital(0, 0, mp.mpf('1.0')))   # bond
    orbs.append(valence_prolate_orbital(1, 0, mp.mpf('1.0')))   # bond xi
    if extended:
        orbs.append(sto_orbital_B(mp.mpf('0.70'), R))    # diffuse H^- lobe
        orbs.append(valence_prolate_orbital(0, 1, mp.mpf('1.0')))  # bond eta (polarization)
        orbs.append(sto_orbital(mp.mpf('1.3'), R))       # inner Li 2s node partner
    return orbs


def h2_valence_orbitals(R, M=4, alpha=1.0):
    """sigma-only H2 valence set (no core): bond ProductFns + centre STOs."""
    orbs = [valence_prolate_orbital(0, 0, mp.mpf(alpha)),
            valence_prolate_orbital(1, 0, mp.mpf(alpha))]
    if M >= 3:
        orbs.append(valence_prolate_orbital(0, 1, mp.mpf(alpha)))   # eta -> u/g mix
    if M >= 4:
        orbs.append(valence_prolate_orbital(2, 0, mp.mpf(alpha)))
    if M >= 5:
        orbs.append(valence_prolate_orbital(1, 1, mp.mpf(alpha)))
    if M >= 6:
        orbs.append(valence_prolate_orbital(0, 2, mp.mpf(alpha)))
    return orbs[:M]


# ==========================================================================
# VALIDATIONS
# ==========================================================================
def _relerr(a, b):
    a, b = float(a), float(b)
    return abs(a - b) / abs(b) if b != 0 else abs(a)


def validate_one_body():
    """(1) one_body_sigma == prolate_recondition single-electron _ov/_kin/_vne
    bit-for-bit at monomial-eta / single-alpha (homonuclear).  (2) C-B core."""
    print("=" * 72)
    print("ONE-BODY VALIDATION")
    print("=" * 72)
    R = 2.0
    alpha = 1.1
    # monomial valence orbitals (single alpha) as Orbital objects
    specs = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1)]
    orbs = [valence_prolate_orbital(j, l, mp.mpf(alpha)) for (j, l) in specs]
    S, H = one_body_sigma(orbs, R, 1.0, 1.0)   # homonuclear Z_A=Z_B=1
    # independent reference from pr single-electron primitives
    with mp.workdps(60):
        A = ngm._mono_moments(mp.mpf(2.0 * alpha), 40)
        hR = mp.mpf(R) / 2
        pref_S = 2 * mp.pi * hR ** 3
        pref_V = 2 * mp.pi * hR ** 2
        pref_T = mp.mpf('0.5') * (4 / mp.mpf(R) ** 2) * hR ** 3 * 2 * mp.pi
        n = len(orbs)
        worstS = worstH = 0.0
        for i in range(n):
            ji, li = specs[i]
            for jj in range(n):
                jjx, ljx = specs[jj]
                ov = pr._ov(ji + jjx, li + ljx, 0, A)
                Sref = float(pref_S * ov)
                # V_ne homonuclear: -1<1/rA> -1<1/rB> = -2 * <xi part> (eta parts cancel)
                #   <1/rA>+<1/rB> = 2 * 2pi hR^2 * A[P+1]<E>  ... use _vne (= A[P+1]<eta^q>)
                vne1 = pr._vne(ji + jjx, li + ljx, 0, A)     # = int xi^{P+1} eta^Q Jac.../rA-form
                Vref = float(-2.0 * pref_V * vne1)
                grad, _ = pr._kin(ji, li, jjx, ljx, 0, alpha, A)
                Tref = float(pref_T * grad)
                Href = Tref + Vref
                worstS = max(worstS, _relerr(S[i, jj], Sref) if Sref else abs(S[i, jj]))
                worstH = max(worstH, _relerr(H[i, jj], Href) if Href else abs(H[i, jj]))
    print(f"  valence monomial single-alpha vs pr._ov/_kin/_vne:")
    print(f"    max rel diff  S={worstS:.2e}  H={worstH:.2e}  "
          f"{'PASS' if max(worstS, worstH) < 1e-12 else 'FAIL'}")
    ok = max(worstS, worstH) < 1e-12

    # C-B: isolated Li core one-body through Route-A, R-INDEPENDENT
    print("\n  C-B  isolated Li core (Route-A Legendre), R-independence + closed form")
    zc = float(ZC_LI)
    print(f"    exact:  norm=1  T=zc^2/2={zc**2/2:.6f}  <1/rA>=zc={zc:.6f}  "
          f"E(Z_A=zc,Z_B=0)=-zc^2/2={-zc**2/2:.6f}")
    Ts, Vs, Ns_, Es = [], [], [], []
    for R in (2.70, 2.85, 3.015, 3.20, 3.45):
        core = [sto_orbital(ZC_LI, R, is_core=True)]
        Snorm, _ = one_body_sigma(core, R, 0.0, 0.0)
        # T alone: Z_A=Z_B=0 -> h1 = T
        _, Hk = one_body_sigma(core, R, 0.0, 0.0)
        Tval = Hk[0, 0]
        # <1/rA>: build V_ne with Z_A=1,Z_B=0 -> h1 = T - <1/rA>; subtract T
        _, Hva = one_body_sigma(core, R, 1.0, 0.0)
        v1rA = Tval - Hva[0, 0]
        # E at Z_A=zc, Z_B=0
        _, Hzc = one_body_sigma(core, R, zc, 0.0)
        Ezc = Hzc[0, 0] / Snorm[0, 0]
        Ts.append(Tval / Snorm[0, 0]); Vs.append(v1rA / Snorm[0, 0])
        Ns_.append(Snorm[0, 0]); Es.append(Ezc)
    print(f"    {'R':>6} {'norm':>12} {'T':>12} {'<1/rA>':>12} {'E(zc,0)':>12}")
    for k, R in enumerate((2.70, 2.85, 3.015, 3.20, 3.45)):
        print(f"    {R:6.3f} {Ns_[k]:12.9f} {Ts[k]:12.8f} {Vs[k]:12.8f} {Es[k]:12.8f}")
    spread = lambda v: max(v) - min(v)
    print(f"    R-SPREAD  norm={spread(Ns_):.2e}  T={spread(Ts):.2e}  "
          f"<1/rA>={spread(Vs):.2e}  E={spread(Es):.2e}")
    okcore = (abs(np.mean(Ns_) - 1) < 1e-10 and abs(np.mean(Ts) - zc**2/2) < 1e-9
              and abs(np.mean(Vs) - zc) < 1e-9 and abs(np.mean(Es) + zc**2/2) < 1e-9
              and spread(Es) < 1e-9)
    print(f"    C-B: {'PASS (R-flat, closed forms hit)' if okcore else 'CHECK'}")
    return ok and okcore


def control_h2(alpha=1.0):
    """C-A: FCI wiring on H2 (2 valence e, no core, Z=1,1, R=1.40)."""
    print("=" * 72)
    print("C-A  H2 wiring control (2 valence e, NO core, Z=1,1, R=1.40)")
    print("=" * 72)
    R = 1.40
    print(f"  known: E_exact=-1.1745  E_HF~-1.128  sigma-only CI is below (no pi)")
    print(f"  {'M':>3} {'kept':>4} {'ndet':>5} {'cond':>9} {'E_tot':>11} {'D_e%':>8}")
    for M in (2, 3, 4, 5, 6):
        E_tot, E_elec, Mk, ndet, cond = assemble_energy(
            h2_valence_orbitals(R, M=M, alpha=alpha), R, 1.0, 1.0, nelec=2, Vnn=1.0/R)
        de = 100.0 * (-1.0 - E_tot) / 0.174475
        print(f"  {M:3d} {Mk:4d} {ndet:5d} {cond:9.1e} {E_tot:11.5f} {de:8.2f}")
    # independent grid cross-check (same fci_energy solver, grid MOs)
    print("  --- independent grid-FCI (build_mo_integrals) cross-check ---")
    for M in (3, 4):
        h1, eri, eps, soff = build_mo_integrals(R, M, 1.0, 1.0, N_xi_solve=5000,
                                                N_grid=44, xi_max=13.0)
        Eg, nd = fci_energy(h1, eri, M, nelec=2)
        print(f"    grid M={M}  E_tot={Eg + 1.0/R:.5f}  (grid sigma MOs)")


def scan_lih(extended=False, wide=False):
    """LiH sigma-only E_tot(R) + R_eq + collapse/bind verdict."""
    tag = "extended" if extended else "minimal"
    print("=" * 72)
    print(f"LiH sigma-only all-electron FCI  ({tag} basis)  V_NN=3/R")
    print("=" * 72)
    Rs = [2.70, 2.85, 3.015, 3.20, 3.45]
    if wide:
        Rs = [2.40, 2.70, 2.85, 3.015, 3.20, 3.45, 5.0, 8.0]
    print(f"  {'R':>6} {'kept':>4} {'ndet':>5} {'cond':>9} {'E_elec':>11} {'E_tot':>11}")
    rows = []
    for R in Rs:
        t = time.time()
        E_tot, E_elec, Mk, ndet, cond = assemble_energy(
            lih_orbitals(R, extended=extended), R, 3.0, 1.0, nelec=4, Vnn=3.0/R)
        rows.append((R, E_tot, E_elec, Mk, ndet, cond))
        print(f"  {R:6.3f} {Mk:4d} {ndet:5d} {cond:9.1e} {E_elec:11.5f} "
              f"{E_tot:11.5f}  [{time.time()-t:.0f}s]", flush=True)
    return rows


def _fit_req(rows):
    """Parabolic fit of E_tot(R) near the minimum (bond-range points only)."""
    band = [(R, E) for (R, E, *_) in rows if 2.4 <= R <= 3.6]
    if len(band) < 3:
        return None
    Rr = np.array([b[0] for b in band]); Ee = np.array([b[1] for b in band])
    imin = int(np.argmin(Ee))
    lo, hi = max(0, imin - 1), min(len(band), imin + 2)
    if hi - lo < 3:
        lo, hi = 0, min(3, len(band))
    c = np.polyfit(Rr[lo:hi], Ee[lo:hi], 2)
    if c[0] <= 0:
        return None
    return -c[1] / (2 * c[0])


def run_lih_report():
    log_path = os.path.join(DATA, "lih_analytic_sigma_req.log")
    os.makedirs(DATA, exist_ok=True)
    import io
    buf = io.StringIO()

    class Tee:
        def write(self, s):
            sys.__stdout__.write(s); buf.write(s)
        def flush(self):
            sys.__stdout__.flush()
    sys.stdout = Tee()
    try:
        print(f"Route C / C3  LiH sigma-only analytic FCI   (mp.dps={mp.mp.dps})")
        print(f"date 2026-09-21   zc={ZC_LI}   LiH ref R_e=3.015 bohr, E~-8.070 Ha\n")
        ok = validate_one_body()
        print(f"\nONE-BODY + C-B: {'PASS' if ok else 'FAIL -- STOP'}\n")
        if not ok:
            return
        control_h2()
        print()
        rows_min = scan_lih(extended=False, wide=True)
        req_min = _fit_req(rows_min)
        print()
        rows_ext = scan_lih(extended=True, wide=False)
        req_ext = _fit_req(rows_ext)
        # verdict
        print("\n" + "=" * 72)
        print("VERDICT")
        print("=" * 72)
        band = [(R, E) for (R, E, *_) in rows_min if R <= 3.45]
        Emin_R = min(band, key=lambda x: x[1])[0]
        e_all = [(R, E) for (R, E, *_) in rows_min]
        variational = all(E > -8.070 for (_, E, *_) in rows_min)
        for label, req in (("minimal", req_min), ("extended", req_ext)):
            if req is not None:
                drift = 100.0 * (req - 3.015) / 3.015
                print(f"  {label:9s} R_eq = {req:.3f} bohr  drift {drift:+.1f}% vs 3.015")
            else:
                print(f"  {label:9s} R_eq = NO INTERIOR MINIMUM (monotone -> collapse)")
        print(f"  lowest-E R in bond band (minimal): {Emin_R:.3f} bohr")
        print(f"  C-C variational (all E_tot > -8.070): {variational}")
        # dissociation sanity
        for (R, E, *_) in rows_min:
            if R >= 8.0:
                print(f"  C-D dissociation  E_tot(R={R})={E:.4f}")
        if req_min is None or Emin_R <= 2.40:
            print("  READOUT: COLLAPSES INWARD (no bound minimum in the bond range).")
        else:
            print("  READOUT: BINDS (interior minimum located).")
    finally:
        with open(log_path, "w") as f:
            f.write(buf.getvalue())
        sys.stdout = sys.__stdout__
        print(f"\n[log written to {log_path}]")


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else "lih"
    if arg == "validate":
        validate_one_body()
    elif arg == "h2":
        control_h2()
    elif arg == "lih":
        run_lih_report()
    else:
        run_lih_report()
