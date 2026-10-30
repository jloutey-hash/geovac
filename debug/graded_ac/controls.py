"""Controls for the graded almost-commutative prototype (no scored quantity).

Run:  python debug/graded_ac/controls.py   -> debug/graded_ac/controls.log

HARD RULE honoured: nothing here evaluates eps2_heat_coefficient /
kinetic_duhamel (the K_A, K_Phi machinery) with the S^3 outer factor.  On
S^3 only exact algebra (cross terms, traces, J signs) and one-form NORMS are
measured.  The heat-coefficient machinery is exercised on the flat-torus toy
only (control (d)).
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import List

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from graded_product import (  # noqa: E402
    QI, SparseQ, FiniteExact, true_window, torus_window,
    dirac_dd, dirac_swap, dirac_sign, cross_term, exact_trace_D2k,
    factorized_trace_D2k, antiunitary_relations, one_form, gauge_one_form,
    higgs_one_form, selfadjoint_fluctuation, eps2_heat_coefficient,
    kinetic_duhamel, eps2_trace_power, torus_continuum_gauge,
    torus_continuum_higgs_kin,
)

LOG: List[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def qstr(q: QI) -> str:
    return f"{q.re}" + (f" + {q.im} i" if q.im != 0 else "")


# Toy Yukawas (rational, NOT a physical assignment; norms / signs only on S^3)
Y_NU, Y_E = Fraction(1, 2), Fraction(1)
H = Fraction(1, 2)


# =============================================================================
# (a) cross term
# =============================================================================


def control_a() -> None:
    log("=" * 78)
    log("(a) CROSS TERM  D^2 - D_GV^2 (x) 1 - 1 (x) D_F^2   (exact, Q(i) arithmetic)")
    log(f"    toy D_F: y_nu = {Y_NU}, y_e = {Y_E}; also a complex Yukawa y_e = 1/3 + i/4")
    for yset in ((Y_NU, Y_E), (Fraction(0), None)):
        pass
    F = FiniteExact(Y_NU, Y_E)
    DF, gF = F.D(), F.gamma()
    # finite sanity: {gamma_F, D_F} = 0 exactly
    log(f"    finite: {{gamma_F, D_F}} == 0 exactly: {((gF @ DF) + (DF @ gF)).is_zero()}")
    for n in (1, 2, 3):
        W = true_window(n)
        Dg = W.D_exact()
        gs = W.swap_grading_exact("none")
        sg = W.sign_exact()
        Ddd, Dsw, Dsg = dirac_dd(Dg, DF, gF), dirac_swap(Dg, gs, DF), dirac_sign(Dg, sg, DF)
        cdd, csw, csg = cross_term(Ddd, Dg, DF), cross_term(Dsw, Dg, DF), cross_term(Dsg, Dg, DF)
        absD = Dg @ sg
        pred = absD.kron(DF).scale(QI(2))
        fro = float(np.sqrt(sum(abs(complex(v)) ** 2 for v in csg.d.values())))
        log(f"  n_max={n} (dim H = {W.dim} x 8 = {8 * W.dim}):")
        log(f"     graded DD  : cross term exactly zero = {cdd.is_zero()}  (nnz {cdd.nnz()})")
        log(f"     graded swap: cross term exactly zero = {csw.is_zero()}  (nnz {csw.nnz()})")
        log(f"     ungraded   : cross term exactly zero = {csg.is_zero()}  nnz {csg.nnz()}, "
            f"||.||_F = {fro:.6f}, == 2|D_GV| (x) D_F exactly: {(csg - pred).is_zero()}")
        for k in (1, 2, 3):
            tdd = exact_trace_D2k(Ddd, k)
            tsg = exact_trace_D2k(Dsg, k)
            tf = factorized_trace_D2k(Dg, DF, k)
            log(f"       Tr D^{2*k}: DD = {qstr(tdd)}  factorized = {qstr(tf)}  equal: {tdd == tf}"
                f" | ungraded = {qstr(tsg)}  equal: {tsg == tf}")
    # complex Yukawa, n_max = 2
    Fc = FiniteExact(Y_NU, Fraction(1, 3))
    DFc = Fc.D()
    # inject an imaginary part in y_e by hand (QI), keeping Hermiticity
    d = dict(DFc.d)
    for off, conjugate in ((0, False), (4, True)):
        y = QI(Fraction(1, 3), Fraction(1, 4))
        yc = y.conj()
        d[(off + 3, off + 1)] = yc if conjugate else y
        d[(off + 1, off + 3)] = y if conjugate else yc
    DFc = SparseQ(8, d)
    herm = (DFc - DFc.dag()).is_zero()
    W = true_window(2)
    Dg = W.D_exact()
    c = cross_term(dirac_dd(Dg, DFc, gF), Dg, DFc)
    log(f"  complex y_e = 1/3 + i/4 at n_max=2: D_F Hermitian {herm}; DD cross term zero = {c.is_zero()}")


# =============================================================================
# (b) J, gamma, KO signs
# =============================================================================


def control_b() -> None:
    log("=" * 78)
    log("(b) REAL STRUCTURE / GRADING SIGNS  (exact zero tests; J = U K)")
    F = FiniteExact(Y_NU, Y_E)
    DF, gF, UF = F.D(), F.gamma(), F.U_J()
    rF = antiunitary_relations(UF, DF, gF)
    log(f"  finite J_F = sigma_x (x) 1 K      : {rF}   (KO-6: (+,+,-))")
    UFg = UF @ gF
    rFg = antiunitary_relations(UFg, DF, gF)
    log(f"  finite J_F gamma_F               : {rFg}   (DD column 6_-: (-,-,-))")
    for n in (1, 2, 3):
        W = true_window(n)
        Dg, UJ = W.D_exact(), W.J_exact()
        gs = W.swap_grading_exact("none")
        gs2 = W.swap_grading_exact("mM")
        log(f"  n_max={n}:")
        log(f"     outer J_GV (exact signed permutation, nnz {UJ.nnz()}): "
            f"{antiunitary_relations(UJ, Dg)}  (KO-3 pair (-,+))")
        I = SparseQ.identity(W.dim)
        log(f"     gamma_s: gamma_s^2 = 1 {((gs @ gs) - I).is_zero()}, {{gamma_s, D_GV}} = 0 "
            f"{((gs @ Dg) + (Dg @ gs)).is_zero()}, J-relations {antiunitary_relations(UJ, Dg, gs)}; "
            f"phase variant (-1)^(m-M+1/2): {antiunitary_relations(UJ, Dg, gs2)}")
        Ddd = dirac_dd(Dg, DF, gF)
        Dsw = dirac_swap(Dg, gs, DF)
        Dsg = dirac_sign(Dg, W.sign_exact(), DF)
        U_plain = UJ.kron(UF)
        U_dd = UJ.kron(UFg)
        g_sw = gs.kron(gF)
        log(f"     DD   with J_GV(x)J_F        : {antiunitary_relations(U_plain, Ddd)}"
            "  (None = no consistent eps'; DD Table 5 row 3_+, col 6_+)")
        log(f"     DD   with J_GV(x)J_F gamma_F: {antiunitary_relations(U_dd, Ddd)}"
            "  (KO-1 = (+,-); 3 + 6 = 1 mod 8)")
        Ifull = SparseQ.identity(8 * W.dim)
        log(f"     DD   grading 1(x)gamma_F anticommutes with D: "
            f"{((SparseQ.identity(W.dim).kron(gF) @ Ddd) + (Ddd @ SparseQ.identity(W.dim).kron(gF))).is_zero()}"
            "  (expected False: product of odd x even is odd)")
        log(f"     swap with J_GV(x)J_F, grading gamma_s(x)gamma_F: "
            f"{antiunitary_relations(U_plain, Dsw, g_sw)}; gamma^2=1 {((g_sw @ g_sw) - Ifull).is_zero()}, "
            f"{{gamma, D}}=0 {((g_sw @ Dsw) + (Dsw @ g_sw)).is_zero()}")
        log(f"     sign with J_GV(x)J_F        : {antiunitary_relations(U_plain, Dsg)}"
            "  (Paper 32: (-,+))")
    # gamma_s against the operator system: the obstruction
    log("  gamma_s vs operator system O^CH (the failure of the swap product):")
    for n in (1, 2, 3):
        W = true_window(n)
        gs = W.swap_grading_exact("none").to_numpy()
        B = W.operator_system_basis()
        rel = max(np.linalg.norm(gs @ b - b @ gs, 2) / max(np.linalg.norm(b, 2), 1e-300) for b in B)
        nz = sum(np.linalg.norm(gs @ b - b @ gs) > 1e-10 for b in B)
        log(f"     n_max={n}: dim O = {len(B)} basis elements; [gamma_s, b] != 0 for {nz}; "
            f"max ||[gamma_s,b]|| / ||b|| = {rel:.4f}")
        if n <= 2:
            d = W.dim
            Id = np.eye(d)
            rows = [np.kron(b, Id) - np.kron(Id, b.T) for b in B]
            Dc = W.D_float()
            sv = np.linalg.svd(np.vstack(rows), compute_uv=False)
            svg = np.linalg.svd(np.vstack(rows + [np.kron(Dc, Id) + np.kron(Id, Dc.T)]), compute_uv=False)
            log(f"        commutant dim = {int((sv < 1e-9 * sv.max()).sum())}, "
                f"grading space dim = {int((svg < 1e-9 * svg.max()).sum())}  (Paper 32 prop:odd(c): 1, 0)")


# =============================================================================
# (c) gauge and Higgs one-form pieces on S^3
# =============================================================================


def _finite_float():
    F = FiniteExact(Y_NU, Y_E)
    DF = F.D().to_numpy()
    gF = F.gamma().to_numpy()
    T = F.action(QI(1), (QI(1), QI(0), QI(0), QI(0))).to_numpy()      # P_matter: commutes with D_F
    beta = F.action(QI(0), (QI(1), QI(0), QI(0), QI(0))).to_numpy()   # P_L^matter: even, [D_F,beta] != 0
    return F, DF, gF, T, beta


def control_c() -> None:
    log("=" * 78)
    log("(c) ONE-FORM PIECES from ONE multiplier pair (a, b) on the true window (norms only)")
    log("    a = P D^{1/2}_{-1/2,+1/2} P, b = P D^{1/2}_{+1/2,+1/2} P (degree-1 functions);")
    log("    gauge: (a(x)T)[D, b(x)1], T = pi_F(1,1) = P_matter;  Higgs: (a(x)1)[D,b(x)beta]-(a(x)beta)[D,b(x)1],")
    log("    beta = pi_F(0,1) = P_L^matter.  Frobenius norms.")
    F, DF, gF, T, beta = _finite_float()
    I8 = np.eye(8)
    for n in (2, 3):
        W = true_window(n)
        Dg = W.D_float()
        a = W.multiplier(H, -H, H)
        b = W.multiplier(H, H, H)
        Ddd = np.kron(Dg, gF) + np.kron(np.eye(W.dim), DF)
        g = gauge_one_form(Ddd, a, b, T)
        h = higgs_one_form(Ddd, a, b, I8, beta)
        g_pred = np.kron(a @ (Dg @ b - b @ Dg), T @ gF)
        h_pred = np.kron(a @ b, DF @ beta - beta @ DF)
        DF0 = FiniteExact().D().to_numpy()
        D0 = np.kron(Dg, gF) + np.kron(np.eye(W.dim), DF0)
        h0 = higgs_one_form(D0, a, b, I8, beta)
        g_triv = gauge_one_form(Ddd, np.eye(W.dim), np.eye(W.dim), T)
        log(f"  n_max={n}  [DD product]")
        log(f"     ||gauge piece|| = {np.linalg.norm(g):.6f};  = a[D_GV,b](x)T gamma_F to "
            f"{np.abs(g - g_pred).max():.1e}")
        log(f"     ||Higgs piece|| = {np.linalg.norm(h):.6f};  = (ab)(x)[D_F,beta] to {np.abs(h - h_pred).max():.1e}")
        log(f"     Y = 0 control: ||Higgs piece|| = {np.linalg.norm(h0):.1e};  b = 1 control: "
            f"||gauge piece|| = {np.linalg.norm(g_triv):.1e}")
        # self-adjoint J-symmetrized fluctuation with the DD real structure (eps' = -1)
        UJ = W.J_unitary_float()
        Udd = np.kron(UJ, F.U_J().to_numpy() @ gF)
        Fg = selfadjoint_fluctuation(g, Udd, eps_prime=-1)
        Fh = selfadjoint_fluctuation(h, Udd, eps_prime=-1)
        for name, X in (("gauge", Fg), ("Higgs", Fh)):
            Dw = Ddd + 0.1 * X
            jd = np.abs(Udd @ np.conj(Dw) + Dw @ Udd).max()
            log(f"     fluctuated D + 0.1*(w + w^* - J w J^-1) [{name}]: Hermitian "
                f"{np.abs(Dw - Dw.conj().T).max():.1e}, J D_w = -D_w J to {jd:.1e}")
        # swap product: the gauge one-form acquires a D_F-carrying piece
        gs = W.swap_grading_exact("none").to_numpy()
        Dsw = np.kron(Dg, I8) + np.kron(gs, DF)
        gsw = gauge_one_form(Dsw, a, b, T)
        clean = np.kron(a @ (Dg @ b - b @ Dg), T)
        extra = np.kron(a @ (gs @ b - b @ gs), T @ DF)
        log(f"  n_max={n}  [swap product, the nearest gamma_GV (x) D_F alternative]")
        log(f"     gauge one-form = a[D_GV,b](x)T + a[gamma_s,b](x)T D_F to {np.abs(gsw - clean - extra).max():.1e};"
            f" ||a[D_GV,b](x)T|| = {np.linalg.norm(clean):.6f}, ||extra D_F piece|| = {np.linalg.norm(extra):.6f}")


# =============================================================================
# (d) toy outer factor: flat 2-torus
# =============================================================================


def _toy_ops(R2, p, q):
    W = torus_window(R2)
    De = W.D_exact()
    a_e = W.shift_exact(*q)
    b_e = W.shift_exact(*p)
    D = De.to_numpy()
    a, b = a_e.to_numpy(), b_e.to_numpy()
    om = a @ (D @ b - b @ D)
    A = om + om.conj().T
    phi = a @ b
    C = D @ phi - phi @ D
    return W, De, a_e, b_e, D, A, C


def control_d() -> None:
    log("=" * 78)
    log("(d) TOY WITH A KNOWN ANSWER: flat 2-torus [0,2pi]^2, antiperiodic spin structure.")
    log("    Window |k| <= R (complete Dirac eigenspaces), multipliers e^{i q.x}, e^{i p.x} compressed.")
    log("    Continuum reference by an independent route: symbolic Gaussian k-integral + mpmath s-integral")
    log("    (Poisson: lattice sum = integral up to O(e^{-pi^2/t})).  Known answers: pure gauge (q = 0)")
    log("    gives continuum k_A = 0 exactly (D + eps i[D,chi] is unitarily equivalent to D); for")
    log("    finite = continuum, the ratio machinery c = 4 (k_A/k_A^c)/(k_Phi/k_Phi^c) must return 4.")
    t = 0.1
    # pure gauge continuum
    for p in ((1, 0), (1, 1), (2, 1)):
        log(f"  continuum pure gauge p={p}, q=(0,0), t={t}: k_A^c = {torus_continuum_gauge(p, (0, 0), t):.3e}")
    p, q = (1, 0), (0, 1)
    kAc = torus_continuum_gauge(p, q, t)
    kPc = torus_continuum_higgs_kin(p, q, t)
    log(f"  scored-type pair p={p}, q={q} (P = (1,1), transverse part nonzero), t={t}:")
    log(f"     continuum k_A^c = {kAc:.12f}, k_Phi^c = {kPc:.12f}")
    # Duhamel vs finite difference, small window
    W, De, a_e, b_e, D, A, C = _toy_ops(12, p, q)
    h = 1e-3
    def S(e):
        ev = np.linalg.eigvalsh(D + e * A)
        return np.sum(np.exp(-t * ev ** 2))
    fd = (S(h) - 2 * S(0) + S(-h)) / (2 * h * h)
    fd2 = (S(2 * h) - 2 * S(0) + S(-2 * h)) / (8 * h * h)
    rich = (4 * fd - fd2) / 3
    du = eps2_heat_coefficient(D, A, t)
    log(f"  Duhamel check (R^2=12, dim {W.dim}): Duhamel {du:.12f}, Richardson FD {rich:.12f}, diff {abs(du-rich):.1e}")
    # exact eps^2 moments vs float
    om_e = a_e @ ((De @ b_e) - (b_e @ De))
    A_e = om_e + om_e.dag()
    for nn in (2, 4, 6):
        ex = eps2_trace_power(De, A_e, nn)
        fl = eps2_trace_power(D, A, nn)
        log(f"  exact eps^2 coefficient of Tr (D+eps A)^{nn} = {qstr(ex)}; float {fl.real:.6f}; "
            f"agree {abs(complex(ex) - fl) < 1e-8}")
    # convergence at fixed t
    log(f"  fixed t = {t}: window R^2 -> infinity")
    log("     R^2    dim    k_A(R)          k_Phi(R)        c(R)            pure-gauge k_A(R)")
    for R2 in (8, 16, 32, 64, 128, 256, 400):
        W, De, a_e, b_e, D, A, C = _toy_ops(R2, p, q)
        kA = eps2_heat_coefficient(D, A, t)
        kP = kinetic_duhamel(D, C, t)
        c = 4 * (kA / kAc) / (kP / kPc)
        _, _, _, _, D0, A0, _ = _toy_ops(R2, (1, 1), (0, 0))
        kpg = eps2_heat_coefficient(D0, A0, t)
        log(f"     {R2:<5}  {W.dim:<5}  {kA:.10f}  {kP:.10f}  {c:.12f}  {kpg:.3e}")
    # non-vacuity: the machinery returns c != 4 at fixed tau = t R^2
    tau = 4.0
    log(f"  fixed tau = t R^2 = {tau} (cutoff tied to window): c(R) on the TOY (shows c != 4 is returnable)")
    for R2 in (16, 36, 64, 144):
        tt = tau / R2
        kAc_t = torus_continuum_gauge(p, q, tt)
        kPc_t = torus_continuum_higgs_kin(p, q, tt)
        W, De, a_e, b_e, D, A, C = _toy_ops(R2, p, q)
        c = 4 * (eps2_heat_coefficient(D, A, tt) / kAc_t) / (kinetic_duhamel(D, C, tt) / kPc_t)
        log(f"     R^2={R2:<4} t={tt:.5f}  c = {c:.8f}")
    # product factorization through the full DD product (toy outer, toy D_F)
    F, DF, gF, T, beta = _finite_float()
    W, De, a_e, b_e, D, A, C = _toy_ops(8, p, q)
    Dprod = np.kron(D, gF) + np.kron(np.eye(W.dim), DF)
    Fg = np.kron(A, gF @ T)
    Kfull = eps2_heat_coefficient(Dprod, Fg, t)
    mu, Q = np.linalg.eigh(DF @ DF)
    inner = float(np.real(np.trace(Q.conj().T @ (T @ T) @ Q @ np.diag(np.exp(-t * mu)))))
    kout = eps2_heat_coefficient(D, A, t)
    log(f"  full DD product, gauge (toy D_F, R^2=8): K_A^full = {Kfull:.12f}; "
        f"Tr_F(T^2 e^(-tD_F^2)) * k_A^out = {inner * kout:.12f}; diff {abs(Kfull - inner*kout):.1e}")
    sig = DF @ beta - beta @ DF
    phi = a_e.to_numpy() @ b_e.to_numpy()
    Xkin = np.kron(D @ phi - phi @ D, gF @ sig)
    Xkin = Xkin + Xkin.conj().T
    Kkin = kinetic_duhamel(Dprod, Xkin, t) * 1.0
    trsig = float(np.real(np.trace(sig.conj().T @ sig)))
    log(f"  full DD product, Higgs kin-kin (R^2=8): (t^2/2)int Tr(X e X e) with X = [D,phi](x)gamma_F sigma + h.c.:"
        f" {Kkin:.10f}; 2*Tr_F(sigma^dag sigma) * k_Phi^out = {2 * trsig * kinetic_duhamel(D, C, t):.10f}"
        " (equal only at D_F -> 0 weights)")
    F0 = FiniteExact(Fraction(0), Fraction(0))
    # D_F-weight-free check: replace D_F by 0 in the heat weights
    Dprod0 = np.kron(D, gF)
    log(f"     same with D_F removed from the weights: {kinetic_duhamel(Dprod0, Xkin, t):.10f} "
        f"vs {2 * trsig * kinetic_duhamel(D, C, t):.10f}")


def main() -> None:
    t0 = time.time()
    log("graded_ac controls  (debug/graded_ac/controls.py)  -- no scored quantity computed")
    control_a()
    control_b()
    control_c()
    control_d()
    log("=" * 78)
    log(f"wall time {time.time() - t0:.1f} s")
    (HERE / "controls.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
