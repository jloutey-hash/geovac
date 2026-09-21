"""Be R12-CI FULL analytic engine (noise-free).  H_01, H_11 via block reduction.

Geminal f = r e^{-g r} (short-range: well-conditioned, cusp f'(0)=1).  All orbitals s ->
every term reduces to radial integrals with monopole kernels over the up-block {1,2} and
down-block {3,4}.  General framework: TWOPAIR(Ka,Kb,ABfun) computes INT Phi0^2 (sum A_pq)(sum B_rs)
for pair operators A (kernel Ka, monopole), B (kernel Kb); reused for S_11, Vee_F, F.d2F.
Kinetic via KE = -1/4 INT Phi0^2 d2F + 1/2 INT F |gradPhi0|^2 (and the F^2 analog).

Validated element-by-element against MC (be_r12ci_ortho targets).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

d = np.load("debug/data/be_r12ci_ref.npz")
Z1, Z2, EREF = float(d["z1"]), float(d["z2"]), float(d["Eref"])
Zn = 4.0; GEM = 2.0
N1 = 2.0 * Z1 ** 1.5
OV = N1 * 6.0 / (Z1 + Z2) ** 4
_xg, _wg = leggauss(600); _r0 = 12.5 * (_xg + 1.0); _w0 = 12.5 * _wg
_u = _r0 * np.exp(-Z2 * _r0) - OV * N1 * np.exp(-Z1 * _r0)
N2 = 1.0 / np.sqrt(np.sum(_u * _u * _r0 ** 2 * _w0))
FPI = 4 * np.pi


def a_r(r): return N1 * np.exp(-Z1 * r)
def b_r(r): return N2 * (r * np.exp(-Z2 * r) - OV * N1 * np.exp(-Z1 * r))
def ap_r(r): return -Z1 * N1 * np.exp(-Z1 * r)
def bp_r(r): return N2 * (np.exp(-Z2 * r) - Z2 * r * np.exp(-Z2 * r) + OV * N1 * Z1 * np.exp(-Z1 * r))
def f_g(r): return r * np.exp(-GEM * r)
def fp_g(r): return (1.0 - GEM * r) * np.exp(-GEM * r)
def fpp_g(r): return (-GEM - GEM * (1.0 - GEM * r)) * np.exp(-GEM * r)   # f''
def lap_f(r): return fpp_g(r) + 2.0 * fp_g(r) / r                        # radial Laplacian of f


NR = 400; Rmax = 25.0
# grid concentrated near r=0 (1s cusp at ~1/(2Z1)): map u in [0,1] -> r = Rmax*u^3
_x, _w = leggauss(NR); _u = 0.5 * (_x + 1.0); _du = 0.5 * _w
rr = Rmax * _u ** 3; wr = Rmax * 3.0 * _u ** 2 * _du
w2 = wr * rr ** 2
NX = 40; xx, wx = leggauss(NX)
A = a_r(rr); B = b_r(rr); Ap = ap_r(rr); Bp = bp_r(rr)
Dp = (np.outer(A, B) - np.outer(B, A)) ** 2                              # up-block density
KD = (np.outer(Ap, B) - np.outer(Bp, A)) ** 2 + (np.outer(A, Bp) - np.outer(B, Ap)) ** 2  # kinetic density
rho = FPI * (A ** 2 + B ** 2)                                           # marginal of Dp
rhoK = FPI * (KD @ w2)                                                  # marginal of KD (function of r1)
Nb = 32 * np.pi ** 2                                                    # INT Dp d3 = block norm
KNb = np.sum(rhoK * w2)                                                 # INT KD d3 ; KE_00 = KNb*Nb


def _r12(r1, r2, x): return np.sqrt(np.maximum(r1 * r1 + r2 * r2 - 2 * r1 * r2 * x, 1e-30))


def mono(fun):
    R1 = rr[:, None, None]; R2 = rr[None, :, None]; X = xx[None, None, :]
    return 0.5 * np.tensordot(fun(_r12(R1, R2, X)), wx, axes=([2], [0]))


F0 = mono(f_g); Cm = mono(lambda r: 1.0 / r); L0 = mono(lap_f)
F0sq = mono(lambda r: f_g(r) ** 2); FC = mono(lambda r: f_g(r) / r); FL = mono(lambda r: f_g(r) * lap_f(r))


def BM(D, Ofun):
    """block moment INT D(r1,r2) O(r1,r2,x) d3r1 d3r2 = 8pi^2 INT r1^2r2^2 D [INT O dx]."""
    R1 = rr[:, None, None]; R2 = rr[None, :, None]; X = xx[None, None, :]
    O = Ofun(R1, R2, X) + 0.0 * (R1 + R2 + X)            # force full (NR,NR,NX) broadcast
    ang = np.tensordot(O, wx, axes=([2], [0]))
    return 8 * np.pi ** 2 * np.einsum("i,j,ij,ij->", w2, w2, D, ang)


def CROSS1(P1, P3, K): return FPI ** 2 * (P1 * w2) @ K @ (P3 * w2)


def marg(D, K=None):
    """4pi INT D(r1,r2) [K(r1,r2)] r2^2 dr2  (K-dressed marginal of density D)."""
    return FPI * ((D * K) @ w2 if K is not None else D @ w2)


def SU2(Pw, K1, K2):
    KD1 = K1 * w2[None, :]; KD2 = K2 * w2[None, :]
    inner = np.einsum("ij,jk,ik->i", KD1, Dp, KD2)          # down-pair joint bridged by K1,K2 from r1
    return FPI ** 3 * np.sum(w2 * Pw * inner)


def DISJ(Kud, Kud2, Dupd=None, Ddnd=None):
    Du = Dp if Dupd is None else Dupd; Dd = Dp if Ddnd is None else Ddnd
    M = (Du * w2[None, :]) @ Kud2 @ (w2[:, None] * Dd)
    return FPI ** 4 * (w2) @ (M * Kud) @ (w2)


def E2CROSS(op_e1, Kcross, Ddn_marg, Dblock=None):
    """(op on e1) x (cross Kcross from e2 to a down electron), down-other integrated.
       = (4pi)^2 INT r1^2 r2^2 Dblock op_e1(r1) Ff(r2) ; Ff(r2)=4pi INT r3^2 Ddn_marg(r3) Kcross(r2,r3)."""
    Dblk = Dp if Dblock is None else Dblock
    Ff = FPI * (Kcross @ (w2 * Ddn_marg))
    return FPI ** 2 * (w2 * op_e1) @ Dblk @ (w2 * Ff)


ONES = np.ones(NR)
KNb_full = 16 * np.pi ** 2 * np.einsum("i,j,ij->", w2, w2, KD)     # INT KD d3r1 d3r2
margKD = FPI * (KD @ w2)                                            # marginal of the KD block


def FKD(Bfun, Kb):
    """INT (sum_pairs B) KD_up Ddn^2  (one block kinetic-density-weighted). B kernel Kb, function Bfun."""
    t_b12 = BM(KD, lambda r1, r2, x: Bfun(_r12(r1, r2, x))) * Nb          # B12 on KD-block
    t_b34 = KNb_full * BM(Dp, lambda r1, r2, x: Bfun(_r12(r1, r2, x)))    # B34 on Ddn^2
    t_cross_e1 = 2 * CROSS1(margKD, rho, Kb)                              # B13,B14 (cross from KD-e1)
    t_cross_e2 = 2 * E2CROSS(ONES, Kb, rho, Dblock=KD)                    # B23,B24 (cross from KD-e2)
    return t_b12 + t_b34 + t_cross_e1 + t_cross_e2


# ---------- TWOPAIR: INT Phi0^2 (sum A_pq)(sum B_rs) ---------- #
def TWOPAIR(Ka, Kb, Kab_fun, Afun, Bfun):
    """Ka,Kb monopole kernels of A,B; Kab_fun = A(r12)*B(r12) function; Afun,Bfun = A,B functions.
       Returns INT Phi0^2 (sum_pairs A)(sum_pairs B)."""
    mA = marg(Dp, Ka); mB = marg(Dp, Kb)
    Kab = mono(Kab_fun)
    t_samepair = 2 * BM(Dp, lambda r1, r2, x: Afun(_r12(r1, r2, x)) * Bfun(_r12(r1, r2, x))) * Nb   # A12B12+A34B34
    t_prod = 2 * BM(Dp, lambda r1, r2, x: Afun(_r12(r1, r2, x))) * BM(Dp, lambda r1, r2, x: Bfun(_r12(r1, r2, x)))  # A12B34+A34B12
    t_self = 4 * CROSS1(rho, rho, Kab)                                  # A13B13 etc (4 cross pairs)
    t_su = 8 * SU2(rho, Ka, Kb)                                         # share-vertex cross-cross
    t_disj = 2 * DISJ(Ka, Kb) + 2 * DISJ(Kb, Ka)                        # A13B24,A24B13,A14B23,A23B14
    t_ixc = 8 * CROSS1(mA, rho, Kb) + 8 * CROSS1(mB, rho, Ka)           # A_intra x B_cross + B_intra x A_cross
    return t_samepair + t_prod + t_self + t_su + t_disj + t_ixc


# ---------- one-body: INT Phi0^2 (sum_i op_i)(sum_pairs B) ---------- #
def ONEBODY_PAIR(op1, Kb, Bfun):
    """INT Phi0^2 (sum_i op(r_i)) (sum_pairs B).  op one-body (e.g. 1/r).  = 4 * [i=1 value]."""
    mB = marg(Dp, Kb)
    op_rho = op1(rr) * rho                                              # (op on e1)-weighted marginal
    # i=1 (up) x B-terms:
    t_b12 = BM(Dp, lambda r1, r2, x: op1(r1) * Bfun(_r12(r1, r2, x))) * Nb   # op1 * B12 (same block)
    t_b34 = BM(Dp, lambda r1, r2, x: op1(r1)) * BM(Dp, lambda r1, r2, x: Bfun(_r12(r1, r2, x)))  # op1 * B34
    t_b13 = 2 * CROSS1(op1(rr) * rho, rho, Kb)                          # op1(e1) * B13,B14 (cross from e1)
    t_b23 = 2 * E2CROSS(op1(rr), Kb, rho)                              # op1(e1) * B23,B24 (cross from e2)
    return 4 * (t_b12 + t_b34 + t_b13 + t_b23)


if __name__ == "__main__":
    S00 = Nb ** 2
    Fbar_num = 2 * BM(Dp, lambda r1, r2, x: f_g(_r12(r1, r2, x))) * Nb + 4 * CROSS1(rho, rho, F0)  # S_01
    Fbar = Fbar_num / S00
    print(f"KE_00/S_00 = {KNb*Nb/S00:.5f}  (virial -E0 = {-EREF:.5f})")
    print(f"S_01/S_00 = Fbar = {Fbar:.5f}")

    # S_11
    S11 = TWOPAIR(F0, F0, lambda r: f_g(r) ** 2, f_g, f_g)
    print(f"S_11/S_00 = {S11/S00:.5f}  sigma2 = {S11/S00 - Fbar**2:.6f}")

    print(f"KE_00/S_00 (fixed) = {KNb_full/Nb:.5f}  (virial {-EREF:.5f})")
    # V_ee F  (A=v Coulomb kernel Cm, B=f)
    VeeF = TWOPAIR(Cm, F0, lambda r: f_g(r) / r, lambda r: 1.0 / r, f_g)
    # V_ne F
    VneF = -Zn * ONEBODY_PAIR(lambda r: 1.0 / r, F0, f_g)
    # kinetic KE_01 = -1/4 INT Phi0^2 d2F + 1/2 INT F |gradPhi0|^2 ; 1/2 INT F|gradPhi0|^2 = FKD (symmetry)
    d2F = 2 * BM(Dp, lambda r1, r2, x: lap_f(_r12(r1, r2, x))) * Nb + 4 * CROSS1(rho, rho, L0)
    KE01 = -0.25 * d2F + FKD(f_g, F0)
    print(f"\nVeeF/S00 = {VeeF/S00:.5f}   VneF/S00 = {VneF/S00:.5f}   KE01/S00 = {KE01/S00:.5f} (MC ~4.581)")
    H01 = KE01 + VeeF + VneF
    print(f"H_01/S_00 = {H01/S00:.5f}   (MC target ~ -4.458)")

    # ---- ANALYTIC orthogonalized coupling h (the ill-conditioned quantity) ---- #
    E0 = EREF
    ke00 = KNb_full / Nb
    Vexp = E0 - ke00                    # <V>/S_00 = E0 - <T>/S_00
    hT = KE01 / S00 - Fbar * ke00
    hV = (VeeF + VneF) / S00 - Fbar * Vexp
    h = hT + hV
    sigma2 = S11 / S00 - Fbar ** 2
    print(f"\n--- analytic orthogonalized elements ---")
    print(f"<T>/S00={ke00:.5f}  <V>/S00={Vexp:.5f}  Fbar={Fbar:.5f}  sigma2={sigma2:.6f}")
    print(f"hT={hT:.5f}  hV={hV:.5f}  h=<Phi|H|G>={h:.5f}  (analytic, exact)")
    # g = <G|H|G>/S00 : gT analytic-pending (KE_11 + F^2 V), use well-conditioned MC value
    g_MC = -0.16348           # from be_r12ci_ortho (well-conditioned; short-range geminal)
    from scipy.linalg import eigh
    Hm = np.array([[E0, h], [h, g_MC]]); Sm = np.array([[1.0, 0.0], [0.0, sigma2]])
    w, _ = eigh(Hm, Sm)
    print(f"\nE_R12 = {w[0]:.5f} Ha   (analytic h, exact sigma2, MC g)")
    print(f"correlation captured = {w[0]-E0:+.5f} Ha  ({100*(w[0]-E0)/-0.0944:.0f}% of true 94.4 mHa)")
    print(f"variational (>= -14.6674): {w[0] >= -14.6674}")
