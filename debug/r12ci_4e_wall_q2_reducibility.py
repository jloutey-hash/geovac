"""N=4 WALL DIAGNOSTIC, Q2 (reducibility / cost).  Follows r12ci_4e_wall_diagnostic.py.

Q1 already decided (that file): the bridging-Coulomb multipole sum in the chain
    T4 = f(r12) f(r34) (1/r13)          [graph 2-1-3-4, four distinct electrons]
TERMINATES at L <= 2*l_bridge (orbital content), so the 4-electron integral is EXACT and
RI-FREE -- finite closed-form angular sum.  The remaining question is COST:

  Q2. Is the connected chain REDUCIBLE to a contraction of vertex kernels (cheap, like the
      N=3 triangle_contract), or an IRREDUCIBLE 4-index angular blow-up?

The Legendre addition theorem makes this decidable exactly:
   P_L(u1.u3) = 4pi/(2L+1) * sum_{M=-L..L} Y_LM(u1) conj(Y_LM(u3))
so the chain FACTORIZES ACROSS THE BRIDGE:
   J_chain(a,c,L; l1234) = 4pi/(2L+1) * sum_M  A_LM  conj(C_LM)
   A_LM = < w1(u1) w2(u2) P_a(u1.u2) Y_LM(u1) >     (vertex-1 half: center e1, leaf e2)
   C_LM = < w3(u3) w4(u4) P_c(u3.u4) Y_LM(u3) >     (vertex-3 half: center e3, leaf e4)
Each half is a THREE-body-ish object (center + leaf + one free bridge index).  If this
identity holds numerically on an ACTIVE (non-null) config, the 4-body chain is REDUCIBLE to
a finite (L,M) contraction of vertex kernels -> the classical wall is SOFT (cost, not RI,
not even a new irreducible kernel).

This corrects Q2 in the first driver, which picked an odd-multipole config that is
identically zero against even diagonal densities |Y_l|^2 (so it divided noise by noise).
Here: even correlation multipoles a=c=2 with l=1 leaves and l=1 bridge -> all legs active.
"""
import numpy as np
from numpy.polynomial.legendre import legval
from scipy.special import lpmv
from math import factorial, sqrt, pi

rng = np.random.default_rng(20260920)


def rand_dirs(m):
    v = rng.normal(size=(m, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def dot(u, v):
    return np.einsum("ij,ij->i", u, v)


def Pleg(L, x):
    cf = np.zeros(L + 1); cf[L] = 1.0
    return legval(x, cf)


def density_weight(u, l, m=0):
    if l == 0:
        return np.ones(u.shape[0])
    ct = u[:, 2]; st = np.sqrt(np.maximum(1 - ct * ct, 0.0))
    if l == 1:
        return ct * ct if m == 0 else 0.5 * st * st
    if l == 2:
        if m == 0:
            return (3 * ct * ct - 1.0) ** 2
        if abs(m) == 1:
            return (st * ct) ** 2
        return st ** 4
    raise ValueError("l<=2")


def Ycomplex(L, M, u):
    """Complex spherical harmonic Y_LM(theta,phi), theta from u_z, phi from (u_x,u_y)."""
    ct = u[:, 2]
    phi = np.arctan2(u[:, 1], u[:, 0])
    am = abs(M)
    norm = sqrt((2 * L + 1) / (4 * pi) * factorial(L - am) / factorial(L + am))
    val = norm * lpmv(am, L, ct) * np.exp(1j * am * phi)
    if M < 0:
        val = (-1) ** am * np.conj(val)
    return val


def J_chain(a, c, L, ls, n=4_000_000):
    l1, l2, l3, l4 = ls
    u1, u2, u3, u4 = rand_dirs(n), rand_dirs(n), rand_dirs(n), rand_dirs(n)
    w = (density_weight(u1, l1) * density_weight(u2, l2)
         * density_weight(u3, l3) * density_weight(u4, l4))
    integ = w * Pleg(a, dot(u1, u2)) * Pleg(c, dot(u3, u4)) * Pleg(L, dot(u1, u3))
    return float(integ.mean()), float(integ.std() / np.sqrt(n))


def vertex_ALM(a, L, M, l_center, l_leaf, n=4_000_000):
    """A_LM = < w_c(uc) w_l(ul) P_a(uc.ul) Y_LM(uc) >  (complex)."""
    uc, ul = rand_dirs(n), rand_dirs(n)
    w = density_weight(uc, l_center) * density_weight(ul, l_leaf)
    integ = w * Pleg(a, dot(uc, ul)) * Ycomplex(L, M, uc)
    re = float(integ.real.mean()); im = float(integ.imag.mean())
    return re + 1j * im


print("=" * 78)
print("Q2.  Bridge factorization test (Legendre addition theorem)")
print("     config: leaves e2,e4 l=1 ; bridge e1,e3 l=1 ; correlation a=c=2 (all active)")
print("=" * 78)
ls = (1, 1, 1, 1)   # (l1,l2,l3,l4)
a = c = 2
L_term = 2 * ls[0]  # Q1 termination bound for this bridge l: L<=2*l_bridge; higher L is 0
print(f"Q1 termination bound here: bridge l={ls[0]} -> nonzero only for L<={L_term} "
      f"(higher L is identically 0; testing factorization there is noise/noise).")
print(f"{'L':>2}  {'J_chain (direct MC)':>22}  {'sum_M ALM conj(CLM) form':>26}  {'rel.diff':>10}  verdict")
allmatch = True
for L in (0, 2, 4):
    jc, ec = J_chain(a, c, L, ls, n=6_000_000)
    # A_LM (vertex 1: center l1=1, leaf l2=1) ; C_LM (vertex 3: center l3=1, leaf l4=1) -- identical here
    recon = 0.0 + 0.0j
    for M in range(-L, L + 1):
        A = vertex_ALM(a, L, M, ls[0], ls[1], n=6_000_000)
        C = vertex_ALM(c, L, M, ls[2], ls[3], n=6_000_000)
        recon += A * np.conj(C)
    recon *= 4 * pi / (2 * L + 1)
    recon = recon.real
    rel = abs(jc - recon) / max(abs(jc), 1e-9)
    signal = L <= L_term and abs(jc) > 10 * ec        # a channel Q1 says is nonzero
    if signal:
        ok = rel < 0.05
        allmatch &= ok
        tag = "FACTORIZES (reducible)" if ok else "does NOT factorize"
    else:
        tag = f"L>{L_term}: ~0 by Q1 (skip; MC noise)"
    print(f"{L:>2}  {jc:>22.6f}  {recon:>26.6f}  {rel:>10.2e}  {tag}")

print()
if allmatch:
    print("VERDICT Q2: the chain FACTORIZES across the bridge -> the connected 4-body angular")
    print("object is a finite (L,M)-contraction of THREE-body vertex kernels A_LM, C_LM.")
    print("It is exact, RI-free (Q1), and reducible: cost is a bridge (L,M) sum of vertex-kernel")
    print("products, structurally the chain analogue of the N=3 triangle_contract -- NOT an")
    print("irreducible 4-index blow-up and NOT a new fundamental kernel class.")
else:
    print("VERDICT Q2: does not factorize -> irreducible 4-vertex object (native 4-body kernel).")

print()
print("-" * 78)
print("RADIAL COST of the reduced form (structural):")
print("  T4 = sum_L (4pi/(2L+1)) sum_M   [ A_LM-radial(r1,r2) ] x g_L(r1,r3) x [ C_LM-radial(r3,r4) ]")
print("  = a separable contraction: build vertex tensors once, contract over the bridge (L,M).")
print("  With n grid points per radial coord this is O(n^2) vertex tensors x O(L_max^2) bridge")
print("  terms x an O(n^2) Coulomb-bridge matmul -> polynomial, no exponential/RI blow-up.")
print("  Contrast the true HARD walls (unchanged, separate axes):")
print("   * KINETIC/vector 4-body terms (grad legs, RULE B) -- NOT tested here; may be harder.")
print("   * QUANTUM ENCODING: an explicit 4-body correlation operator is a 4-body Pauli string")
print("     (ledger: TC 2nd-quant plateau 3.4%, TC angular-gradient 2.66x Pauli) -- REAL wall")
print("     for the quantum-simulation product, independent of the classical closed form.")
