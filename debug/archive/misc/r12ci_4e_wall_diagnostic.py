"""N=4 WALL DIAGNOSTIC — is the disjoint-pair 4-body term RI-forced, or exact-no-RI-but-costly?

Context (2026-09-20). The exact-algebraic explicit-r12 (James-Coolidge style, r12 in the
BASIS, NO strong-orthogonality projector -> NO resolution-of-identity) closes with no RI
through N=3 (memo debug/sprint_neumann_r12_build_memo.md sec.9b): the term inventory of
<Phi| F H F |Phi> (F = sum_{p<q} f_pq, multiplicative) reduces via three angular rules --
  RULE A  <P_a(1.2) P_b(1.3)>          = d_a0 d_b0     (shared vertex, scalar -> L=0)
  RULE B  <(r12^.r13^) h g>            = factorizes    (shared vertex, vec-vec)
  TRIANGLE <P_a(12) P_b(13) P_c(23)>   = d_abc/(2a+1)^2 (closed loop, 3 shared vertices)
The TRIANGLE closes the N=3 case because every pair shares a vertex with every other pair,
so the delta collapses the triple multipole sum to ONE finite index.

At N=4 a genuinely new object appears in <Phi| F H F |Phi>:
    T4  =  f(r12) * f(r34) * (1/r13)      [left F -> f12, H -> 1/r13, right F -> f34]
Four DISTINCT electrons; the graph is a CHAIN 2-1-3-4 (edge f on 12, bridging Coulomb on
13, edge f on 34), NOT a triangle. Vertices 2 and 4 are LEAVES; vertices 1 and 3 each carry
two legs. The memo calls this "the wall -- first case that genuinely needs 4-body operators,
no <=3-body reduction."

WHAT THIS DIAGNOSTIC DECIDES (not what it asserts):
  Q1 TERMINATION. Expand 1/r13 = sum_L g_L(r1,r3) P_L(u1.u3). Is the angular factor
     J(a,c,L; l1..l4) = < w1 w2 w3 w4  P_a(u1.u2) P_c(u3.u4) P_L(u1.u3) >  ZERO for L beyond
     a finite bound set by the orbital angular momenta l_i and the correlation multipoles a,c?
       terminates at finite L  -> the infinite Coulomb sum TRUNCATES EXACTLY -> exact, no RI.
       nonzero for unbounded L  -> RI / truncation forced.
     (w_i = |Y_{l_i m_i}(u_i)|^2, the orbital density's angular weight.)
  Q2 REDUCIBILITY. Is the connected 4-body kernel J expressible as a product/contraction of
     the existing <=3-body kernels (RULE A / TRIANGLE), or is it an IRREDUCIBLE 4-vertex
     (9j-type) object -- i.e. does the corpus's "needs 4-body operators" mean "new native
     kernel" (still exact, no RI) rather than "RI-forced"?

HOW THIS DIFFERS FROM THE DOCUMENTED DEAD END (ledger 2026-08-23, "TC three-body operator
collapse via Gaunt/6j"): that tried to COLLAPSE a 3-body operator to 2-body via an abelian
plane-wave momentum trick and FAILED (GeoVac's Y.Y=sum_Lambda is non-abelian). This does the
OPPOSITE -- it accepts no collapse and measures whether the genuine 4-body integral is still
finite/closed-form (exact, no RI) or truly RI-forced. Hermitian variational scalar object,
not the non-Hermitian TC commutator with vector legs.

Pure diagnostic: builds and validates angular machinery only, no physics/energies.
"""
import numpy as np
from numpy.polynomial.legendre import legval

rng = np.random.default_rng(20260920)


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def rand_dirs(m):
    v = rng.normal(size=(m, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def P(L, x):
    cf = np.zeros(L + 1); cf[L] = 1.0
    return legval(x, cf)


def dot(u, v):
    return np.einsum("ij,ij->i", u, v)


# real-spherical-harmonic DENSITY weight |Y_lm|^2 as a function of a direction u,
# built from the polar angle cos(theta) = u_z.  (m choice only sets which axis;
# the multipole CONTENT -- P_0..P_2l -- is what matters for termination.)
def density_weight(u, l, m=0):
    if l == 0:
        return np.ones(u.shape[0])
    ct = u[:, 2]
    st = np.sqrt(np.maximum(1.0 - ct * ct, 0.0))
    if l == 1:
        if m == 0:                    # |Y10|^2 ~ cos^2 : content P0 + P2
            return ct * ct
        return 0.5 * st * st          # |Y1,+-1|^2 ~ sin^2 : content P0 + P2
    if l == 2:
        if m == 0:                    # |Y20|^2 ~ (3cos^2-1)^2 : content P0+P2+P4
            return (3 * ct * ct - 1.0) ** 2
        if abs(m) == 1:
            return (st * ct) ** 2
        return st ** 4
    raise ValueError("l<=2 only")


# --------------------------------------------------------------------------- #
# the core object:  J(a,c,L; l1,l2,l3,l4)  by Monte-Carlo orientation average
#   chain 2-1-3-4 :  P_a on (1,2) , P_c on (3,4) , P_L on (1,3) ; densities on all four
# --------------------------------------------------------------------------- #
def J_chain(a, c, L, ls, ms=(0, 0, 0, 0), n=2_000_000):
    """< w1 w2 w3 w4  P_a(u1.u2) P_c(u3.u4) P_L(u1.u3) >  (isotropic 4-direction average)."""
    l1, l2, l3, l4 = ls
    m1, m2, m3, m4 = ms
    u1, u2, u3, u4 = rand_dirs(n), rand_dirs(n), rand_dirs(n), rand_dirs(n)
    w = (density_weight(u1, l1, m1) * density_weight(u2, l2, m2)
         * density_weight(u3, l3, m3) * density_weight(u4, l4, m4))
    integ = w * P(a, dot(u1, u2)) * P(c, dot(u3, u4)) * P(L, dot(u1, u3))
    val = float(integ.mean())
    err = float(integ.std() / np.sqrt(n))
    return val, err


# --------------------------------------------------------------------------- #
# PART 1 : term inventory extension to N=4  (analytic, printed)
# --------------------------------------------------------------------------- #
print("=" * 78)
print("PART 1.  Term inventory of <Phi| F H F |Phi> at N=4  (F = sum f_pq, s-focus)")
print("=" * 78)
rows = [
    ("indices within a 3-subset", "reduces to the N=3 inventory", "RULE A / TRIANGLE (done)"),
    ("f_ij f_ij (1/r_kl)  {ij}!={kl}", "same-pair f^2 x disjoint Coulomb", "factorizes: 2-body x 2-body"),
    ("f_ij (1/r_ij) f_kl  {ij},{kl} disj", "Coulomb ON an f-edge + disjoint f", "factorizes: 2-body x 2-body"),
    ("f_12 f_34 (1/r_13)  4 DISTINCT   ", "CHAIN 2-1-3-4 (Coulomb BRIDGES)  ", ">>> THE NEW 4-BODY OBJECT <<<"),
]
print(f"{'term shape':<38}{'graph':<34}{'status'}")
for t, g, s in rows:
    print(f"{t:<38}{g:<34}{s}")
print("\nThe bridging Coulomb (1/r13, or 1/r24) is the only way to connect the two disjoint")
print("f-edges into ONE 4-electron object. That chain is what N<=3 never produced.")


# --------------------------------------------------------------------------- #
# PART 2 : Q1 TERMINATION -- does the Coulomb multipole sum truncate at finite L?
# --------------------------------------------------------------------------- #
print("\n" + "=" * 78)
print("PART 2 (Q1).  Does J(a,c,L) vanish for L beyond a finite, l-set bound?")
print("=" * 78)
print("For each orbital-l configuration we scan L and report the largest L with |J|>>noise.")
print("Leaves (electrons 2,4) carry the f-multipoles a,c; if a leaf orbital is s (l=0),")
print("RULE A forces its f-multipole to 0.  Correlation content taken as a,c in {0,1}.\n")


def largest_nonzero_L(a, c, ls, ms=(0, 0, 0, 0), Lmax=10, n=2_000_000, k=6.0):
    """Return (bound, table) : largest L with |J|>k*mc_err, and the per-L values."""
    tab = []
    bound = -1
    for L in range(Lmax + 1):
        v, e = J_chain(a, c, L, ls, ms, n=n)
        sig = abs(v) > k * e
        tab.append((L, v, e, sig))
        if sig:
            bound = L
    return bound, tab


configs = [
    # (label, ls, ms, a, c)   -- a,c are the f-edge multipoles on the LEAF pairs
    ("all-s  (l=0000), a=c=0            ", (0, 0, 0, 0), (0, 0, 0, 0), 0, 0),
    ("all-s  (l=0000), a=1,c=0          ", (0, 0, 0, 0), (0, 0, 0, 0), 1, 0),
    ("bridge p on e1,e3 (l=1010), a=c=0 ", (1, 0, 1, 0), (0, 0, 0, 0), 0, 0),
    ("bridge d on e1,e3 (l=2020), a=c=0 ", (2, 0, 2, 0), (0, 0, 0, 0), 0, 0),
    ("leaf p on e2,e4  (l=0101), a=c=1  ", (0, 1, 0, 1), (0, 0, 0, 0), 1, 1),
    ("all-p (l=1111), a=c=1             ", (1, 1, 1, 1), (0, 0, 0, 0), 1, 1),
]
print(f"{'configuration':<36}{'a':>2}{'c':>2}  {'largest-L nonzero':>18}   {'2*l1+2*l3(+a) bound'}")
results = {}
for label, ls, ms, a, c in configs:
    bound, tab = largest_nonzero_L(a, c, ls, ms)
    # heuristic termination bound: Coulomb vertex-1 couples P_L to density-l1 and P_a;
    # vertex-3 couples P_L to density-l3 and P_c.  So L <= 2*l1 + a  AND  L <= 2*l3 + c.
    l1, l2, l3, l4 = ls
    pred = min(2 * l1 + a, 2 * l3 + c)
    results[label] = (bound, pred, tab)
    print(f"{label:<36}{a:>2}{c:>2}  {bound:>18}   L<=min(2l1+a,2l3+c)={pred}")

print("\nPer-L detail (value +/- mc_err ; * = significant):")
for label, ls, ms, a, c in configs:
    bound, pred, tab = results[label]
    cells = "  ".join(f"L{L}:{v:+.4f}{'*' if s else ' '}" for L, v, e, s in tab if L <= max(pred + 3, 5))
    print(f"  {label.strip():<34} {cells}")


# --------------------------------------------------------------------------- #
# PART 3 : Q2 REDUCIBILITY -- is the connected chain a PRODUCT of <=3-body kernels?
# --------------------------------------------------------------------------- #
print("\n" + "=" * 78)
print("PART 3 (Q2).  Is the connected chain reducible to <=3-body kernels, or 4-body-native?")
print("=" * 78)
print("Hypothesis H0 (reducible): the chain factorizes across the bridge as")
print("   J_chain(a,c,L; l1234) =?  Jvertex1(a,L; l1,l2) * Jvertex3(c,L; l3,l4) / norm_L")
print("i.e. the two ends are independent given the shared multipole L on the bridge.")
print("Test: compare J_chain to that product for a case where all legs are ACTIVE (l>0).\n")


def J_vertex(mult_leaf, L, l_center, l_leaf, m=(0, 0), n=2_000_000):
    """< w_center w_leaf  P_{mult_leaf}(uc.ul) P_L(uc.uref) >  with uref the bridge partner
    (isotropic).  This is the 'half chain' at one vertex: leaf-leg P_mult and bridge-leg P_L
    meeting at the center, averaged with the bridge partner isotropic."""
    uc, ul, uref = rand_dirs(n), rand_dirs(n), rand_dirs(n)
    w = density_weight(uc, l_center, m[0]) * density_weight(ul, l_leaf, m[1])
    integ = w * P(mult_leaf, dot(uc, ul)) * P(L, dot(uc, uref))
    return float(integ.mean()), float(integ.std() / np.sqrt(n))


# active test case: bridge carries l on e1,e3; leaves carry l on e2,e4; a=c=1, L on bridge.
ls_test = (1, 1, 1, 1)
for L in (0, 2):
    jc, ec = J_chain(1, 1, L, ls_test, n=4_000_000)
    j1, e1 = J_vertex(1, L, 1, 1, n=4_000_000)   # vertex 1: center e1(l1) leaf e2(l2), a=1
    j3, e3 = J_vertex(1, L, 1, 1, n=4_000_000)   # vertex 3: center e3(l3) leaf e4(l4), c=1
    # candidate reducible reconstruction (independent-ends product, normalized by the
    # bridge self-overlap <w1 w3 P_L(u1.u3)> so units match a factorization):
    ub1, ub3 = rand_dirs(4_000_000), rand_dirs(4_000_000)
    br = density_weight(ub1, 1) * density_weight(ub3, 1) * P(L, dot(ub1, ub3))
    norm_L = float(br.mean())
    recon = (j1 * j3 / norm_L) if abs(norm_L) > 1e-9 else np.nan
    rel = abs(jc - recon) / max(abs(jc), 1e-12)
    print(f"L={L}:  chain={jc:+.6f}   product/normL={recon:+.6f}   rel.diff={rel:.2e}"
          f"   {'REDUCIBLE' if rel < 0.03 else 'IRREDUCIBLE (4-body native)'}")

print("\n(If rel.diff is large at L>0, the chain does NOT factor across the bridge -> the")
print(" angular kernel is an irreducible 4-vertex object, i.e. a native 4-body kernel is")
print(" required.  It is still FINITE/closed-form iff Part 2 shows termination.)")


# --------------------------------------------------------------------------- #
# PART 4 : verdict scaffold (printed; interpretation in the session write-up)
# --------------------------------------------------------------------------- #
print("\n" + "=" * 78)
print("PART 4.  Verdict scaffold")
print("=" * 78)
print("Q1 (termination) decides RI-forced vs exact-no-RI:")
print("   finite L-bound tracking 2*l_max  -> exact, no RI (soft wall)")
print("   unbounded L                       -> RI / truncation forced (hard wall)")
print("Q2 (reducibility) decides the COST of the soft wall:")
print("   product of <=3-body kernels       -> cheap, no new machinery")
print("   irreducible 4-vertex (9j) kernel  -> new native 4-body kernel + 4-index radial")
print("The quantum-encoding axis is separate: an explicit 4-body correlation operator is a")
print("4-body Pauli string regardless of classical closed-form (cf. ledger: TC 2nd-quant")
print("plateau 3.4%, TC angular-gradient 2.66x Pauli) -- that wall is real for the product.")
