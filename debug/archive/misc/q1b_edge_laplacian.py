"""
Track Q1-B: Edge Laplacian of the GeoVac S^3 graph at n_max=3.

Computes the edge (1-form) Laplacian L_edge = B^T B where B is the signed
incidence matrix.  Checks spectral invariants against Paper 2's
K = pi(B + F - Delta) ingredients.

Author: GeoVac sub-agent (Q1-B)
Date: 2026-04-15
"""

import json
import math
from pathlib import Path
from fractions import Fraction

import numpy as np

# ---------------------------------------------------------------------------
# 1.  Build the n_max=3 graph explicitly
# ---------------------------------------------------------------------------

# Nodes: all (n, l, m) with 1 <= n <= 3, 0 <= l < n, -l <= m <= l
nodes = []
for n in range(1, 4):
    for l in range(n):
        for m in range(-l, l + 1):
            nodes.append((n, l, m))

node_index = {s: i for i, s in enumerate(nodes)}
N = len(nodes)
print(f"Nodes ({N}):")
for i, s in enumerate(nodes):
    print(f"  {i:2d}: {s}")

# Edges: angular (same n, l; m <-> m+1) and radial (same l, m; n <-> n+1)
edges = []
for i, (n, l, m) in enumerate(nodes):
    # Angular: m -> m+1
    nb_ang = (n, l, m + 1)
    if nb_ang in node_index:
        j = node_index[nb_ang]
        edge = (min(i, j), max(i, j))
        if edge not in edges:
            edges.append(edge)
    # Radial: n -> n+1 (same l, m, and l < n+1 which is always true)
    nb_rad = (n + 1, l, m)
    if nb_rad in node_index:
        j = node_index[nb_rad]
        edge = (min(i, j), max(i, j))
        if edge not in edges:
            edges.append(edge)

M = len(edges)
print(f"\nEdges ({M}):")
for k, (i, j) in enumerate(edges):
    print(f"  e{k:2d}: {nodes[i]} -- {nodes[j]}  (idx {i}-{j})")

print(f"\nGraph summary: {N} nodes, {M} edges")
print(f"  beta_1 = M - N + 1 = {M - N + 1}")
is_tree = (M == N - 1)
print(f"  Tree: {is_tree}")

# ---------------------------------------------------------------------------
# 2.  Signed incidence matrix B  (N x M)
# ---------------------------------------------------------------------------
# Convention: for edge k = (i, j) with i < j, B[i, k] = +1, B[j, k] = -1.

B = np.zeros((N, M), dtype=int)
for k, (i, j) in enumerate(edges):
    B[i, k] = +1
    B[j, k] = -1

print("\nIncidence matrix B (N x M):")
print(f"  shape = {B.shape}")

# ---------------------------------------------------------------------------
# 3.  Edge Laplacian  L_edge = B^T B  (M x M)
# ---------------------------------------------------------------------------
L_edge = B.T @ B
print("\nEdge Laplacian L_edge = B^T B:")
print(L_edge)

# ---------------------------------------------------------------------------
# 4.  Node Laplacian  L_node = B B^T  (N x N)
# ---------------------------------------------------------------------------
L_node = B @ B.T

# Sanity: compare with D - A
A = np.zeros((N, N), dtype=int)
for (i, j) in edges:
    A[i, j] = 1
    A[j, i] = 1
D = np.diag(A.sum(axis=1))
L_check = D - A
assert np.allclose(L_node, L_check), "L_node != D - A"
print("\nNode Laplacian L_node = B B^T matches D - A: OK")

# ---------------------------------------------------------------------------
# 5.  Eigenvalues — exact rational via sympy
# ---------------------------------------------------------------------------
try:
    import sympy as sp

    L_edge_sym = sp.Matrix(L_edge.tolist())
    L_node_sym = sp.Matrix(L_node.tolist())

    # Characteristic polynomials
    lam = sp.Symbol('lambda')
    char_edge = L_edge_sym.charpoly(lam)
    char_node = L_node_sym.charpoly(lam)

    print("\nEdge Laplacian characteristic polynomial:")
    print(f"  {char_edge.as_expr()}")

    evals_edge_sym = L_edge_sym.eigenvals()   # {eigenvalue: multiplicity}
    evals_node_sym = L_node_sym.eigenvals()

    print("\nEdge Laplacian eigenvalues (sympy exact):")
    evals_edge_list = []
    for ev, mult in sorted(evals_edge_sym.items(), key=lambda x: float(x[0])):
        print(f"  {ev}  (mult {mult})  ≈ {float(ev):.10f}")
        for _ in range(mult):
            evals_edge_list.append(ev)

    print("\nNode Laplacian eigenvalues (sympy exact):")
    evals_node_list = []
    for ev, mult in sorted(evals_node_sym.items(), key=lambda x: float(x[0])):
        print(f"  {ev}  (mult {mult})  ≈ {float(ev):.10f}")
        for _ in range(mult):
            evals_node_list.append(ev)

    use_sympy = True
except Exception as e:
    print(f"\nsympy eigenvalue computation failed ({e}), using numpy")
    evals_edge_list = sorted(np.linalg.eigvalsh(L_edge.astype(float)))
    evals_node_list = sorted(np.linalg.eigvalsh(L_node.astype(float)))
    use_sympy = False

# ---------------------------------------------------------------------------
# 6.  Node / edge spectral relationship  (SVD theorem)
# ---------------------------------------------------------------------------
# Nonzero eigenvalues of BB^T and B^T B are the same.
evals_edge_float = sorted([float(e) for e in evals_edge_list])
evals_node_float = sorted([float(e) for e in evals_node_list])

node_nonzero = sorted([e for e in evals_node_float if e > 1e-12])
edge_nonzero = sorted([e for e in evals_edge_float if e > 1e-12])

print(f"\nNode nonzero eigenvalues ({len(node_nonzero)}): {[f'{e:.6f}' for e in node_nonzero]}")
print(f"Edge nonzero eigenvalues ({len(edge_nonzero)}): {[f'{e:.6f}' for e in edge_nonzero]}")

# For a tree (beta_1=0), the edge Laplacian has no zero eigenvalues,
# but only N-1 = 13 eigenvalues total (= M). The node Laplacian has one zero eigenvalue.
# So node has N-1 nonzero eigenvalues and edge has M = N-1 eigenvalues, ALL nonzero.
# These N-1 nonzero sets should match.
if len(node_nonzero) == len(edge_nonzero):
    max_diff = max(abs(a - b) for a, b in zip(node_nonzero, edge_nonzero))
    print(f"SVD theorem check: max |node_nz - edge_nz| = {max_diff:.2e}  {'OK' if max_diff < 1e-10 else 'FAIL'}")
else:
    print(f"SVD theorem: different counts ({len(node_nonzero)} vs {len(edge_nonzero)})")

# ---------------------------------------------------------------------------
# 7.  Spectral invariants of L_edge vs K ingredients
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("SPECTRAL INVARIANTS vs K = pi(B + F - Delta) INGREDIENTS")
print("=" * 70)

# Target values
B_target = 42
F_target = math.pi**2 / 6
Delta_target = 1 / 40
BpF = B_target + F_target
BpFmD = B_target + F_target - Delta_target
K_target = math.pi * BpFmD
alpha_inv = 137.035999084

targets = {
    "B = 42": 42.0,
    "F = pi^2/6": F_target,
    "Delta = 1/40": Delta_target,
    "B + F": BpF,
    "B + F - Delta": BpFmD,
    "K/pi": BpFmD,
    "K = pi(B+F-Delta)": K_target,
    "1/alpha": alpha_inv,
}

# Compute invariants from the edge spectrum
e_nz = edge_nonzero  # all 13 for a tree
e_all = evals_edge_float

trace_edge = sum(e_all)
det_edge = 1.0
for e in e_all:
    det_edge *= e

# Spectral zeta (on nonzero eigenvalues)
def spectral_zeta(evals, s):
    return sum(e**(-s) for e in evals if e > 1e-12)

zeta1 = spectral_zeta(e_all, 1)
zeta2 = spectral_zeta(e_all, 2)
zeta3 = spectral_zeta(e_all, 3)
zeta4 = spectral_zeta(e_all, 4)

# log-determinant
logdet = sum(math.log(e) for e in e_all if e > 1e-12)

# von Neumann entropy
total = sum(e_all)
probs = [e / total for e in e_all if e > 1e-12]
vn_entropy = -sum(p * math.log(p) for p in probs if p > 0)

# Cheeger-like: algebraic connectivity (smallest nonzero eigenvalue)
alg_conn = min(e for e in e_all if e > 1e-12)

# Largest eigenvalue
largest = max(e_all)

# Sum of squares
sum_sq = sum(e**2 for e in e_all)

# Product of nonzero eigenvalues (tree number for node Laplacian = Kirchhoff)
det_nz = 1.0
for e in e_all:
    if e > 1e-12:
        det_nz *= e

invariants = {
    "Tr(L_edge)": trace_edge,
    "det(L_edge)": det_edge,
    "det'(L_edge) [product of nonzero]": det_nz,
    "log det'(L_edge)": logdet,
    "zeta_edge(1)": zeta1,
    "zeta_edge(2)": zeta2,
    "zeta_edge(3)": zeta3,
    "zeta_edge(4)": zeta4,
    "von Neumann entropy": vn_entropy,
    "algebraic connectivity": alg_conn,
    "largest eigenvalue": largest,
    "sum of squares": sum_sq,
    "M (edge count)": float(M),
    "N (node count)": float(N),
}

# Additional combinations
extra = {
    "Tr/M": trace_edge / M,
    "Tr/N": trace_edge / N,
    "Tr * pi": trace_edge * math.pi,
    "det' * pi": det_nz * math.pi,
    "zeta1 * pi": zeta1 * math.pi,
    "zeta2 * pi^2": zeta2 * math.pi**2,
    "logdet / pi": logdet / math.pi,
    "Tr + 42": trace_edge + 42,
    "Tr - 42": trace_edge - 42,
    "42 / Tr": 42.0 / trace_edge if trace_edge != 0 else float('inf'),
    "Tr / zeta1": trace_edge / zeta1 if zeta1 != 0 else float('inf'),
    "det'^(1/M)": det_nz**(1.0 / M),
    "M * zeta1": M * zeta1,
    "N * zeta1": N * zeta1,
    "pi * zeta1": math.pi * zeta1,
    "pi * zeta2": math.pi * zeta2,
    "pi^2 * zeta1": math.pi**2 * zeta1,
}

all_candidates = {**invariants, **extra}

print("\nEdge Laplacian spectral invariants:")
for name, val in invariants.items():
    print(f"  {name:40s} = {val:.10f}")

print("\nAdditional combinations:")
for name, val in extra.items():
    print(f"  {name:40s} = {val:.10f}")

# Check against targets
print("\n--- Target matching ---")
hits = []
near_misses = []
for cname, cval in all_candidates.items():
    for tname, tval in targets.items():
        if tval == 0:
            continue
        rel_err = abs(cval - tval) / abs(tval)
        if rel_err < 1e-6:
            hits.append((cname, cval, tname, tval, rel_err))
        elif rel_err < 0.05:
            near_misses.append((cname, cval, tname, tval, rel_err))

if hits:
    print("\n  EXACT HITS (rel_err < 1e-6):")
    for cname, cval, tname, tval, rel_err in hits:
        print(f"    {cname} = {cval:.10f}  matches  {tname} = {tval:.10f}  (rel_err {rel_err:.2e})")
else:
    print("\n  No exact hits found.")

if near_misses:
    print(f"\n  NEAR MISSES (rel_err < 5%):")
    near_misses.sort(key=lambda x: x[4])
    for cname, cval, tname, tval, rel_err in near_misses[:15]:
        print(f"    {cname} = {cval:.10f}  vs  {tname} = {tval:.10f}  (rel_err {rel_err:.4f})")
else:
    print("\n  No near misses found.")

# ---------------------------------------------------------------------------
# 8.  Weighted edge Laplacians (Casimir, degeneracy weighting)
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("WEIGHTED EDGE LAPLACIANS")
print("=" * 70)

# Edge Casimir weight: for each edge, weight by l(l+1) of one or both endpoints
# or by (2l+1) degeneracy, etc.

# Weight scheme 1: each edge weighted by average l(l+1) of endpoints
def edge_weighted_laplacian(B_mat, edge_weights):
    """Build B^T W B where W = diag(edge_weights) -- NO, this weights edges in B^T B.
    Actually: L_edge = B^T B for unweighted. For edge-weighted: L_edge = W^{1/2} B^T B W^{1/2}?
    More standard: just weight the adjacency and rebuild."""
    W = np.diag(edge_weights)
    return B_mat.T @ B_mat  # unweighted

# Instead: node-weighted edge Laplacian. Weight each NODE by some function,
# then the edge Laplacian inherits.

# Casimir weight on nodes: w_i = l_i(l_i + 1) for node i = (n_i, l_i, m_i)
casimir_node_weights = np.array([l * (l + 1) for (n, l, m) in nodes], dtype=float)
# Degeneracy weight: w_i = (2l_i + 1)
degen_node_weights = np.array([2 * l + 1 for (n, l, m) in nodes], dtype=float)
# n^2 weight (Fock degeneracy)
fock_node_weights = np.array([n**2 for (n, l, m) in nodes], dtype=float)

# For edge weighting: assign each edge the average (or product, or sum) of endpoint weights
for wname, wvec in [("Casimir l(l+1)", casimir_node_weights),
                      ("Degeneracy (2l+1)", degen_node_weights),
                      ("Fock n^2", fock_node_weights)]:
    # Edge weight = average of endpoint node weights
    ew = np.array([(wvec[i] + wvec[j]) / 2.0 for (i, j) in edges])

    # Weighted incidence: B_w[i, k] = sqrt(ew[k]) * B[i, k]
    sqw = np.sqrt(np.maximum(ew, 0))
    B_w = B * sqw[np.newaxis, :]
    L_edge_w = B_w.T @ B_w

    evals_w = sorted(np.linalg.eigvalsh(L_edge_w))
    tr_w = sum(evals_w)
    nz_w = [e for e in evals_w if e > 1e-12]
    det_w = 1.0
    for e in nz_w:
        det_w *= e

    print(f"\n  {wname}-weighted edge Laplacian:")
    print(f"    Eigenvalues: {[f'{e:.6f}' for e in evals_w]}")
    print(f"    Trace = {tr_w:.6f}")
    print(f"    det' = {det_w:.6f}")

    # Check trace against targets
    for tname, tval in targets.items():
        if tval == 0:
            continue
        rel = abs(tr_w - tval) / abs(tval)
        if rel < 0.05:
            print(f"    ** Trace near {tname}: rel_err = {rel:.4f}")
    for tname, tval in targets.items():
        if tval == 0:
            continue
        rel = abs(det_w - tval) / abs(tval)
        if rel < 0.05:
            print(f"    ** det' near {tname}: rel_err = {rel:.4f}")

# ---------------------------------------------------------------------------
# 9.  Edge-Casimir trace (analog of node Casimir trace = 42)
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("EDGE-CASIMIR ANALOG OF B = 42")
print("=" * 70)

# On nodes: B = sum_{(n,l)} (2l+1) * l(l+1) = 42
# On edges: can we define an analogous sum?

# Classify edges by type
angular_edges = []
radial_edges = []
for k, (i, j) in enumerate(edges):
    ni, li, mi = nodes[i]
    nj, lj, mj = nodes[j]
    if ni == nj and li == lj:  # angular: same n, l
        angular_edges.append(k)
    else:  # radial: same l, m
        radial_edges.append(k)

print(f"\n  Angular edges: {len(angular_edges)}")
print(f"  Radial edges:  {len(radial_edges)}")
print(f"  Total:         {len(angular_edges) + len(radial_edges)}")

# Edge "quantum numbers": for angular edge (n,l,m)-(n,l,m+1), assign (n, l, m+1/2)
# For radial edge (n,l,m)-(n+1,l,m), assign (n+1/2, l, m)
print("\n  Edge classification:")
for k, (i, j) in enumerate(edges):
    ni, li, mi = nodes[i]
    nj, lj, mj = nodes[j]
    if k in angular_edges:
        etype = "angular"
        # Casimir of the edge = l(l+1) of the shell
        edge_casimir = li * (li + 1)
    else:
        etype = "radial"
        # Casimir of the edge: average of endpoint n-shells?
        edge_casimir = li * (li + 1)  # same l for both endpoints
    print(f"    e{k:2d}: {nodes[i]}-{nodes[j]}  {etype:8s}  l={li}  l(l+1)={edge_casimir}")

# Edge Casimir trace: sum over edges of l(l+1)
edge_casimir_sum = sum(li * (li + 1) for k, (i, j) in enumerate(edges)
                       for li in [nodes[i][1]])  # l from first endpoint (same for both on radial)
print(f"\n  Sum over edges of l(l+1) = {edge_casimir_sum}")

# Weighted edge Casimir trace: weight by (2l+1)
edge_weighted_casimir = sum((2 * nodes[i][1] + 1) * nodes[i][1] * (nodes[i][1] + 1)
                            for k, (i, j) in enumerate(edges))
print(f"  Sum over edges of (2l+1)*l(l+1) = {edge_weighted_casimir}")

# Check: node Casimir trace should be 42
node_casimir = sum((2 * l + 1) * l * (l + 1) for (n, l, m) in nodes)
# Actually Paper 2 sums over (n,l) SECTORS, not individual (n,l,m) nodes
sector_casimir = 0
for n in range(1, 4):
    for l in range(n):
        sector_casimir += (2 * l + 1) * l * (l + 1)
print(f"  Node sector Casimir = {sector_casimir}  (should be 42)")
print(f"  Node (n,l,m) Casimir = {node_casimir}")

# ---------------------------------------------------------------------------
# 10.  Exact sympy eigenvalues — display
# ---------------------------------------------------------------------------
if use_sympy:
    print("\n" + "=" * 70)
    print("EXACT EDGE LAPLACIAN SPECTRUM (sympy)")
    print("=" * 70)
    for ev, mult in sorted(evals_edge_sym.items(), key=lambda x: float(x[0])):
        print(f"  {ev}  (mult {mult})")

    # Exact trace
    exact_trace = sum(ev * mult for ev, mult in evals_edge_sym.items())
    print(f"\n  Exact trace = {exact_trace} = {sp.simplify(exact_trace)}")

    # Exact determinant
    exact_det = sp.prod(ev**mult for ev, mult in evals_edge_sym.items())
    print(f"  Exact det = {sp.simplify(exact_det)}")

    # Exact spectral zeta at s=1 (sum of reciprocals)
    zeta1_exact = sum(mult / ev for ev, mult in evals_edge_sym.items() if ev != 0)
    print(f"  Exact zeta(1) = {sp.simplify(zeta1_exact)} ≈ {float(zeta1_exact):.10f}")

    zeta2_exact = sum(mult / ev**2 for ev, mult in evals_edge_sym.items() if ev != 0)
    print(f"  Exact zeta(2) = {sp.simplify(zeta2_exact)} ≈ {float(zeta2_exact):.10f}")

# ---------------------------------------------------------------------------
# 11.  Kirchhoff's theorem check
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("KIRCHHOFF'S THEOREM CHECK")
print("=" * 70)
# For a connected graph, det(L_node with any row/col deleted) = number of spanning trees.
# For a tree, this should equal 1.
# Also: product of nonzero eigenvalues of L_node / N = number of spanning trees
node_nz_prod = 1.0
for e in evals_node_float:
    if e > 1e-12:
        node_nz_prod *= e
kirchhoff = node_nz_prod / N
print(f"  Product of nonzero node Laplacian eigenvalues / N = {kirchhoff:.6f}")
print(f"  (For a tree, should be 1.0)")

# For edge Laplacian of a tree: det(L_edge) = product of ALL eigenvalues (no zeros)
print(f"  det(L_edge) = {det_edge:.6f}")
print(f"  This equals the product of all edge eigenvalues (tree has no zero edge eigenvalues)")

# ---------------------------------------------------------------------------
# 12.  Save results
# ---------------------------------------------------------------------------
results = {
    "meta": {
        "track": "Q1-B",
        "description": "Edge Laplacian of GeoVac S^3 graph at n_max=3",
        "n_nodes": N,
        "n_edges": M,
        "beta_1": M - N + 1,
        "is_tree": is_tree,
    },
    "nodes": [list(s) for s in nodes],
    "edges": [list(e) for e in edges],
    "edge_types": {
        "angular_count": len(angular_edges),
        "radial_count": len(radial_edges),
        "angular_indices": angular_edges,
        "radial_indices": radial_edges,
    },
    "edge_laplacian_spectrum": evals_edge_float,
    "node_laplacian_spectrum": evals_node_float,
    "svd_theorem_max_diff": max(abs(a - b) for a, b in zip(node_nonzero, edge_nonzero)),
    "edge_invariants": {k: v for k, v in invariants.items()},
    "extra_combinations": {k: v for k, v in extra.items()},
    "exact_eigenvalues_str": [str(ev) for ev in evals_edge_list] if use_sympy else None,
    "edge_casimir_sum": edge_casimir_sum,
    "edge_weighted_casimir_sum": edge_weighted_casimir,
    "node_sector_casimir": sector_casimir,
    "targets": {k: v for k, v in targets.items()},
    "exact_hits": [(cn, cv, tn, tv, re) for cn, cv, tn, tv, re in hits] if hits else [],
    "near_misses": [(cn, cv, tn, tv, re) for cn, cv, tn, tv, re in near_misses[:15]] if near_misses else [],
    "verdict": "POSITIVE" if hits else "NEGATIVE",
}

outpath = Path(__file__).parent / "data" / "q1b_edge_laplacian.json"
outpath.parent.mkdir(parents=True, exist_ok=True)
with open(outpath, "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nResults saved to {outpath}")

# ---------------------------------------------------------------------------
# 13.  Summary verdict
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("VERDICT")
print("=" * 70)
if hits:
    print("POSITIVE: The following edge Laplacian invariants match K ingredients:")
    for cn, cv, tn, tv, re in hits:
        print(f"  {cn} = {cv:.10f}  matches  {tn}")
else:
    print("NEGATIVE: No edge Laplacian spectral invariant matches any K ingredient")
    print("within 1e-6 relative error.")
    if near_misses:
        print(f"\nClosest near-misses (within 5%):")
        for cn, cv, tn, tv, re in near_misses[:5]:
            print(f"  {cn} = {cv:.10f}  vs  {tn} = {tv:.10f}  (rel_err {re:.4f})")
