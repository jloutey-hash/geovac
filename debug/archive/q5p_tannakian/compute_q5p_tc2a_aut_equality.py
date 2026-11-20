r"""
debug/compute_q5p_tc2a_aut_equality.py

Sprint Q5'-Tannakian-Closure TC-2a:\ finite-cutoff reconstruction test.

Computes $\dim \mathrm{Aut}^\otimes(\omega)$ on
$\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}}(n_{\max}=2))$
directly as the dimension of the bit-exact naturality variety on the
faithful witness panel, and compares to the classical Tannakian
prediction $\dim = 3 N(2) = 15$.

This is the converse direction of TC-1e/TC-1f. TC-1e/TC-1f showed
$U^*_{\mathrm{Levi}} \subseteq \mathrm{Aut}^\otimes(\omega)$ (the
inclusion). TC-2a asks the reverse:\ at finite cutoff, is the
inclusion an equality?

Strategy
--------
1. Build the faithful witness panel:\ 15 reps $V_g$, one per primitive
   generator $g = (n, l, k)$, each 2-dim with $X_g^{V_g} = E_{12}$ and
   $X_h^{V_g} = 0$ for $h \ne g$. Plus the unit object $T = \mathbf{1}$.
2. Parameterize $\eta_{V_g}$ on each panel rep by 4 symbolic
   $\mathbb{Q}$-indeterminates $(p_g, q_g, r_g, s_g)$. Plus $\eta_T$ for
   the unit. Total:\ $15 \cdot 4 + 1 = 61$ symbols.
3. For each ordered pair $(V_a, V_b)$ from the panel (including $T$),
   compute the $\mathbb{Q}$-basis of $\mathrm{Hom}_{\mathrm{Rep}_{\mathrm{fin}}}(V_a, V_b)$
   by solving the intertwining system bit-exactly.
4. For each morphism $f \in \mathrm{Hom}(V_a, V_b)$, impose naturality
   $\eta_{V_b} f = f \eta_{V_a}$:\ four linear equations per morphism in
   the 61 symbols.
5. Impose unit normalization $\eta_T = 1$.
6. Build the constraint matrix over $\mathbb{Q}$, compute rank and
   nullity bit-exactly.
7. Report $\dim$ of the solution variety; compare to predicted 15.

If $\dim = 15$:\ Aut$^\otimes(\omega) = \mathbb{G}_a^{15}$ confirmed
bit-exactly at finite cutoff. Inclusion (TC-1e/TC-1f) $\to$ equality
on the $n_{\max}$-axis substrate. ($SL_2$ on the PW axis is a separate
layer;\ deferred to TC-2b.)

If $\dim > 15$:\ extra automorphisms detected. Identify them.

If $\dim < 15$:\ inclusion proof has a hole. Audit.

Discipline
----------
Bit-exact ``sympy.Rational`` / ``sympy.Integer`` throughout. No floats.
No PSLQ. Single thread.

References
----------
- Deligne, P.; Milne, J. S. ``Tannakian categories'' (1982) Theorem 2.11.
- Sprint Q5'-Tannakian-Closure TC-1e memo
  ``debug/sprint_q5p_tc1e_aut_inclusion_memo.md``.
- Sprint Q5'-Tannakian-Closure TC-1f memo
  ``debug/sprint_q5p_tc1f_sl2_inclusion_memo.md``.
- Paper 56 (math.OA standalone v3.74.0).
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from sympy import (
    Integer, Matrix, Symbol, Rational,
    eye as sp_eye, zeros as sp_zeros,
)

from geovac.pro_system import primitive_generators
from geovac.tannakian import FinDimRep, trivial_rep


N_MAX = 2


# ---------------------------------------------------------------------
# Witness panel
# ---------------------------------------------------------------------


def build_witness_panel(n_max):
    r"""Build the faithful witness panel:\ $V_g$ for each primitive
    generator $g$.

    $V_g$ is 2-dim, $X_g^{V_g} = E_{12}$, $X_h^{V_g} = 0$ for $h \ne g$.
    The panel is *faithful*:\ every primitive generator is activated by
    exactly one panel rep, distinguishable.

    Returns
    -------
    gens : list of (n, l, k)
        Canonical ordering of primitive generators at $n_{\max}$.
    panel : dict
        $g \mapsto V_g$ (a `FinDimRep`).
    """
    gens = primitive_generators(n_max)
    E12 = Matrix([[Integer(0), Integer(1)], [Integer(0), Integer(0)]])
    panel = {}
    for g in gens:
        V_g = FinDimRep(
            n_max=n_max, dim=2, endos={g: E12},
            label=f"V_{g[0]}_{g[1]}_{g[2]}",
        )
        panel[g] = V_g
    return gens, panel


# ---------------------------------------------------------------------
# Hom-space basis
# ---------------------------------------------------------------------


def _kron_q(A, B):
    r"""Kronecker product over $\mathbb{Q}$ in column-major lex convention
    compatible with the standard ``vec`` operator.

    For a $d_W \times d_V$ matrix $F$ with column-major vectorisation,
    $\mathrm{vec}(A F B) = (B^T \otimes A) \mathrm{vec}(F)$.
    """
    ra, ca = A.rows, A.cols
    rb, cb = B.rows, B.cols
    out = sp_zeros(ra * rb, ca * cb)
    for i in range(ra):
        for j in range(ca):
            aij = A[i, j]
            if aij == 0:
                continue
            for k in range(rb):
                for l in range(cb):
                    out[i * rb + k, j * cb + l] = aij * B[k, l]
    return out


def _vec_to_matrix(v, d_W, d_V):
    r"""Reconstruct a $d_W \times d_V$ matrix from a column-major
    ``vec`` column vector of length $d_W \cdot d_V$.
    """
    F = sp_zeros(d_W, d_V)
    for idx in range(d_W * d_V):
        i = idx % d_W
        j = idx // d_W
        F[i, j] = v[idx]
    return F


def hom_basis(V, W):
    r"""Bit-exact $\mathbb{Q}$-basis of $\mathrm{Hom}_{\mathrm{Rep}_{\mathrm{fin}}}(V, W)$.

    A morphism $f:\ V \to W$ is a $\dim(W) \times \dim(V)$ matrix
    satisfying $f X_g^V = X_g^W f$ for every primitive generator $g$.
    Vectorising column-major,
    $(X_V^T \otimes I_W - I_V \otimes X_W) \mathrm{vec}(f) = 0$
    for each $g$.

    Stack one block per generator with non-trivial action on $V$ or
    $W$, compute the right null space bit-exactly, reshape back to
    matrices.
    """
    d_V = V.dim
    d_W = W.dim
    I_V = sp_eye(d_V) if d_V > 0 else sp_zeros(0, 0)
    I_W = sp_eye(d_W) if d_W > 0 else sp_zeros(0, 0)
    relevant_gens = set(V.non_zero_endos().keys()) | set(W.non_zero_endos().keys())
    n_entries = d_W * d_V
    if n_entries == 0:
        return []
    blocks = []
    for g in relevant_gens:
        X_V = V.X(g)
        X_W = W.X(g)
        M_g = _kron_q(X_V.T, I_W) - _kron_q(I_V, X_W)
        blocks.append(M_g)
    if not blocks:
        # No active generators; every linear map is a morphism.
        # Basis = standard E_{ij} for all (i, j).
        morphisms = []
        for j in range(d_V):
            for i in range(d_W):
                E = sp_zeros(d_W, d_V)
                E[i, j] = Integer(1)
                morphisms.append(E)
        return morphisms
    M = blocks[0]
    for B in blocks[1:]:
        M = M.col_join(B)
    null_vectors = M.nullspace()
    morphisms = [_vec_to_matrix(v, d_W, d_V) for v in null_vectors]
    return morphisms


# ---------------------------------------------------------------------
# Symbolic parameterisation of $\eta$ on the witness panel
# ---------------------------------------------------------------------


def parameterize_eta(gens):
    r"""Parameterize $\eta_{V_g}$ on each panel rep by 4 indeterminates,
    and $\eta_T$ by a single indeterminate.

    Returns
    -------
    etas : dict $g \mapsto 2 \times 2$ symbolic matrix.
    eta_T : 1x1 symbolic matrix.
    all_symbols : ordered list of all 61 indeterminates.
    """
    etas = {}
    all_symbols = []
    for g in gens:
        n, l, k = g
        p = Symbol(f"p_{n}_{l}_{k}")
        q = Symbol(f"q_{n}_{l}_{k}")
        r = Symbol(f"r_{n}_{l}_{k}")
        s = Symbol(f"s_{n}_{l}_{k}")
        etas[g] = Matrix([[p, q], [r, s]])
        all_symbols.extend([p, q, r, s])
    eta_T_sym = Symbol("eta_T")
    eta_T_mat = Matrix([[eta_T_sym]])
    all_symbols.append(eta_T_sym)
    return etas, eta_T_mat, all_symbols


# ---------------------------------------------------------------------
# Naturality constraint collection
# ---------------------------------------------------------------------


def collect_naturality_constraints(gens, panel, etas, eta_T_mat, n_max):
    r"""Build all naturality constraints $\eta_b \cdot f = f \cdot \eta_a$
    for $f \in \mathrm{Hom}(V_a, V_b)$ on the panel including $T$.

    Constraints are linear in the symbolic entries of $\eta$ (matrices
    over $\mathbb{Q}[symbols]$);\ each ``diff`` matrix entry yields one
    linear equation.

    Returns a list of sympy expressions, each linear in
    `all_symbols`. Trivial (identically-zero) entries are skipped.
    """
    T = trivial_rep(n_max, dim=1)
    panel_with_T = dict(panel)
    panel_with_T['_T'] = T
    eta_with_T = dict(etas)
    eta_with_T['_T'] = eta_T_mat
    labels = list(panel.keys()) + ['_T']
    constraints = []
    n_morphisms_per_pair = {}
    for a_label in labels:
        for b_label in labels:
            V_a = panel_with_T[a_label]
            V_b = panel_with_T[b_label]
            eta_a = eta_with_T[a_label]
            eta_b = eta_with_T[b_label]
            morphs = hom_basis(V_a, V_b)
            n_morphisms_per_pair[(str(a_label), str(b_label))] = len(morphs)
            for f_mat in morphs:
                lhs = eta_b * f_mat
                rhs = f_mat * eta_a
                diff = lhs - rhs
                for i in range(diff.rows):
                    for j in range(diff.cols):
                        entry = diff[i, j].expand()
                        if entry != 0:
                            constraints.append(entry)
    return constraints, n_morphisms_per_pair


def add_unit_normalization(constraints, eta_T_mat):
    """Add the constraint $\\eta_T = 1$."""
    constraints.append(eta_T_mat[0, 0] - Integer(1))


# ---------------------------------------------------------------------
# Linear-system solver
# ---------------------------------------------------------------------


def solve_linear_system(constraints, all_symbols):
    r"""Solve $A v = b$ bit-exactly where each constraint is linear in
    ``all_symbols``.

    Returns rank, augmented-rank, nullity, and consistency flag.
    """
    from sympy import linear_eq_to_matrix
    A, b = linear_eq_to_matrix(constraints, all_symbols)
    rank_A = A.rank()
    aug = A.row_join(b)
    rank_aug = aug.rank()
    consistent = (rank_A == rank_aug)
    n_vars = len(all_symbols)
    nullity = n_vars - rank_A
    return {
        'n_constraints_kept': len(constraints),
        'A_shape': [A.rows, A.cols],
        'rank_A': rank_A,
        'rank_aug': rank_aug,
        'consistent': consistent,
        'n_vars': n_vars,
        'nullity': nullity,
        'dim_solution_variety': nullity if consistent else None,
    }


def extract_solution_structure(constraints, all_symbols, gens, etas, eta_T_mat):
    r"""Solve the system in closed form and report which symbols are free
    and which are determined.

    Uses `sympy.solve` (linear; bit-exact). Returns a dict mapping each
    symbol to its value (in terms of free symbols).
    """
    from sympy import solve
    sol = solve(constraints, all_symbols, dict=True)
    if not sol:
        return {'error': 'no solution found'}
    sol_dict = sol[0]
    # Identify free variables: symbols not appearing as keys, OR appearing as values
    # with free RHS. sympy.solve typically eliminates determined symbols, leaving
    # the free ones as the parametrising set.
    free_vars = set()
    for sym in all_symbols:
        if sym not in sol_dict:
            free_vars.add(sym)
        else:
            rhs = sol_dict[sym]
            for s in rhs.free_symbols:
                free_vars.add(s)
    free_vars_ordered = [s for s in all_symbols if s in free_vars]
    # Classify each panel rep's η structure
    rep_structure = {}
    for g in gens:
        eta_g = etas[g]
        eta_g_sub = eta_g.subs(sol_dict)
        rep_structure[str(g)] = {
            'eta_solved': [
                [str(eta_g_sub[0, 0]), str(eta_g_sub[0, 1])],
                [str(eta_g_sub[1, 0]), str(eta_g_sub[1, 1])],
            ],
        }
    eta_T_solved = eta_T_mat.subs(sol_dict)[0, 0]
    return {
        'n_free_vars': len(free_vars),
        'free_vars': [str(s) for s in free_vars_ordered],
        'rep_structure': rep_structure,
        'eta_T_solved': str(eta_T_solved),
    }


# ---------------------------------------------------------------------
# Recovery of $\Phi$ from TC-1e: verify each free $q_g$ matches
# the unipotent matrix exponential parameter
# ---------------------------------------------------------------------


def verify_phi_recovery(gens, sol_struct):
    r"""For each $V_g$, check that the solved $\eta_{V_g}$ matches
    the unipotent matrix exponential form
    $\begin{pmatrix} 1 & q_g \\ 0 & 1 \end{pmatrix}$ — which is exactly
    $\exp(q_g E_{12}) = \Phi((0, \ldots, q_g, \ldots, 0))(V_g)$ from TC-1e.
    """
    matches = []
    expected_q_per_g = {g: f"q_{g[0]}_{g[1]}_{g[2]}" for g in gens}
    for g in gens:
        eta_solved = sol_struct['rep_structure'][str(g)]['eta_solved']
        # Expected: [[1, q_g], [0, 1]] for q_g = free param
        ok = (
            eta_solved[0][0] == '1'
            and eta_solved[1][0] == '0'
            and eta_solved[1][1] == '1'
            and eta_solved[0][1] == expected_q_per_g[g]
        )
        matches.append({
            'g': list(g),
            'eta': eta_solved,
            'expected_q': expected_q_per_g[g],
            'match': ok,
        })
    all_match = all(m['match'] for m in matches)
    return {'all_match': all_match, 'per_rep': matches}


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------


def main():
    print("=" * 72)
    print("Sprint Q5'-Tannakian-Closure TC-2a")
    print(f"Finite-cutoff Aut^otimes(omega) computation at n_max = {N_MAX}")
    print("=" * 72)

    t_start = time.time()

    gens, panel = build_witness_panel(N_MAX)
    n_gens = len(gens)
    predicted_dim = n_gens  # 3 * N(n_max)
    print(f"\nPrimitive generators at n_max = {N_MAX}: {n_gens}")
    print(f"Predicted dim Aut^otimes(omega) = 3 * N({N_MAX}) = {predicted_dim}")
    print("  (classical Tannakian reconstruction for Sym(V) primitive Hopf algebra)")

    etas, eta_T_mat, all_symbols = parameterize_eta(gens)
    print(f"Symbolic indeterminates: {len(all_symbols)} = {n_gens} * 4 + 1")

    print("\nCollecting naturality constraints (this is the bulk of the work)...")
    t_constr = time.time()
    constraints, n_morphisms_per_pair = collect_naturality_constraints(
        gens, panel, etas, eta_T_mat, N_MAX
    )
    n_total_morphisms = sum(n_morphisms_per_pair.values())
    add_unit_normalization(constraints, eta_T_mat)
    print(f"  Total morphisms in panel: {n_total_morphisms}")
    print(f"  Total constraints (incl. unit): {len(constraints)}")
    print(f"  Wall time so far: {time.time() - t_constr:.2f} s")

    print("\nSolving linear system bit-exactly...")
    t_solve = time.time()
    result = solve_linear_system(constraints, all_symbols)
    print(f"  Solve time: {time.time() - t_solve:.2f} s")
    for k, v in result.items():
        print(f"  {k}: {v}")

    dim_solved = result['dim_solution_variety']
    print(f"\n  Predicted: {predicted_dim}")
    print(f"  Computed:  {dim_solved}")

    # Verdict
    if dim_solved == predicted_dim:
        verdict = "POSITIVE — equality at finite cutoff"
    elif dim_solved is None:
        verdict = "INCONSISTENT — system has no solution (bug)"
    elif dim_solved > predicted_dim:
        verdict = f"SURPRISE — found {dim_solved - predicted_dim} extra automorphisms"
    else:
        verdict = f"NEGATIVE — short by {predicted_dim - dim_solved} (audit needed)"
    print(f"\nVerdict: {verdict}")

    # Closed-form solution structure
    print("\nExtracting closed-form solution structure...")
    t_struct = time.time()
    sol_struct = extract_solution_structure(constraints, all_symbols, gens, etas, eta_T_mat)
    print(f"  Solve time: {time.time() - t_struct:.2f} s")
    if 'error' not in sol_struct:
        print(f"  Free variables: {sol_struct['n_free_vars']}")
        print(f"  eta_T solved: {sol_struct['eta_T_solved']}")

    # Phi recovery
    print("\nVerifying recovery of TC-1e Phi structure...")
    phi_check = verify_phi_recovery(gens, sol_struct)
    print(f"  All reps match unipotent exp form: {phi_check['all_match']}")
    if not phi_check['all_match']:
        for m in phi_check['per_rep']:
            if not m['match']:
                print(f"    MISMATCH at {m['g']}: eta = {m['eta']}, expected q = {m['expected_q']}")

    elapsed = time.time() - t_start
    print(f"\nTotal wall time: {elapsed:.2f} s")
    print("=" * 72)

    # Save
    out_data = {
        'sprint': 'Q5-prime-Tannakian-Closure-TC-2a',
        'n_max': N_MAX,
        'predicted_dim': predicted_dim,
        'verdict': verdict,
        **result,
        'n_total_morphisms_in_panel': n_total_morphisms,
        'n_free_vars': sol_struct.get('n_free_vars'),
        'eta_T_solved': sol_struct.get('eta_T_solved'),
        'phi_recovery_all_match': phi_check['all_match'],
        'wall_time_seconds': float(elapsed),
    }
    data_path = Path(__file__).parent / 'data' / 'sprint_q5p_tc2a_aut_equality.json'
    data_path.parent.mkdir(parents=True, exist_ok=True)
    with open(data_path, 'w') as f:
        json.dump(out_data, f, indent=2)
    print(f"\nData saved to {data_path}")

    return result, sol_struct, phi_check


if __name__ == "__main__":
    main()
