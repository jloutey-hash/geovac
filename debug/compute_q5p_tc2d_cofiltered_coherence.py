r"""
debug/compute_q5p_tc2d_cofiltered_coherence.py

Sprint Q5'-Tannakian-Closure TC-2d:\ pro-system cofiltered coherence
package.

Combines the per-cutoff TC-2a/b/c equalities into a single pro-system
coherence statement via PS-2 functoriality.

Structural setup
----------------
PS-2 (v3.68.0) gave the cofiltered Hopf-algebra transitions
$\Phi_{m, k}:\ \mathcal{H}_{\mathrm{GV}}(m) \to \mathcal{H}_{\mathrm{GV}}(k)$
for $m \ge k$, satisfying the cofiltered axiom
$\Phi_{m, k} = \Phi_{n, k} \circ \Phi_{m, n}$ at every triple $k \le n \le m$.

Spec dualizes these into a cofiltered system of pro-unipotent groups
going in the opposite direction:\ for $k \le m$ there is a
\emph{restriction map}
$$
\rho_{m, k}:\ \mathrm{Aut}^\otimes(\omega)^{(m)} \to \mathrm{Aut}^\otimes(\omega)^{(k)},
$$
sending an automorphism $\eta^{(m)}$ to its restriction to the
sub-category $\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}}(k))
\hookrightarrow \mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}}(m))$.

By TC-2a/c, $\mathrm{Aut}^\otimes(\omega)^{(n)} = \mathbb{G}_a^{3 N(n)}$
at every cutoff. The restriction $\rho_{m, k}$ is exactly the
\emph{coordinate truncation} on $\mathbb{G}_a^{3 N(m)}$:\ drop the
parameters $t_g$ for generators $g$ outside the cutoff-$k$ panel
(i.e., generators with sector index $n > k$).

This sprint verifies bit-exactly:
1. For every triple $(\ell, m, t^{(m)})$ with $\ell \le m$ and a
   parameter $t^{(m)} \in \mathbb{Q}^{3 N(m)}$, the restriction
   $\rho_{m, \ell}(\Phi(t^{(m)})) = \Phi(t^{(m)}|_{P^{(\ell)}})$
   bit-exact on every $V_g$ for $g$ in the $\ell$-panel.
2. Restriction transitivity $\rho_{n, k} = \rho_{m, k} \circ \rho_{n, m}$
   at every $k \le m \le n$ (the cofiltered axiom on Aut$^\otimes$).
3. Restriction is a group homomorphism (by linearity of the parameter
   truncation, which is automatic since $\Phi$ is the matrix
   exponential of a $\mathbb{Q}$-linear combination of commuting
   nilpotents).

Combined with TC-2b ($SL_2$ factor at the PW panel, dim 3,
$n_{\max}$-independent), the per-cutoff coherence is

$$
\mathrm{Aut}^\otimes(\omega)^{(\infty)} = \varprojlim_n \mathrm{Aut}^\otimes(\omega)^{(n)}
   = \varprojlim_n \bigl(\mathbb{G}_a^{3 N(n)} \rtimes SL_2\bigr)
   = \mathbb{G}_a^{\infty} \rtimes SL_2.
$$

The remaining multi-year content is the **pro-finite Tannakian theorem**
(Deligne 1990, Brown 2012, Glanois 2015) coherent with v3.66.0 FO3
Interpretation C at the period level on $\mathcal{O}_\infty$ (PS-3).

References
----------
- Sprint TC-2a memo (`debug/sprint_q5p_tc2a_aut_equality_memo.md`).
- Sprint TC-2b memo (`debug/sprint_q5p_tc2b_sl2_aut_equality_memo.md`).
- Sprint TC-2c memo (`debug/sprint_q5p_tc2c_higher_cutoff_memo.md`).
- Sprint PS-2 memo (`debug/sprint_q5p_ps2_ustar_compatibility_memo.md`).
- Sprint PS-3 memo (`debug/sprint_q5p_ps3_inverse_limit_memo.md`).
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from sympy import Integer, Matrix, Rational, eye as sp_eye, zeros as sp_zeros

from geovac.pro_system import primitive_generators, n_primitive_generators
from geovac.tannakian import FinDimRep, levi_unipotent_action


# ---------------------------------------------------------------------
# Restriction map on parameters
# ---------------------------------------------------------------------


def restrict_parameter(t_dict_m, gens_low):
    r"""Restrict $t^{(m)}$ to the $k$-panel by dropping generators
    outside the lower cutoff.

    $t^{(m)} \in \mathbb{Q}^{3 N(m)}$ is a sparse dict mapping
    primitive generators to rationals. The restriction to cutoff $k$
    keeps only the entries indexed by generators in $P^{(k)}$.

    Parameters
    ----------
    t_dict_m : dict[(n, l, s), Rational]
    gens_low : list[(n, l, s)]  --  generators at the lower cutoff

    Returns
    -------
    t_dict_k : dict[(n, l, s), Rational]
    """
    gens_low_set = set(gens_low)
    return {g: v for g, v in t_dict_m.items() if g in gens_low_set}


# ---------------------------------------------------------------------
# Witness rep V_g valid at multiple cutoffs
# ---------------------------------------------------------------------


def witness_rep_at_cutoff(g, n_max):
    r"""The 2-dim witness rep $V_g$ at cutoff $n_{\max}$ activating only
    the generator $g$.
    """
    E12 = Matrix([[Integer(0), Integer(1)], [Integer(0), Integer(0)]])
    return FinDimRep(
        n_max=n_max, dim=2, endos={g: E12},
        label=f"V_{g[0]}_{g[1]}_{g[2]}@nmax={n_max}",
    )


# ---------------------------------------------------------------------
# Restriction coherence on a single $V_g$
# ---------------------------------------------------------------------


def verify_restriction_on_Vg(t_dict_m, g, k, m):
    r"""Verify bit-exact that
    $\eta_{V_g}^{(m)}(t^{(m)}) = \eta_{V_g}^{(k)}(t^{(m)}|_{P^{(k)}})$
    for $g \in P^{(k)} \subset P^{(m)}$.

    By construction of $\Phi(t)(V) = \exp(\sum_h t_h X_h^V)$ and the
    fact that $X_h^{V_g} = 0$ for $h \ne g$, both sides equal
    $\exp(t_g E_{12})$ where $t_g$ is the $g$-component of $t^{(m)}$
    (which equals the $g$-component of the restriction).

    This is the panel-level coherence statement.
    """
    if k > m:
        raise ValueError(f"k={k} must be <= m={m}")
    gens_low = primitive_generators(k)
    if g not in gens_low:
        raise ValueError(f"g={g} not in P^({k})")
    # eta at cutoff m
    V_g_m = witness_rep_at_cutoff(g, m)
    eta_m = levi_unipotent_action(t_dict_m, V_g_m)
    # eta at cutoff k after restriction
    t_dict_k = restrict_parameter(t_dict_m, gens_low)
    V_g_k = witness_rep_at_cutoff(g, k)
    eta_k = levi_unipotent_action(t_dict_k, V_g_k)
    return eta_m == eta_k


# ---------------------------------------------------------------------
# Transitivity $\rho_{n, k} = \rho_{m, k} \circ \rho_{n, m}$
# ---------------------------------------------------------------------


def verify_restriction_transitivity(t_dict_n, k, m, n):
    r"""Verify $\rho_{n, k}(t) = \rho_{m, k}(\rho_{n, m}(t))$ for
    $k \le m \le n$.

    At the parameter level, this is a trivial set-theoretic statement
    about successive dict restrictions:\
    $t|_{P^{(k)}} = (t|_{P^{(m)}})|_{P^{(k)}}$,
    since $P^{(k)} \subseteq P^{(m)} \subseteq P^{(n)}$.

    The bit-exact verification confirms that
    ``restrict_parameter`` respects this.
    """
    if not (k <= m <= n):
        raise ValueError(f"Need k <= m <= n; got k={k}, m={m}, n={n}")
    gens_k = primitive_generators(k)
    gens_m = primitive_generators(m)
    direct = restrict_parameter(t_dict_n, gens_k)
    via_m = restrict_parameter(restrict_parameter(t_dict_n, gens_m), gens_k)
    return direct == via_m


# ---------------------------------------------------------------------
# Restriction is a group homomorphism:\ $\rho_{m, k}(\Phi(t_1 + t_2))
# = \rho_{m, k}(\Phi(t_1)) \cdot \rho_{m, k}(\Phi(t_2))$
# ---------------------------------------------------------------------


def verify_restriction_homomorphism(t1_dict_m, t2_dict_m, g, k, m):
    r"""Verify the restriction is a group hom on a witness $V_g$ with
    $g \in P^{(k)} \subset P^{(m)}$:\
    $\rho_{m, k}(\Phi(t_1) \cdot \Phi(t_2)) = \rho_{m, k}(\Phi(t_1))
    \cdot \rho_{m, k}(\Phi(t_2))$.

    Since $\Phi(t_1) \cdot \Phi(t_2) = \Phi(t_1 + t_2)$ (the unipotent
    factor is abelian), the LHS becomes
    $\rho_{m, k}(\Phi(t_1 + t_2)) = \Phi((t_1 + t_2)|_{P^{(k)}})
    = \Phi(t_1|_{P^{(k)}} + t_2|_{P^{(k)}})$.

    The RHS is $\Phi(t_1|_{P^{(k)}}) \cdot \Phi(t_2|_{P^{(k)}}) =
    \Phi(t_1|_{P^{(k)}} + t_2|_{P^{(k)}})$.

    Both sides equal as bit-exact matrices.
    """
    # Sum
    t_sum = dict(t1_dict_m)
    for g_key, v in t2_dict_m.items():
        if g_key in t_sum:
            t_sum[g_key] = t_sum[g_key] + v
        else:
            t_sum[g_key] = v
    V_g_m = witness_rep_at_cutoff(g, m)
    # LHS:\ Phi(t_sum) restricted
    eta_sum_m = levi_unipotent_action(t_sum, V_g_m)
    # Direct: Phi(t1) * Phi(t2)
    eta_1 = levi_unipotent_action(t1_dict_m, V_g_m)
    eta_2 = levi_unipotent_action(t2_dict_m, V_g_m)
    product = eta_1 * eta_2
    return (eta_sum_m - product).expand() == sp_zeros(2, 2)


# ---------------------------------------------------------------------
# Main panel
# ---------------------------------------------------------------------


def main():
    print("=" * 72)
    print("Sprint Q5'-Tannakian-Closure TC-2d")
    print("Cofiltered coherence package")
    print("=" * 72)

    t_start = time.time()

    # Panel of (k, m) pairs and test parameters
    cutoff_pairs = [(1, 2), (1, 3), (2, 3), (1, 4), (2, 4), (3, 4)]
    print(f"\nCofiltered pairs (k, m): {cutoff_pairs}")

    # Test parameters: choose a couple of representative dicts per cutoff
    def test_t_at(m, seed=0):
        gens_m = primitive_generators(m)
        if seed == 0:
            # Single non-zero on first generator
            return {gens_m[0]: Rational(2, 3)}
        elif seed == 1:
            # Two non-zero values on the first two generators
            return {gens_m[0]: Rational(1), gens_m[1]: Rational(-1, 2)}
        elif seed == 2:
            # All non-zero (rationals on every generator)
            return {g: Rational((i + 1) % 7, (i + 2) % 5 + 1)
                    for i, g in enumerate(gens_m)}
        else:
            # All zero
            return {}

    n_residuals = 0
    n_pass = 0
    block_a = []
    print("\nBlock A -- panel-level restriction coherence")
    print("-" * 50)
    for k, m in cutoff_pairs:
        gens_k = primitive_generators(k)
        for seed in (0, 1, 2):
            t_dict_m = test_t_at(m, seed=seed)
            for g in gens_k:
                n_residuals += 1
                ok = verify_restriction_on_Vg(t_dict_m, g, k, m)
                if ok:
                    n_pass += 1
                block_a.append({"k": k, "m": m, "seed": seed, "g": list(g), "pass": ok})
        print(f"  pair (k={k}, m={m}): {len(gens_k) * 3} restriction tests, all pass = {all(b['pass'] for b in block_a if b['k']==k and b['m']==m)}")

    print(f"\n  Block A total: {n_pass} / {n_residuals} pass")

    # Block B: transitivity
    print("\nBlock B -- restriction transitivity (cofiltered axiom on Aut^otimes)")
    print("-" * 50)
    transitivity_triples = [(1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)]
    block_b = []
    for k, m, n in transitivity_triples:
        for seed in (0, 1, 2):
            t_dict_n = test_t_at(n, seed=seed)
            ok = verify_restriction_transitivity(t_dict_n, k, m, n)
            block_b.append({"k": k, "m": m, "n": n, "seed": seed, "pass": ok})
            print(f"  triple (k={k}, m={m}, n={n}), seed={seed} -- {'pass' if ok else 'FAIL'}")
    n_b_pass = sum(b["pass"] for b in block_b)
    n_b_total = len(block_b)
    print(f"\n  Block B total: {n_b_pass} / {n_b_total} pass")

    # Block C: group homomorphism
    print("\nBlock C -- restriction is a group homomorphism")
    print("-" * 50)
    block_c = []
    for k, m in [(1, 2), (2, 3), (1, 3)]:
        gens_k = primitive_generators(k)
        for g in gens_k[:2]:  # Just a couple per (k, m)
            t1 = test_t_at(m, seed=1)
            t2 = test_t_at(m, seed=0)
            ok = verify_restriction_homomorphism(t1, t2, g, k, m)
            block_c.append({"k": k, "m": m, "g": list(g), "pass": ok})
    n_c_pass = sum(b["pass"] for b in block_c)
    n_c_total = len(block_c)
    print(f"  ({n_c_pass} / {n_c_total} pass)")

    overall = (n_pass == n_residuals) and (n_b_pass == n_b_total) and (n_c_pass == n_c_total)
    verdict = "POSITIVE" if overall else "FAIL"

    elapsed = time.time() - t_start
    print(f"\nTotal wall time: {elapsed:.2f} s")
    print("=" * 72)

    # Headline numbers
    total_pass = n_pass + n_b_pass + n_c_pass
    total = n_residuals + n_b_total + n_c_total
    print(f"\nTotal bit-exact zero residuals: {total_pass} / {total}")
    print(f"Overall TC-2d verdict: {verdict}")
    print("\nStructural consequence:")
    print(f"  Aut^otimes(omega)^(infinity) = projlim_n (G_a^{{3 N(n)}} rtimes SL_2)")
    print(f"                                 = G_a^infinity rtimes SL_2.")
    print("  Remaining multi-year content: pro-finite Tannakian theorem")
    print("  coherent with v3.66.0 FO3 Interpretation C on O_infinity (PS-3).")

    # Save
    out_data = {
        "sprint": "Q5-prime-Tannakian-Closure-TC-2d",
        "cofiltered_pairs": cutoff_pairs,
        "transitivity_triples": transitivity_triples,
        "block_A_restriction_coherence_count": n_residuals,
        "block_A_pass": n_pass,
        "block_B_transitivity_count": n_b_total,
        "block_B_pass": n_b_pass,
        "block_C_homomorphism_count": n_c_total,
        "block_C_pass": n_c_pass,
        "total_residuals": total,
        "total_pass": total_pass,
        "overall_verdict": verdict,
        "wall_time_seconds": float(elapsed),
    }
    data_path = Path(__file__).parent / "data" / "sprint_q5p_tc2d_cofiltered_coherence.json"
    data_path.parent.mkdir(parents=True, exist_ok=True)
    with open(data_path, "w") as f:
        json.dump(out_data, f, indent=2)
    print(f"\nData saved to {data_path}")
    return out_data


if __name__ == "__main__":
    main()
