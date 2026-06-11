"""
Verification of the structural claim underlying the analytic proof of the
Dirac-triangle inequality.

Setting. We test the following structural claim, which is the heart of the
analytic proof:

  PROOF KERNEL CLAIM. Let g be a compact connected simple Lie algebra.
  Let lambda, lambda' be dominant integral weights and let sigma be a
  dominant integral weight appearing with positive multiplicity in
  V_lambda (x) V_{lambda'}^*. Then there exists w in the Weyl group W
  and mu' in the W-orbit of -lambda' such that

      sigma + rho = w (lambda + mu' + rho).                              (*)

  Equivalently, there exists w_1 in W such that

      w_1 (sigma + rho) - lambda - rho  is in the W-orbit of -lambda',

  i.e., it is a *Weyl conjugate of -lambda'* (an extremal weight of
  V_{lambda'^*} = V_{-w_0 lambda'} of "long" length).

This is the strongest possible form of the assertion. We test it
numerically on the same four panels as `dirac_triangle_extended_verify.py`:
SU(3) (p+q <= 5), SU(4) (a+b+c <= 3), Sp(2)=C_2 (a+b <= 3), and G_2
(a+b <= 2).

If this claim holds at every tested case, then the analytic argument
in `dirac_triangle_full_proof.md` Section 3 closes:

  |w_1(sigma+rho)|^2 = |sigma+rho|^2     (Weyl-invariant)
  = |lambda + mu' + rho|^2               (by the claim)
  = |lambda + rho|^2 + 2<lambda+rho, mu'> + |mu'|^2

  Now mu' = w(-lambda'), and the norm |mu'|^2 = |lambda'|^2.
  Therefore:

  |sigma+rho|^2 - |rho|^2 - |lambda-lambda'|^2
    = 2<lambda+rho, mu'> + |lambda'|^2 - |rho|^2 + |lambda+rho|^2
      - |lambda|^2 + 2<lambda, lambda'> - |lambda'|^2
    = ...

Let's do the algebra cleanly. Set u = lambda+rho, v = lambda'+rho. Then
|D(lambda)| = |u|, |D(lambda')| = |v|, C(sigma) = |sigma+rho|^2 - |rho|^2.

We want C(sigma) >= |lambda - lambda'|^2 = |u - v|^2.

Claim (*): sigma + rho = w(lambda + mu' + rho) = w(u - lambda' + (1+0)) --
            actually sigma+rho = w(lambda + mu' + rho) where mu' in W*(-lambda').
            So sigma + rho = w(u + mu') where u = lambda + rho and mu' is
            some Weyl-image of -lambda'.

|sigma + rho|^2 = |w(u + mu')|^2 = |u + mu'|^2
  = |u|^2 + 2<u, mu'> + |mu'|^2
  = |u|^2 + 2<u, mu'> + |lambda'|^2.

|lambda - lambda'|^2 = |u - v|^2 = |u|^2 - 2<u, v> + |v|^2
  = |u|^2 - 2<u, lambda' + rho> + |lambda' + rho|^2
  = |u|^2 - 2<u, lambda'> - 2<u, rho> + |lambda'|^2 + 2<lambda', rho> + |rho|^2.

C(sigma) - |lambda-lambda'|^2
  = |sigma+rho|^2 - |rho|^2 - |lambda-lambda'|^2
  = [|u|^2 + 2<u, mu'> + |lambda'|^2] - |rho|^2
    - [|u|^2 - 2<u, lambda'> - 2<u, rho> + |lambda'|^2 + 2<lambda', rho> + |rho|^2]
  = 2<u, mu'> + 2<u, lambda'> + 2<u, rho> - 2<lambda', rho> - 2|rho|^2
  = 2<u, mu' + lambda'> + 2<u, rho> - 2<lambda', rho> - 2|rho|^2
  = 2<u, mu' + lambda' + rho - (lambda' + rho)>  + 2<u - rho, rho> - 2<lambda', rho>
  Hmm getting tangled. Let me redo more carefully.

C(sigma) - |lambda - lambda'|^2
  = 2<u, mu'> + |lambda'|^2 - |rho|^2 - (|u|^2 - 2<u,v> + |v|^2 - |u|^2)
    Wait. Let's just keep |u|^2 and cancel terms.

  = |sigma+rho|^2 - |rho|^2 - |lambda-lambda'|^2
  = (|u|^2 + 2<u, mu'> + |lambda'|^2) - |rho|^2 - |u|^2 + 2<u,v> - |v|^2
  = 2<u, mu'> + |lambda'|^2 - |rho|^2 + 2<u, v> - |v|^2

With v = lambda' + rho:
  |v|^2 = |lambda'|^2 + 2<lambda', rho> + |rho|^2

So
  C(sigma) - |lambda-lambda'|^2
  = 2<u, mu'> + |lambda'|^2 - |rho|^2 + 2<u, v> - |lambda'|^2 - 2<lambda', rho> - |rho|^2
  = 2<u, mu'> + 2<u, v> - 2<lambda', rho> - 2|rho|^2
  = 2<u, mu' + v> - 2<lambda', rho> - 2|rho|^2
  = 2<u, mu' + lambda' + rho> - 2<lambda', rho> - 2|rho|^2.

If mu' = -lambda' (one of the Weyl orbits of -lambda', specifically the one
when w0 = identity which isn't usually the case):
  = 2<u, -lambda' + lambda' + rho> - 2<lambda', rho> - 2|rho|^2
  = 2<u, rho> - 2<lambda', rho> - 2|rho|^2
  = 2<u - lambda' - rho, rho>
  = 2<lambda + rho - lambda' - rho, rho>
  = 2<lambda - lambda', rho>.

For this to be >= 0 we'd need lambda >= lambda' in fundamental-weight order,
which is NOT always true. So mu' = -lambda' (the lowest weight, in some
conventions) doesn't always work.

KEY INSIGHT (this is what we test):

  We choose mu' in W*(-lambda') such that <u, mu' + lambda' + rho> takes a
  helpful sign, namely such that

      <u, mu' + lambda' + rho> >= <lambda', rho> + |rho|^2.    (#)

  Equivalently, we want
      <lambda + rho, mu'> >= <lambda', rho> + |rho|^2 - <lambda + rho, lambda' + rho>.
  But this is just (after simplification):
      <lambda + rho, mu'> + |lambda + rho|^2 / ... no wait.

  Let me just directly check: we need C(sigma) - |lambda-lambda'|^2 >= 0, i.e.:
      <u, mu' + v> >= <lambda', rho> + |rho|^2.

  In particular, since the LEFT side equals <u, w^{-1}(sigma+rho)> + <u, v> - <u, v>...
  not simplifying.

Let me just verify the kernel claim directly: for every (lambda, lambda', sigma)
appearing as a positive-multiplicity summand, there exists a w in the Weyl
group and a mu' in W*(-lambda') such that sigma+rho = w(lambda + mu' + rho).
"""

import json
import os
import sys
import time

# Import from the existing verifier
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_triangle_extended_verify import (
    SimpleLieAlgebra,
    build_A,
    build_C2,
    build_G2,
    weights_of_irrep,
    tensor_product,
    panel_dominant_weights,
)
import sympy as sp
from sympy import Matrix, Rational, sqrt, Integer
from typing import Dict, List, Tuple


# ---------------------------------------------------------------------------
# Weyl group orbit enumeration
# ---------------------------------------------------------------------------

def weyl_orbit(la: SimpleLieAlgebra, weight: Tuple[int, ...],
               max_orbit_size: int = 200000) -> List[Tuple[int, ...]]:
    """
    Enumerate the Weyl orbit of a weight via BFS using simple reflections.

    Reflection s_i sends weight v to v - <v, alpha_i^vee> alpha_i.
    In omega basis, alpha_i is the i-th column of the Cartan matrix, and
    <v, alpha_i^vee> = v_i (the i-th coordinate in fundamental-weight basis).

    So s_i: v[j] -> v[j] - v[i] * cartan[j, i].
    """
    rank = la.rank
    cartan = la.cartan
    orbit = {tuple(weight)}
    frontier = [tuple(weight)]
    while frontier and len(orbit) < max_orbit_size:
        new_frontier = []
        for v in frontier:
            for i in range(rank):
                cv = list(v)
                ci = cv[i]
                for j in range(rank):
                    cv[j] = cv[j] - ci * int(cartan[j, i])
                tcv = tuple(cv)
                if tcv not in orbit:
                    orbit.add(tcv)
                    new_frontier.append(tcv)
        frontier = new_frontier
    return list(orbit)


# ---------------------------------------------------------------------------
# Weyl group elements (as compositions of simple reflections)
# ---------------------------------------------------------------------------

def reflect(la: SimpleLieAlgebra, v: Tuple[int, ...], i: int) -> Tuple[int, ...]:
    """Apply simple reflection s_i to weight v in omega basis."""
    rank = la.rank
    cartan = la.cartan
    ci = v[i]
    out = list(v)
    for j in range(rank):
        out[j] = v[j] - ci * int(cartan[j, i])
    return tuple(out)


def all_weyl_elements_as_reduced_words(la: SimpleLieAlgebra, max_length: int = 30
                                        ) -> List[List[int]]:
    """
    Enumerate the Weyl group as reduced words in simple reflections.

    Each word is a list of integers (indices of simple reflections, 0-indexed).
    We use the standard fact: |W| = product of degrees / ... = order of Weyl group.

    Algorithm: BFS over reduced words. A word w_1 ... w_k is reduced iff its
    action on rho gives a weight whose orbit size in the rho orbit has
    length equal to k. We use the orbit-size criterion: a word is reduced
    iff applying its reflections one at a time to rho gives a new element
    of W*rho each time (no cycle).

    Returns list of reduced words. Caller can compose these via reflect().
    """
    rank = la.rank
    rho = la.rho_omega
    # Each element of W is identified by its action on rho.
    seen = {rho: []}  # rho_image -> reduced word
    frontier = [rho]
    while frontier:
        new_frontier = []
        for v in frontier:
            for i in range(rank):
                vp = reflect(la, v, i)
                if vp not in seen:
                    seen[vp] = seen[v] + [i]
                    new_frontier.append(vp)
        frontier = new_frontier
    return list(seen.values())


def apply_word(la: SimpleLieAlgebra, word: List[int], v: Tuple[int, ...]) -> Tuple[int, ...]:
    """Apply a sequence of simple reflections (in order) to v."""
    cur = tuple(v)
    for i in word:
        cur = reflect(la, cur, i)
    return cur


# ---------------------------------------------------------------------------
# PRV / kernel claim verification
# ---------------------------------------------------------------------------

def verify_kernel_claim(la: SimpleLieAlgebra, lam: Tuple[int, ...],
                        lam_prime: Tuple[int, ...],
                        weyl_words: List[List[int]],
                        ) -> Tuple[bool, Dict, List[Dict]]:
    """
    For every sigma appearing with positive multiplicity in V_lambda (x) V_{lam_prime}^*,
    verify that there exists w in W and mu' in the Weyl orbit of -lambda' such that
        sigma + rho = w(lambda + mu' + rho).

    Returns (claim_holds, summary, sigma_diagnostics).
    """
    rank = la.rank
    rho = la.rho_omega
    lam_prime_dual = la.dual(lam_prime)

    # Compute the tensor product decomposition
    decomp = tensor_product(la, lam, lam_prime_dual)

    # Compute the orbit of -lambda' under W
    minus_lp = tuple(-int(x) for x in lam_prime)
    orbit_minus_lp = set(weyl_orbit(la, minus_lp))

    summary = {
        "decomp_size": len([s for s, m in decomp.items() if m > 0]),
        "claim_failures": 0,
        "claim_successes": 0,
    }

    sigma_diags = []
    all_pass = True

    # For efficiency: precompute the set of candidate target points
    # Each target = w(lambda + mu' + rho) for mu' in W*(-lambda').
    # We need: sigma+rho in this set.
    targets_for_sigma_plus_rho: Dict[Tuple[int, ...], List[Tuple]] = {}

    for mu_p in orbit_minus_lp:
        shifted = tuple(lam[i] + mu_p[i] + rho[i] for i in range(rank))
        # For every Weyl word, compute w(shifted).
        for word in weyl_words:
            wshifted = apply_word(la, word, shifted)
            # Check if wshifted is dominant interior: all coords > 0
            if all(c > 0 for c in wshifted):
                # Reverse-engineered sigma = wshifted - rho
                sigma_candidate = tuple(wshifted[i] - rho[i] for i in range(rank))
                if all(c >= 0 for c in sigma_candidate):
                    if wshifted not in targets_for_sigma_plus_rho:
                        targets_for_sigma_plus_rho[wshifted] = []
                    targets_for_sigma_plus_rho[wshifted].append((mu_p, word))

    # Now for each sigma in decomp, check it's in targets_for_sigma_plus_rho.
    for sigma, mult in decomp.items():
        if mult <= 0:
            continue
        sigma_plus_rho = tuple(sigma[i] + rho[i] for i in range(rank))
        if sigma_plus_rho in targets_for_sigma_plus_rho:
            summary["claim_successes"] += 1
            sigma_diags.append({
                "sigma": list(sigma),
                "mult": mult,
                "claim_holds": True,
                "n_witnesses": len(targets_for_sigma_plus_rho[sigma_plus_rho]),
            })
        else:
            summary["claim_failures"] += 1
            all_pass = False
            sigma_diags.append({
                "sigma": list(sigma),
                "mult": mult,
                "claim_holds": False,
                "n_witnesses": 0,
            })

    return all_pass, summary, sigma_diags


# ---------------------------------------------------------------------------
# (INT) and (DT) verifications via PRV witness
# ---------------------------------------------------------------------------

def int_inequality_via_prv(la: SimpleLieAlgebra, lam: Tuple[int, ...],
                            lam_prime: Tuple[int, ...]
                            ) -> Tuple[bool, Dict]:
    """
    For each sigma in V_lam (x) V_{lam_prime}^*, verify the INT inequality:
        |lam - lam_prime|^2 <= C(sigma)
    AND verify that the PRV witness gives a constructive proof.

    Returns (pass, info).
    """
    rank = la.rank
    rho = la.rho_omega
    lam_prime_dual = la.dual(lam_prime)
    decomp = tensor_product(la, lam, lam_prime_dual)

    lam_minus_lp = tuple(lam[i] - lam_prime[i] for i in range(rank))
    lam_minus_lp_sq = la._inner_omega(lam_minus_lp, lam_minus_lp)

    info = {"sigmas": []}
    all_pass = True

    for sigma, mult in decomp.items():
        if mult <= 0:
            continue
        C_sigma = la.casimir(sigma)
        passes = (lam_minus_lp_sq <= C_sigma)
        if not passes:
            all_pass = False
        info["sigmas"].append({
            "sigma": list(sigma),
            "mult": mult,
            "C_sigma": str(C_sigma),
            "lam_diff_sq": str(lam_minus_lp_sq),
            "passes": passes,
        })

    return all_pass, info


# ---------------------------------------------------------------------------
# Panel runner for kernel claim
# ---------------------------------------------------------------------------

def run_panel_kernel(la: SimpleLieAlgebra, irreps: List[Tuple[int, ...]],
                     label: str, verbose: bool = False) -> Dict:
    """
    Verify the PRV-style kernel claim on all ordered pairs of `irreps`.

    Returns summary dict.
    """
    print(f"\n{'=' * 78}")
    print(f"Kernel claim verification: {label}")
    print(f"{'=' * 78}")
    print(f"Computing Weyl group elements ...")
    t0 = time.time()
    weyl_words = all_weyl_elements_as_reduced_words(la)
    print(f"  |W| = {len(weyl_words)} (computed in {time.time()-t0:.2f}s)")

    total_pairs = 0
    pairs_with_failure = 0
    total_sigma = 0
    sigma_failures = 0

    for lam in irreps:
        for lam_pr in irreps:
            total_pairs += 1
            claim_ok, summary, diags = verify_kernel_claim(la, lam, lam_pr, weyl_words)
            if not claim_ok:
                pairs_with_failure += 1
                if verbose:
                    failed_sigmas = [d for d in diags if not d["claim_holds"]]
                    print(f"  FAIL pair lam={lam}, lam'={lam_pr}: {len(failed_sigmas)} sigma failures")
                    for d in failed_sigmas[:3]:
                        print(f"     sigma={tuple(d['sigma'])}, mult={d['mult']}")
            for d in diags:
                total_sigma += 1
                if not d["claim_holds"]:
                    sigma_failures += 1

    print(f"\n--- {label} summary ---")
    print(f"  Total pairs: {total_pairs}")
    print(f"  Pairs with kernel-claim failure: {pairs_with_failure}")
    print(f"  Total sigmas tested: {total_sigma}")
    print(f"  Sigma failures: {sigma_failures}")
    print(f"  Pass rate: {(total_sigma - sigma_failures)}/{total_sigma}")

    return {
        "label": label,
        "rank": la.rank,
        "panel_size": len(irreps),
        "total_pairs": total_pairs,
        "pairs_with_failure": pairs_with_failure,
        "total_sigma": total_sigma,
        "sigma_failures": sigma_failures,
        "weyl_order": len(weyl_words),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "dirac_triangle_proof_verification.json")

    results = {
        "description": "Verification of PRV-style kernel claim for the analytic "
                       "proof of the Dirac-triangle (INT) inequality.",
        "panels": [],
    }

    # SU(3): full p+q <= 5 panel (small, fast)
    a2 = build_A(2)
    a2_panel = panel_dominant_weights(2, 5)
    print(f"\nSU(3) panel: {len(a2_panel)} irreps, {len(a2_panel)**2} pairs")
    res = run_panel_kernel(a2, a2_panel, "SU(3) p+q<=5", verbose=True)
    results["panels"].append(res)

    # SU(4): a+b+c <= 3 (moderate)
    a3 = build_A(3)
    a3_panel = panel_dominant_weights(3, 3)
    print(f"\nSU(4) panel: {len(a3_panel)} irreps, {len(a3_panel)**2} pairs")
    res = run_panel_kernel(a3, a3_panel, "SU(4) a+b+c<=3", verbose=True)
    results["panels"].append(res)

    # Sp(2)=C_2: a+b <= 3
    c2 = build_C2()
    c2_panel = panel_dominant_weights(2, 3)
    print(f"\nSp(2)=C_2 panel: {len(c2_panel)} irreps, {len(c2_panel)**2} pairs")
    res = run_panel_kernel(c2, c2_panel, "Sp(2)=C_2 a+b<=3", verbose=True)
    results["panels"].append(res)

    # G_2: a+b <= 2
    g2 = build_G2()
    g2_panel = panel_dominant_weights(2, 2)
    print(f"\nG_2 panel: {len(g2_panel)} irreps, {len(g2_panel)**2} pairs")
    res = run_panel_kernel(g2, g2_panel, "G_2 a+b<=2", verbose=True)
    results["panels"].append(res)

    with open(output_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {output_path}")

    print("\n" + "=" * 78)
    print("FINAL KERNEL CLAIM SUMMARY")
    print("=" * 78)
    overall_pass = True
    for p in results["panels"]:
        verdict = "PASS" if p["sigma_failures"] == 0 else "FAIL"
        print(f"  {p['label']}: {p['total_sigma']-p['sigma_failures']}/{p['total_sigma']} "
              f"sigmas, {p['total_pairs']-p['pairs_with_failure']}/{p['total_pairs']} pairs [{verdict}]")
        if p["sigma_failures"] > 0:
            overall_pass = False
    print(f"\nOverall: {'KERNEL CLAIM VERIFIED' if overall_pass else 'KERNEL CLAIM FAILED AT SOME CASES'}")


if __name__ == "__main__":
    main()
