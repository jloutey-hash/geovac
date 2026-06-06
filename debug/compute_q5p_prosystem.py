r"""
Sprint Q5'-Stage1-Prosystem driver — pro-system functoriality of the
JLO HP^even and CM-eta cocycle classes across truncations
$\mathcal{T}_{n_{\max}}$ for $n_{\max} \in \{1, 2, 3, 4\}$ on the
truncated Camporesi--Higuchi spectral triple.

Goal
----
Close the "pro-system functoriality" named follow-on from Sprint
Q5'-Stage1-Followon (v3.59.0). Build truncation
$P_{n+1 \to n}: \mathcal{O}_{n+1} \to \mathcal{O}_n$ and Berezin
$B_{n \to n+1}: \mathcal{O}_n \to \mathcal{O}_{n+1}$ maps between the
truncated spectral triples, then verify that the JLO HP^even cocycle
class and CM-eta residue class form a coherent pro-system under these
maps.

Setup
-----
Sectors at cutoff n_max are (n, l) with 1 <= n <= n_max, 0 <= l <= n.
Number of sectors: N(n_max) = n_max*(n_max + 3)/2.
- n_max=1: 2 sectors {(1,0),(1,1)}
- n_max=2: 5 sectors (existing Sub-Sprint 2c baseline)
- n_max=3: 9 sectors
- n_max=4: 14 sectors

Sectors are NESTED: sectors at cutoff n+1 = sectors at cutoff n union
{(n+1, l) : 0 <= l <= n+1}. The increment is N(n+1) - N(n) = n+2.

Truncation map P_{n+1 -> n}
---------------------------
Algebra-level: P sends a sector idempotent at level n+1 to the same
sector idempotent at level n if that sector exists at level n, zero
otherwise. This is an algebra homomorphism C^{N(n+1)} -> C^{N(n)}.

On the class level (Q^{N(n+1)} -> Q^{N(n)}), the dual P^* is sector-
projection: drop the coordinates corresponding to sectors (n+1, l)
that appeared at the higher cutoff.

Berezin map B_{n -> n+1}
------------------------
The Berezin map goes FORWARD (coarser -> finer). On the algebra level,
B sends e_s^{(n)} to e_s^{(n+1)} for sectors s present at both levels
(sector inclusion). New sectors (n+1, l) at level n+1 have no preimage,
so B is NOT surjective (consistent with finer truncation having more
degrees of freedom).

Pull-back of cocycle classes
----------------------------
Given a cocycle class psi^{(n+1)} at level n+1 (as a Q^{N(n+1)} vector
of per-sector values), the pull-back via P^* is its projection to the
N(n)-dim subspace corresponding to old sectors. The pull-back identity
to verify:
    P^*_{n+1 -> n} [psi^{(n+1)}] = [psi^{(n)}]
sector-by-sector, where [psi^{(n)}] is the cocycle class computed
INTRINSICALLY at cutoff n.

For JLO HP^even: per-sector value is chi_s = Tr(gamma e_s).
For CM-eta: per-sector value is eta_s = Tr(gamma D e_s) = Tr(gamma Lambda e_s).

Both are LOCAL quantities (sector-supported), and chi_s, eta_s depend
ONLY on the local Dirac structure at sector s — which is sector-
invariant across n_max (since adding a higher shell does not modify
the Dirac structure of lower shells). Hence we EXPECT bit-exact pull-
back identity, with sector-resolved polynomial closed forms in n_max
being TRIVIAL (each sector value is n_max-INDEPENDENT).

Closed-form sector evolution
----------------------------
JLO per-sector: chi_{(n,l)} = dim(gamma_+ e_{(n,l)}) - dim(gamma_- e_{(n,l)}).
CM-eta per-sector: eta_{(n,l)} = dim_{(n,l)} * (n + 1/2).

Both depend ONLY on the sector label (n, l), not on the cutoff n_max.
This is the key structural finding: the cocycle classes are sector-
LOCAL and therefore form a strictly compatible pro-system.

Polynomial closed forms for total invariants
--------------------------------------------
Summing per-sector values across all sectors at cutoff n_max gives:
  M_1(n_max) = dim H = sum_{n,l} dim_{(n,l)} = polynomial cubic in n_max.
  Sum of JLO class = sum_{n,l} chi_{(n,l)} = 0 (McKean-Singer, all n_max).
  Sum of CM-eta class = M_3(n_max) = quartic polynomial in n_max.

The McKean-Singer sum is constant (zero) across cutoffs because each
shell contributes a chirality-balanced pair. The CM-eta sum is the
polynomial growth M_3(n) = n(n+1)^2(n+2)/2 from Sub-Sprint 2a.

Pull-back validation across n_max in {1, 2, 3, 4}
--------------------------------------------------
For each pair (n, n+1) in {(1,2), (2,3), (3,4)}, verify:
1. Sector inclusion: sectors at n_max = n+1 contain all sectors at n_max = n.
2. JLO class pull-back: chi_{(n,l)}^{(n+1)} = chi_{(n,l)}^{(n)} bit-exact.
3. CM-eta class pull-back: eta_{(n,l)}^{(n+1)} = eta_{(n,l)}^{(n)} bit-exact.
4. New-sector class values match the expected polynomial extrapolation.

Output
------
- debug/data/sprint_q5p_prosystem.json — exact rational data
- prints to stdout, summary at end.

Discipline
----------
Bit-exact sympy.Rational throughout. No floats. No PSLQ.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

import sympy as sp
from sympy import Integer, Matrix, Rational, eye as sp_eye, zeros as sp_zeros

from geovac.spectral_triple import FockSpectralTriple


# =====================================================================
# Helpers
# =====================================================================


def sector_idempotent(st: FockSpectralTriple, sector_idx: int) -> Matrix:
    """Return the sector idempotent e_s as dim_H x dim_H diagonal."""
    N = st.dim_H
    M = sp_zeros(N, N)
    for i in range(N):
        if st._state_to_sector[i] == sector_idx:
            M[i, i] = Integer(1)
    return M


def N_sectors_closed_form(n_max: int) -> int:
    """Closed form for the number of (n, l) sectors at cutoff n_max."""
    return n_max * (n_max + 3) // 2


# =====================================================================
# Per-sector cocycle class values at cutoff n_max
# =====================================================================


def jlo_class_per_sector(st: FockSpectralTriple) -> Dict[Tuple[int, int], int]:
    """JLO HP^even class leading vector at this cutoff.

    chi_s = Tr(gamma e_s)
    """
    gamma = st.grading
    out = {}
    for s_idx, (n, l) in enumerate(st.sectors):
        e_s = sector_idempotent(st, s_idx)
        chi_s = (gamma * e_s).trace()
        out[(n, l)] = int(chi_s)
    return out


def cm_eta_class_per_sector(st: FockSpectralTriple) -> Dict[Tuple[int, int], int]:
    """CM-eta residue class leading vector at this cutoff.

    eta_s = Tr(gamma D e_s) -- on the Lambda-only part (off-diagonal kappa*A
    contribution to Tr(gamma A e_s) vanishes per-sector because e_s is
    diagonal idempotent and A has zero diagonal).
    """
    gamma = st.grading
    D = st.dirac_operator
    out = {}
    for s_idx, (n, l) in enumerate(st.sectors):
        e_s = sector_idempotent(st, s_idx)
        eta_s = (gamma * D * e_s).trace()
        out[(n, l)] = int(eta_s)
    return out


def sector_dim_per_sector(st: FockSpectralTriple) -> Dict[Tuple[int, int], int]:
    """Dimension of each sector at this cutoff."""
    out = {}
    for s_idx, (n, l) in enumerate(st.sectors):
        e_s = sector_idempotent(st, s_idx)
        out[(n, l)] = int(e_s.trace())
    return out


# =====================================================================
# Pull-back verification
# =====================================================================


def verify_pullback_jlo(
    classes_by_nmax: Dict[int, Dict[Tuple[int, int], int]],
) -> Dict:
    """Verify P^*_{n+1 -> n} [chi^{(n+1)}] = [chi^{(n)}] sector-by-sector.

    Returns a dict with per-pair (n, n+1) bit-exact equality checks.
    """
    results = {}
    nmax_levels = sorted(classes_by_nmax.keys())
    for i in range(len(nmax_levels) - 1):
        n_lo = nmax_levels[i]
        n_hi = nmax_levels[i + 1]
        cls_lo = classes_by_nmax[n_lo]
        cls_hi = classes_by_nmax[n_hi]
        sectors_lo = set(cls_lo.keys())
        sectors_hi = set(cls_hi.keys())
        # Inclusion: sectors at n_lo should be subset of sectors at n_hi
        inclusion_ok = sectors_lo.issubset(sectors_hi)
        # Pull-back identity: for each sector at n_lo, the value at n_hi must match
        mismatches = []
        for s in sorted(sectors_lo):
            if cls_hi[s] != cls_lo[s]:
                mismatches.append({
                    "sector": list(s),
                    "value_n_lo": cls_lo[s],
                    "value_n_hi": cls_hi[s],
                })
        new_sectors = sorted(sectors_hi - sectors_lo)
        new_sector_values = {str(list(s)): cls_hi[s] for s in new_sectors}
        results[f"P_{n_hi}_to_{n_lo}"] = {
            "n_lo": n_lo,
            "n_hi": n_hi,
            "inclusion_ok": inclusion_ok,
            "pullback_bit_exact": (len(mismatches) == 0),
            "mismatches": mismatches,
            "new_sectors_at_n_hi": [list(s) for s in new_sectors],
            "new_sector_values": new_sector_values,
        }
    return results


# =====================================================================
# Sums (total invariants) and polynomial closed forms
# =====================================================================


def total_class_sums(
    classes_by_nmax: Dict[int, Dict[Tuple[int, int], int]],
) -> Dict[int, int]:
    """Sum of per-sector class values at each cutoff."""
    return {nm: sum(cls.values()) for nm, cls in classes_by_nmax.items()}


def jlo_per_sector_closed_form(n: int, l: int) -> int:
    """Closed form for chi_{(n, l)} = Tr(gamma e_{(n,l)}).

    Each (n, l) sector contains Dirac states with kappa in {-(l+1), +l}
    (the two allowed kappa for orbital angular momentum l), and chirality
    chi = sign(-kappa) = +1 for kappa < 0, -1 for kappa > 0.

    State count per kappa is 2*j+1 = 2*|kappa| (since j = |kappa| - 1/2).

    For l in {0, ..., n-1}: kappa = -(l+1) (chi=+1) gives 2(l+1) states,
                            kappa = +l (chi=-1) gives 2l states (if l > 0,
                            else 0 — l=0 has no kappa=+0).
        chi_{(n,l)} = 2(l+1) - 2l = 2 for l > 0
                    = 2(l+1) - 0 = 2 for l = 0
    Both give +2 for l < n.

    For l = n (only one kappa = +l, since -(l+1) = -(n+1) does not appear):
        Actually the CH labeling has kappa range {-n, ..., -1, +1, ..., +(n-1)}
        for n_max shell, so at orbital l = n-1 we get the boundary.

    The empirical pattern from n_max=2: (1,0)+2, (1,1)-2, (2,0)+2, (2,1)+2, (2,2)-4.
    The l = n sector has chi = -2n (all states at chi = -1, count = 2n).

    Derivation from the data:
    - (1, 0): chi = +2 (2 states all chi=+1)
    - (1, 1): chi = -2 (2 states all chi=-1)
    - (2, 0): chi = +2
    - (2, 1): chi = +2 (6 states, 4 chi=+1 + 2 chi=-1 = +2)
    - (2, 2): chi = -4 (4 states all chi=-1)

    Pattern: for l = n (top angular momentum at shell n), chi = -2n.
             for l < n, chi = +2.
    """
    if l == n:
        return -2 * n
    else:
        return 2


def cm_eta_per_sector_closed_form(n: int, l: int) -> int:
    """Closed form for eta_{(n, l)} = Tr(gamma D e_{(n,l)}) = dim_{(n,l)} * (n + 1/2).

    From the CM-eta memo: each sector (n, l) contributes
    dim_s * (n_s + 1/2) to Tr(gamma Lambda).

    Sector dimensions per the empirical n_max=2 data:
    - (1, 0): dim 2, eta = 2 * 3/2 = 3
    - (1, 1): dim 2, eta = 2 * 3/2 = 3
    - (2, 0): dim 2, eta = 2 * 5/2 = 5
    - (2, 1): dim 6, eta = 6 * 5/2 = 15
    - (2, 2): dim 4, eta = 4 * 5/2 = 10
    Sum = 36 = M_3(2). ✓

    Pattern for dim_{(n, l)}: for l < n: dim = 2(2l + 1); for l = n: dim = 2n.
    Hence:
    eta_{(n, l)} = (2l + 1)(2n + 1) for l < n
                 = n(2n + 1)        for l = n
    """
    if l == n:
        return n * (2 * n + 1)
    else:
        return (2 * l + 1) * (2 * n + 1)


def verify_closed_forms(
    classes_jlo: Dict[Tuple[int, int], int],
    classes_eta: Dict[Tuple[int, int], int],
) -> Dict:
    """Verify the per-sector closed forms against computed data."""
    jlo_check = []
    eta_check = []
    for (n, l), val in classes_jlo.items():
        predicted = jlo_per_sector_closed_form(n, l)
        jlo_check.append({
            "sector": [n, l],
            "computed": val,
            "predicted": predicted,
            "match": val == predicted,
        })
    for (n, l), val in classes_eta.items():
        predicted = cm_eta_per_sector_closed_form(n, l)
        eta_check.append({
            "sector": [n, l],
            "computed": val,
            "predicted": predicted,
            "match": val == predicted,
        })
    return {
        "jlo_per_sector_closed_form": jlo_check,
        "jlo_all_match": all(c["match"] for c in jlo_check),
        "cm_eta_per_sector_closed_form": eta_check,
        "cm_eta_all_match": all(c["match"] for c in eta_check),
    }


# =====================================================================
# Berezin-direction check (forward inclusion)
# =====================================================================


def berezin_forward_check(
    classes_by_nmax: Dict[int, Dict[Tuple[int, int], int]],
) -> Dict:
    """Berezin direction (coarser -> finer): is the class at level n
    obtained by SUBSET selection from the class at level n+1?

    B_{n -> n+1} maps idempotent e_s^{(n)} -> e_s^{(n+1)}; on classes
    this is just sector inclusion. The Berezin pull-back consistency is
    that the cocycle class at coarser cutoff equals the restriction of
    the cocycle class at finer cutoff to old sectors. Same condition as
    the truncation pull-back (just verifying the other direction).
    """
    out = {}
    nmax_levels = sorted(classes_by_nmax.keys())
    for i in range(len(nmax_levels) - 1):
        n_lo = nmax_levels[i]
        n_hi = nmax_levels[i + 1]
        cls_lo = classes_by_nmax[n_lo]
        cls_hi = classes_by_nmax[n_hi]
        # Berezin condition: for every sector in coarser class, its value
        # at finer class equals its value at coarser class.
        mismatches = []
        for s, v_lo in cls_lo.items():
            v_hi = cls_hi.get(s, None)
            if v_hi != v_lo:
                mismatches.append({"sector": list(s), "v_lo": v_lo, "v_hi": v_hi})
        out[f"B_{n_lo}_to_{n_hi}"] = {
            "n_lo": n_lo,
            "n_hi": n_hi,
            "compatibility_bit_exact": (len(mismatches) == 0),
            "mismatches": mismatches,
        }
    return out


# =====================================================================
# Polynomial closed forms for total invariants (sums)
# =====================================================================


def total_polynomial_predictions(n_max: int) -> Dict[str, int]:
    """Closed forms from Sub-Sprint 2a re-derived:
    M_1(n) = 2 n (n+1)(n+2) / 3  (dim H)
    M_3(n) = n (n+1)^2 (n+2) / 2  (Tr(gamma Lambda))
    Sum of JLO class = 0 (McKean-Singer = 0 for all n_max).
    Sum of CM-eta class = M_3(n_max).
    """
    M1 = 2 * n_max * (n_max + 1) * (n_max + 2) // 3
    M3 = n_max * (n_max + 1) ** 2 * (n_max + 2) // 2
    return {
        "M_1_predicted": M1,
        "M_3_predicted": M3,
        "jlo_class_sum_predicted": 0,
        "cm_eta_class_sum_predicted": M3,
    }


# =====================================================================
# Continuum-limit growth rate cross-check
# =====================================================================


def growth_rate_check(
    cm_eta_sums: Dict[int, int],
    dim_H_by_nmax: Dict[int, int],
) -> Dict:
    """Verify the polynomial growth rate of M_3 sum vs n_max.

    Expected: M_3(n) ~ n^4 / 2 leading order.
    M_3(1)=6, M_3(2)=36, M_3(3)=120, M_3(4)=300.

    Ratios: 36/6=6, 120/36~3.33, 300/120=2.5.
    Leading-order ratio for n^4: (n+1)^4 / n^4 -> 1 as n -> oo.
    """
    nmax_levels = sorted(cm_eta_sums.keys())
    M_3_values = [cm_eta_sums[n] for n in nmax_levels]
    M_3_predicted = [
        n * (n + 1) ** 2 * (n + 2) // 2 for n in nmax_levels
    ]
    ratios = []
    for i in range(1, len(M_3_values)):
        n_lo = nmax_levels[i - 1]
        n_hi = nmax_levels[i]
        if M_3_values[i - 1] > 0:
            ratios.append({
                "from_n": n_lo,
                "to_n": n_hi,
                "ratio": Rational(M_3_values[i], M_3_values[i - 1]),
            })
    return {
        "M_3_computed": M_3_values,
        "M_3_predicted": M_3_predicted,
        "match": M_3_values == M_3_predicted,
        "consecutive_ratios": [
            {"from_n": r["from_n"], "to_n": r["to_n"], "ratio": str(r["ratio"])}
            for r in ratios
        ],
        "asymptotic_leading_order": "n_max^4 / 2 (M_3); n_max^3 * 2/3 (M_1); cubic growth in N(n_max) ~ n_max^2 / 2",
    }


# =====================================================================
# Main driver
# =====================================================================


def main() -> None:
    print("=" * 72)
    print("Sprint Q5'-Stage1-Prosystem")
    print("Pro-system functoriality of JLO HP^even and CM-eta classes")
    print("across n_max in {1, 2, 3, 4}")
    print("=" * 72)

    t_global = time.time()

    cutoffs = [1, 2, 3, 4]
    triples: Dict[int, FockSpectralTriple] = {}
    sector_lists: Dict[int, List[Tuple[int, int]]] = {}
    N_actual: Dict[int, int] = {}
    N_closed_form: Dict[int, int] = {}
    dim_H: Dict[int, int] = {}
    jlo_classes: Dict[int, Dict[Tuple[int, int], int]] = {}
    eta_classes: Dict[int, Dict[Tuple[int, int], int]] = {}
    sector_dims: Dict[int, Dict[Tuple[int, int], int]] = {}

    print("\n[1] Constructing spectral triples and per-sector class data...")
    for nm in cutoffs:
        t0 = time.time()
        st = FockSpectralTriple(n_max=nm)
        triples[nm] = st
        sector_lists[nm] = list(st.sectors)
        N_actual[nm] = st.n_sectors
        N_closed_form[nm] = N_sectors_closed_form(nm)
        dim_H[nm] = st.dim_H
        jlo_classes[nm] = jlo_class_per_sector(st)
        eta_classes[nm] = cm_eta_class_per_sector(st)
        sector_dims[nm] = sector_dim_per_sector(st)
        elapsed = time.time() - t0
        print(f"    n_max={nm}: dim_H={st.dim_H}, n_sectors={st.n_sectors}"
              f" (closed form {N_closed_form[nm]}), {elapsed:.2f}s")

    # Verify sector count closed form
    print("\n[2] Sector count closed form N(n_max) = n_max(n_max+3)/2:")
    sector_count_match_all = True
    for nm in cutoffs:
        match = N_actual[nm] == N_closed_form[nm]
        sector_count_match_all = sector_count_match_all and match
        print(f"    n_max={nm}: N_actual={N_actual[nm]}, N_closed_form={N_closed_form[nm]}, {'MATCH' if match else 'MISMATCH'}")
    print(f"    All cutoffs match: {sector_count_match_all}")

    # Display per-sector classes
    print("\n[3] Per-sector class data at each cutoff:")
    print("\n    JLO HP^even class chi_s = Tr(gamma e_s):")
    for nm in cutoffs:
        cls = jlo_classes[nm]
        vec = [cls[s] for s in sector_lists[nm]]
        print(f"      n_max={nm}: {vec}, sum={sum(vec)}")

    print("\n    CM-eta class eta_s = Tr(gamma D e_s):")
    for nm in cutoffs:
        cls = eta_classes[nm]
        vec = [cls[s] for s in sector_lists[nm]]
        print(f"      n_max={nm}: {vec}, sum={sum(vec)}")

    # Per-sector closed-form verification at each cutoff
    print("\n[4] Per-sector closed-form verification:")
    closed_form_checks = {}
    for nm in cutoffs:
        cf = verify_closed_forms(jlo_classes[nm], eta_classes[nm])
        closed_form_checks[nm] = cf
        print(f"    n_max={nm}: JLO all-match={cf['jlo_all_match']}, CM-eta all-match={cf['cm_eta_all_match']}")

    # Pull-back verification (truncation P^*)
    print("\n[5] Pull-back P^*_{n+1 -> n} on JLO HP^even class:")
    jlo_pullback = verify_pullback_jlo(jlo_classes)
    for key, res in jlo_pullback.items():
        print(f"    {key}: inclusion={res['inclusion_ok']}, "
              f"pullback_bit_exact={res['pullback_bit_exact']}, "
              f"new_sectors={res['new_sectors_at_n_hi']}")

    print("\n[6] Pull-back P^*_{n+1 -> n} on CM-eta class:")
    eta_pullback = verify_pullback_jlo(eta_classes)
    for key, res in eta_pullback.items():
        print(f"    {key}: inclusion={res['inclusion_ok']}, "
              f"pullback_bit_exact={res['pullback_bit_exact']}, "
              f"new_sectors={res['new_sectors_at_n_hi']}")

    # Berezin forward direction
    print("\n[7] Berezin direction B_{n -> n+1} compatibility (coarser -> finer):")
    jlo_berezin = berezin_forward_check(jlo_classes)
    eta_berezin = berezin_forward_check(eta_classes)
    for key in jlo_berezin:
        jlo_ok = jlo_berezin[key]['compatibility_bit_exact']
        eta_ok = eta_berezin[key]['compatibility_bit_exact']
        print(f"    {key}: JLO bit-exact={jlo_ok}, CM-eta bit-exact={eta_ok}")

    # Total sums and polynomial closed forms
    print("\n[8] Total class sums vs polynomial closed forms:")
    jlo_sums = total_class_sums(jlo_classes)
    eta_sums = total_class_sums(eta_classes)
    total_predictions = {nm: total_polynomial_predictions(nm) for nm in cutoffs}
    for nm in cutoffs:
        pred = total_predictions[nm]
        print(f"    n_max={nm}: dim_H={dim_H[nm]} (M_1 pred {pred['M_1_predicted']}), "
              f"JLO_sum={jlo_sums[nm]} (pred 0), "
              f"eta_sum={eta_sums[nm]} (M_3 pred {pred['M_3_predicted']})")

    # Growth rate check
    print("\n[9] M_3 growth rate (continuum-limit consistency with Track 2):")
    growth = growth_rate_check(eta_sums, dim_H)
    print(f"    M_3 computed: {growth['M_3_computed']}")
    print(f"    M_3 predicted (Sub-Sprint 2a closed form): {growth['M_3_predicted']}")
    print(f"    Bit-exact match: {growth['match']}")
    for r in growth['consecutive_ratios']:
        print(f"    M_3(n={r['to_n']}) / M_3(n={r['from_n']}) = {r['ratio']}")
    print(f"    {growth['asymptotic_leading_order']}")

    # Pro-system formalization summary
    print("\n" + "=" * 72)
    print("PRO-SYSTEM FORMALIZATION SUMMARY")
    print("=" * 72)
    all_jlo_pullback = all(r['pullback_bit_exact'] for r in jlo_pullback.values())
    all_eta_pullback = all(r['pullback_bit_exact'] for r in eta_pullback.values())
    all_jlo_berezin = all(r['compatibility_bit_exact'] for r in jlo_berezin.values())
    all_eta_berezin = all(r['compatibility_bit_exact'] for r in eta_berezin.values())
    all_jlo_cf = all(closed_form_checks[nm]['jlo_all_match'] for nm in cutoffs)
    all_eta_cf = all(closed_form_checks[nm]['cm_eta_all_match'] for nm in cutoffs)

    print(f"\n  JLO HP^even pull-back P^* bit-exact across (1->2, 2->3, 3->4):  "
          f"{all_jlo_pullback}")
    print(f"  CM-eta pull-back P^* bit-exact across (1->2, 2->3, 3->4):       "
          f"{all_eta_pullback}")
    print(f"  JLO HP^even Berezin B bit-exact (coarser embeds in finer):     "
          f"{all_jlo_berezin}")
    print(f"  CM-eta Berezin B bit-exact (coarser embeds in finer):          "
          f"{all_eta_berezin}")
    print(f"  JLO per-sector closed form chi_{{(n,l)}} matches all cutoffs:  "
          f"{all_jlo_cf}")
    print(f"  CM-eta per-sector closed form eta_{{(n,l)}} matches all cutoffs: "
          f"{all_eta_cf}")
    print(f"  M_3 polynomial closed form n(n+1)^2(n+2)/2 matches all cutoffs:  "
          f"{growth['match']}")

    verdict = (all_jlo_pullback and all_eta_pullback and
               all_jlo_berezin and all_eta_berezin and
               all_jlo_cf and all_eta_cf and growth['match'])
    print(f"\n  VERDICT: {'POSITIVE (bit-exact pro-system functoriality)' if verdict else 'NOT-POSITIVE'}")

    elapsed = time.time() - t_global
    print(f"\nTotal wall time: {elapsed:.2f}s")

    # Save data
    data_out = {
        "sprint": "Q5'-Stage1-Prosystem",
        "cutoffs": cutoffs,
        "wall_time_seconds": elapsed,
        "sector_count_closed_form": "N(n_max) = n_max(n_max+3)/2",
        "sector_counts": {nm: {"N_actual": N_actual[nm], "N_closed_form": N_closed_form[nm]} for nm in cutoffs},
        "sectors_by_nmax": {nm: [list(s) for s in sector_lists[nm]] for nm in cutoffs},
        "dim_H_by_nmax": dim_H,
        "jlo_classes_by_nmax": {
            nm: [{"sector": list(s), "chi_s": v} for s, v in jlo_classes[nm].items()]
            for nm in cutoffs
        },
        "cm_eta_classes_by_nmax": {
            nm: [{"sector": list(s), "eta_s": v} for s, v in eta_classes[nm].items()]
            for nm in cutoffs
        },
        "sector_dims_by_nmax": {
            nm: [{"sector": list(s), "dim_s": v} for s, v in sector_dims[nm].items()]
            for nm in cutoffs
        },
        "closed_form_checks_per_sector": {
            nm: {
                "jlo_all_match": closed_form_checks[nm]["jlo_all_match"],
                "cm_eta_all_match": closed_form_checks[nm]["cm_eta_all_match"],
                "jlo_per_sector": closed_form_checks[nm]["jlo_per_sector_closed_form"],
                "cm_eta_per_sector": closed_form_checks[nm]["cm_eta_per_sector_closed_form"],
            }
            for nm in cutoffs
        },
        "jlo_pullback": jlo_pullback,
        "cm_eta_pullback": eta_pullback,
        "jlo_berezin": jlo_berezin,
        "cm_eta_berezin": eta_berezin,
        "total_sums": {
            "jlo_class_sums": jlo_sums,
            "cm_eta_class_sums": eta_sums,
            "predictions": {str(nm): total_predictions[nm] for nm in cutoffs},
        },
        "growth_rate_check": growth,
        "per_sector_closed_forms": {
            "JLO": "chi_{(n,l)} = 2 for l < n; chi_{(n,n)} = -2n; n_max-INDEPENDENT",
            "CM_eta": "eta_{(n,l)} = (2l+1)(2n+1) for l < n; eta_{(n,n)} = n(2n+1); n_max-INDEPENDENT",
        },
        "verdict": "POSITIVE" if verdict else "NOT-POSITIVE",
        "pro_system_statement": (
            "Both JLO HP^even and CM-eta classes are sector-LOCAL: per-sector "
            "value depends ONLY on sector label (n, l), NOT on cutoff n_max. "
            "Truncation pull-back P^*_{n+1 -> n} is exact sector projection. "
            "Berezin B_{n -> n+1} is exact sector inclusion. The cocycle classes "
            "form a strictly compatible inverse system in {Q^{N(n)}} with sector "
            "addition along (n+1, l), 0 <= l <= n+1, and the cohomological "
            "growth is the polynomial M_3(n) = n(n+1)^2(n+2)/2 from sub-sprint 2a."
        ),
    }

    out_path = Path("debug/data/sprint_q5p_prosystem.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(data_out, f, indent=2, default=str)
    print(f"\nData saved: {out_path}")


if __name__ == "__main__":
    main()
