r"""
Sprint Q5'-Stage1-StrictStrong driver — strict-strong-form pro-system
functoriality at the COCHAIN-MORPHISM level of the JLO entire-cyclic
bicomplex of the truncated CH spectral triple, computed at $n_{\max} = 3$
with a $P_{3 \to 2}^*$ pull-back diagnostic.

Background
----------
- Sub-Sprint 2c at $n_{\max} = 2$ (`compute_jlo_bicomplex.py`,
  `sprint_q5p_2c_bicomplex_memo.md`) found:
    * Degree-1 closure $(b\phi_0 + B\phi_2) = 0$ bit-exactly on the full
      idempotent panel.
    * Degree-3 residual $\pm 1/(3 \cdot 2^{16})$ on the palindromes
      $(e_2, e_3, e_3, e_2)$ and $(e_3, e_2, e_2, e_3)$;
      $-1/196608$ on $(e_2, e_2, e_3, e_3)$; bit-exact zero on
      $(e_2, e_3, e_2, e_3)$ and $(e_0, e_1, e_1, e_0)$.

- v3.60.0 pro-system (`compute_q5p_prosystem.py`) closed pro-system
  functoriality at the level of cocycle CLASSES.

The strict-strong-form question
-------------------------------
Does the $\pm 1/196608$ degree-3 residual persist at $n_{\max} = 3$
when one
  (a) computes $(b\phi_2^{(3)} + B\phi_4^{(3)})$ directly at $n_{\max} = 3$
      on the same palindromic 4-tuples (now lifted to the $n_{\max} = 3$
      ambient Hilbert space, but inputs still drawn from sectors
      $e_2 = (2,0)$ and $e_3 = (2,1)$), or
  (b) computes $P_{3 \to 2}^* [\phi_4^{(3)}]$ and checks the corrected
      degree-3 closure
        $b\phi_2^{(2)} + B(\phi_4^{(2)} + \Delta\phi_4) = 0$
      where $\Delta\phi_4 = P_{3 \to 2}^*[\phi_4^{(3)}] - \phi_4^{(2)}$.

If (a) gives the same $\pm 1/196608$ bit-exact: STRUCTURAL artifact
(persists at every cutoff; closes the Connes 1994 / Loday published
gap empirically).
If (b) gives bit-exact zero: FINITE-CUTOFF artifact (strict-strong
holds in the pro-limit; resolves the published gap, structurally
unexpected).

Diagnostic sub-panel
--------------------
At $n_{\max} = 2$: 5 sectors $e_0 \ldots e_4$ with
  $e_0 = (1,0), e_1 = (1,1), e_2 = (2,0), e_3 = (2,1), e_4 = (2,2)$.

At $n_{\max} = 3$: 9 sectors $e_0 \ldots e_8$ with
  $e_5 = (3,0), e_6 = (3,1), e_7 = (3,2), e_8 = (3,3)$.

OLD palindromic 4-tuples (involving only $n \le 2$ sectors):
  (e_2, e_3, e_3, e_2), (e_3, e_2, e_2, e_3),
  (e_2, e_2, e_3, e_3), (e_2, e_3, e_2, e_3),
  (e_0, e_1, e_1, e_0).

NEW palindromic 4-tuples involving $n = 3$ sectors:
  (e_5, e_6, e_6, e_5)  -- the structural analog (l<n adjacency)
  (e_5, e_5, e_6, e_6)
  (e_6, e_5, e_5, e_6)
  (e_5, e_6, e_5, e_6)
  (e_2, e_5, e_5, e_2)  -- cross-shell n=2 vs n=3
  (e_3, e_6, e_6, e_3)
  (e_3, e_5, e_5, e_3)  -- adjacent (l=1 vs l=0)

CM-$\eta$ cross-check
---------------------
Compute the analog of Steps 1-2 on the CM-$\eta$ cochain tower
($\gamma \to \gamma D$ in the leftmost slot). Track 1 found a clean
degree-1 closure at $n_{\max} = 2$ without the JLO artifact; the
question is whether the CM-$\eta$ degree-3 closure is clean at
$n_{\max} = 3$ on the same palindromic structures.

Output
------
- `debug/data/sprint_q5p_strict_strong.json` -- panel of bit-exact
  residuals at $n_{\max} \in \{2, 3\}$ for JLO and CM-$\eta$.
- `debug/sprint_q5p_strict_strong_memo.md` -- this sprint's memo.

Discipline
----------
- Bit-exact `sympy.Rational` throughout. No PSLQ. No floats.
- The diagnostic sub-panel is small (~14 four-tuples × 2 cocycle towers
  = ~28 panel cells per cutoff at $t^0$).
- $\phi_4$ on a single 4-tuple at $n_{\max} = 3$ (40x40 matrices)
  requires 5 matrix products + traces -- manageable.
"""

from __future__ import annotations

import json
import time
from itertools import product
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import sympy as sp
from sympy import Integer, Matrix, Rational, eye as sp_eye, factorial, zeros as sp_zeros

from geovac.spectral_triple import FockSpectralTriple


# =====================================================================
# Algebra helpers
# =====================================================================


def sector_idempotent(st: FockSpectralTriple, sector_idx: int) -> Matrix:
    N = st.dim_H
    M = sp_zeros(N, N)
    for i in range(N):
        if st._state_to_sector[i] == sector_idx:
            M[i, i] = Integer(1)
    return M


def commutator(D: Matrix, a: Matrix) -> Matrix:
    return D * a - a * D


# =====================================================================
# Unified cochain coefficient: gamma_left * a_0 * D^{2m_0} * [D,a_1] ...
# - JLO: gamma_left = gamma (even) or I (odd)
# - CM-eta: gamma_left = gamma * D (even) or D (odd)
# =====================================================================


def _partitions(total: int, n_parts: int) -> List[Tuple[int, ...]]:
    if n_parts == 1:
        return [(total,)]
    out = []
    for m0 in range(total + 1):
        for sub in _partitions(total - m0, n_parts - 1):
            out.append((m0,) + sub)
    return out


def cochain_n_coeff(
    D: Matrix,
    D2: Matrix,
    gamma_left: Matrix,
    a_list: List[Matrix],
    m_target: int,
) -> sp.Expr:
    """Compute the t^{m_target} coefficient of the cochain
       phi_n(a_0, ..., a_n; t) := gamma_left-prefactored
       JLO simplex/moment expansion.

    n = len(a_list) - 1.

    Formula:
      c_m = (-1)^m / (m + n)! *
            sum_{(m_0,...,m_n) | sum = m}
              Tr(gamma_left a_0 D^{2m_0} [D, a_1] D^{2m_1} ... [D, a_n] D^{2m_n})

    For JLO: gamma_left = gamma  (even) or I (odd).
    For CM-eta: gamma_left = gamma * D (even) or D (odd).
    """
    n = len(a_list) - 1
    a_0 = a_list[0]
    a_rest = a_list[1:]
    comm_list = [commutator(D, a) for a in a_rest]

    D2_powers: List[Matrix] = [sp_eye(D.shape[0])]
    for _ in range(m_target):
        D2_powers.append(D2_powers[-1] * D2)

    partitions = _partitions(m_target, n + 1)
    s_sum = Integer(0)
    for part in partitions:
        prod = gamma_left * a_0 * D2_powers[part[0]]
        for j in range(n):
            prod = prod * comm_list[j] * D2_powers[part[j + 1]]
        s_sum = s_sum + prod.trace()
    prefac = Rational((-1) ** m_target, sp.factorial(m_target + n))
    return sp.simplify(prefac * s_sum)


# =====================================================================
# Hochschild b and Connes B (same as Sub-Sprint 2c / Track 1)
# =====================================================================


def hochschild_b_on_phi_n(
    phi_n_eval,
    inputs: List[Matrix],
    n: int,
) -> sp.Expr:
    if len(inputs) != n + 2:
        raise ValueError(f"need n+2 = {n+2} inputs, got {len(inputs)}")
    total = Integer(0)
    for i in range(n + 1):
        new_args = inputs[:i] + [inputs[i] * inputs[i + 1]] + inputs[i + 2:]
        assert len(new_args) == n + 1
        sign = (-1) ** i
        total = total + sign * phi_n_eval(*new_args)
    cyclic_args = [inputs[n + 1] * inputs[0]] + inputs[1: n + 1]
    sign = (-1) ** (n + 1)
    total = total + sign * phi_n_eval(*cyclic_args)
    return sp.simplify(total)


def connes_B_on_phi_n(
    phi_n_eval,
    inputs: List[Matrix],
    n: int,
    unit: Matrix,
) -> sp.Expr:
    if len(inputs) != n:
        raise ValueError(f"need n = {n} inputs, got {len(inputs)}")
    total = Integer(0)
    for j in range(n):
        permuted = [inputs[(j + k) % n] for k in range(n)]
        args = [unit] + permuted
        sign = (-1) ** ((n - 1) * j)
        total = total + sign * phi_n_eval(*args)
    return sp.simplify(total)


# =====================================================================
# Make t^m cochain evaluators with gamma_left baked in
# =====================================================================


def make_cochain_eval(
    D: Matrix, D2: Matrix, gamma_left: Matrix, n: int, m: int,
):
    """Closure returning t^m coefficient of cochain (defined by gamma_left)
    at given (a_0, ..., a_n)."""
    def evaluator(*args):
        a_list = list(args)
        assert len(a_list) == n + 1
        return cochain_n_coeff(D, D2, gamma_left, a_list, m)
    return evaluator


# =====================================================================
# Degree-3 residual at a single 4-tuple
# =====================================================================


def degree3_residual_single(
    D: Matrix,
    D2: Matrix,
    gamma_left: Matrix,
    inputs: List[Matrix],
    m: int,
    unit: Matrix,
) -> Dict[str, sp.Expr]:
    """Compute (b phi_2 + B phi_4) at t^m on a single 4-tuple."""
    phi_2_eval = make_cochain_eval(D, D2, gamma_left, n=2, m=m)
    phi_4_eval = make_cochain_eval(D, D2, gamma_left, n=4, m=m)

    b_term = hochschild_b_on_phi_n(phi_2_eval, inputs, n=2)
    B_term = connes_B_on_phi_n(phi_4_eval, inputs, n=4, unit=unit)
    residual = sp.simplify(b_term + B_term)
    return {"b_phi_2": b_term, "B_phi_4": B_term, "residual": residual}


# =====================================================================
# Degree-1 residual (b phi_0 + B phi_2) at a single pair
# =====================================================================


def degree1_residual_single(
    D: Matrix,
    D2: Matrix,
    gamma_left: Matrix,
    inputs: List[Matrix],
    m: int,
    unit: Matrix,
) -> Dict[str, sp.Expr]:
    phi_0_eval = make_cochain_eval(D, D2, gamma_left, n=0, m=m)
    phi_2_eval = make_cochain_eval(D, D2, gamma_left, n=2, m=m)
    b_term = hochschild_b_on_phi_n(phi_0_eval, inputs, n=0)
    B_term = connes_B_on_phi_n(phi_2_eval, inputs, n=2, unit=unit)
    residual = sp.simplify(b_term + B_term)
    return {"b_phi_0": b_term, "B_phi_2": B_term, "residual": residual}


# =====================================================================
# Diagnostic sub-panel definitions
# =====================================================================


def build_diagnostic_panel_nmax2(idem):
    """Sub-Sprint 2c original palindromic 4-tuples at n_max=2."""
    e0, e1, e2, e3, e4 = idem[0], idem[1], idem[2], idem[3], idem[4]
    return [
        ("e2_e3_e3_e2_palindrome", [e2, e3, e3, e2]),
        ("e3_e2_e2_e3_palindrome", [e3, e2, e2, e3]),
        ("e2_e2_e3_e3", [e2, e2, e3, e3]),
        ("e2_e3_e2_e3", [e2, e3, e2, e3]),
        ("e0_e1_e1_e0_palindrome", [e0, e1, e1, e0]),
    ]


def build_diagnostic_panel_nmax3(idem):
    """OLD + NEW palindromic 4-tuples at n_max=3.

    OLD = same indices as in n_max=2 panel (still sectors 0..4 because
    sector ordering is canonical). NEW = palindromes involving n=3
    sectors (idem[5..8]).
    """
    e0, e1, e2, e3, e4 = idem[0], idem[1], idem[2], idem[3], idem[4]
    e5, e6, e7, e8 = idem[5], idem[6], idem[7], idem[8]
    panel = [
        # OLD palindromes (n<=2 sectors only) -- should be bit-exact analogs
        # of Sub-Sprint 2c if STRUCTURAL.
        ("e2_e3_e3_e2_palindrome_OLD", [e2, e3, e3, e2]),
        ("e3_e2_e2_e3_palindrome_OLD", [e3, e2, e2, e3]),
        ("e2_e2_e3_e3_OLD", [e2, e2, e3, e3]),
        ("e2_e3_e2_e3_OLD", [e2, e3, e2, e3]),
        ("e0_e1_e1_e0_palindrome_OLD", [e0, e1, e1, e0]),
        # NEW palindromes (involving n=3 sectors)
        # Analog of (e_2, e_3, ...) at next shell: (e_5, e_6, ...)
        ("e5_e6_e6_e5_palindrome_NEW", [e5, e6, e6, e5]),
        ("e6_e5_e5_e6_palindrome_NEW", [e6, e5, e5, e6]),
        ("e5_e5_e6_e6_NEW", [e5, e5, e6, e6]),
        ("e5_e6_e5_e6_NEW", [e5, e6, e5, e6]),
        # Cross-shell (n=2 vs n=3) palindromes
        ("e2_e5_e5_e2_palindrome_CROSS", [e2, e5, e5, e2]),
        ("e3_e6_e6_e3_palindrome_CROSS", [e3, e6, e6, e3]),
        ("e3_e5_e5_e3_palindrome_CROSS", [e3, e5, e5, e3]),
        # Higher-l NEW palindrome (top angular momentum at shell 3)
        ("e7_e8_e8_e7_palindrome_NEW_high_l", [e7, e8, e8, e7]),
        ("e7_e7_e8_e8_NEW_high_l", [e7, e7, e8, e8]),
    ]
    return panel


# =====================================================================
# Pull-back P^*_{3 -> 2} on cochains: restrict input matrices to
# sectors 0..4 at n_max = 3, evaluate the cochain in the n_max=3
# ambient. The "pull-back of phi_4^{(3)} restricted to A^{(2)} inputs"
# is just phi_4^{(3)} evaluated on idempotents from sectors 0..4.
# =====================================================================


def pull_back_cochain_eval(
    D3: Matrix, D2_3: Matrix, gamma_left_3: Matrix,
    n: int, m: int,
):
    """Returns evaluator of phi_n^{(3)} on idempotent inputs (which must
    be n_max=3 idempotents living in sectors 0..4 for proper pull-back).

    The pull-back P^*_{3 -> 2}[phi_n^{(3)}] is the restriction of the
    n_max=3 cochain to inputs in the n_max=2 sub-algebra.
    """
    return make_cochain_eval(D3, D2_3, gamma_left_3, n, m)


# =====================================================================
# Convert n_max=2 idempotents to n_max=3 idempotents (same sector labels)
# =====================================================================


def get_idempotents_for(st: FockSpectralTriple) -> List[Matrix]:
    return [sector_idempotent(st, s) for s in range(st.n_sectors)]


# =====================================================================
# Main verification
# =====================================================================


def verify_at_cutoff(
    st: FockSpectralTriple,
    panel: List[Tuple[str, List[Matrix]]],
    cocycle_kind: str,  # "JLO" or "CM_eta"
    flavor: str,  # "even" or "odd"
    m: int = 0,
) -> Dict:
    """Compute degree-3 residual on the panel at the given t-order."""
    N = st.dim_H
    unit = sp_eye(N)
    D = st.dirac_operator
    D2 = D * D
    gamma = st.grading

    if flavor == "even":
        if cocycle_kind == "JLO":
            gamma_left = gamma
        elif cocycle_kind == "CM_eta":
            gamma_left = gamma * D
        else:
            raise ValueError(f"unknown cocycle_kind: {cocycle_kind}")
    elif flavor == "odd":
        if cocycle_kind == "JLO":
            gamma_left = sp_eye(N)
        elif cocycle_kind == "CM_eta":
            gamma_left = D
        else:
            raise ValueError(f"unknown cocycle_kind: {cocycle_kind}")
    else:
        raise ValueError(f"unknown flavor: {flavor}")

    results = {}
    for label, inputs in panel:
        try:
            res = degree3_residual_single(
                D, D2, gamma_left, inputs, m=m, unit=unit
            )
            results[label] = {
                "b_phi_2": str(res["b_phi_2"]),
                "B_phi_4": str(res["B_phi_4"]),
                "residual": str(res["residual"]),
            }
        except Exception as e:
            results[label] = {"error": str(e)}
    return results


def verify_degree1_at_cutoff(
    st: FockSpectralTriple,
    cocycle_kind: str,
    flavor: str,
    m: int = 0,
) -> Dict:
    """Sanity check: degree-1 closure on a small pair panel."""
    N = st.dim_H
    unit = sp_eye(N)
    D = st.dirac_operator
    D2 = D * D
    gamma = st.grading

    if flavor == "even":
        gamma_left = gamma if cocycle_kind == "JLO" else gamma * D
    else:
        gamma_left = sp_eye(N) if cocycle_kind == "JLO" else D

    idem = get_idempotents_for(st)
    # Test pairs same as Sub-Sprint 2c / Track 1 representative subset
    pairs = []
    n_sec = st.n_sectors
    # All idempotent pairs on sectors {2, 3} (relevant to JLO artifact)
    for s in [2, 3]:
        for t in [2, 3]:
            pairs.append((f"e{s}_e{t}", [idem[s], idem[t]]))
    # And one (1, e_s) pair
    pairs.append(("1_e2", [unit, idem[2]]))
    pairs.append(("1_e3", [unit, idem[3]]))

    out = {}
    for label, inputs in pairs:
        res = degree1_residual_single(D, D2, gamma_left, inputs, m=m, unit=unit)
        out[label] = {
            "b_phi_0": str(res["b_phi_0"]),
            "B_phi_2": str(res["B_phi_2"]),
            "residual": str(res["residual"]),
        }
    return out


# =====================================================================
# Pull-back analysis: compare phi_4^{(3)} restricted to n<=2 sectors
# against phi_4^{(2)} directly.
# =====================================================================


def pull_back_delta_phi4_analysis(
    st2: FockSpectralTriple,
    st3: FockSpectralTriple,
    panel_inputs_indices: List[Tuple[str, List[int]]],
    cocycle_kind: str,
    flavor: str,
) -> Dict:
    """For each OLD palindromic 4-tuple (indices in {0,1,2,3,4}), compute
       phi_4^{(2)}(...) (using n_max=2 triple) and phi_4^{(3)}(...) (using
       n_max=3 triple, with idempotents lifted to dim_H=40 matrices).

       Then compute Delta_phi_4 = phi_4^{(3)} - phi_4^{(2)}.

       This is the cochain-level pull-back diagnostic at degree 4.
       If Delta_phi_4 == 0 bit-exactly: strict cochain-morphism on
         this input (no correction needed).
       If Delta_phi_4 != 0: the n_max=3 cochain provides additional
         content on the input.

       Then check whether (b phi_2^{(2)} + B[phi_4^{(2)} + Delta_phi_4])
       = (b phi_2^{(2)} + B phi_4^{(3)}_restricted) at the n_max=2 ambient.

       Here B is the Connes operator on the (n_max=2) algebra; we use
       the n_max=3 cochain values evaluated at n_max=3 idempotents
       which restrict to sectors 0..4.
    """
    N2 = st2.dim_H
    N3 = st3.dim_H
    D2_op = st2.dirac_operator
    D2_sq = D2_op * D2_op
    D3_op = st3.dirac_operator
    D3_sq = D3_op * D3_op
    gamma_2 = st2.grading
    gamma_3 = st3.grading

    if flavor == "even":
        gleft_2 = gamma_2 if cocycle_kind == "JLO" else gamma_2 * D2_op
        gleft_3 = gamma_3 if cocycle_kind == "JLO" else gamma_3 * D3_op
    else:
        gleft_2 = sp_eye(N2) if cocycle_kind == "JLO" else D2_op
        gleft_3 = sp_eye(N3) if cocycle_kind == "JLO" else D3_op

    idem2 = get_idempotents_for(st2)
    idem3 = get_idempotents_for(st3)

    # phi_4 evaluators at t^0
    phi4_2_eval = make_cochain_eval(D2_op, D2_sq, gleft_2, n=4, m=0)
    phi4_3_eval = make_cochain_eval(D3_op, D3_sq, gleft_3, n=4, m=0)

    # Also need phi_2 evaluators for the corrected residual check.
    phi2_2_eval = make_cochain_eval(D2_op, D2_sq, gleft_2, n=2, m=0)
    phi2_3_eval = make_cochain_eval(D3_op, D3_sq, gleft_3, n=2, m=0)

    unit2 = sp_eye(N2)

    results = {}
    for label, idx_4tuple in panel_inputs_indices:
        # Build inputs at n_max=2 and n_max=3 with same sector indices
        # phi_4 takes 5 args (a_0,...,a_4); we evaluate on (1, e_a, e_b, e_c, e_d)
        # for the Connes B sub-computation.
        # The diagnostic is at the cochain level: phi_4(1, e_a, e_b, e_c, e_d)
        # for each cyclic shift.
        # But the natural pull-back diagnostic is at the FULL B-input level:
        # compute B phi_4 on the 4-tuple (e_a, e_b, e_c, e_d) directly,
        # which already involves the cyclic-shift sum over phi_4 inputs.

        inputs2 = [idem2[i] for i in idx_4tuple]
        inputs3 = [idem3[i] for i in idx_4tuple]

        # Direct B phi_4 at n_max=2 and n_max=3
        B_phi4_2 = connes_B_on_phi_n(phi4_2_eval, inputs2, n=4, unit=unit2)
        unit3 = sp_eye(N3)
        B_phi4_3 = connes_B_on_phi_n(phi4_3_eval, inputs3, n=4, unit=unit3)

        # b phi_2 at n_max=2 and n_max=3
        b_phi2_2 = hochschild_b_on_phi_n(phi2_2_eval, inputs2, n=2)
        b_phi2_3 = hochschild_b_on_phi_n(phi2_3_eval, inputs3, n=2)

        residual_2 = sp.simplify(b_phi2_2 + B_phi4_2)
        residual_3 = sp.simplify(b_phi2_3 + B_phi4_3)

        delta_B_phi4 = sp.simplify(B_phi4_3 - B_phi4_2)
        delta_b_phi2 = sp.simplify(b_phi2_3 - b_phi2_2)

        results[label] = {
            "input_idx_4tuple": idx_4tuple,
            "b_phi_2_nmax2": str(b_phi2_2),
            "b_phi_2_nmax3": str(b_phi2_3),
            "delta_b_phi_2_3_minus_2": str(delta_b_phi2),
            "B_phi_4_nmax2": str(B_phi4_2),
            "B_phi_4_nmax3": str(B_phi4_3),
            "delta_B_phi_4_3_minus_2": str(delta_B_phi4),
            "residual_nmax2": str(residual_2),
            "residual_nmax3": str(residual_3),
            "delta_residual_3_minus_2": str(sp.simplify(residual_3 - residual_2)),
        }
    return results


# =====================================================================
# Main entry
# =====================================================================


def main() -> None:
    print("=" * 78)
    print("Sprint Q5'-Stage1-StrictStrong")
    print("Strict-strong-form pro-system functoriality at the cochain-morphism")
    print("level of the JLO entire-cyclic bicomplex")
    print("=" * 78)
    t_global = time.time()

    print("\n[1] Building FockSpectralTriple at n_max = 2 and n_max = 3...")
    st2 = FockSpectralTriple(n_max=2)
    st3 = FockSpectralTriple(n_max=3)
    print(f"    n_max=2: dim_H = {st2.dim_H}, n_sectors = {st2.n_sectors}")
    print(f"             sectors = {st2.sectors}")
    print(f"    n_max=3: dim_H = {st3.dim_H}, n_sectors = {st3.n_sectors}")
    print(f"             sectors = {st3.sectors}")

    # Sanity check: sectors 0..4 must agree between n_max=2 and n_max=3
    assert list(st2.sectors) == list(st3.sectors[:5]), \
        "Sector ordering mismatch between cutoffs"

    idem2 = get_idempotents_for(st2)
    idem3 = get_idempotents_for(st3)

    panel2 = build_diagnostic_panel_nmax2(idem2)
    panel3 = build_diagnostic_panel_nmax3(idem3)

    print(f"\n[2] Diagnostic panel:")
    print(f"    n_max=2: {len(panel2)} 4-tuples (the original Sub-Sprint 2c panel)")
    print(f"    n_max=3: {len(panel3)} 4-tuples (5 OLD + 9 NEW/CROSS/high-l)")

    out = {
        "sprint": "Q5'-Stage1-StrictStrong",
        "n_max": [2, 3],
        "dim_H": {"n_max=2": st2.dim_H, "n_max=3": st3.dim_H},
        "n_sectors": {"n_max=2": st2.n_sectors, "n_max=3": st3.n_sectors},
        "sectors": {
            "n_max=2": [list(s) for s in st2.sectors],
            "n_max=3": [list(s) for s in st3.sectors],
        },
        "panel_labels": {
            "n_max=2": [lbl for lbl, _ in panel2],
            "n_max=3": [lbl for lbl, _ in panel3],
        },
    }

    # =================================================================
    # Step 1: JLO at n_max=2 (re-confirm baseline)
    # =================================================================
    print("\n[3] JLO degree-3 panel at n_max=2 (re-confirm Sub-Sprint 2c baseline)...")
    t0 = time.time()
    jlo_n2_even = verify_at_cutoff(st2, panel2, "JLO", "even", m=0)
    print(f"    JLO even at n_max=2: done in {time.time() - t0:.1f}s")
    for lbl, vals in jlo_n2_even.items():
        if "residual" in vals:
            print(f"      {lbl}: {vals['residual']}")
    out["jlo_nmax2_degree3_even_t0"] = jlo_n2_even

    # =================================================================
    # Step 2a: JLO at n_max=3 on OLD palindromes
    # =================================================================
    print("\n[4] JLO degree-3 panel at n_max=3...")
    print("    This is the LOAD-BEARING computation:")
    print("    if OLD palindromic residuals are bit-exact ± 1/196608 -> STRUCTURAL.")
    print("    if zero or different scaling -> RESOLVES or FINITE-CUTOFF.")
    t0 = time.time()
    jlo_n3_even = verify_at_cutoff(st3, panel3, "JLO", "even", m=0)
    elapsed = time.time() - t0
    print(f"    JLO even at n_max=3: done in {elapsed:.1f}s")
    for lbl, vals in jlo_n3_even.items():
        if "residual" in vals:
            print(f"      {lbl}: {vals['residual']}")
    out["jlo_nmax3_degree3_even_t0"] = jlo_n3_even
    out["jlo_nmax3_wall_seconds"] = elapsed

    # =================================================================
    # Step 3: Pull-back P^*_{3 -> 2} analysis on OLD palindromes
    # =================================================================
    print("\n[5] Pull-back P^*_{3 -> 2} analysis on OLD palindromes...")
    old_idx = [
        ("e2_e3_e3_e2_palindrome", [2, 3, 3, 2]),
        ("e3_e2_e2_e3_palindrome", [3, 2, 2, 3]),
        ("e2_e2_e3_e3", [2, 2, 3, 3]),
        ("e2_e3_e2_e3", [2, 3, 2, 3]),
        ("e0_e1_e1_e0_palindrome", [0, 1, 1, 0]),
    ]
    t0 = time.time()
    pullback_jlo = pull_back_delta_phi4_analysis(
        st2, st3, old_idx, "JLO", "even"
    )
    elapsed = time.time() - t0
    print(f"    Pull-back analysis done in {elapsed:.1f}s")
    for lbl, vals in pullback_jlo.items():
        print(f"      {lbl}:")
        print(f"        residual nmax=2: {vals['residual_nmax2']}")
        print(f"        residual nmax=3: {vals['residual_nmax3']}")
        print(f"        delta(B phi_4):  {vals['delta_B_phi_4_3_minus_2']}")
    out["pullback_jlo_even"] = pullback_jlo

    # =================================================================
    # Step 4: CM-eta cross-check at n_max=2 (Track 1 follow-on: degree-3)
    # =================================================================
    print("\n[6] CM-eta degree-3 panel at n_max=2 (extending Track 1 to degree 3)...")
    t0 = time.time()
    cm_n2_even = verify_at_cutoff(st2, panel2, "CM_eta", "even", m=0)
    print(f"    CM-eta even at n_max=2 done in {time.time() - t0:.1f}s")
    for lbl, vals in cm_n2_even.items():
        if "residual" in vals:
            print(f"      {lbl}: {vals['residual']}")
    out["cm_eta_nmax2_degree3_even_t0"] = cm_n2_even

    # =================================================================
    # Step 5: CM-eta at n_max=3 (smaller panel: OLD palindromes only)
    # =================================================================
    print("\n[7] CM-eta degree-3 panel at n_max=3 (OLD palindromes only)...")
    cm_panel3_subset = build_diagnostic_panel_nmax3(idem3)[:5]  # OLD only
    t0 = time.time()
    cm_n3_even = verify_at_cutoff(st3, cm_panel3_subset, "CM_eta", "even", m=0)
    elapsed = time.time() - t0
    print(f"    CM-eta even at n_max=3 (OLD only) done in {elapsed:.1f}s")
    for lbl, vals in cm_n3_even.items():
        if "residual" in vals:
            print(f"      {lbl}: {vals['residual']}")
    out["cm_eta_nmax3_degree3_even_t0_OLD_only"] = cm_n3_even
    out["cm_eta_nmax3_wall_seconds"] = elapsed

    # =================================================================
    # Step 6: Sanity -- degree-1 closure at n_max=3 for JLO
    # =================================================================
    print("\n[8] Sanity: JLO degree-1 closure at n_max=3...")
    t0 = time.time()
    jlo_n3_deg1 = verify_degree1_at_cutoff(st3, "JLO", "even", m=0)
    print(f"    done in {time.time() - t0:.1f}s")
    for lbl, vals in jlo_n3_deg1.items():
        if "residual" in vals:
            print(f"      {lbl}: {vals['residual']}")
    out["jlo_nmax3_degree1_even_t0_sanity"] = jlo_n3_deg1

    # =================================================================
    # Verdict
    # =================================================================
    print("\n" + "=" * 78)
    print("Verdict analysis")
    print("=" * 78)

    # Check OLD palindromic residuals at n_max=3
    old_palindrome_labels = [
        "e2_e3_e3_e2_palindrome_OLD",
        "e3_e2_e2_e3_palindrome_OLD",
        "e2_e2_e3_e3_OLD",
    ]
    expected_value_strs = [
        "1/196608", "1/196608", "-1/196608",
    ]
    old_residuals_n3 = [jlo_n3_even[lbl]["residual"] for lbl in old_palindrome_labels]
    print(f"\nOLD palindrome residuals at n_max=3:")
    for lbl, val, exp in zip(old_palindrome_labels, old_residuals_n3, expected_value_strs):
        match = "MATCH (structural)" if val == exp else "DIFFERENT"
        print(f"  {lbl}: {val}  vs expected {exp}  --> {match}")

    matches_structural = all(
        old_residuals_n3[i] == expected_value_strs[i] for i in range(3)
    )
    all_resolve = all(r == "0" for r in old_residuals_n3)

    if matches_structural:
        verdict = "POSITIVE-STRUCTURAL"
        verdict_text = (
            "The ± 1/196608 degree-3 residual on OLD palindromes persists "
            "bit-exactly at n_max=3, confirming the truncation artifact is "
            "STRUCTURAL (not finite-cutoff). Strict-strong-form pro-system "
            "functoriality at the cochain-morphism level FAILS bit-exactly "
            "on commutative algebra A at every cutoff. Empirical confirmation "
            "of the Connes 1994 / Loday published-open gap."
        )
    elif all_resolve:
        verdict = "POSITIVE-RESOLVES"
        verdict_text = (
            "The degree-3 residual at n_max=3 is bit-exact zero on the same "
            "palindromic 4-tuples. The artifact is finite-cutoff and the "
            "strict-strong-form holds in the pro-limit. Structurally "
            "unexpected; substantial Stage-2 advance."
        )
    else:
        verdict = "BORDERLINE"
        verdict_text = "Mixed pattern; see panel for details."

    print(f"\nVerdict: {verdict}")
    print(f"\n{verdict_text}")

    out["verdict"] = verdict
    out["verdict_text"] = verdict_text
    out["old_palindrome_residuals_nmax3"] = dict(zip(
        old_palindrome_labels, old_residuals_n3
    ))
    out["wall_seconds_total"] = time.time() - t_global

    # Save data
    out_path = Path("debug/data/sprint_q5p_strict_strong.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\n[OUT] Data written to {out_path}")
    print(f"[TIME] Total wall: {time.time() - t_global:.1f}s")


if __name__ == "__main__":
    main()
