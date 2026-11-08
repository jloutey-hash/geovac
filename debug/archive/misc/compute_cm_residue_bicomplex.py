r"""
Sprint Q5'-Stage1-CM-Bicomplex driver — explicit Connes--Moscovici (CM)
residue cyclic $(b, B)$ bicomplex on the truncated CH spectral triple
at $n_{\max} = 2$ (and $n_{\max} = 3$ for cross-cutoff sanity), housing
the M3 mechanism $\mathrm{Tr}(\gamma D\,e^{-t D^2})\big|_{t^0} = 36$
at $n_{\max} = 2$, with cocycle verification on the full idempotent
panel and identification of the HP class.

Goal
----
Stage 1 Sub-Sprint CM (the named follow-on from the Q5'-Stage1-Arc
umbrella memo: "CM-residue bicomplex: construct the full
Connes-Moscovici cyclic complex housing M3 explicitly; verify the
$\eta$-pairing cocycle condition; identify the analog HP^odd or
HP^even class on that complex"). Companion to Sub-Sprint 2c (JLO
bicomplex; `debug/compute_jlo_bicomplex.py`).

Setup recap (same as Sub-Sprint 1 / 2c)
---------------------------------------
- Algebra A = C^5: 5 sector idempotents e_0..e_4 on Fock sectors
  (1,0), (1,1), (2,0), (2,1), (2,2). Commutative.
- Hilbert space H = C^16, dim 16, balanced chirality 8+/8-.
- Dirac D = Lambda + kappa A with kappa = -1/16.
- Grading gamma diagonal.
- Bit-exact sympy.Rational throughout.

Connes--Moscovici residue cocycle: the structural shape
------------------------------------------------------
Reference: Connes, A.; Moscovici, H. "The local index formula in
noncommutative geometry." GAFA 5 (1995), 174-243 (arXiv:math/9806109
gives the cleaner formulation).

For a $p$-summable spectral triple $(A, H, D)$ with dimension spectrum
in arithmetic progression and isolated simple poles, the
CM-residue cocycle is the entire-cyclic cocycle $\{\psi_n\}_{n \ge 0}$
(odd or even parity) defined by

  $\psi_n(a^0, \ldots, a^n)
     = \sum_{k \ge 0} c_{n, k}\,
       \mathrm{Res}_{s = 0}\,\mathrm{Tr}\!\left(\gamma\,a^0\,
         [D, a^1]^{(k_1)}\,\cdots\,[D, a^n]^{(k_n)}\,
         |D|^{-2(n + 2|k| + s)}\right)$

with $X^{(k)} = \nabla^k(X)$, $\nabla = [D^2, \cdot]$, $|k| = \sum k_i$,
and $c_{n, k}$ explicit Connes--Moscovici coefficients. The CM cocycle
is in the SAME periodic-cyclic-cohomology class as JLO (CM 1995
Theorem II.3 = "local index theorem" equivalence).

KEY structural difference vs JLO at finite cutoff
--------------------------------------------------
On a FINITE-DIMENSIONAL truncated spectral triple, the spectral zeta
$\zeta_{D^2}(s) = \mathrm{Tr}(|D|^{-2s})$ is an ENTIRE function (no
poles; the truncated dimension is "effectively zero" in the CM sense).
Therefore $\mathrm{Res}_{s = 0}$ of any analytic continuation REDUCES
to the constant term (coefficient of $s^0$) of the meromorphic
extension --- which on a finite-dim spectral triple is simply the
finite trace evaluated at $s = 0$.

The CM residue cocycle at finite cutoff THEREFORE collapses to:

  $\psi_n^{\mathrm{CM}, \mathrm{finite}}(a^0, \ldots, a^n)
     = \sum_{k \ge 0} c_{n, k}\,
       \mathrm{Tr}\!\left(\gamma\,a^0\,[D, a^1]^{(k_1)}\,\cdots\,
         [D, a^n]^{(k_n)}\,(D^2)^{-(n + |k|)}\right)$

(formal expression; when $D^2$ has zero eigenvalues this needs care,
but on the truncated CH triple $D^2$ is bounded below away from zero
in any non-trivial sector --- $\Lambda$ has minimum $\pm 3/2$ from
$n_{\mathrm{Fock}} = 1$).

For the LEADING M3 mechanism at degree 0, the structure simplifies:

  $\psi_0^{\eta}(a) = \mathrm{Tr}(\gamma\,D\,a\,e^{-t D^2})\big|_{t^0}
                   = \mathrm{Tr}(\gamma\,D\,a)$

at the t^0 coefficient. This is the "$\eta$-density" analog at finite
cutoff: an OFF-DIAGONAL pairing between the chirality $\gamma$ and the
Dirac $D$, with $a$ inserted as algebra element.

The CM bicomplex we build
-------------------------
We construct the CM residue cochains $\{\psi_n^{\eta, \mathrm{even}},
\psi_n^{\eta, \mathrm{odd}}\}_{n \in \{0, 1, 2\}}$ on the truncated CH
triple at $n_{\max} = 2$, where:

  $\psi_n^{\eta, \mathrm{even}}(a^0, \ldots, a^n; t)$
     := t^0-coefficient of
        $\mathrm{Tr}\!\left(\gamma\,a^0\,D\,e^{-s_0 t D^2}\,[D, a^1]\,
         e^{-s_1 t D^2}\,\cdots\,[D, a^n]\,e^{-s_n t D^2}\right)$
        integrated over $\Delta_n$

  (and similarly $\psi_n^{\eta, \mathrm{odd}}$ without $\gamma$).

Equivalently --- using the same simplex/moment expansion as JLO ---
each $\psi_n^{\eta}$ is the JLO cochain with the $\gamma$ replaced by
$\gamma D$ (factoring D out of the leftmost slot). This is the
"$\eta$-style trace" identified in Sub-Sprint 1 as the M3 host.

The cocycle condition for the CM-residue bicomplex
---------------------------------------------------
On the parity-graded CM residue cochain tower
$\psi^{\eta, \mathrm{even}} = (\psi_0^{\eta}, \psi_2^{\eta}, \ldots)$,
the entire-cyclic-cocycle condition is the SAME $(b + B)\psi = 0$
as for JLO --- but applied to the $\eta$-modified cochain tower. We
verify:

  $b\psi_0^{\eta} + B\psi_2^{\eta} = 0$ bit-exactly on the full
  idempotent panel, at $t$-orders $\{t^0, t^1, t^2\}$, both flavors.

We then identify the leading $\mathrm{HP}^{??}$ class
$[\psi_0^{\eta}(e_s)]_{s = 0..4}$ in $\mathbb{Q}^5$ (the Morita-trivial
baseline; cross-checked against R3 closure).

We CROSS-CHECK the M3 sub-sprint 1 ground truth:
$\psi_0^{\eta, \mathrm{even}}(1) = \mathrm{Tr}(\gamma D) = 0$ on the
unit (the global M3 leading is structurally zero by chirality balance
combined with the $\Lambda$-spectrum being $\pm$-symmetric), and the
substantive M3 content lives in the SECTOR-RESOLVED data
$\psi_0^{\eta, \mathrm{even}}(e_s)$, with
$\sum_s \psi_0^{\eta, \mathrm{even}}(e_s) = 0$ (global) and
$\sum_s \psi_0^{\eta, \mathrm{odd}}(e_s) = \mathrm{Tr}(D)$.

Sub-Sprint 1 ground truth: $\mathrm{Tr}(\gamma D)|_{\text{full panel}}
= ???$ --- let's see. CH-1 reported $\mathrm{Tr}(\gamma \Lambda) = 36$
on the t^0 coefficient of $\mathrm{Tr}(\gamma D e^{-tD^2})$. The unit
trace $\mathrm{Tr}(\gamma\Lambda) = \sum_i \chi_i^2 (n_i + 1/2)
= \sum_i (n_i + 1/2)$ on the chirality-signed sum, and the
SIGNED structure is $\sum_i \chi_i (n_i + 1/2) \cdot \chi_i =
\sum_i (n_i + 1/2)$, which gives the sum of |eigenvalues| weighted by
chirality. At $n_{\max} = 2$, this is $4 \cdot 3/2 + 12 \cdot 5/2 =
6 + 30 = 36$. ✓

We expect the CM residue cocycle leading vector to have the structure
$[\mathrm{Tr}(\gamma D\,e_s)]_{s = 0..4}$ which, on the diagonal $\Lambda$,
is $[2 \cdot 3/2, -2 \cdot 3/2, 2 \cdot 5/2, ?, ?]$ where ? involves
sector dimensions and chirality balance.

Output
------
- `debug/data/sprint_q5p_cm_bicomplex.json` --- bit-exact rational data
  for $\psi_n^{\eta}$ at $n_{\max} = 2$, panel residuals, HP class
  cross-checks vs JLO HP^even.
- `debug/sprint_q5p_cm_bicomplex_memo.md` --- this sprint's memo (PI
  to apply Paper 55 edit).

Discipline
----------
- Bit-exact sympy.Rational throughout.
- Reuse JLO infrastructure from `compute_jlo_bicomplex.py` --- the
  CM cochain at degree n shares the simplex/moment expansion with
  the JLO cochain, modulo a $D$ insertion on the leftmost slot
  (factor D into a_0 -> D * a_0).
- No PSLQ, no floats, no transcendentals.
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
# Algebra helpers --- same as JLO bicomplex driver
# =====================================================================


def sector_idempotent(st: FockSpectralTriple, sector_idx: int) -> Matrix:
    """Return the sector idempotent e_s as a dim_H x dim_H diagonal matrix."""
    N = st.dim_H
    M = sp_zeros(N, N)
    for i in range(N):
        if st._state_to_sector[i] == sector_idx:
            M[i, i] = Integer(1)
    return M


def algebra_unit(st: FockSpectralTriple) -> Matrix:
    return sp_eye(st.dim_H)


def commutator(D: Matrix, a: Matrix) -> Matrix:
    return D * a - a * D


# =====================================================================
# CM-residue cochain coefficient extraction (t^m order)
# =====================================================================
#
# The CM cochain $\psi_n^{\eta}$ at t^m order is the SAME simplex/moment
# expansion as the JLO cochain, but with $\gamma$ replaced by $\gamma D$
# at the leftmost slot. Equivalently, we treat the leftmost slot as
# "(gamma * D) * a_0" instead of "gamma * a_0" --- this is the
# $\eta$-style trace inserting $D$ at the leftmost position.
#
# At order t^m, the CM cochain coefficient is:
#
#   psi_n^{eta, m}(a_0, ..., a_n)
#     = (-1)^m / (m + n)!
#       * sum_{(m_0, ..., m_n) | sum = m}
#         Tr(gamma D a_0 D^{2 m_0} [D, a_1] D^{2 m_1} ... [D, a_n] D^{2 m_n})
#
# The only structural difference from JLO is the (gamma D) insertion.
# =====================================================================


def _partitions(total: int, n_parts: int) -> List[Tuple[int, ...]]:
    if n_parts == 1:
        return [(total,)]
    out = []
    for m0 in range(total + 1):
        for sub in _partitions(total - m0, n_parts - 1):
            out.append((m0,) + sub)
    return out


def cm_psi_n_coeff(
    D: Matrix,
    D2: Matrix,
    gamma: Matrix,
    a_list: List[Matrix],
    m_target: int,
    use_gamma: bool,
) -> sp.Expr:
    """Compute the t^{m_target} coefficient of the CM-residue $\eta$-style
    cochain $\psi_n^{\eta}(a_0, ..., a_n; t)$ where n = len(a_list) - 1.

    Same simplex/moment expansion as JLO, but with leftmost insertion
    being $\gamma D \cdot a_0$ rather than $\gamma \cdot a_0$.

    use_gamma = True  --> "even" CM eta cochain (with gamma)
    use_gamma = False --> "odd" CM eta cochain (no gamma)
    """
    n = len(a_list) - 1
    a_0 = a_list[0]
    a_rest = a_list[1:]
    comm_list = [commutator(D, a) for a in a_rest]

    # D^{2k} powers up to k = m_target
    D2_powers: List[Matrix] = [sp_eye(D.shape[0])]
    for _ in range(m_target):
        D2_powers.append(D2_powers[-1] * D2)

    sign_gamma = gamma if use_gamma else sp_eye(D.shape[0])
    # CM-residue / eta insertion: leftmost factor is (gamma * D) instead of gamma
    leftmost = sign_gamma * D

    partitions = _partitions(m_target, n + 1)
    s_sum = Integer(0)
    for part in partitions:
        prod = leftmost * a_0 * D2_powers[part[0]]
        for j in range(n):
            prod = prod * comm_list[j] * D2_powers[part[j + 1]]
        s_sum = s_sum + prod.trace()
    prefac = Rational((-1) ** m_target, sp.factorial(m_target + n))
    return sp.simplify(prefac * s_sum)


# =====================================================================
# Hochschild b and Connes B on the CM-residue cochains
# (formulas IDENTICAL to JLO --- the bicomplex structure does NOT
#  depend on whether the cochain is JLO or CM-residue; it depends on
#  the parity-graded multilinear cochain tower having "(b + B) phi = 0".)
# =====================================================================


def hochschild_b_on_phi_n(
    phi_n_eval,
    inputs: List[Matrix],
    n: int,
) -> sp.Expr:
    r"""(b phi_n)(a_0, ..., a_{n+1})
        = sum_{i=0}^{n} (-1)^i phi_n(a_0, ..., a_i a_{i+1}, ..., a_{n+1})
        + (-1)^{n+1} phi_n(a_{n+1} a_0, a_1, ..., a_n)
    """
    if len(inputs) != n + 2:
        raise ValueError(f"need n+2 = {n+2} inputs, got {len(inputs)}")

    total = Integer(0)
    for i in range(n + 1):
        new_args = inputs[:i] + [inputs[i] * inputs[i + 1]] + inputs[i + 2:]
        assert len(new_args) == n + 1
        sign = (-1) ** i
        total = total + sign * phi_n_eval(*new_args)

    cyclic_args = [inputs[n + 1] * inputs[0]] + inputs[1: n + 1]
    assert len(cyclic_args) == n + 1
    sign = (-1) ** (n + 1)
    total = total + sign * phi_n_eval(*cyclic_args)

    return sp.simplify(total)


def connes_B_on_phi_n(
    phi_n_eval,
    inputs: List[Matrix],
    n: int,
    unit: Matrix,
) -> sp.Expr:
    """(B phi_n)(a_0, ..., a_{n-1}) for an n-cochain phi_n.

    Formula:
        (B phi_n)(a_0, ..., a_{n-1})
          = sum_{j=0}^{n-1} (-1)^{(n-1) j}
              phi_n(1, a_j, a_{j+1 mod n}, ..., a_{j-1 mod n})

    Output is an (n-1)-cochain taking n arguments.
    """
    if len(inputs) != n:
        raise ValueError(f"need n = {n} inputs, got {len(inputs)}")

    total = Integer(0)
    for j in range(n):
        permuted = [inputs[(j + k) % n] for k in range(n)]
        args = [unit] + permuted
        assert len(args) == n + 1
        sign = (-1) ** ((n - 1) * j)
        total = total + sign * phi_n_eval(*args)

    return sp.simplify(total)


# =====================================================================
# Bicomplex assembly: $(b + B)\psi$ on the CM-residue cochain tower
# =====================================================================


def make_psi_n_t_coeff(
    D: Matrix, D2: Matrix, gamma: Matrix, n: int, m: int, use_gamma: bool
):
    """Closure returning t^m coefficient of CM-residue eta cochain
    psi_n^eta at given (a_0, ..., a_n)."""
    def psi_n_eval(*args):
        a_list = list(args)
        assert len(a_list) == n + 1
        return cm_psi_n_coeff(D, D2, gamma, a_list, m, use_gamma)
    return psi_n_eval


def bb_residual_n_to_np1(
    psi_n_eval,
    psi_np2_eval,
    inputs: List[Matrix],
    n: int,
    unit: Matrix,
) -> Dict[str, sp.Expr]:
    """(b psi_n + B psi_{n+2})(a_0, ..., a_{n+1})."""
    b_term = hochschild_b_on_phi_n(psi_n_eval, inputs, n)
    B_term = connes_B_on_phi_n(psi_np2_eval, inputs, n + 2, unit)
    residual = sp.simplify(b_term + B_term)
    return {"b_psi_n": b_term, "B_psi_np2": B_term, "residual": residual}


# =====================================================================
# Main verification panel
# =====================================================================


def verify_cm_bicomplex_cocycle(
    st: FockSpectralTriple,
    D: Matrix,
    D2: Matrix,
    gamma: Matrix,
    M_max: int = 2,
    use_gamma: bool = True,
) -> Dict:
    """Verify the CM-residue $\eta$-cochain bicomplex cocycle condition:

        $(b\,\psi_0^{\eta} + B\,\psi_2^{\eta})(a_0, a_1) = 0$

    at each $t$-order $m = 0, ..., M_{\max}$ on every idempotent pair
    (e_s, e_t), plus unit-mixed pairs.
    """
    N = st.dim_H
    unit = sp_eye(N)
    n_sectors = st.n_sectors
    idempotents = [sector_idempotent(st, s) for s in range(n_sectors)]

    panel_residuals = {}
    all_zero_at_order = {m: True for m in range(M_max + 1)}

    test_inputs_label_pairs: List[Tuple[Tuple[int, int], Tuple[Matrix, Matrix]]] = []
    for s in range(n_sectors):
        for t_idx in range(n_sectors):
            test_inputs_label_pairs.append(((s, t_idx), (idempotents[s], idempotents[t_idx])))
    for s in range(n_sectors):
        test_inputs_label_pairs.append(((-1, s), (unit, idempotents[s])))
        test_inputs_label_pairs.append(((s, -1), (idempotents[s], unit)))
    test_inputs_label_pairs.append(((-1, -1), (unit, unit)))

    for m in range(M_max + 1):
        psi_0_eval = make_psi_n_t_coeff(D, D2, gamma, n=0, m=m, use_gamma=use_gamma)
        psi_2_eval = make_psi_n_t_coeff(D, D2, gamma, n=2, m=m, use_gamma=use_gamma)

        order_residuals = {}
        for ((s, t_idx), (a_0, a_1)) in test_inputs_label_pairs:
            res = bb_residual_n_to_np1(
                psi_0_eval, psi_2_eval, [a_0, a_1], n=0, unit=unit
            )
            label = f"({s},{t_idx})"
            order_residuals[label] = {
                "b_psi_0": str(res["b_psi_n"]),
                "B_psi_2": str(res["B_psi_np2"]),
                "residual": str(res["residual"]),
            }
            if res["residual"] != Integer(0):
                all_zero_at_order[m] = False
        panel_residuals[f"t^{m}"] = order_residuals

    return {
        "panel_residuals": panel_residuals,
        "all_zero_at_each_order": all_zero_at_order,
        "M_max": M_max,
    }


# =====================================================================
# HP class identification for the CM-residue cochain
# =====================================================================


def identify_cm_hp_class(
    st: FockSpectralTriple,
    D: Matrix,
    D2: Matrix,
    gamma: Matrix,
) -> Dict:
    """Identify the leading periodic-cyclic-cohomology class for the
    CM-residue $\eta$-cochain.

    On commutative A = C^5, K_0(A) = Z^5 with basis {[e_s]}.
    The CM-eta-cochain pairing at leading t^0 is:

       psi_0^{eta, even}(e_s) = Tr(gamma D e_s)

    Sector resolution of the eta vector:
       eta_s := Tr(gamma D e_s) on the sector idempotent e_s.

    Cross-check vs Sub-Sprint 1:
       sum_s eta_s = Tr(gamma D) on the full Hilbert space.
       On the diagonal Lambda, Tr(gamma Lambda) = sum_i chi_i^2 (n_i + 1/2)
                                              = sum_i (n_i + 1/2) = 36
                                              at n_max = 2.
       On the full D = Lambda + kappa A:
       Tr(gamma D) = Tr(gamma Lambda) + kappa Tr(gamma A)
                   = 36 + kappa Tr(gamma A).
       The off-diagonal A connects sectors of opposite chirality
       (Delta l = +/-1 flips chirality on Dirac states), so
       Tr(gamma A) should typically be 0 by chirality block structure.

    HP class identification:
      Leading vector [psi_0^{eta, even}(e_s)]_{s = 0..4} in Q^5.
      sum_s = sub-sprint-1 ground truth = 36 (or kappa-corrected for D).

    Compare to JLO HP^even class: (+2, -2, +2, +2, -4) sum = 0 = chi(S^3).
    The CM eta class has DIFFERENT sector resolution and sums to
    a NON-zero value (the M3 leading coefficient itself).
    """
    N = st.dim_H
    n_sectors = st.n_sectors
    unit = sp_eye(N)
    Lam = st.diagonal_part

    cm_eta_per_sector_LambdaOnly = []
    cm_eta_per_sector_fullD = []
    for s in range(n_sectors):
        e_s = sector_idempotent(st, s)
        # Lambda-only (cleanest match to CH-1 ground truth)
        val_Lam = sp.simplify((gamma * Lam * e_s).trace())
        cm_eta_per_sector_LambdaOnly.append(val_Lam)
        # Full D (= Lambda + kappa A)
        val_D = sp.simplify((gamma * D * e_s).trace())
        cm_eta_per_sector_fullD.append(val_D)

    sum_Lam = sum(cm_eta_per_sector_LambdaOnly, Integer(0))
    sum_D = sum(cm_eta_per_sector_fullD, Integer(0))

    # Cross-check: JLO HP^even class (from Sub-Sprint 2c)
    jlo_hp_even_class = []
    for s in range(n_sectors):
        e_s = sector_idempotent(st, s)
        val = sp.simplify((gamma * e_s).trace())
        jlo_hp_even_class.append(int(val))

    # ODD CM eta = Tr(D e_s) per sector
    cm_eta_odd_per_sector_LambdaOnly = []
    for s in range(n_sectors):
        e_s = sector_idempotent(st, s)
        val_Lam = sp.simplify((Lam * e_s).trace())
        cm_eta_odd_per_sector_LambdaOnly.append(val_Lam)

    return {
        "cm_eta_even_LambdaOnly_per_sector": [
            str(v) for v in cm_eta_per_sector_LambdaOnly
        ],
        "cm_eta_even_LambdaOnly_sum": str(sum_Lam),
        "cm_eta_even_fullD_per_sector": [
            str(v) for v in cm_eta_per_sector_fullD
        ],
        "cm_eta_even_fullD_sum": str(sum_D),
        "cm_eta_odd_LambdaOnly_per_sector": [
            str(v) for v in cm_eta_odd_per_sector_LambdaOnly
        ],
        "cm_eta_odd_LambdaOnly_sum": str(sum(cm_eta_odd_per_sector_LambdaOnly, Integer(0))),
        "jlo_hp_even_class_cross_check_Q5": jlo_hp_even_class,
        "jlo_hp_even_class_sum": int(sum(jlo_hp_even_class)),
        "_note": (
            "JLO HP^even class sums to 0 (McKean-Singer index = chi(S^3) = 0). "
            "CM-eta-even class sums to 36 (Lambda-only) "
            "= Sub-Sprint 1 M3 ground truth Tr(gamma Lambda)."
        ),
    }


# =====================================================================
# Sub-Sprint 1 cross-check: confirm Tr(gamma D e_s) ground truth
# =====================================================================


def cross_check_subsprint1_ground_truth(
    st: FockSpectralTriple,
    D: Matrix,
    gamma: Matrix,
) -> Dict:
    """Confirm: M_3(n_max=2) = Tr(gamma Lambda) = 36 (Stage 1 Sub-Sprint 1).

    Also confirm: on the FULL D = Lambda + kappa A, what does Tr(gamma D)
    give? The kappa A part should integrate out by chirality block structure
    (A is parity-respecting in (n,l), but chirality of Dirac states depends
    on kappa-sign which flips with l-parity --- net effect on Tr(gamma A)
    needs explicit check).
    """
    N = st.dim_H
    Lam = st.diagonal_part
    tr_g_Lam = sp.simplify((gamma * Lam).trace())
    tr_g_D = sp.simplify((gamma * D).trace())

    # Sub-Sprint 1 reported eta sums:
    # t^0: 36 (Lambda-only)
    # t^1: -201
    # t^2: 4809/8
    # Let's reproduce t^1: -Tr(gamma D^3) on Lambda would give -(36 * 3 * (1/2)^2)?
    # Actually t^1 of Tr(gamma D e^{-tD^2}) = -Tr(gamma D^3)
    Lam3 = Lam * Lam * Lam
    tr_g_Lam3 = sp.simplify((gamma * Lam3).trace())
    t1_from_eta = -tr_g_Lam3  # matches Sub-Sprint 1 t^1 = -201

    return {
        "M3_Tr_gamma_Lambda_t0_n_max_2": str(tr_g_Lam),
        "M3_Tr_gamma_D_t0_n_max_2": str(tr_g_D),
        "M3_eta_t1_minus_Tr_gamma_Lambda_3": str(t1_from_eta),
        "_ground_truth_match_t0": (str(tr_g_Lam) == "36"),
        "_ground_truth_match_t1_LambdaOnly": (str(t1_from_eta) == "-201"),
    }


# =====================================================================
# Cross-cutoff sanity check at n_max = 3
# =====================================================================


def cross_check_n_max_3(M_max: int = 0) -> Dict:
    """Quick cross-check: verify M3(3) = 120 (Sub-Sprint 2a closed form)
    on the truncated CH triple at n_max = 3, and verify degree-1 cocycle
    condition on a SMALL panel (computational cost is much higher at
    n_max = 3: dim_H = 40, 14 sectors).

    Only minimal coverage: full diagonal panel at t^0 only.
    """
    st3 = FockSpectralTriple(n_max=3)
    D3 = st3.dirac_operator
    D3sq = D3 * D3
    gamma3 = st3.grading
    Lam3 = st3.diagonal_part
    N3 = st3.dim_H
    n_sec3 = st3.n_sectors
    unit3 = sp_eye(N3)

    # Cross-check: M3(3) = Tr(gamma Lambda) on n_max=3 truncation
    tr_g_Lam3 = sp.simplify((gamma3 * Lam3).trace())

    # Small panel: only diagonal idempotent pairs (e_s, e_s) at t^0
    # Then cocycle condition (b psi_0 + B psi_2)(e_s, e_s) at t^0.
    diag_residuals_t0 = {}
    if M_max >= 0:
        psi_0_eval = make_psi_n_t_coeff(D3, D3sq, gamma3, n=0, m=0, use_gamma=True)
        psi_2_eval = make_psi_n_t_coeff(D3, D3sq, gamma3, n=2, m=0, use_gamma=True)
        for s in range(n_sec3):
            e_s = sector_idempotent(st3, s)
            res = bb_residual_n_to_np1(
                psi_0_eval, psi_2_eval, [e_s, e_s], n=0, unit=unit3
            )
            diag_residuals_t0[f"({s},{s})"] = str(res["residual"])

    return {
        "n_max": 3,
        "dim_H": N3,
        "n_sectors": n_sec3,
        "M3_Tr_gamma_Lambda_3": str(tr_g_Lam3),
        "ground_truth_M3_3_eq_120": (str(tr_g_Lam3) == "120"),
        "diag_panel_residuals_t0_(degree_1_eta_cocycle)": diag_residuals_t0,
        "all_diag_zero": all(v == "0" for v in diag_residuals_t0.values()),
    }


# =====================================================================
# Main
# =====================================================================


def main() -> None:
    print("=" * 70)
    print("Sprint Q5'-Stage1-CM-Bicomplex")
    print("Connes-Moscovici residue $(b, B)$ bicomplex on the truncated CH triple")
    print("Housing M3 = Tr(gamma D e^{-tD^2})|_{t^0} = 36 at n_max = 2")
    print("=" * 70)

    t_global = time.time()

    print("\n[1] Building truncated CH spectral triple at n_max = 2...")
    st = FockSpectralTriple(n_max=2)
    D = st.dirac_operator
    gamma = st.grading
    Lam = st.diagonal_part
    N = st.dim_H
    print(f"    dim_H = {N}, n_sectors = {st.n_sectors}, sectors = {st.sectors}")
    print(f"    kappa = {st._kappa}")
    D2 = D * D

    print("\n[2] Sub-Sprint 1 ground truth cross-check...")
    gt = cross_check_subsprint1_ground_truth(st, D, gamma)
    print(f"    Tr(gamma Lambda) = {gt['M3_Tr_gamma_Lambda_t0_n_max_2']}  (expected 36)")
    print(f"    Tr(gamma D)      = {gt['M3_Tr_gamma_D_t0_n_max_2']}")
    print(f"    M3 eta t^1 (Lambda-only) = {gt['M3_eta_t1_minus_Tr_gamma_Lambda_3']}  (expected -201)")
    print(f"    Match t^0: {gt['_ground_truth_match_t0']}")
    print(f"    Match t^1 (Lambda): {gt['_ground_truth_match_t1_LambdaOnly']}")

    print("\n[3] Verifying EVEN CM-eta bicomplex cocycle:")
    print("    b psi_0^eta_even + B psi_2^eta_even = 0")
    print("    on all idempotent pairs at t-orders m = 0, 1, 2...")
    t0 = time.time()
    even_check = verify_cm_bicomplex_cocycle(st, D, D2, gamma, M_max=2, use_gamma=True)
    elapsed = time.time() - t0
    print(f"    done in {elapsed:.1f}s")
    print(f"    EVEN cocycle all-zero per order: {even_check['all_zero_at_each_order']}")
    nonzero_even = []
    for order_label, order_dict in even_check["panel_residuals"].items():
        for pair_label, vals in order_dict.items():
            if vals["residual"] != "0":
                nonzero_even.append((order_label, pair_label, vals["residual"]))
    print(f"    EVEN: {len(nonzero_even)} non-zero residual(s) across panel")
    if nonzero_even and len(nonzero_even) <= 10:
        for nz in nonzero_even:
            print(f"      {nz[0]} pair {nz[1]}: residual = {nz[2]}")

    print("\n[4] Verifying ODD CM-eta bicomplex cocycle:")
    print("    b psi_0^eta_odd + B psi_2^eta_odd = 0")
    t0 = time.time()
    odd_check = verify_cm_bicomplex_cocycle(st, D, D2, gamma, M_max=2, use_gamma=False)
    elapsed = time.time() - t0
    print(f"    done in {elapsed:.1f}s")
    print(f"    ODD cocycle all-zero per order: {odd_check['all_zero_at_each_order']}")
    nonzero_odd = []
    for order_label, order_dict in odd_check["panel_residuals"].items():
        for pair_label, vals in order_dict.items():
            if vals["residual"] != "0":
                nonzero_odd.append((order_label, pair_label, vals["residual"]))
    print(f"    ODD: {len(nonzero_odd)} non-zero residual(s) across panel")
    if nonzero_odd and len(nonzero_odd) <= 10:
        for nz in nonzero_odd:
            print(f"      {nz[0]} pair {nz[1]}: residual = {nz[2]}")

    print("\n[5] Identifying CM-residue HP class...")
    t0 = time.time()
    hp_data = identify_cm_hp_class(st, D, D2, gamma)
    print(f"    done in {time.time() - t0:.1f}s")
    print(f"    CM eta-even class (Lambda-only) per sector: {hp_data['cm_eta_even_LambdaOnly_per_sector']}")
    print(f"    CM eta-even class (Lambda-only) sum: {hp_data['cm_eta_even_LambdaOnly_sum']} (expected 36)")
    print(f"    CM eta-even class (full D)    per sector: {hp_data['cm_eta_even_fullD_per_sector']}")
    print(f"    CM eta-even class (full D)    sum: {hp_data['cm_eta_even_fullD_sum']}")
    print(f"    CM eta-odd  class (Lambda)    per sector: {hp_data['cm_eta_odd_LambdaOnly_per_sector']}")
    print(f"    CM eta-odd  class (Lambda)    sum: {hp_data['cm_eta_odd_LambdaOnly_sum']}")
    print(f"    Cross-check vs JLO HP^even class (Q5): {hp_data['jlo_hp_even_class_cross_check_Q5']}")
    print(f"    JLO HP^even class sum: {hp_data['jlo_hp_even_class_sum']} (expected 0 = chi(S^3))")

    print("\n[6] Cross-cutoff sanity at n_max = 3 (small panel)...")
    t0 = time.time()
    try:
        n3_check = cross_check_n_max_3()
        elapsed = time.time() - t0
        print(f"    done in {elapsed:.1f}s")
        print(f"    dim_H_3 = {n3_check['dim_H']}, n_sectors_3 = {n3_check['n_sectors']}")
        print(f"    M3(3) = Tr(gamma Lambda)_3 = {n3_check['M3_Tr_gamma_Lambda_3']}  (expected 120)")
        print(f"    Match: {n3_check['ground_truth_M3_3_eq_120']}")
        print(f"    Diagonal panel cocycle (b psi_0 + B psi_2)(e_s, e_s) at t^0 all-zero: {n3_check['all_diag_zero']}")
    except Exception as exc:
        elapsed = time.time() - t0
        print(f"    SKIPPED ({elapsed:.1f}s): {exc}")
        n3_check = {"error": str(exc), "skipped": True}

    # ===== Output =====
    out_dir = Path(__file__).parent / "data"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "sprint_q5p_cm_bicomplex.json"
    payload = {
        "sprint": "Q5'-Stage1-CM-Bicomplex",
        "date": "2026-06-05",
        "n_max": 2,
        "dim_H": N,
        "n_sectors": st.n_sectors,
        "kappa": str(st._kappa),
        "subsprint1_ground_truth_cross_check": gt,
        "even_cm_eta_bicomplex_check": {
            "all_zero_at_each_order": even_check["all_zero_at_each_order"],
            "nonzero_residuals": nonzero_even,
            "panel_residuals_t0_full": even_check["panel_residuals"]["t^0"],
        },
        "odd_cm_eta_bicomplex_check": {
            "all_zero_at_each_order": odd_check["all_zero_at_each_order"],
            "nonzero_residuals": nonzero_odd,
            "panel_residuals_t0_full": odd_check["panel_residuals"]["t^0"],
        },
        "hp_class_identification": hp_data,
        "n_max_3_cross_check": n3_check,
        "total_wall_time_s": time.time() - t_global,
    }
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"\nWrote data: {out_path}")

    print("\n" + "=" * 70)
    print(f"Total wall time: {time.time() - t_global:.1f}s")
    print("=" * 70)


if __name__ == "__main__":
    main()
