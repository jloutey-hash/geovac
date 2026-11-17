# Dirac-seminorm item vs the SPECTRAL-WINDOW-COMPRESSION edge-crossing term (2026-10-31)

**Verdict: NOT-A-NEGATIVE.** On the true Camporesi-Higuchi window there is no concrete failing quantity
attached to the Dirac seminorm. The one measured fact is a success (the kernel condition holds), and the
open item is an unproved estimate. It cannot be the cluster's third wall. Driver:
`debug/dirac_seminorm_edge_crossing.py` (ALL CHECKS PASS, exit 0).

## 1. What the item is, from the owning papers

- Paper 38, `rem:dirac_degeneracy` (and the abstract, and the opening of `sec:named_gaps`): on the true window
  the kernel of a -> [D, a] on O_{n_max} is **exactly the scalars** at n_max = 1, 2, 3 (dim O = 5, 30, 91), so
  the Dirac-commutator seminorm *is* a norm on O/C1 there. That is the kernel condition a compact quantum
  metric space requires: it **holds**. The kernel that was too large (10/14, 26/55) belongs to the doubled
  Fock spin-orbital model, which the paper says "is not a property of the truncated spectral triple".
- Paper 38 lines 330-339, and Paper 32 item (iv) near line 2900, both say it in so many words: the
  translation seminorm "is used because its compression and lifted-state estimates close, **not because the
  Dirac seminorm fails**; whether the convergence can be proved in the Dirac seminorm is [OPEN]".
- No test, memo or register entry records a measured failure of a Dirac-seminorm quantity on the true window.
  The history (`debug/sprint_p38_g1g2_phaseA_memo.md`) shows the "degeneracy" lived on the doubled model only.

So the register's 2026-10-30 wording "its kernel is the scalars at the tested cutoffs" describes a property
that holds, and the item is an open proof, as the register itself already flags ("OPEN item ... not a filed
negative").

## 2. Why the edge-crossing sum cannot appear in it (structural, checked)

The two walls are product identities: the order-one defect [[D_W,a_W],b_W] and the gauge defect (its
a = b = chi case with an f' weight). Each has an intermediate index, and compression removes the k outside W,
leaving Sum_{k not in W} with a Dirac-gap weight.

The Dirac seminorm is **linear** in a. Because P commutes with D, [D_W, a_W] = P[D, a]P exactly: a single
commutator with diagonal D, ([D,a])_il = (lambda_i - lambda_l) a_il, has no intermediate index and so no sum
over k. The compressed seminorm has no defect against the continuum commutator; it is the continuum
commutator, compressed.

| n_max | dim W | linear defect | order-one defect | order-one minus edge sum | rotated-window control |
|:-----:|:-----:|:-------------:|:----------------:|:------------------------:|:----------------------:|
| 1 | 4  | 0 (exact) | 0.750 | 4.4e-16 | 0.182 |
| 2 | 16 | 0 (exact) | 0.778 | 8.9e-16 | 0.210 |
| 3 | 40 | 0 (exact) | 0.812 | 1.1e-15 | 0.223 |

Functions tested: the four spin-1/2 coefficients D^{1/2}_{MM'} and their ten products (band 1), ambient
cutoff 2j <= 2 n_max + 4. The order-one column reproduces the existing wall (guard
`tests/test_paper32_window_order_one.py`) as a positive control of the setup.

**Control (non-tautology).** The zero is the spectral-projection property, not an identity of every
compression: replacing P by P' = U P U^dag, with U a 0.3 rad rotation mixing the top window state with the
lowest outside state (||[D, P']|| = 0.28), gives a linear defect of 0.18 to 0.22. (|D| in place of D is not a
useful control here, because P commutes with |D| too and the linear identity survives it; the rotated window
is the control that can fail, and does.)

## 3. Where the same matrix elements do sit: a norm loss in the safe direction

The crossing block Q[D,a]P has entries (lambda_k - lambda_i) a_ki, i in W, k not in W: the same matrix
elements, with the Dirac-gap weight, but at **first order** and as a discarded block, not as the residue of a
broken identity. Its effect is that compression is a contraction:

||[D_W, a_W]|| <= ||[D, a] P|| <= ||[D, a]|| = Lip(a).

Measured ratio ||[D_W,a_W]|| / ||[D,a]P|| for the four generators: 0.949, 0.992, 0.995 at n_max = 1, 2, 3
(rising to 1); crossing block / column norm 0.447, 0.458, 0.464. The forward estimate a convergence proof
needs from compression therefore holds with no boundary term. What is open is the reverse (lifted-state)
estimate, controlling the continuum Lipschitz constant of a lift by ||[D_W, T]|| for T in O_W; for the
translation seminorm that step is the Fejer smoothing identity, and for the Dirac seminorm it is simply not
proved. No measurement of it failing exists, and this sprint did not construct one.

## 4. Consequence for the register

- The Dirac-seminorm item is not a third negative of SPECTRAL-WINDOW-COMPRESSION, and not a wall. It fits the
  sharpened mechanism statement on the positive side: compression keeps what commutes with D, and [D, .] is
  linear, so it survives compression unbroken. Only the product identities leak.
- The cluster stays at two walls, FORMING. A third independent negative has to be another *product* identity
  that compression breaks (an order condition, a derivation/Leibniz property, a multiplicative-domain
  property), not a linear functional of a.
- If the reverse estimate is ever attempted and fails, that would be a new candidate; whether it would be the
  same boundary term is not predicted here and should not be assumed.

## Files

- `debug/dirac_seminorm_edge_crossing.py` (driver; reuses the S^3 construction in
  `tests/test_paper32_window_gauge_breaking.py`)
- `tests/test_p38_true_ch_window.py::test_dirac_commutator_kernel_is_the_scalars` re-run: 3 passed.
