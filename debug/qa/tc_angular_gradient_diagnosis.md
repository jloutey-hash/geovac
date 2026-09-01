# TC angular-gradient assembly: diagnosis (not a fix)

**2026-08-31.** The last open item before group4 is reviewable. The bug is
now located precisely. It is **not** repaired — the repair is a re-derivation,
not a patch, and an attempted patch is recorded below as a dead end so nobody
walks it twice.

Baseline unchanged: He n_max=2 radial-only gives 107 nonzeros with 0 L_z
violations; turning the angular term on adds 26 entries, **26/26 violating**.

---

## 1. What is actually wrong

### (a) Two conventions, mixed

`_gaunt_integral` is the **complex**-harmonic Gaunt integral: it conjugates its
third argument, carries the `(-1)^{m3}` phase, and selects on `m1 + m2 = m3`.

The assembly around it reasons with **real** harmonics — the comments say
"for real Y, `Y* = Y`" and drop conjugates on that basis, in two places: the
Neumann pairing `Y*_{KQ}(O1) Y_{KQ}(O2)`, and the gradient identity's
`Y*_{1,q}(O1) Y_{1,q}(O2)`.

Dropping a conjugate flips the sign of that factor's m. Doing it on both sides
makes the two electrons' m-shifts add instead of cancel:
`ma - mc = +(q+Q)` and `mb - md = +(q+Q)`, summing to `2(q+Q)` rather than
zero. That is the L_z violation, in one line.

### (b) The structural error: `a` is delta-selected, so the Neumann sum collapses

This is the load-bearing one. The block first selects the bra orbital by

```python
a_indices = lm_to_idx.get((L_eff, m_eff), [])     # ma == m_eff, exactly
```

which comes from an early assumption — written in the comments — that the O1
integral is a simple delta. The same block then works out, correctly and at
length, that O1 is a *triple* product because the Neumann `Y_{KQ}(O1)` also
lives there ("This is NOT simply a delta function! I was wrong earlier").

**The realisation never reached the code.** The delta-selection stayed, so
`ma == m_eff` identically, and therefore

```python
Q_val = ma - m_eff     # == 0, always
```

The sum over Neumann orders `Q` is silently truncated to the single term
`Q = 0`. The expansion is gutted, and no test notices because the only
observable anyone checked was the L_z violation count.

### (c) A second, unpatched block with the same lineage

The emitted entries include rows with `lc = 0` — an s orbital, whose angular
gradient `_angular_gradient_coefficients` returns as empty. Those cannot come
from the electron-1 block at all. They come from the separate electron-2
gradient block further down the same function, which repeats the construction
with the roles swapped. **Any fix must cover both**; patching one and measuring
the total is how a partial fix looks like a regression.

### (d) Dead code that encodes the confusion

- `m_eff = mc - q` (top of `_angular_gradient_coefficients`) is superseded by
  `m_eff_correct = q + mc` — but the stale value is still used by the
  `if abs(m_eff) > L: continue` guard a few lines down, which therefore drops
  valid terms on a quantity the function has already abandoned.
- `Q_val = m_eff - ma` is computed, used in a guard, and then overwritten ~40
  lines later by `Q_val = ma - m_eff  # CORRECTED`. (The guard is
  `abs()`-symmetric, so this one is harmless — but it reads as two competing
  derivations left in place.)
- The module carries the unfinished derivation verbatim: *"Wait -- I need to
  reconsider"*, *"I was wrong earlier"*, *"That can't be right in general!"*.
  It was committed mid-derivation.

---

## 2. Dead end: restoring the two conjugates does NOT fix it

Tried, measured, reverted. Restoring both conjugates —
`m_eff = mc - q` with `C_L = (-1)^q G(1,-q; lc,mc; L,m_eff)`, and
`Q = m_eff - ma` with `gaunt_O1 = (-1)^Q G(L,m_eff; K,-Q; la,ma)` — makes the
m-algebra close **on paper** (`ma - mc = -(q+Q)`, `mb - md = +(q+Q)`).

Measured result: added entries **26 → 46, still 46/46 violating.** Worse.

Why: the delta-selection in (b) sits *upstream* of both conjugates. With
`ma == m_eff` forced, `Q ≡ 0`, so the conjugate correction on the O1 Neumann
factor has nothing to act on, and the electron-2 block (c) is untouched.
Fixing the signs without fixing the selection changes which wrong entries are
emitted, not whether they are wrong.

---

## 3. What a real fix requires

1. Remove the delta-selection of `a`. The O1 integral is
   `int Y*_{la,ma} Y_{L,m_eff} Y*_{K,Q} dO1` — a genuine triple product — so
   `a` ranges over all states satisfying the Gaunt selection, and `Q` becomes
   a real summation index rather than an identically-zero quantity.
2. Fix the conventions once, explicitly: either move the whole assembly to
   complex harmonics (matching `_gaunt_integral`) or write a real-harmonic
   Gaunt and use it consistently. Do not mix.
3. Apply both to the electron-2 block as well.
4. **Validate against an independent route** — direct numerical quadrature of
   the two-electron angular integral on a product Lebedev grid, at a handful
   of `(la ma | lb mb | lc mc | ld md)` tuples. L_z conservation is a
   necessary condition, not a sufficient one: an assembly can conserve m and
   still have wrong magnitudes, and this operator has now been wrong twice in
   ways that a selection-rule check alone would not separate.

That is a sprint with its own validation gate, not a patch.

---

## 4. The interim position is already correct

Nothing needs to change while this is open:

- Paper 14's angular-gradient resource multipliers (2.66×/2.31×/2.67×, total
  4.49×) are **withdrawn, not re-measured** — correct, since no valid input
  survives to rebuild them.
- The paper's *verdict* — angular gradient is a net negative for quantum
  efficiency — stands and is over-determined by the radial-only result.
- The violations are deliberately **not filtered**. They are the only
  detector, and masking them would hide whatever the same mis-assignment does
  to entries that happen to land in a valid M_L sector. That reasoning is
  reinforced by (b): the operator is not merely mis-indexed, it is also
  truncated, and the violation count would not reveal the truncation.

**Recommendation for the sweep:** this does not block `/qa group4`. The
group4 papers make no live claim that depends on the angular-gradient
operator; the one claim that did is withdrawn, and the withdrawal is
registered. Fixing it is worth doing on its own terms — the cusp/TC arc would
regain an axis — but it is not a certification prerequisite.

Reproduce the baseline:

```
python -c "
from geovac.tc_integrals import compute_tc_integrals_block
s=[(1,0,0),(2,0,0),(2,1,-1),(2,1,0),(2,1,1)]
r=compute_tc_integrals_block(2.0,s,n_grid=400,include_angular=False)
a=compute_tc_integrals_block(2.0,s,n_grid=400,include_angular=True)
m=lambda i:s[i][2]; add=set(a)-set(r)
print(len(r), len(add), sum(1 for k in add if m(k[0])+m(k[1])!=m(k[2])+m(k[3])))"
```

Existing backing: `tests/test_tc_angular.py`, incl.
`test_composed_angular_adds_only_lz_violating_entries`.
