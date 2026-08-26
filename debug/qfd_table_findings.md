# QFD table extension -- certified minimal diatomics beyond H2 and LiH

**Task.** Extend the quadrature-free diatomic (QFD) build from the two systems
it was demonstrated on (H2, LiH) to a certified minimal-diatomic TABLE, and add
the certified rows to the reference artefact.

**Machinery (unchanged, used as promoted):** `geovac/qfd_core.py` +
`geovac/qfd_assemble.py` (v4.105.0; backing tests `tests/test_paper58_qfd.py`).
**Driver:** `debug/qfd_table_ext.py`.
**Artefact module:** `benchmarks/certified_reference/entries_qfd_diatomic.py`
(new category `qfd_diatomic`).
**Raw records:** `debug/data/qfd_table_ext_certified.json`,
`debug/data/qfd_ext_light.log`, `debug/data/qfd_ext_behplus.log`.

---

## 1. Verdict

**GO.** Five certifications across four systems, every one at or above the
25-digit gate, all closed-form-vs-quadrature validations inside tolerance, and
the reference artefact regenerated with the test suite green.

| system | R (bohr) | N_e | basis (hydrogenic, a = Z_orb/n) | tau series | digits | E_total (Ha) |
|:--|--:|--:|:--|:--|--:|:--|
| H2+ | 2.0 | 1 | 1s_A a=1, 1s_B a=1 | terminates | 40 | `-0.5537714953184827365067633613319649091848` |
| H2+ | 1.4 | 1 | 1s_A a=1, 1s_B a=1 | terminates | 40 | `-0.4713457017330697952317131946245614997863` |
| HeH+ | 1.46 | 2 | He 1s a=2, H 1s a=1 | infinite, tau_max 16 | 40 | `-2.895950902302325174310463610145323586925` |
| He2^2+ | 1.3 | 2 | 1s a=2, 1s a=2 | terminates | 40 | `-3.586503080128554308997786255254462141274` |
| BeH+ | 2.5 | 4 | Be 1s a=4, Be 2s a=2, H 1s a=1 | infinite, tau_max 22/18/16 | 31 | `-14.69882770595311489573042272632` |

LiH (30 digits) and H2 (84 digits) were already certified and are unchanged.

---

## 2. What was certified, and what was not

Every number above is the exact ground-state energy of a **completely
specified** model: s-type hydrogenic orbitals with decay rate `a = Z_orbital/n`,
a stated internuclear separation, and full CI in that basis. Anyone who
implements the same model must reproduce these digits. That is the claim.

The claim is **not** accuracy. These are two- and three-function s-only bases
with unoptimised exponents, and they are far from the exact energies. The
comparisons kept in the artefact are ordering checks, not accuracy claims:

* H2+ at R = 2.0: the model gives -0.55377 Ha. The exact Born-Oppenheimer value
  at that separation is -0.6026 Ha, and the zeta = 1 LCAO curve this model *is*
  has its own minimum at R = 2.49 bohr, E = -0.5648 Ha. So the model number is
  where it should be, and says nothing about H2+ accuracy.
* HeH+ at R = 1.46: -2.8960 Ha against an exact value near -2.978 Ha. Above the
  exact energy, as a variational minimal basis must be.
* He2^2+ at R = 1.3: -3.5865 Ha, which is **above** the He+ + He+ dissociation
  limit that this same basis reproduces exactly (-4 Ha = 2 x -Z^2/2 at
  Z_orbital = 2). The two-function model does not bind the real
  barrier-protected minimum. A fact about the basis, not about the number.

---

## 3. The structural axis: which systems have a truncation at all

The exchange class is a Neumann (Legendre) expansion in an index `tau`. It
**terminates** exactly when the two centres of a charge density carry the same
orbital exponent -- the Phase 0-e criterion `q = (alpha - beta) R / 2 = 0`.
That single condition splits the table in two, and it is why the digit counts
have different *kinds* of ceiling.

**Homonuclear at equal exponent (H2+, He2^2+; and the companion H2).** The tau
sum is finite. This was verified, not assumed: `exchange_closed_form` was asked
for every term through tau = 8 and each one above the last nonzero term is a
**symbolic zero**.

| system | last nonzero tau | terms tau+1 .. 8 |
|:--|--:|:--|
| H2+ (R = 2.0) | 2 | all symbolic zero |
| H2+ (R = 1.4) | 2 | all symbolic zero |
| He2^2+ (R = 1.3) | 2 | all symbolic zero |

There is therefore **no truncation error at all** in these rows, and the only
ceiling is arithmetic: the whole pipeline was re-run at working precisions 40
and 60, agreeing to ~1e-65 relative, and the claim is capped at 40 (a claim may
never exceed the weaker run).

**Heteronuclear (HeH+, BeH+; and the companion LiH).** The sum is infinite and
each quartet is truncated at its own `tau_max`. The per-tau magnitudes were
measured to decrease monotonically through the tail in every case, so the
geometric bound `|a_taumax| * r/(1-r)` with `r` the last observed ratio
dominates the true tail. Every amplification factor below was **recomputed for
its own system** -- two-particle density-matrix weight `N(N-1)/2` times
`||S^{-1/2}||^4 = lambda_min(S)^{-2}`, rounded up to the next power of two at
>= 4x margin -- and none is inherited from LiH.

| system | quartet (n labels) | tau_max | last ratio r | integral tail bound |
|:--|:--|--:|--:|:--|
| HeH+ | (1,1\|1,1) | 16 | 6.01e-4 | 1.19e-42 |
| BeH+ | (1,1\|1,1) | 22 | 7.95e-3 | 2.64e-35 |
| BeH+ | (1,1\|2,1) | 18 | 4.36e-3 | 6.46e-33 |
| BeH+ | (2,1\|2,1) | 16 | 2.06e-3 | 3.79e-33 |

| system | lambda_min(S) | raw amplification | used | energy tail bound (Ha) | digits from tail |
|:--|--:|--:|--:|:--|--:|
| HeH+ | 0.515423 | 3.764 | 16 | 1.90e-41 | 41 |
| BeH+ | 0.579724 | 17.853 | 128 | 1.31e-30 | 31 |

The exchange accumulator's own arithmetic was separately checked by re-running
it at working precisions 30 and 50 at fixed tau: HeH+ 1.23e-44 (43 digits),
BeH+ 1.49e-43 (40 digits).

Each heteronuclear claim is then
`min(linear algebra, tau tail, exchange arithmetic)`.

---

## 4. Independent-route validation

Nothing on the production path calls a quadrature routine, so every closed form
was falsified against `debug/qfd_quad.py`, which shares no code with it. Worst
deviations per system (absolute, mpmath at dps 20 unless noted):

| system | one-electron (S / T / V) | one-centre, (AA\|BB), hybrid | exchange (tau-matched) |
|:--|:--|:--|:--|
| H2+ (2.0) | 5.8e-21 / 1.3e-23 / 1.7e-22 | 4.2e-21 / 1.9e-22 / 1.7e-22 | 1.9e-22 |
| H2+ (1.4) | 1.0e-21 / 9.7e-23 / 1.0e-21 | 4.2e-21 / 3.4e-22 / 3.4e-22 | 1.1e-21 |
| HeH+ | 1.7e-22 / 3.9e-22 / 5.8e-22 | 1.7e-21 / 4.7e-22 / 6.3e-22 | 3.4e-21 |
| He2^2+ | 4.9e-22 / 6.3e-22 / 1.0e-21 | 1.7e-21 / 7.8e-22 / 2.8e-22 | 1.7e-22 |
| BeH+ | 1.1e-23 / 2.0e-24 / 1.5e-22 | 7.0e-21 / 1.6e-22 / 4.3e-23 | 4.6e-20 |

All are 10 to 13 orders of magnitude inside the 1e-10 gate.

Two further checks that are stronger than quadrature:

* **Two independent closed-form routes to `h`.** The explicit radial-Laplacian
  route and the hydrogenic eigen-trick route share no code, and their difference
  is **exactly zero, symbolically**, for every system in the table.
* **Literature closed forms** (available only for the zeta = 1 H2+ pair): the
  classical 1s two-centre expressions for `S_AB`, `T_AB`, `V^B_AA` and `V^A_AB`
  agree with the built closed forms to 4.8e-42, 1.6e-42, 4.5e-43 and 6.9e-42 --
  external to the corpus entirely, and the strongest available check on the
  one-electron layer.
* **Cross-consistency with the already-published H2 row.** H2+ at R = 1.4 shares
  its basis with H2 at R = 1.4, so its two-electron tensor must be bit-identical
  to the published H2 table. It is: `(00|00) = 0.625`,
  `(00|01) = 0.42588266110507069323712000634`,
  `(00|11) = 0.50352093294397668656160463148`,
  `(01|01) = 0.32329114155307318238571155893`, digit for digit.

---

## 5. Two findings worth recording

**(a) The LiH driver's exchange-vs-quadrature line was mislabelled.** In
`debug/qfd_lih.py` the exchange row is printed as "tau-matched at tau <= 4", but
the comparison is between the *full* closed form (summed to tau_max = 20) and a
numeric Neumann reference truncated at tau <= 4. The reported deviations
(2.18e-6, 8.56e-7, 3.76e-7) are therefore dominated by the tau terms the
reference simply omits, not by quadrature error -- they are consistent with the
LiH per-tau table's own tau = 5 relative magnitudes. The extension driver does
the comparison properly: it rebuilds the closed form *at the same truncation*
before differencing, and the exchange class then validates at 1.9e-22 (H2+),
1.7e-22 (He2^2+), 3.4e-21 (HeH+) and 2.7e-20 / 4.6e-20 (BeH+) -- i.e. at the
accuracy of the quadrature route itself, 14 orders of magnitude better than the
LiH line suggests. **The LiH energy and its certified digit count are
unaffected** (the tail bound, not this comparison, is what nets out the LiH
truncation); only the diagnostic line in that driver overstates the residual.

**(b) No engine domain gap was hit.** The known `l > 0` hybrid restriction
(`Z_B < Z_A` strictly, via the shell route) is not reachable from an s-only
basis: every hybrid here goes through the direct `V_L` route and stays
elementary. All four classes evaluated cleanly at every rate combination
attempted -- `(Z_A, Z_B)` in `{(1,1), (2,1), (2,2), (4,1)}` with `a` in
`{1, 2, 4}` -- including the equal-rate `(2,2)` case that makes the `l > 0`
route return NaN. Nothing was worked around.

---

## 6. Cost, for whoever runs this next

The exchange class dominates everything. Rough wall times on one core:

| stage | H2+ | HeH+ | He2^2+ | BeH+ |
|:--|--:|--:|--:|--:|
| two-electron closed forms | 1 s | 195 s | 1 s | 2203 s |
| independent quadratures | ~380 s | ~550 s | ~330 s | ~1500 s |
| total per system | ~460 s | ~880 s | ~400 s | 4492 s |

The symbolic route is used wherever the tau sum terminates (it is then both
faster and exact); the numeric accumulator `qfd_core.exchange_hp` is used only
for the heteronuclear systems, where the symbolic tree at tau ~ 20 is what makes
the symbolic route slow.

---

## 7. Artefact integration

New module `benchmarks/certified_reference/entries_qfd_diatomic.py`, category
`qfd_diatomic` (section 4 of the table), five rows. Registration was kept to the
four shared points the generator requires -- the `CATEGORY_TITLES` line, the
`CATEGORY_BLURBS` entry and the `build_all` import/append in
`generate_table.py`, plus the category-order line in `_common.py` -- because a
second agent was adding a `helium_ci` category to the same files concurrently.
The two sets merged without conflict; the table now carries 61 entries across
seven categories.

Rows for the homonuclear systems are **recomputed live** by the generator (they
are cheap and exact, ~17 s for all three). The two heteronuclear rows are
**quoted** from this campaign with their full accounting, exactly as the T2 row
is, because their exchange assembly runs to tau = 22 and cannot sit inside a
table generator.

Regenerated with `python -m benchmarks.certified_reference.generate_table`;
`tests/test_certified_reference_values.py` (12) and
`tests/test_paper58_qfd.py` (4 + 2 skipped) both green, as are the 18 symbolic
S3 proofs and `tests/test_two_center_eri_aabb.py` (113 passed, 4 skipped).
