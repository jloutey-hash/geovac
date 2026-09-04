# Slow-test audit: is there a closed form? (2026-09-04)

PI question: *"audit our slow tests. the slow ones that are marked slow and
the ones that aren't. maybe there's a closed form or something that can
compute them faster."*

The 2026-08-31 cost memo (`test_suite_cost_memo.md`) measured the suite and
concluded the cost is **inherent** — 6.7 h serial, 2.4 h with `-n auto`, 98 of
361 files over budget, dominated by "everything that builds a Hamiltonian or
solves an eigenproblem." It tested parallelism, thread-pinning and process
isolation. **It never asked whether the computation could be replaced.**

That is the algebraic-first question (§4), applied to the test suite rather
than to the physics, and it has an answer.

---

## 1. The finding

**The production lattice's Laplacian spectrum is closed form, and it has been
since Paper 0 §VI / Paper 1 §III — but the closed form lived as a private
helper inside the test that proves it**, so nothing else could use it. Fourteen
other test files were still calling `numpy.linalg.eigh` on the dense operator.

Measured, n_max = 30 (9455 × 9455):

| route | time | agreement vs dense | speedup |
|:--|--:|:--|--:|
| dense `eigh` | 69.2 s | — | 1× |
| blockwise `eigh` (eigenpairs) | 1.05 s | 1.4e-14 | **66×** |
| closed form (eigenvalues) | 1.0 ms | 1.4e-13 | **43,800×** |

Both are **exact**, not approximations. `L` is block diagonal in `ℓ` (no edge
changes `ℓ`) and each block is a Cartesian product of two path graphs, so

```
spec(L_ℓ) = { 2 − 2cos(jπ/(n_max−ℓ)) + 2 − 2cos(kπ/(2ℓ+1)) }
```

Eigenvectors are the tensor products of path eigenvectors
`v_j(i) = cos(jπ(i+½)/m)` — residual 3e-15.

Promoted to `geovac/lattice_spectrum.py`.

## 2. What it bought, measured end to end

| test | before | after | note |
|:--|--:|--:|:--|
| `test_sp_splitting_identification_is_ill_conditioned` | 69.70 s | **0.94 s** | also now basis-free |
| `test_lambda_max_deficit_at_nmax_70` | ~180 s, `@slow` | **3.17 s** | marker removed; now runs by default |

Both keep their assertions **verbatim**; only how the quantity is obtained
changed. Both still diagonalise the operator that was actually built, so
neither now *rests on* the closed form — they remain independent checks **on**
it. Fire-tested: both FIRE on a planted perturbation.

**A capability, not just a speedup.** `lambda_max(n)` is O(n_max):

| n_max | nodes | time |
|--:|--:|--:|
| 70 | 116,795 | 0.03 ms |
| 5,000 | 4.2e10 | 1.9 ms |
| 40,000 | 2.1e13 | 16 ms |

Claims at cutoffs the dense route could never reach are now testable. (The
first version of `lambda_max` materialised every block and could not finish at
n_max = 5000 — the profile caught it; the top of a Cartesian product is the sum
of the factor tops, so it is O(n).)

## 3. The marker does not track cost

Of the **six costliest** tests measured in this slice, **zero** carry
`@pytest.mark.slow`:

| test | s | marked |
|:--|--:|:--|
| `test_sp_splitting_identification_is_ill_conditioned` | 69.7 | no |
| `test_paper27_ep2l_nmax5_below_two` | 62.8 | no |
| `test_spectrum_confined_and_bottom_dense_at_nmax_30` | 11.5 | no |
| `test_paper27_ep2b_ho_gs_is_not_a_single_determinant` | 9.6 | no |
| `test_paper27_proposition_nondegeneracy_qualifier` | 9.5 | no |
| `test_cg_construction_splitting_does_not_decay` | 2.4 | no |

Meanwhile the one marked test in scope (n_max = 70) was the *only* one with a
closed-form escape, and is now 3.2 s. In this sample the marker is
**anti-correlated** with cost. 307 markers across 139 files are therefore not a
cost map, and should not be read as one.

This does not overturn the 2026-08-31 conclusion "do not mark the 98 slow" —
that conclusion was about coverage, and it stands. It adds: the markers already
present do not identify the expensive tests either.

## 4. Where there is no closed form

Checked and negative — recorded so the next audit does not re-derive them:

- **`test_paper27_*` (62.8 s + 9.6 + 9.5 + …)** — many-electron CI +
  entanglement. Genuine physics; no closed form. *Separately:* this file
  imports its builder from `debug/archive/misc/`, i.e. the prune-by-design
  directory, and doubly so from `archive/`. One of C22 check-D's three
  baselined offenders.
- **Dirac lattice** — already closed form (|λ| = n + 3/2) and small (n_max ≤ 7);
  not a cost centre.
- **Remaining dense-`eigh` sites** — `test_paper25_s2_quotient` (n_max = 3),
  `test_graph_qed_vertex` (n_max = 2), `test_nuclear_lattice`: all small. The
  lattice opportunity was concentrated in the two tests above and is harvested.
- **`test_composed_*`, `test_balanced_*`, `test_prolate_*`, `test_level3/4_*`,
  `test_n_electron_*`** (the 98-file bulk) — Hamiltonian construction. The
  2026-08-31 verdict stands: inherent.

## 5. Found while auditing, not caused by it

`test_paper27_entropy.py` has **3 failures that predate this session.**
Established by dates, not by inference — every input last changed on or before
2026-08-30:

| file | last change |
|:--|:--|
| `geovac/lattice.py` | 2026-04-12 |
| `geovac/nuclear/ho_two_fermion.py` | 2026-04-15 |
| `geovac/molecular_spec.py` | 2026-06-07 |
| `debug/archive/misc/energy_entanglement_decoupling.py` | 2026-06-11 |
| `tests/test_paper27_entropy.py` | 2026-08-30 |

Two distinct assertions fail:

```
V_diagonal_fraction_in_H1_eigenbasis
  obtained 0.9202400368022998
  expected 0.9202580563910042 ± 9.2e-09
```

and

```
test_paper27_ep2n_be_analytical_degenerate_pt
  obtained 1.8986013727937123
  expected 1.924 ± 0.002
```

The first is a **diagonal fraction in an eigenbasis** pinned to nine
significant figures — a basis-dependent quantity of exactly the class
/qa DELTA #5 flagged, where a no-op relabelling moved a reported 2.70% to
4421%. The tolerance asserts reproducibility the quantity does not have. Not
fixed here; logged.

### 5a. The process finding underneath

**A paper-backing test file has been red since at least 2026-08-30 and no gate
noticed.** That is not an accident of this file — it is a consequence of the
cadence decision made in the 2026-08-31 cost memo:

> *the full scope is a **scheduled baseline**, not a close gate.*

That was the right call on cost. But **there is no scheduler.** `/regression
touched` derives its selection from the diff and the import graph, so a file
whose inputs have not changed is never selected — which is exactly the file
that can rot unobserved. The cadence has a gap where the schedule was assumed.

Cheap closure, not implemented here: run the full suite on a fixed cadence with
`-n auto` (≈2.4 h) and record the red list, or add a nightly/weekly job whose
only output is *which files are red*. Either makes "scheduled baseline" real
rather than nominal.

## 6. Owed (guard pass, per CLAUDE.md §9)

1. **`tests/test_lattice_spectrum.py` does not exist.** The module is verified
   against dense `eigh` at n_max ∈ {8, 20, 30} and exercised by two converted
   tests, but has no dedicated file.
2. **The degeneracy-pooling property is unverified.** Fire-test evidence:
   planting `while False:` into `eigenspace_overlap`'s pooling loop — i.e.
   disabling the very thing that makes the result basis-free — leaves
   `test_sp_splitting_identification_is_ill_conditioned` **green**. The
   converted test checks the near-tie, not the invariance. That guard is owed
   and must name pooled-vs-unpooled as the wrong answer it rejects.
3. `tests/_durations.json` is still the empty 2-byte file from 2026-06-07
   (carried from the 2026-08-31 memo), so `/regression fast` still selects
   nothing.
