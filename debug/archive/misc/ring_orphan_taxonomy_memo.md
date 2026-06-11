# Ring-orphan taxonomy — does the orphan club look like decomposed alpha? (2026-05-30)

**Origin:** Josh — the ring-orphan constants (correlation entropy S_full(GS), K, the
L2 constant c, the Wolfenstein parameters) sit at one table; does that table look like
the decomposed alpha K = pi(B+F-Delta)? "Maybe nothing there, but look and think."

**Receipts:** debug/ring_orphan_probe.py + debug/data/ring_orphan_probe.json.
S_full = -sum p_i log p_i reproduced from the 1-RDM occupations to 2e-17 (He n3) and
5e-16 (Li+ n3); occupations are eigenvalues of an EXACT-RATIONAL 1-RDM (algebraic);
dominant occupation p_max is algebraic of degree > 6.

## Verdict: same CLUB, different TIER. The parallel has teeth — and isolates K.

The instinct is right at the **club** level: every member is "a simple functional of
COMPUTABLE pieces, whose WHOLE is ring-orphan." K = pi*(B + F - Delta), pieces each
computable, sum orphan. S_full = -sum p_i log p_i, occupations each algebraic, entropy
orphan. Same shape.

But the **mechanism** of orphan-hood differs, and naming it isolates alpha as the one
hard case:

| orphan | what it is | WHY it is orphan | tier |
|:--|:--|:--|:--|
| **K = pi(B+F-Delta)** | a graph combination | it **equals an external physical constant** (alpha^-1, to 8.8e-8) AND the combination rule resists any common generator (12 mechanisms eliminated) | **DEEP** (coincidence-with-physics) |
| **c (L2 next-order)** | the propinquity convergence **rate** constant | it is a Stein-Weiss IBP **derivative of M1**, one step downstream of the ring (MR-C) | intermediate (derivative) |
| **S_full(GS)** | von Neumann entropy of the GS | it is a **Baker-class transcendental**: a Q-combination of log(algebraic occupations). Generic fate of ANY entropy of an algebraic state. Not even a universal constant — system-specific. | **TAME** (derived functional) |
| **Wolfenstein** | CKM mixing params | one period in the ratio is **off-framework** (the Higgs vev, AC inner factor) | off-framework |

So three of the four orphans have a **mechanism** explanation. Only K/alpha is a genuine
coincidence-with-physics. Putting the table together does NOT unify alpha with the others
— it **sharpens what the alpha mystery is**: not "why is this combination uncomputable"
(lots of natural numbers are orphan-in-the-ring, for boring reasons), but "why does THIS
particular orphan equal an external physical constant."

## The confinement reading (the unifying map — honest status: retrodiction)

The Mellin ring (M1/M2/M3) is the **interior** = the fully-confined skeleton (pi-clean,
forced). The orphans are **boundary observables**:

- alpha / K  = boundary **permeability** (coupling).
- c          = boundary-**crossing rate** (how fast discrete -> continuum; literally the
               propinquity rate constant).
- S_full     = boundary **position** (how far a state has slid from the closed-shell,
               S=0, single-determinant confined point toward the open correlated regime).
- Wolfenstein = a boundary the framework **does not contain**.

The interior is computable-in-the-ring; the seam is not, by its nature. One genuine
**retrodiction**: the reframe predicts the propinquity RATE constant c should be a boundary
observable (orphan) — and it is (MR-C null). That is the reframe earning slightly more than
relabeling.

## Honest status (charter sec 6 bar)

This is a **taxonomy of the orphan club**, not a unification of alpha with the rest. It
retro-explains (why each orphan is orphan) and isolates K. It does NOT yet predict past the
known. The earnable next step is **classify-then-compute**: take an UNclassified quantity,
decide interior-vs-boundary from its DEFINITION, predict ring-membership BEFORE computing,
then compute. Candidate targets: a thermal entropy at a non-BW point (predict: interior /
bridged / modular); a Hawking-side crossing-rate analog of c (predict: boundary / orphan).

## One-liner for charter sec 7 (proposed, not yet applied)

> The ring-orphan club (K, c, S_full, Wolfenstein) is the set of BOUNDARY observables —
> permeability, crossing-rate, position, off-framework. The Mellin ring is the confined
> interior. Of the four, only alpha/K is orphan by coincidence-with-physics; the rest are
> mechanism-explained. Retrodiction, not yet prediction.
