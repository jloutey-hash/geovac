# Is chemistry disjoint from the ELLIPTIC period class?

**2026-08-22, exploratory (PI-directed). Drivers:
`debug/chemistry_elliptic_period_class.py`, `debug/chemistry_elliptic_hiprec.py`.
No paper claim.**

## The gap this closes

The W1e period-class sprint (2026-06-04) tested 11 chemistry correction terms
against the outer-factor classes **M1/M2/M3** and got 0/11, concluding chemistry
is calibration tier, "categorically disjoint from outer-factor periods."

But M1/M2/M3 are the **pure-Tate** classes. Paper 59's **elliptic Γ(2)** periods
are a different class that did not exist when that sprint ran. Whether chemistry
is disjoint from *that* class had never been asked.

## Targets

The certified closed-form H₂ PES constants — the only chemistry quantities in
the corpus with PSLQ-grade precision, because E(R) is one symbolic expression
over {exp, E₁, log, γ} and closed-form Newton on its exact derivatives certifies
its critical point.

The stored 49 digits turned out to be a **formatting cap** (`mp.nstr(..., 50)`
inside `h2_pes_certify`), not a computational limit. Recomputed at dps 150/170
with cross-precision agreement:

| constant | agreement |
|---|---|
| R_eq | ~150 digits |
| D_e | ~151 digits |
| k | ~149 digits |

Working precision used: **130 dps** (certified, with margin).

## Result: DECISIVE-NEG at weight ≤ 2, height ≤ 10⁸

| ring | \|ring\| | ceiling | R_eq | D_e | k |
|---|---|---|---|---|---|
| wt≤1 {π, ϖ} | 4 | 10^11.5 | **NEG** (h≤10⁶) | **NEG** | **NEG** |
| wt≤1 {π, ϖ, G} | 4 | 10^11.5 | **NEG** (h≤10⁶) | **NEG** | **NEG** |
| wt≤2 {π, ϖ} | 9 | 10^14.4 | **NEG** (h≤10⁸) | **NEG** | **NEG** |
| wt≤2 {π, ϖ, G} | 10 | 10^13.0 | **NEG** (h≤10⁸) | **NEG** | **NEG** |
| wt≤2 {π, ϖ, G, P8} | 20 | 10^6.5 | INCONCLUSIVE | INCONCLUSIVE | INCONCLUSIVE |
| wt≤3 {π, ϖ, G} | 20 | 10^6.5 | INCONCLUSIVE | INCONCLUSIVE | INCONCLUSIVE |

ϖ = K(1/2) = Γ(¼)²/(4√π), the disc-4 CM period, signed (negative powers admitted,
so quasiperiod directions are in the ring). G = Catalan = β(2), the Eisenstein
L-value. P8 = the disc-8 period.

**So: the H₂ chemistry constants are not low-height elements of the disc-4
elliptic period ring at weight ≤ 2.** The W1e finding extends to the elliptic
class — chemistry is disjoint from that too, at the heights reachable.

## Why the INCONCLUSIVE rows are honest, not evasive

They are *decoy-matched*: PSLQ found a "relation" for the real target AND for a
random decoy of the same magnitude, at comparable height. That is the signature
of a search finding noise. With 20 basis elements at 130 digits the honest
detectable height is only 10^6.5, and the search ran right at it.

Closing wt≤3 and the disc-8 ring needs roughly **250–300 digits** (20 terms ×
8 decades + margin). The Newton runs comfortably at 170 dps, so this is
reachable — it is a compute question, not a structural one.

## Gates

* **G1 decoy calibration** — every target paired with a decoy; a real hit is
  only reportable if the decoy does *not* match. This is what converted four
  apparent "relations" into INCONCLUSIVE rather than into a false discovery.
* **G2 positive control** — PSLQ must recover a planted element built from each
  ring's *own* keys. PASS at every ring size (height ≤ 40). *(First version
  planted π·ϖ, a weight-2 element, which is absent from a weight-1 ring; the
  control "failed" for a reason unrelated to the ring's power. Fixed.)*
* **G3 honest ceiling** — 10^(D/n) reported for every ring, and no negative
  claimed above it. The search coefficient cap is set *at* the ceiling so that
  "no relation found" is a statement rather than an artifact.

## The Yukawa question: structurally blocked, and not by the ring

The same extension was proposed for the Yukawa couplings (the H1 non-selection
theorem + the 162-cell PSLQ sweep also covered only M1 ∪ M2). It cannot be run.

PSLQ resolves a relation among n terms only if precision D exceeds roughly
n·log₁₀(H). **Yukawas are measured**, so D is capped by experiment at ~8 digits
for the charged leptons (worse for light quarks). That gives:

| ring size | detectable height at 8 digits |
|---|---|
| 4 | ~10² |
| 9 | ~8 |
| 20 | ~2.5 |

Even the smallest meaningful ring reaches only height 100. The corpus's own
sweep memo already recorded this — "charged-lepton precision (8 digits at M_Z)
honest at M=10."

**The structural point:** this test worked for chemistry because the target is a
*mathematically defined* constant we can compute to 150 digits. A Yukawa is a
*measured* number. PSLQ needs the former. No ring refinement fixes a
precision-starved target, and running it anyway across 9 fermions would
manufacture coincidences rather than test for them — the Paper-2 failure mode.

## Standing

Chemistry now measured disjoint from **both** period families at low weight:
pure-Tate M1/M2/M3 (W1e sprint, 2026-06-04) and elliptic Γ(2) disc-4 (here).
That strengthens the calibration-tier placement of chemistry rather than
overturning it.
