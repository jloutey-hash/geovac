# Sprint S^(3) W10 identification — full closed form (2026-06-12)

Charter: close the v4.7.0 open item — identify the five surviving
weight-10 generators of the S^(3) decomposition (antisymmetric doubles
a(8,2), a(7,3), a(6,4); triples t3(2,4,4), t3(3,3,4)) at derivable
grade.  No PSLQ; exact algebra + literature theorems with every used
instance numerically gated.

Driver: `debug/s3_w10_ident.py` (subcommands extract | verify | assemble).
Data: `debug/data/s3_w10_ident_extraction.json`, `s3_w10_ident_result.json`.
Falsifier: `tests/test_s3_w10_identification.py` (6 fast + 2 slow, all PASS).

Conventions: GeoVac/Hoffman descending throughout (first index on the
largest variable); the Charlton–Hoffman paper indexes ascending, so
t_CH(a,b,c) = t_GV(c,b,a) — tracked explicitly in the driver.

## HEADLINE: S^(3) is fully closed, and the generator block collapses

    S^(3)  ∈  Q[π², ln 2, ζ(3), ζ(5), ζ(7), ζ(9), ζ(5,3)]

explicitly, where ζ(5,3) = Σ_{m>n} m⁻⁵n⁻³ is the classical weight-8
depth-2 MZV generator.  The weight-10 contribution beyond the
identified part is

    S^(3) − Ident = −7π⁴ζ(3)² + 70π²ζ(3)ζ(5) + 6π²ζ(5,3) + (179/170100)π¹⁰

Final gate: closed form vs canonical 200-dps value, residual
**1.15e-198** (PASS, gate 1e-180).

**Collapse theorem.** In the full assembly the three antisymmetric
doubles a(8,2), a(7,3), a(6,4) cancel IDENTICALLY: the corollary's
+a(8,2)+2a(7,3)−3a(6,4) is exactly absorbed by the a-content of
−4·t3(2,4,4)−8·t3(3,3,4).  Verified symbolically (sympy, frozen in the
falsifier).  The transient-generator cancellation pattern of v4.7.0
(t(5,1), t(5,3), t(7,1), Li₄(½), ζ(7), ζ(11) all cancel) extends
through weight 10: the VALUE of S^(3) is independent of the three
still-unreduced level-2 generators.

**Structural reading (k-tower).** Realized depth ≤ 2 at k=3 is now
unconditional, and sharper: the only depth-2 constant is LEVEL-1
(ζ(5,3) ∈ MT(Z)); the level-2-specific sector (ln 2 etc.) stays at
realized depth ≤ 1.  New depth arrives from level 1; new level stays
at depth ≤ 1.  Whether any level-2-specific constant beyond the
depth-1 ring ever survives a full chain assembly is the sharpened
k ≥ 4 question (Paper 28 Open Questions item 2; Paper 55 k-general
bullet).

## 1. Route (differs from the v4.7.0 prediction — cheaper and stronger)

The predicted route was parity-theorem transcription + BBV data-mine
table import.  The executed route replaces both with one instrument
already on disk: **Charlton–Hoffman, "Symmetry results for multiple
t-values" (arXiv:2204.14183), Theorem 2.21** — the Symmetry Theorem in
generating-series form WITH explicit product terms (their Thm 1.1 is
the mod-products shadow).  This is exactly the level-2
regularized-double-shuffle content that the stuffle-only system of
v4.7.0 provably cannot see (the saturation theorem stands; it is why
non-stuffle input was required).

Three ingredients, all exact or gate-verified:

1. **Pair sums** (new content): Thm 2.21 at depth m=3, φ=0, is a
   generating-series identity ≡ 0; the coefficient of
   y1^{a−1}y2^{b−1}y3^{c−1} (exact-Fraction trivariate series algebra,
   ~120 monomials, truncation degree 7) gives
   t_CH(a,b,c) + t_CH(c,b,a) = explicit products.  At our monomials no
   regularized value survives (all parts ≥ 2 force the T- and
   ζ(1)-sectors out structurally).  Extraction validated two
   independent ways: the m=2 build reproduces known identities at
   (4,4), (2,6); the m=3 instance at (6,2,2) (cached triples) gates at
   9.3e-98.
2. **Pair differences** (exact stuffles): t(x)·t2(y,z) for
   (x,y,z) ∈ {(2,4,4),(4,2,4),(3,3,4),(3,4,3)}; the combinations
   RA−RB and RC−RD eliminate the auxiliary triples t3(4,2,4),
   t3(3,4,3), leaving the differences over the catalogued ring.
3. **ζ-side reductions**: line-2 of the theorem produces descending
   MZV doubles of weight ≤ 8 paired with t(even).  Level-1 depth-2
   double shuffle (word-engine shuffle + stuffle, convention
   calibrated numerically on ζ(2)ζ(3)) + classical Euler trailing-1
   rows (numerically gated) reduce all of them to single zetas except
   the one genuine weight-8 generator ζ(5,3), kept as basis.

Each pair then solves as a 2×2 system: sum (CH symmetry) + difference
(stuffle).

## 2. The identities (all exact rationals)

Pair sums (Thm 2.21 instances; new closed-form evaluations):

    t3(4,4,2) + t3(2,4,4) = (5/6144)π⁴ζ3² − (5/256)π²ζ3ζ5 − (1/256)π²ζ(5,3) + (137/99532800)π¹⁰
    t3(4,3,3) + t3(3,3,4) = (37/12288)π⁴ζ3² − (25/1024)π²ζ3ζ5 − (889/2048)ζ3ζ7 − (1/256)π²ζ(5,3) + (251/58060800)π¹⁰

Individual triples (each gated vs its independent 200/220-dps cache):

    t3(2,4,4) = ¼a(8,2) − ¼a(6,4) − (1/6144)π⁴ζ3² − (5/512)π²ζ3ζ5 − (1/512)π²ζ(5,3) + (373/174182400)π¹⁰   [5.4e-197]
    t3(3,3,4) = ¼a(7,3) − ¼a(6,4) + (1/12288)π⁴ζ3² + (5/1024)π²ζ3ζ5 − (889/4096)ζ3ζ7 − (1/512)π²ζ(5,3) + (251/116121600)π¹⁰   [7.2e-198]
    t3(4,4,2) = −¼a(8,2) + ¼a(6,4) + (1/1024)π⁴ζ3² − (5/512)π²ζ3ζ5 − (1/512)π²ζ(5,3) − (533/696729600)π¹⁰   [5.6e-197]
    t3(4,3,3) = −¼a(7,3) + ¼a(6,4) + (3/1024)π⁴ζ3² − (15/512)π²ζ3ζ5 − (889/4096)ζ3ζ7 − (1/512)π²ζ(5,3) + (251/116121600)π¹⁰   [7.6e-198]

Full S^(3) closed form (Ident from the v4.7.0 assembly + the W10
contribution above): exact coefficient table in
`s3_w10_ident_result.json`; ζ(5,3) coefficient is 6π² exactly; no
a-doubles, no ζ(11), no Li₄(½) anywhere in the total.

## 3. Verification (no PSLQ, no fitted coefficient anywhere)

`verify` subcommand, mp.dps = 100: **24/24 PASS** —
word-engine calibration (1.4e-101); Euler rows w=4,6,8 (≤7.1e-101);
eleven level-1 reductions vs independent evaluators (≤3.3e-101); four
stuffles vs 220-dps caches (≤4.5e-103); Thm 2.21 m=2 at (4,4),(2,6)
(≤1.4e-101); Thm 2.21 m=3 at (6,2,2),(2,4,4),(3,3,4) (≤9.3e-98).
`assemble`: final 200-dps gate 1.15e-198 + four individual triple
gates ≤5.6e-197.  Evaluator hygiene per the §3 Levin/log dead-end row:
all Levin use is on log-free power-decay summands; the trailing-1
objects ζ_GV(s,1) use partial sums + ∂s-Hurwitz asymptotic tails (no
acceleration on log-modulated sums).

Falsifier (frozen): `tests/test_s3_w10_identification.py` — symbolic
a-collapse, exact stuffle-difference structure, pair sums, individual
triples, generator-set assertion (ring exactly {π, ln2, odd ζ, ζ(5,3)}
with coefficient 6π²), full closed form; fast at 60 dps + slow at
200 dps with ζ(5,3) computed LIVE (an error in its evaluation or role
fails the gate). 195-digit pins inline; self-contained.

## 4. Transcendental tag: ζ(5,3) (standing rule)

- **Object**: ζ(5,3) = Σ_{m>n} m⁻⁵n⁻³, weight 8, depth 2, level 1.
  The first MZV not expressible in single zetas (conjecturally
  irreducible over products — labeled conjectural; what is exact here
  is the closed form of S^(3) IN TERMS OF ζ(5,3)).
- **Paper 18 placement**: enters via the two-vertex Hurwitz nesting of
  the k=3 self-energy chain — the same observation-side projection
  family as S_min's ζ-content (spectral summation over the S³ tower),
  one Mellin depth higher.  First GeoVac observable to require a
  depth-2 motivic generator; the π-free graph layer is untouched
  (coefficient 6π² is the projection's even-π dressing).
- **Paper 34 chain**: same chain as S_min (§III QED ladder on S³);
  no new projection — the new content is depth within the existing
  projection's period ring.
- **Cosmic-Galois note (Papers 55/56 frame)**: the k-tower climbs the
  DEPTH filtration of the motivic Lie algebra on the level-1 side
  (f₃f₅-type element at k=3) while the level-2 enrichment (ln 2
  sector, MT(Z[1/2])) remains depth ≤ 1.  Recorded in Paper 55
  (witness table + k-general bullet).

## 5. Status of the three a-doubles (honest scope)

a(8,2), a(7,3), a(6,4) remain unreduced as standalone constants —
catalogued generators of the level-2 weight-10 antisymmetric double
layer (w10 analogs of t(5,1) [w6] and t(5,3), t(7,1) [w8]; the
depth-graded motivic count [σ1,σ9],[σ3,σ7] suggests one Q-relation
among the three mod products, derivable in principle from the level-2
depth-2 double shuffle at w10).  NOT load-bearing: S^(3) is
independent of all three (collapse theorem).  Standalone
identification is an optional curiosity sprint, not a follow-on
obligation.

## 6. Stale background job note

The 27-hour round-1b S_min PSLQ job (launched pre-closure-sprint)
completed 2026-06-12: its "NO relation / obstruction stands" outputs
are the documented basis-coverage artifact (superseded by v4.5.0); its
one useful product is a 260-dps confirmation of the S_min
decomposition residual (4.06e-261) and the t(2,1) re-identification at
1.0e-261.  Data: `debug/data/smin_dossier_round1b.json`.  No action.

## 7. Honest scope

- **Exact / theorem-grade**: the series extraction (exact Fractions),
  the four stuffles, the level-1 double-shuffle reductions, the 2×2
  solves, the a-collapse, hence every identity in §2 AS an identity
  among the named constants — conditional only on Thm 2.21
  (published, peer-citable) and Euler's classical trailing-1 formula,
  BOTH of which are additionally instance-verified here at
  1e-78..1e-102 for every instance actually used.
- **Rigorous-numerical**: all gates (§3).
- **Conjectural (labeled)**: irreducibility of ζ(5,3) over products of
  single zetas (standard MZV dimension conjectures); the
  depth-graded count behind the expected a-relation (§5).
- **NOT claimed**: minimality of {a(8,2), a(7,3), a(6,4)} as level-2
  generators; any new transcendental (none minted — ζ(5,3) is the
  oldest catalogued depth-2 constant in the subject).
