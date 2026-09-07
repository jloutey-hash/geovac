# Adversarial external-literature scan: WH7 / WH8 vs the information-theoretic and reconstruction literature

**Date:** 2026-09-06. **Question:** where does the reconstruction / information-theoretic
literature say DISCRETENESS comes from, and does anyone tie it to COMPACTNESS or to a
finite observation window? **Task:** novelty only, not truth.

**Method:** ~90 web queries plus fetch-verification of every identifier below. Unfetched items
are labelled UNVERIFIED and are not load-bearing.

---

## HEADLINE

1. **WH7 is anticipated; the closest claimant is Rovelli - twice, plus Hoehn with theorems.**
   "Compact => finite information => discreteness" is RQM Postulate 1 (1996); "the observer's
   finite lifetime supplies the KMS beta" is Martinetti-Rovelli (2003). WH7's temporal
   conjunction is a one-step composition of the two, and Elze (2003) already performs it.
2. **WH8's position (Born imported, not derived) is QBism's, 2013**, stated more sharply by
   Zhang (2026). Its uniqueness leg is Gleason 1957 / Busch 2003 verbatim.
3. **A 2026 preprint derives Born from a graded spectral triple plus an observer bottleneck
   plus Gleason** (Fall & Kondo). Same tools as GeoVac, opposite verdict - the most dangerous
   single item in the scan.

---

## RANKED FINDINGS

### COLLISION

| # | Verified identifier | What it actually says | Hits |
|---|---|---|---|
| 1 | Rovelli, `quant-ph/9609002`, Int. J. Theor. Phys. **35** (1996) 1637; postulates quoted from the SEP entry "Relational Quantum Mechanics" (Rovelli, rev. 2025-02-04) | Postulate 1 verbatim: *"relevant information is finite for a system with **compact phase space**"*; and *"first postulate implies the characteristic discreteness of quantum theory."* | WH7 A/C; WH8(i) |
| 2 | Martinetti & Rovelli, `gr-qc/0212074`, Class. Quantum Grav. **20** (2003) 4919, DOI 10.1088/0264-9381/20/22/015 | An observer *with a finite lifetime T* sees only the diamond algebra; the vacuum's modular flow over it gives temperature 2 hbar / (pi k_B T). The **observer's compact window supplies the KMS parameter**, and pi is in the answer. | WH7 B+C |
| 3 | Chataignier, Hoehn, Lock & Mele, `arXiv:2409.06479`, New J. Phys. **28** (2026) 034504, DOI 10.1088/1367-2630/ae46d0 | Relative to a **periodic (compact) clock** the dynamics is necessarily periodic; a system periodic w.r.t. a periodic clock can be monotonic w.r.t. an aperiodic one. Temporal structure is a property of the chosen clock, not the system. | WH7 C, with theorems |
| 4 | Fall & Kondo, `arXiv:2604.27125` (29 Apr 2026) | Real **graded spectral triples**; boundary reduction to the commutative centre = the "Aperture", *"a permanent information bottleneck for any embedded observer"*; coherence plus **Gleason** then *uniquely determine the Born rule*. Concludes DERIVED. | WH8 (adversarial) |
| 5 | Fuchs & Schack, `arXiv:0906.2187`; Rev. Mod. Phys. **85**, 1693 (2013) | QBism: the Born rule *"should be seen as an empirical addition to Bayesian reasoning itself ... a normative rule in addition to usual Dutch-book coherence."* Imported by construction. | WH8 (closest) |
| 6 | Zhang, `arXiv:2603.06211` (6 Mar 2026) | Additivity proved indispensable across all existing Born derivations: *"the Born rule cannot be derived solely from other non-probabilistic quantum or additional postulates."* Names the imported ingredient more precisely than WH8 does. | WH8 |
| 7 | Bartlett, Rudolph & Spekkens, `quant-ph/0610030`, Rev. Mod. Phys. **79**, 555 (2007) | Lacking a shared frame = G-twirling. Verified in text: *"we will restrict our attention to Lie groups that (i) are compact, so that they possess a group-invariant (Haar) measure"*; the twirl decomposes states into a **discrete sum over irrep labels q** (block-diagonal). WH8(ii) is already a theorem, with compactness doing exactly the work WH8 assigns it. | WH8(ii) |
| 8 | Elze, `gr-qc/0301109`, Phys. Lett. A **310** (2003) 110 | Two internal dimensions **compactified to a torus**; a **discrete physical time** built from a *quasi-local* (observer-restricted) invariant. Compactness plus observer-restriction => discrete time. | WH7 |
| 9 | Favalli & Smerzi, Quantum **4**, 354 (2020), DOI 10.22331/q-2020-10-29-354 | A **bounded** clock Hamiltonian is what makes a Hermitian time operator exist at all inside Page-Wootters. | WH7 |

### ADJACENT

| # | Verified identifier | Note |
|---|---|---|
| 10 | Dolce, `arXiv:0903.3680`, Found. Phys. **41** (2011) 178 | *"this compactification naturally leads to a quantized energy spectrum."* WH7's mechanism as the whole thesis - but the compactness is the **particle's** de Broglie period, explicitly not the observer's. |
| 11 | Peres, Am. J. Phys. **48**, 552 (1980), DOI 10.1119/1.12061 | A finite (N = 2j+1, compact-SU(2)) clock gives **discrete** time readings. WH7 at model level, in 1980. |
| 12 | Pegg & Barnett, Phys. Rev. A **39**, 1665 (1989), DOI 10.1103/PhysRevA.39.1665 | Truncate to (s+1) dimensions and the phase spectrum is discrete at spacing **2 pi/(s+1)**: "2 pi enters at compactification", worked. |
| 13 | Carroll, `arXiv:2307.11927` | Finite-dim Hilbert space + spectrum conditions => periodic evolution => discrete time. Finite-dimensionality **assumed**. |
| 14 | Bekenstein, Phys. Rev. D **23**, 287 (1981), DOI 10.1103/PhysRevD.23.287 | S/E <= **2 pi** R / hbar c for **bounded** systems: boundedness and 2 pi both present, but it bounds an entropy ratio, with no discreteness claim. |
| 15 | Kempf, `arXiv:1010.4354`, NJP **12** (2010) 115001 | Bandlimitation + Shannon sampling: continuous and discrete descriptions are **equivalent**, not rival. **WH7's "honest cap" (empirical undecidability) is this paper, 2010.** |
| 16 | Zeh, `arXiv:0809.2904`, Found. Phys. **40** (2010) 1476 | *"Quantum discreteness is an illusion"* - discreteness as interpretive artifact, via oscillator spectra + decoherence. WH7's spirit, other mechanism. |

### BACKGROUND / TRIVIALIZERS (do not claim as novelty)

- **Rellich / compact-resolvent => purely discrete spectrum**, and **Peter-Weyl**. "Records are
  bound/compact => spectra discrete" is a graduate-textbook theorem; WH7's *automatically* is
  carrying a 1930s functional-analysis lemma.
- Chirvasitu, `arXiv:2310.15139` (to appear J. Noncommut. Geom.): *discreteness <=> spectrum
  finiteness <=> uniform continuity* for compact (quantum) group actions - WH8(i)'s
  identification is a **proved equivalence family in GeoVac's own field**.
- **GPT axiomatics already build compactness in.** Plavala, `arXiv:2103.07469` (Phys. Rep., DOI
  10.1016/j.physrep.2023.09.001), Prop. 2.2: every state space is **compact** convex in a
  finite-dim real vector space. Barnum & Hilgert, `arXiv:1904.03753`: "spectrality" (every state
  decomposes into finitely many perfectly distinguishable pure states) is an **axiom** on compact
  convex bodies, not a consequence of compactness.
- **Gleason**, J. Math. Mech. **6** (1957) 885, extended to effects by **Busch**,
  `quant-ph/9909073`, PRL **91**, 120403 (2003) - WH8's "unique given the projection lattice".

### Does ANY reconstruction DERIVE discrete outcome menus? **No.**

Verified individually: **Hardy** `quant-ph/0101012` (N finite by fiat);
**Chiribella-D'Ariano-Perinotti** `arXiv:1011.6451` (finite dimensions and finite-outcome tests
restricted at the outset; probability pairing primitive); **Masanes & Mueller** `arXiv:1004.1483`
(state space *taken* finite-dim compact convex); **Dakic & Brukner** `arXiv:0911.0695` (Axiom 1 =
capacity of one bit); **Hoehn** `arXiv:1412.8323` (Rule 1 = a limit on the observer's
information; questions binary by construction); **Zeilinger**, Found. Phys. **29**, 631 (1999).
Finiteness is always an axiom. Rovelli alone *reads discreteness off it* - WH8(i)'s move, thirty
years earlier.

---

## OVERALL VERDICT - WH7 (time is discrete because the observer's window is compact)

**ANTICIPATED in substance; not identical in statement.** Leg (A) compact => discrete and leg
(B) compactification => 2 pi (Matsubara 2 pi n / beta; KMS beta = 2 pi) are textbook. Leg (C) -
that the compactness is the **observer's** - has four prior claimants: RQM Postulate 1 (1996),
Martinetti-Rovelli's diamond temperature (2003), Elze (2003), Chataignier-Hoehn-Lock-Mele (2026,
with theorems). **Closest existing position: Martinetti & Rovelli 2003 read with RQM Postulate
1** - Rovelli publishes both halves; WH7's distinctive act is joining them on the temporal axis,
and Elze already performed a version of that join. The corpus's defensible contribution is the
Toeplitz falsifier programme, not the thesis.

## OVERALL VERDICT - WH8 (the Born measure is the exchange constant of the observation projection)

**ANTICIPATED as a position; contested as a claim.** "Born is an irreducible import, not
derived" is QBism (Fuchs & Schack 2013) and, sharper, Zhang (2026). WH8(i) (records compact =>
menus discrete) is Rovelli's Postulate 1 plus Rellich. WH8(ii) (diagonality of effective observer
states) is the compact-group G-twirl superselection theorem (Bartlett-Rudolph-Spekkens 2007),
where compactness is exactly the condition that makes the Haar average exist. The uniqueness leg
is Gleason. **Closest existing position: QBism**; runner-up Zhang `arXiv:2603.06211`. No prior
use of "exchange constant" for the Born measure was found, and no prior placement of it inside a
named projection taxonomy - **the taxonomic slot is unoccupied; the content that would fill it is
not.**

---

## FALSIFIERS / TRIVIALIZERS - flag these

1. **WH7's sharpest hole: a compact *interval* is not a compact *circle*.** A finite real-time
   window with no BC gives a broadened but still **continuous** spectrum; spacing 2 pi / beta
   needs periodic BCs, i.e. S^1. And the Matsubara/KMS circle is the **imaginary**-time circle
   while an observer's lifetime is a **real**-time interval. WH7 identifies the two silently.
   Aim the falsifier here.
2. **WH7 counterexample candidate (my inference - verify first):** local algebras of bounded
   regions, diamonds included, are Type III_1, with **continuous** modular spectrum. In the one
   rigorous computation of "observer with a compact window" - Martinetti-Rovelli - the window
   yields a **temperature but no discrete spectrum**. If it holds, fatal to strong WH7.
3. **WH7's honest cap is already published:** Kempf, NJP **12** (2010) 115001.
4. **WH8's sharpest threat is self-inflicted:** Masanes, Galley & Mueller, `arXiv:1811.11060`,
   Nat. Commun. **10**, 1361 (2019) derive Born + state-update from unitary QM **plus** "the
   assumption that ensembles on finite-dimensional Hilbert spaces are characterised by finitely
   many parameters" - i.e. **WH8's own finiteness premise**. Grant WH8(i) and MGM says Born
   *follows*, contradicting "external calibration input". Contested, not settled: Kent,
   `arXiv:2307.06191`, Quantum **9**, 1749 (2025) refutes MGM with explicit alternative rules;
   Galley & Masanes, Quantum **1**, 15 (2017) show Born-alternatives need extra axioms to kill.
5. **WH8 attacked in GeoVac's own idiom:** Fall & Kondo, `arXiv:2604.27125`. If it stands it
   contradicts WH8's classification directly. Also Torres Alegre, `arXiv:2512.12636`
   (UNVERIFIED): no-signaling plus purification selects Born.
6. **Vocabulary collision:** Lax, `arXiv:2604.27339` (30 Apr 2026), "The Born Rule ... from
   Metric Non-Expansion and **Calibration**" - "calibration pins the vertex values" - DERIVED.

## UNVERIFIED (do not cite without checking)

Pegg PRA **58**, 4307 (1998); Salecker-Wigner Phys. Rev. **109**, 571 (1958); Page & Wootters
PRD **27**, 2885 (1983); Sewell Ann. Phys. **141**, 201 (1982); Hayden & Wang Quantum **9**,
1664 (2025); Torres Alegre `arXiv:2512.12636`; Banks `arXiv:2509.17856`.

**Withdrawn here:** Goheer-Kleban-Susskind `hep-th/0212209` (JHEP 0307:056) does **not** state in
its abstract that spectral discreteness follows from finite entropy.
