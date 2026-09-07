# Adversarial literature scan: "Discreteness is compactness"

**Date:** 2026-09-06 | **Target:** CLAUDE.md S1.7 organizing observation + Paper 18 SIII
("Compactness as the source of discreteness", `papers/group3_foundations/paper_18_exchange_constants.tex`
L707-762) | **Scope:** who states the thesis, how generally, and who states the converse.

All identifiers below were resolved against arXiv/publisher pages. Anything not
directly opened is marked **UNVERIFIED**.

---

## 1. COLLISION tier

### C1. Liu & Noui, *Gravity as an SU(1,1) gauge theory in four dimensions* -- arXiv:1702.06793
Class. Quantum Grav. 34 (2017) 135008. **Verbatim (full text, verified):** "The discreteness of
the quantum geometry at the Planck scale predicted in Loop Quantum Gravity can be interpreted as
a direct consequence of the compactness (via Harmonic analysis) of the residual symmetry group
SU(2)." Abstract (verified): "space-like areas have discrete spectra ... whereas time-like areas
have continuous spectra" -- the space-like sector is compact SU(2), the time-like sector
non-compact SU(1,1). **This is the thesis AND its converse in one paper**, stated causally, with
"harmonic analysis" (= Peter-Weyl) named as the mechanism. Scope: LQG's gauge group only.

### C2. Pontryagin duality: G compact <=> G-hat discrete
The textbook biconditional for locally compact abelian groups; the literature calls it "the
duality between compactness and discreteness". Named-principle evidence: *Duality between
compactness and discreteness beyond Pontryagin duality*, Proc. Steklov Inst. Math. 271 (2010),
DOI 10.1134/S0081543810040164 (title + DOI verified via publisher listing; **abstract
UNVERIFIED**, paywalled). This is the strongest evidence that the slogan is *already named* -- but
it is named in harmonic analysis, not physics, and carries no transcendental content.

### C3. Connes' compact-resolvent axiom (spectral triples)
A spectral triple (A, H, D) requires (D - lambda)^-1 compact; this is *equivalent to* D having
discrete spectrum with finite multiplicities. Connes, *Noncommutative Geometry* (Academic Press,
1994); nLab "spectral triple" (verified). Here "compactness => discreteness" is not a thesis but
a **definition**: NCG builds the equivalence into the axiom and never argues for it. Directly
relevant to WH1 -- the corpus's own keystone framework already presupposes the thesis.

### C4. Dolce, Elementary Cycles Theory -- the WH7 collision
*Introduction to the Quantum Theory of Elementary Cycles*, arXiv:1707.00677, ch. 4 in Licata &
't Hooft (eds.), *Beyond Peaceful Coexistence* (World Scientific 2016), 93-135. **Verbatim
abstract (verified):** "the unification of quantum and relativistic physics is fully achieved by
imposing an intrinsically cyclic (or compact) nature for relativistic space-time coordinates. In
particular the Minkowskian time must be cyclic." Companion: *Gauge Interaction as Periodicity
Modulation*, arXiv:1110.0315 = Ann. Phys. 327 (2012) 1562-1592 (verified) -- "periodic conditions
at the boundary ... as semi-classical quantization condition". Also *Compact Time and Determinism
for Bosons: Foundations*, Found. Phys. 41 (2011), DOI 10.1007/s10701-010-9485-4 (**abstract
UNVERIFIED**). Body-level (search snippet, **not** abstract-verified): "In analogy with finite
temperature field theory and with extra-dimensional field theories, this compactification
naturally leads to a quantized energy spectrum" -- i.e. Dolce explicitly unifies Matsubara + KK +
particle-in-a-box under compactness. **WH7's "time is discrete because the window is compact" is
Dolce's founding postulate**, though he takes compact time as ontological, not observer-relative.

---

## 2. ADJACENT tier (one instance, or the theorem without the thesis)

- **Rovelli & Smolin**, *Discreteness of area and volume in quantum gravity*, gr-qc/9411005 =
  Nucl. Phys. B442 (1995) 593-622, Erratum B456 (1995) 753 (verified). The founding discreteness
  result; the abstract makes **no** compactness claim. The thesis is attributed retroactively.
- **Ben Achour, Geiller, Noui, Yu**, *Spectra of geometric operators in 3D LQG: from discrete to
  continuous*, arXiv:1306.3246 = PRD 89, 064064 (2014) (verified). SU(2) gives discrete, self-dual
  variables give continuous spectra -- **the converse as an outcome**, attributed to the choice of
  variables rather than to compactness. Cf. Frodden-Geiller-Noui-Perez, arXiv:1212.4060 =
  EPL 107 (2014) 10005 (verified).
- **Maz'ya & Shubin**, arXiv:math/0305278 = Annals of Math. 162 (2005) 919-942 (verified). Sharp
  necessary-and-sufficient conditions for discreteness of Schrodinger spectra (sharpening
  Molchanov 1953, solving a Gelfand 1953 problem). **This is the strongest attack on the thesis:**
  discreteness on a *non-compact* R^n is achievable by a confining potential alone -- the harmonic
  oscillator is the standing counterexample (and is GeoVac's own Paper 24 object). The honest
  statement is "compact *resolvent* <=> discrete spectrum", a property of the operator, not the space.
- **Chirvasitu**, *(Quantum) discreteness, spectrum compactness and uniform continuity*,
  arXiv:2310.15139, to appear J. Noncommut. Geom. (verified; recency guard 2023-2025). Compact
  quantum group finite <=> discrete-type conditions. Conditional equivalences, explicitly **not** a
  blanket principle.
- **Luscher**, Commun. Math. Phys. 104 (1986) 177-206 and 105 (1986) 153-188 (verified). Finite
  volume => discrete spectrum, with the quantization condition as the exact discrete/continuum
  dictionary. The cleanest physics instance of "an exchange constant at the compactness boundary".
- **Rovelli & Vidotto**, *Compact phase space, cosmological constant, discrete time*,
  arXiv:1502.00278 = PRD 91, 084037 (2015) (verified). Compact phase space => finite-dimensional
  Hilbert space => discrete spectra, including discrete time. WH7-adjacent, reported as an
  outcome of the quantization rather than argued as a principle.
- **Sornette**, *Discrete-scale invariance and complex dimensions*, Phys. Rep. 297 (1998) 239-270
  (= cond-mat/9707012) (verified). The **inverse direction**: continuous scale invariance broken
  to a discrete subgroup, i.e. compactification of the log-radial axis; the exchange constant is
  literally 2*pi/ln(lambda) (complex exponents). Efimov is the physical instance; see also
  arXiv:1909.05505 (title verified via arXiv listing, **authors UNVERIFIED**).
- **Kontsevich-Zagier periods / Viu-Sos**: every non-zero real period is the volume of a *compact*
  semi-algebraic set -- Cresson & Viu-Sos, arXiv:1912.01751 = J. Theor. Nombres Bordeaux (verified
  listing); Viu-Sos, arXiv:1509.01097 = Int. J. Number Theory (verified listing). The
  number-theory converse: the transcendentals at issue (pi, log) *are* compact volumes. Never
  joined to the physics statement anywhere found.
- **Poisson summation**: the actual 2*pi-carrying identity underlying Matsubara, KK,
  decompactification limits, Selberg and Weyl. Standard; nLab "Poisson summation formula". No
  source found that names its 2*pi as the price of crossing the discrete/continuous boundary.

## 3. BACKGROUND
Peter-Weyl (1927); Weyl's law (1911) and its omega_d/(2*pi)^d constant; Selberg trace formula
(1956) and the compact/non-compact (cusp form vs Eisenstein) dichotomy; Rellich-Kondrachov;
Molchanov (1953). All are the theorem, none is the thesis.

---

## 4. OVERALL VERDICT

**"Discreteness is compactness" is a named principle in *mathematics* and an unnamed folk theorem
in *physics*.** In mathematics it is owned three times over -- Pontryagin duality
("compactness-discreteness duality"), the compact-resolvent characterisation, and Peter-Weyl -- but
always as a theorem about groups or operators, with no physical or transcendental content. In
physics it is **assembled per-instance**: every textbook re-derives it for the ring, the box, the
KK circle, the Matsubara circle, the finite-volume lattice, and none states the common cause. Two
exceptions were found, both narrow: **LQG states it for its gauge group** (C1), and **Dolce states
it for time** (C4).

**What no one found states.** (i) The *accounting* -- "each released axis costs one transcendental,
and the exchange constant is identifiable". Weyl's omega_d/(2*pi)^d is universally known and never
framed as a toll. (ii) The *corpus-wide* application -- a claim that the exact/rational skeleton of
a quantum system as a whole is its compact regime. (iii) The join between the physics statement and
the number-theoretic converse (periods = compact volumes).

**Converse-staters.** Only two, both partial: LQG's self-dual/SU(1,1) sector (continuous spectrum
from a non-compact group -- an outcome, not a law: C1, and arXiv:1306.3246), and
Kontsevich-Zagier/Viu-Sos (periods are compact volumes -- number theory, no physics attached).

**LQG collision risk: REAL, and it must be cited.** Liu-Noui (arXiv:1702.06793) states the thesis
as a one-sentence causal claim, names harmonic analysis as the mechanism, and exhibits the converse
in the same paper. Any GeoVac text presenting "compactness is the source of discreteness" as its
own observation without citing this is an overclaim. What survives as GeoVac's: the generalisation
beyond one gauge group, and the transcendental-toll corollary. What must be conceded: the slogan
itself, in physics, is owned by LQG for its own domain.

**Standing defect in the thesis (flag).** Maz'ya-Shubin makes "discreteness IS compactness" false
as a biconditional over spaces: a confining potential on non-compact R^n gives purely discrete
spectrum. Paper 18 SIII should say *compact resolvent*, not *compact manifold*. The corpus's own
Paper 24 (HO on S^5) is exactly this counterexample.

## 5. Searches run (support for "nobody states it")
Peter-Weyl + "discreteness of the spectrum" as physical principle; LQG compactness-origin-of-
discreteness; self-dual/SU(1,1) continuous area spectrum; Maz'ya-Shubin/Molchanov; Efimov/DSI/
Sornette; KK + Matsubara + finite-volume "unified statement"; Casimir/zeta where pi enters;
Connes compact resolvent; Pontryagin compact-vs-discrete; Kontsevich-Zagier periods as compact
volumes; Poisson summation as discrete/continuous bridge; thermal time/KMS; philosophy-of-physics
"why are quantum numbers integers"; three separate 2022-2026 recency sweeps for a unifying
principle or slogan. The recency sweeps returned **no** 2022-2026 statement of the general thesis.
