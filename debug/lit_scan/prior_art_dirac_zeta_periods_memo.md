# Prior-art scan: period classification of S^3 Dirac/Laplace spectral zeta values and Mellin moments

Date 2026-11-01. Scope: Paper 18 §III.7 tie-break ("a value reached from two slots is placed by the period ring it lands in"), Paper 55 placements, Paper 28 T9. Read-only scan (22 web searches). Downloads are in `debug/lit_scan/dl_prior_dirac_zeta/` (one subfolder per item, read with pdftotext). Not re-reported: works the corpus already cites (Fathizadeh-Marcolli 2016, Deligne 2010, Glanois 2015, Brown 2012, Kontsevich-Zagier, Camporesi-Higuchi, Chamseddine-Connes 1997/2010, Klebanov-Pufu-Safdi, Dowker 1977/1994). Chamseddine-Connes *Uncanny* is cited by Paper 51 only, so it appears below where other papers owe it.

Confidence labels: VERIFIED = read in the downloaded text; ABSTRACT ONLY; SECONDARY.

## Headline

1. **The tie-break is not in the literature.** No source I found classifies Dirac or Laplace zeta values on S^3 as periods by the ring they land in, or gives a rule for a value that two different Mellin moments reach. The nearest work (Fathizadeh-Marcolli; Fathizadeh-Ghorbanpour-Khalkhali) is about heat-coefficient rationality and motives, not zeta values at regular points. Nothing I found contradicts the rule.
2. **Two closed forms are already published, without credit in the corpus:** the closed form of D(s), and the two-term heat trace on S^3.
3. **Published material supports both stated reasons.** The spectral-action expansion uses only the residues and ζ_D(0), never regular positive values. Brown's motivic coaction is trivial on even zeta values and non-trivial on odd ones.

## Q1. Dirac zeta and eta on spheres and lens spaces

**F1. Ball & Marcolli, "Spectral action models of gravity on packed Swiss cheese cosmology", Class. Quantum Grav. 33 (2016) 115018; arXiv:1506.01401.** VERIFIED. Prop./Cor. 2.4-2.5 gives ζ_{D_{S^3}}(s) = 2ζ(s−2,3/2) − ½ζ(s,3/2). Verbatim: "The Hurwitz zeta function ζ(s, q) has a simple pole at s = 1 with residue one, hence ζ_{D_{S^3}}(s) has simple poles at s = 1 and s = 3, respectively with residues ... −1/2 and ... 2 ... ζ_{D_{S^3}}(0) = 2ζ(−2, 3/2) − ζ(0, 3/2)/2 = 0". The same formula is repeated in Guicardi-Marcolli, arXiv:2110.08878 (VERIFIED). *Bearing:* this is the corpus's D(s) exactly, and it was already published. **Credit owed** at Paper 55 eq. `eq:dirac_dirichlet` (line ~280, "The associated spectral zeta function is D(s) := ... = 2ζ(s−2,3/2) − ½ζ(s,3/2)") and at Paper 55's proof sketch (~L803). Neither paper does any period classification.

**F2. Eckstein & Iochum, *Spectral Action in Noncommutative Geometry*, SpringerBriefs Math. Phys. 27 (2018); arXiv:1902.05306, Example 3.15 (and Example 2.41, Remark 3.16).** VERIFIED (π factors were lost in the text extraction). Verbatim: "ζ_{D²}(s) = ... = 2ζ(2s − 2, 3/2) − ½ζ(2s, 3/2)"; the residues at 3/2 and 1/2 are given, ζ(−2k,1/2) = −B_{2k+1}(1/2)/(2k+1) = 0, so all residues at negative integers vanish; "the asymptotic expansion of the heat trace associated with the square of the standard Dirac operator on any odd dimensional sphere is convergent for all t > 0, but is not exact." The remainder is tied to the Jacobi theta identity. *Bearing:* this is prior art for:
- Paper 18 ~L1338: "the continuum spectral zeta ζ_{D^2}^cont(s) = 2ζ(2s−2,3/2) − ½ζ(2s,3/2) has a simple pole at s = 3/2 with ... residue 1". Watch the normalization: Eckstein-Iochum quote residue 2 for their convention.
- Paper 28 §"SD two-term exactness for the Dirac heat kernel" (eq. `eq:sd_two_term`, proved by the θ₂/θ₄ modular identity).
- Paper 51's G1 two-term exactness.

**Credit owed:** Eckstein-Iochum (the same theta route) and Chamseddine-Connes, *The uncanny precision of the spectral action*, CMP 293 (2010) 867, arXiv:0812.0165 (VERIFIED: "all higher order terms a_{2n} vanish"). Paper 51 cites *Uncanny*; Papers 28 and 18 do not.

**F3. Bär**, "The Dirac operator on space forms of positive curvature", J. Math. Soc. Japan 48 (1996); **Bär**, "Dependence of the Dirac spectrum on the spin structure", Sémin. Congr. 4 (2000), Thm 5.2: η-invariant of spherical space forms. **Bär & Schopka**, "The Dirac determinant of spherical space forms", in *Geometric Analysis and Nonlinear PDEs* (Springer 2003). ABSTRACT ONLY / SECONDARY: the Potsdam server refused the connection. These give spectra, η and determinants. I found no period classification in them. On lens spaces the spectrum splits into arithmetic progressions, so lens-space Dirac zeta values are Hurwitz functions at rational shifts a/N. That is cyclotomic level N, which extends beyond the N ∈ {2,3,4,6,8} that Deligne and Glanois cover. This last point is my inference and is not stated in any source.

**F4. Dowker, "Spherical Dirac GJMS operator determinants", arXiv:1310.5563 (J. Phys. A 2015).** VERIFIED (abstract and conclusion). The odd-dimensional spinor log-determinant is "the traditional combination of Dirichlet η-functions and log 2"; for the other factorisation it is "a sum of Riemann, not Dirichlet, zetas. Also there is no log 2". **Dowker, arXiv:2403.09576 (2024).** VERIFIED: "In the simplest case of the 3-sphere, C = ζ_R′(−2) = −ζ_R(3)/(4π²)", which is the constant term of the high-temperature expansion on S^1×S^3. *Bearing:* both concern derivative values (ζ′(0)) rather than values. Paper 18 ~L2146 already has the matching caveat: "functional determinants on round spheres generally carry ζ(odd) and log 2". Citing Dowker there would be optional.

## Q2. Arithmetic of spectral-action coefficients

**F5. Fathizadeh, Ghorbanpour & Khalkhali, "Rationality of spectral action for Robertson-Walker metrics", JHEP 12 (2014) 064; arXiv:1407.5972.** ABSTRACT ONLY. They prove Chamseddine-Connes' conjecture (Chamseddine-Connes, "Spectral action for Robertson-Walker metrics", JHEP 10 (2012) 101; arXiv:1105.4637) that the a_{2n} are polynomials in a(t) and its derivatives with rational coefficients. *Bearing:* this is the published root of "M2 = π^{2k}·Q after dividing out the volume normalisation". It is the pre-motivic statement that Fathizadeh-Marcolli refine. It would be a natural extra cite in Paper 18 ~L1365, next to Fathizadeh-Marcolli. It says nothing about zeta values.

**F6. Fathizadeh-Marcolli (already cited), re-read.** VERIFIED: the coefficient prefactors are "in Q[(2πi)^{−1}]". The paper makes no statement about zeta values, odd zeta, or two-moment ambiguity.

**F7. Iochum, Levy & Vassilevich, "Global and local aspects of spectral actions", J. Phys. A 45 (2012) 374020; arXiv:1201.6637.** ABSTRACT ONLY. Local (heat-kernel) versus global (non-local) parts of the spectral action. It has no arithmetic content and no tie-break. Dabrowski-Sitarz: I found no S^3 spectral-action arithmetic paper. Eckstein-Iochum-Sitarz on the Podleś sphere (CMP 2014, arXiv:1307.4188) has an exact heat trace but no period content (ABSTRACT ONLY).

**Q2 conclusion:** I found no published rule for a value that two Mellin moments reach.

## Q3. Build-on: separating heat coefficients from regular values

**F8. The spectral-action expansion sees only residues and ζ_D(0).** Ball-Marcolli eq. (1.2), VERIFIED: "Tr(f(D/Λ)) ∼ Σ_{β∈Σ⁺} f_β Λ^β ⨍|D|^{−β} + f(0) ζ_D(0)". With the standard expansion over the dimension spectrum, heat coefficients are residues and values at non-positive points. On S^3 these are Bernoulli-rational times the volume factor. Values at regular positive points never enter. *Bearing:* this **supports reason 2** of the tie-break ("ζ_{D²} at half-integer s ≥ 5/2 is not a Seeley-DeWitt coefficient") with a standard citation. It is a local/global split, not a motivic one.

**F9. Brown, "Motivic periods and P¹∖{0,1,∞}", Proc. ICM 2014; arXiv:1407.5165, §2.3.** VERIFIED. Verbatim: "It is important to note that ζ^m(2) is non-zero, whereas ζ^u(2) = 0"; "D ζ^m(n) = [ζ^u(n)] ⊗ 1 ... where we recall that ζ^u(2n) = 0. One easily shows that ζ^u(2n+1) ≠ 0 for n ≥ 1". The "≠" sign was lost in extraction; the meaning is clear from context. *Bearing:* this is the precise motivic reason for **reason 1**. The coaction kills π^{2k} (pure Tate) and acts non-trivially on ζ(odd) (genuine extensions Ext¹(Q(0),Q(2n+1)) = Q; also Dupont, "Odd zeta motive and linear forms in odd zeta values", Compositio 154 (2018), arXiv:1601.00950, ABSTRACT ONLY). It is a better citation than Brown 2012 for the line "odd zeta values violate M2's own signature". **Caveat for the PI:** π^{odd} (M1's ring Q[π,1/π]) is also a pure-Tate period up to i. So "pure Tate" separates M2/M1 from M3, not M2 from M1. The tie-break wording "M2's pure-Tate ring ⊕π^{2k}Q" is correct but should not be read as "only M2 is pure Tate". The rule's placement of rationals under M2 is a convention. Q is the period ring of Q(0), so the literature neither supports nor contradicts it.

**F10. Chamseddine-Connes, *Uncanny*, §2.3 (VERIFIED, arXiv:0812.0165).** Verbatim: "one cannot expect to use the Poisson summation in the same way for the 4 sphere. The problem ... is that, e.g. for the heat kernel, one is dealing with a function like |x|e^{−tx²} which is not smooth and whose Fourier transform does not have rapid decay". *Bearing (build-on, my inference):* on S^3 the CH multiplicity 2λ² − ½ is even in λ. So slots k = 0, 2 have even weight, Poisson summation is exact up to e^{−π²/t}, and the values are π-powers. Slot k = 1 (|D| e^{−tD²}) has weight 2λ³ − λ/2, which is odd in λ. This is exactly Chamseddine-Connes' non-smooth case, and its values are odd zeta. In general D(m) = Σ g(λ) λ^{−m}, so the tie-break "even m → M2, odd m → M3" is equivalent to **the parity in λ of the summand of the Mellin integrand**. That property belongs to the value, not to the slot. It gives an analytic restatement of the rule that is independent of which slot reached the value, with a published anchor.

## Credits owed (summary)

| Locus | Item |
|---|---|
| Paper 55 `eq:dirac_dirichlet` (~L280) and proof sketch (~L803) | Ball-Marcolli 2016 (arXiv:1506.01401); Eckstein-Iochum 2018 Ex. 3.15 |
| Paper 18 ~L1338 (continuum ζ_{D²}, residue at 3/2) | Eckstein-Iochum Ex. 3.15 (check residue normalization) |
| Paper 28 §SD two-term exactness (`eq:sd_two_term`, theta proof) | Chamseddine-Connes 2010 *Uncanny*; Eckstein-Iochum Ex. 2.41/3.15 |
| Paper 18 ~L1365 (F-M paragraph) | optionally Fathizadeh-Ghorbanpour-Khalkhali 2014 + Chamseddine-Connes RW 2012 |
| Paper 18 tie-break reason 1 (~L1425) | Brown ICM 2014 (ζ^u(2n)=0, ζ^u(2n+1)≠0) |
| Paper 18 tie-break reason 2 | standard expansion (Ball-Marcolli (1.2) / Chamseddine-Connes 1997, already cited) |

## Top build-on items

1. **The λ-parity restatement (F10):** classify D(m) by the parity of g(λ)λ^{−m}. This removes the slot ambiguity outright, and Chamseddine-Connes' smooth/non-smooth Poisson obstruction is the published analytic mechanism.
2. **Brown's coaction (F9) as the motivic justification** for odd zeta ∉ M2. Add the caveat that π^{odd} is also pure Tate.
3. **Lens spaces (F3, inference):** Dirac zeta on L(N;…) gives Hurwitz values at shifts a/N, which is cyclotomic level N. This tests whether M3's "level ≤ 4" is an S^3 artefact, and N ∉ {2,3,4,6,8} exits the Deligne-Glanois coverage.

Process note: one download command ran with the download folder as the working directory, for curl and pdftotext only. No interpreter ran there; all Python reads used `python -I` from outside.
