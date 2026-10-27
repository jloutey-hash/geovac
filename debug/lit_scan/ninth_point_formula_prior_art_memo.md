# Literature scan: prior art for the kite's ninth singular point, searched by its closed form

Date: 2026-10-27. Extends `debug/lit_scan/topology_rung_ninth_point_memo.md` (Q2). Inputs read: the NP record
`docs/predictions/2026-10-26_kite_ninth_point.md` (incl. Disposition), the PK record (labels, cuts, adjacency),
Paper 59 around `eq:kite_ninth_point`. Read-only plus two downloads (`debug/lit_scan/dl_ninth_formula/*/`) and two
sympy identity checks of a *published* formula against ours (scratchpad, no repo file). No derivation of the open
questions was attempted. Nothing committed; no paper or record edited.

Key: **[V]** = primary text fetched and read this session (page image where noted); **[V-abs]** = abstract page only;
**[UNV]** = search-snippet or memory only.

---

## 1. Headline: the formula is published, exactly

**Broadhurst, "Multivariate elliptic kites and tetrahedral tadpoles", arXiv:2212.01962 (v2, 14 Apr 2024)** [V; PDF
downloaded, Sec. 3 read from the rendered page image `dl_ninth_formula/broadhurst_2212.01962/pg-07.png`].

His kite: D = 4 Minkowski, five arbitrary masses, momenta (p1..p5) = (l, l-q, l-k, k, k-q). So edges 1,2 meet at one
external vertex, 4,5 at the other, edge 3 is the rung, and the internal vertices are {1,3,4} and {2,3,5} ("edge 1 is
not adjacent to edge 5. Neither is edge 2 adjacent to edge 4"). Thresholds at s = (m_j + m_k)^2.

Sec. 3 ("Zero-mass limits"), as m3 -> 0 "with m1 ≠ m4 and m2 ≠ m5", the two-particle-cut part of the derivative of
the discontinuity becomes (Eq. 24)

  Delta_{1,2}(s) sigma'_{1,2}(s) -> Re( (2s - s3) D_{4,5}(s) + L^_{4,5} + sum_{i=0,3} C_i (D_{4,5}(s) - D_{4,5}(s_i))/(s - s_i) ),

with (Eq. 25)

  **s3 = (m2^2 m4^2 - m1^2 m5^2)(m1^2 - m2^2 - m4^2 + m5^2) / M,   M = (m1^2 - m4^2)(m2^2 - m5^2),**

and s0 = 0 (Eq. 12, C0 = (m2^2 - m1^2)(m4^2 - m5^2)).

**Dictionary.** Our lines 1 (a), 2 (b) meet at V_in; 3 (g1), 4 (g2) at V_out; 5 (massless) joins the internal vertices
{1,3,5}, {2,4,5}. So Broadhurst 1,2,3,4,5 = our a, b, Coulomb, g1, g2 (his rung m3 is our massless line; his adjacent
pairs (1,4), (2,5) are our (a,g1), (b,g2)). Our z is Euclidean: thresholds at z = -(m_i + m_j)^2, so z = -s. Then

  -s3 = (A G2 - B G1)(A - B - G1 + G2) / ((A - G1)(B - G2)) = z9.

**Exact match, factor by factor, including the denominator.** Checked by hand and symbolically (identity = 0); at
P2 it gives -111/1120.

**What s3 is, in his framework.** Eq. 9 gives, for a massive rung, two "(possibly complex) constants"
s± = [m1^2 + m2^2 - 2m3^2 + m4^2 + m5^2 - alpha]/2 ± Delta_{1,3,4} Delta_{2,3,5}/(2 m3^2), with
alpha = (m1^2 - m4^2)(m2^2 - m5^2)/m3^2 - m3^2 and beta = (m1^2 m5^2 - m2^2 m4^2)(m1^2 - m2^2 - m4^2 + m5^2)/m3^2 (Eq. 8).
In his words they "locate leading Landau singularities [2] of triangles that form the kite" ([2] = Amati-Fubini,
Ann. Rev. Nucl. Sci. 12 (1962) 359). I checked symbolically that s± are exactly the two roots in s of the Cayley
determinant of the one-loop triangle {1,2,rung} with its two other legs on shell at m4 and m5 (equivalently, by the
kite's symmetry, the triangle {4,5,rung} with legs on shell at m1, m2; both give the same quadratic). The quadratic's
s^2 coefficient is proportional to m3^2. As m3 -> 0 one root runs to infinity like M/m3^2 and the other tends to the
root of the now-linear equation, which is s3. In our notation: the Cayley matrix

  [[2A, A+B-s, A-G1], [A+B-s, 2B, B-G2], [A-G1, B-G2, 0]]

has a determinant that is linear in s, with unique root s = -z9. So **z9 is the anomalous threshold (leading Landau
singularity) of the triangle formed by the two lines of one two-particle cut and the massless line, with the
other two lines on shell as its external legs.**

**Why the kite's own Landau analysis misses it, per Broadhurst.** It enters only through the dispersion relation of
the triangle that sits inside a two-particle cut: across the {3,4} cut the triangle's legs are pinned on shell, so the
discontinuity carries the triangle's singularity. Sec. 5, p. 12-13: "The kite has no anomalous thresholds, yet the
triangles that it contains may have dispersion relations with anomalous thresholds [2] and complex branch points."
For m3 = 0 he notes "There is no anomalous term, since both Delta_{1,3,4} and Delta_{2,3,5} are real", i.e. on the
physical sheet s3 does nothing. But the term C3 (D_{4,5}(s) - D_{4,5}(s3))/(s - s3) in Eq. 24 is regular at s3 only on
the principal branch of D_{4,5}. Continued round the {4,5} threshold (where D_{4,5} jumps by -2πi/Delta_{4,5}, his
Eq. 7), the numerator no longer vanishes and s3 becomes a pole of sigma'. That is a singular point of any
differential equation the integral satisfies, which is what the NP operator saw. This last step is my reading of his
Eq. 24, not a sentence of his. **[my inference]**

**Poles at equal adjacent masses (task item 4).** Directly in his formula: M = (m1^2 - m4^2)(m2^2 - m5^2). In the
Cayley matrix the entry for (massive internal line, massless line) is m_i^2 + 0 - p_leg^2 = A - G1 (or B - G2); it
vanishes when an internal mass equals the on-shell mass of the adjacent external leg at the vertex where the massless
line attaches. That is the soft (IR) configuration of a vertex with massless exchange between equal masses, the
on-shell QED vertex. Broadhurst treats this case separately (Eq. 27, "As m3 -> 0 with m1 = m4"): the s3 term is gone
and only the s = 0 term remains, i.e. the point has left to infinity. The IR reading is mine; the exceptional status
of m1 = m4 and the denominator M are his.

**Zeros (A G2 = B G1, or A + G2 = B + G1)** are the two numerator factors of his Eq. 25. When s3 = 0 it lands on his
s0 = 0, the other pole of Eq. 24 (coefficient C0 = (m2^2 - m1^2)(m4^2 - m5^2)). This is the base-point case
(a g2 = b g1): the z9 that NP1(c) found sitting on z = 0 is s3 merging with s0. **[my inference from Eqs. 12, 24]**

**Caveats on transfer.** His kite is D = 4 with unit indices; ours is D = 3 with indices (3/2, 3/2, 2, 2, 1). Landau
loci depend only on masses and topology, not on D or the indices (the location transfers), but the local exponent
(our -7/2) and whether a given sheet's singularity survives in a given D are not transferred. He studies the
discontinuity sigma', not a differential equation, and never calls s3 a singular point of an ODE.

**Closeness: EXACT FORMULA, same mechanism stated in the source.**

---

## 2. Same formula family, earlier and less explicit

**Bauberger and Böhm, "Simple one-dimensional integral representations for two-loop self-energies: the master
diagram", hep-ph/9501201** [V, html text via fetch]. The master diagram, arbitrary masses, via the dispersion relation
for the vertex function C0 that sits inside it. Sec. 3.2: "The Landau equations show that a leading singularity and
thus an anomalous threshold occurs in C0" when m3^2 < m3thr^2, with branch point
s1 = (p2^2 + p3^2 + m1^2 + m2^2 - m3^2)/2 - [(p2^2 - m2^2)(p3^2 - m1^2) + l(p2^2,m2^2,m3^2) l(p3^2,m1^2,m3^2)]/(2 m3^2)
(their labels; l = product of the two Källén square roots), with six cases for real/complex s1 relative to the normal
threshold. This is the same triangle Landau point, written for the triangle inside the master diagram. Their m3 -> 0
limit is not given in what I read. **Closeness: same object (triangle anomalous threshold inside the master
diagram), formula in a different parametrization; massless limit not stated.**

Other names in the brief, not matched: Broadhurst-Fleischer-Tarasov, Davydychev-Tausk, Berends-Davydychev-Smirnov
(small/zero-threshold expansions; OSTI listings only [UNV]), Remiddi-Tancredi arXiv:1602.01481 (equal-mass QED kite;
listing only [UNV]), Adams-Bogner-Schweitzer-Weinzierl 1607.01571 (equal-mass; [V-abs] in the earlier memo). None
was needed once Broadhurst matched, and I did not read them for this point.

---

## 3. Non-leading / missed singularities (task item 2): what applies, what does not

- **Not a second-type singularity.** z9 is a first-type (anomalous-threshold) Landau point of a one-loop triangle,
  appearing on a non-physical sheet of the kite because the triangle lives inside a two-particle cut. The earlier
  memo's "mixed/second type" guess is superseded. The Euler-characteristic drop and the PLD miss are consistent with
  this: the critical-point count sees all sheets; the kite's physical-sheet leading Landau equations do not have this
  solution (with α3 = α4 = 0 the q-loop equation forces q = 0, which needs A = G1).
- **General framework for singularities of discontinuities on other sheets.** Hannesdottir-McLeod-Schwartz-Vergu,
  arXiv:2211.07633 [V-abs]: Pham's on-shell-space picture, sequential discontinuities and compatibility of Landau
  singularities; proofs at generic masses. Berghoff-Panzer, arXiv:2212.06661 [V-abs via listing]: relative
  Picard-Lefschetz, hierarchy constraints, degenerate cases "e.g. the massless triangle". These are the language in
  which "a triangle anomalous threshold inside a cut" is a singularity of a sequential discontinuity. Neither is
  specific to the kite or contains the formula.
- **Why PLD/face scans miss it, and Whitney/Euler catches it**: as in the earlier memo (Fevola-Mizera-Telen
  2311.16219, 2311.14669; Helmer-Papathanasiou-Tellander 2402.14787; all [V-abs] there). Not re-fetched. Their
  examples (parachute, slashed box) are a **same-class phenomenon**, not this formula.
- **Massless line.** The massless line matters only by making the triangle's Cayley equation linear (one finite
  root instead of two). With mass m5 > 0 the point is the pair s±(m5), so it moves; that alone explains why the
  Euler count at m5^2 = 1/100 no longer drops at z9 (earlier memo, Q2(c)).

---

## 4. Atomic three-centre / Hylleraas integrals (task item 3)

- Pachucki-Puchalski-Remiddi, "Recursion relations for Hylleraas three-electron integral", physics/0405057 [V, full
  text]. IBP in momentum space for the six-exponent K4 integral (Fromm-Hill generic integral). Mentions that Harris's
  form of Fromm-Hill removed "many spurious singularities" and thanks Korobov for "two-electron recursions in the
  singular case". No singular-point formula in the exponents.
- Fromm-Hill, PRA 36, 1013 (1987); Harris, PRA 55, 1820 (1997): APS pages returned 403; content known only from
  search snippets (multi-valued dilogarithms, branch tracking, "formal cancellation of the individual-term
  singularities") [UNV].
- **No match found** for a non-threshold singular point in the exponents. **Same-class at most** (spurious or
  branch-tracking singularities of the generating integral), and that is unverified.

---

## 5. Ranking

| rank | source | closeness | what it gives |
|:--|:--|:--|:--|
| 1 | Broadhurst arXiv:2212.01962, Eq. 25 (+ Eqs. 8-9, 12, 24, 27) | **exact formula** (z9 = -s3), mechanism stated | s3 as the m3 -> 0 limit of the triangle Landau points s±; M in the denominator; m1 = m4 case separate |
| 2 | Bauberger-Böhm hep-ph/9501201, Sec. 3.2 | same object, other parametrization | triangle anomalous threshold s1 inside the master diagram's dispersion relation, massive rung |
| 3 | Hannesdottir et al. 2211.07633; Berghoff-Panzer 2212.06661 | framework | singularities of (sequential) discontinuities; degenerate massless-triangle cases |
| 4 | FMT / HPT (earlier memo) | same class of phenomenon (missed by PLD, caught by Euler drop) | the diagnostic, not the formula |
| 5 | Fromm-Hill / Harris / PPR | related at most | atomic K4 integral; no matching point found |

## 6. Suggestions for the derivation (not done here)

1. **Derive z9 as the Cayley determinant of the triangle {1,2,5} with legs 3,4 on shell** (or {3,4,5} with legs 1,2),
   via the {3,4}-cut discontinuity: in D = 3 the two-body phase space factorizes and the triangle depends only on
   (s, G1, G2), so disc = phase space × triangle. Cite Broadhurst Eq. 25 as [PRIOR ART] for the location.
2. **NP2 is answerable from Eq. 9 without new reductions**: with m5^2 = eps the point splits into s±(eps), one tending
   to s3, the other to infinity like (A - G1)(B - G2)/eps. Prediction to check against the Stage B unscored data: at
   z = z9, the eps-roots of the triangle quadratic should be eps = 0 and one other value; compare it with the
   eps-root -1982123/372960 seen in 120 entries. I did not compute this.
3. **Exponent**: compare -7/2 with the triangle's leading-singularity exponent in D = 3 with indices (3/2, 3/2, 1),
   composed with the {3,4} phase space and the operator's derivative shifts.

## Sources

Verified this session: arXiv:2212.01962 (full text, Sec. 2-3, 5; page 7 image); hep-ph/9501201 (html text);
physics/0405057 (full text); arXiv:2211.07633 (abstract); arXiv:2212.06661 (abstract via listing). From the earlier
memo, not re-fetched: 2311.16219, 2311.14669, 2402.14787, 1607.01571. Unverified: Fromm-Hill PRA 36 1013, Harris PRA
55 1820 (403), Amati-Fubini 1962, Remiddi-Tancredi 1602.01481 content, Berends-Davydychev-Smirnov papers.
Downloads: `debug/lit_scan/dl_ninth_formula/broadhurst_2212.01962/`, `debug/lit_scan/dl_ninth_formula/ppr_physics0405057/`.
