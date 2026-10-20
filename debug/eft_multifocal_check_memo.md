# EFT / precision-literature check of the finite-nucleus doors (J14, J16 / HF-4, Dirac)

**Date:** 2026-10-19. Read-only research; this file is the only output. Not a /walls or /qa pass.
**Question:** how much of the multi-focal wall's priced finite-nucleus mechanism ("labels survive; one
consumed moment per order; then a universal ln(2ZR) + 2 gamma_E") is known physics, and does the corpus
credit it.

**Short verdict.** Almost all of it is known. The mechanism *is* the standard EFT (NRQED) statement:
nuclear structure enters only through matching coefficients (moments), and logs come from iterating
lower-order contacts, which fix their size. The corpus credits the leading terms but no EFT source, and
not Karshenboim, who published the J16 psi(n) form for arbitrary n in 1997/2001. The check also finds
three overstatements (section 4).

Sources downloaded and read for this memo (primary texts, arXiv versions):
- [PPY25] Pachucki, Patkos, Yerokhin, "Recoil corrections with finite nuclear size in hydrogenic
  systems", arXiv:2502.09455v2 (2025); re-derives Friar's (Z alpha)^6 nonrecoil result.
- [K97] Karshenboim, Phys. Lett. A 225, 97 (1997), arXiv:hep-ph/9608484.
- [K01] Karshenboim, "2s hyperfine structure in hydrogen atom and helium-3 ion", Springer (2001) p. 335,
  arXiv:physics/0102085.
- [KI02] Karshenboim, Ivanov, Phys. Lett. B 524, 259 (2002), arXiv:hep-ph/0109128.
- [EGS01] Eides, Grotch, Shelyuto, Phys. Rep. 342, 63 (2001), arXiv:hep-ph/0002158.
- [HP11] Hill, Paz, Phys. Rev. Lett. 107, 160402 (2011), arXiv:1103.4617.
- [FP05] Friar, Payne, Phys. Rev. C 72, 014002 (2005), arXiv:nucl-th/0504015.
- [Sh05] Shabaev et al., arXiv:physics/0510083 (quotes Shabaev, J. Phys. B 26, 1103 (1993)).
- [Vol06] Volotka, TU Dresden thesis (2006), App. A, reproducing Volotka-Shabaev-Plunien-Soff, Eur.
  Phys. J. D 23, 51 (2003), the EPJ D DOI Paper 23 cites (resolved via Crossref).
- [KUP12] Kotov, Uchoa, Pereira et al., Rev. Mod. Phys. 84, 1067 (2012), arXiv:1012.3484.

Not opened (cited only through the works above): Friar 1979, Zemach 1956, Bohr-Weisskopf 1950,
Caswell-Lepage 1986, Pieper-Greiner 1969, Zel'dovich-Popov 1972. Gershtein-Zeldovich, Pineda and
Jentschura were not checked and are not cited.

---

## 1. Item-by-item status

| # | Item (corpus) | Status | Source + exact location | Verified in primary? |
|:-|:-|:-|:-|:-|
| 1 | Labels l, m survive; n released; E = -Z^2/(2(n-delta)^2) (quantum-defect reading) | KNOWN in substance (standard QDT; Seaton) | Corpus already cites Seaton 1983. As a reading of the finite-size shift it rewrites -Z^2 delta/n^3 | n/a (rewriting) |
| 2 | J14 leading s defect -(2/3)Z^2<r^2> = -(2/5)x^2, same every n | KNOWN | PPY25 Eq. (3) and (72): (2 pi/3) Z alpha phi^2(0) r_C^2; EGS01 Eq. (158), (186); Sh05 Eq. (3) (Shabaev 1993), whose non-relativistic limit (Z alpha)^2/(10n) (2 Z alpha R/n)^2 gives exactly 2/5 | yes (PPY25, EGS01, Sh05); Friar 1979 not opened |
| 3 | J14 l = 1 defect c_1 = 1/105, factor (1 - 1/n^2) | KNOWN (l = 1); general-l c_l product NOT FOUND | PPY25 Eq. (77): non-relativistic p term (1/9)(1/n^3 - 1/n^5)(1/5)<r^4>. For the sphere <r^4> = 3R^4/7 gives (1/105)(n^2-1)/n^5, exactly the corpus's value. Note that the l = 1 defect consumes <r^4>, not R | yes (PPY25) |
| 4 | x^3 term = Friar (1/3)Z^3<r^3>_(2) | KNOWN, credited | PPY25 Eq. (6)-(7); EGS01 Eq. (168)-(169); Friar 1979 | yes (PPY25, EGS01) |
| 5 | x^4: -delta_0^2 [ln(2ZR) + 2 gamma_E] + rational, "universal" | PARTLY KNOWN; the "universal 2 gamma_E" part is contradicted (see 4.1) | PPY25 Eq. (28) and (30), second (non-relativistic) line: (4/(9n^3)) r_C^4 [1 - 1/n - ln n + gamma + psi(n) + kappa_2 + ln(m r_C Z alpha)] + r_CC^4/(15 n^5), with kappa_2 model dependent (Eq. 32, exponential model). The log coefficient matches the corpus's -delta_0^2. The gamma coefficient is 1 there, not 2, and EGS01 Eq. (188) writes the relativistic analogue with <ln(2 m r Z alpha/n)> averaged over the nuclear density, so the constant beside the log depends on the profile | yes (PPY25, EGS01) |
| 6 | Other profiles (shell exact; Gaussian, a pi^(-1/2) from the profile) | KNOWN in substance | PPY25 Table II (exponential versus Gaussian model values); the pi^(-1/2) is the Gaussian's Friar moment | yes |
| 7 | Continuation across threshold with J13's factor; Seaton limit | KNOWN, credited (Seaton) | Standard QDT | not re-checked |
| 8 | Dirac 1s shift scales as R^(2 gamma) | KNOWN | Sh05 Eq. (3), quoting Shabaev, J. Phys. B 26, 1103 (1993): Delta E_ns = (alpha Z)^2/(10n) [1 + (alpha Z)^2 f_ns] (2 alpha Z R/(n lambda_C))^(2 gamma) mc^2 | via a citing work (same author group) |
| 9 | At Z alpha = 1, W -> 0 slower than any power of R | NOT FOUND (not searched deeply; Zel'dovich-Popov is the likely home) | - | no |
| 10 | Z_cr = 169.7 for R = 1.2(2.5Z)^(1/3) fm | KNOWN, credited | Near 170: Pieper-Greiner 1969; Zel'dovich-Popov 1972 (corpus cites both); about 170 for a ball and 173 for a shell is quoted in later supercritical-QED papers | via citing works only |
| 11 | Diving radii in ratio e^(pi/sigma) (discrete scale invariance) | PARTLY KNOWN | Log-periodic supercritical spectrum is standard: KUP12 Eq. (4.8), E_n proportional to (1/R) e^(-pi n/sqrt(g^2 - g_c^2)) (2D massless Dirac, citing Shytov et al. 2007). The 3D statement (fall to the centre) is in the Zel'dovich-Popov tradition, not verified. The diving-radius ratio follows directly | 2D form yes (KUP12); 3D no |
| 12 | J16 leading hyperfine change -2 Z r_Z, same every n | KNOWN, credited | Zemach 1956; EGS01 Eq. (347)-(348); Vol06 Eq. (2.16), (7.14), (7.26) (moment form for the sphere charge) | via EGS01 and Vol06 |
| 13 | Point charge: entire function of x, coefficients polynomial in 1/n^2; x^2 coefficient 1 + 1/(5n^2) (ball) | PARTLY KNOWN | The n-dependent magnetic-radius term: K01 Eq. (26), -(4/3)(Z alpha)^2 (R_M/R_E)^2 (n^2-1)/(4n^2) (m R_E)^2 E_F, i.e. +(1/3) Z^2 <r^2>_M / n^2, which is 1/(5n^2) for the ball; also K97 below Eq. (7). That the series has no log at any order is not stated, but it is elementary | yes (K01, K97) |
| 14 | J16 (iii): (4/5) x^2 ln x = -2 x J14's leading defect, set by the charge alone | KNOWN (coefficient); the "-2 times J14" link is implicit | K97 Sec. 2.5, "structure-charge" term: nu_F x [-(2/3)(Z alpha m_e)^2 r_p^2 log(1/(Z alpha)^2)], at leading-log accuracy (the log argument uses m_e r_p ~ alpha for hydrogen), i.e. (4/3) Z^2<r^2> ln x = (4/5) x^2 ln x for the sphere. The link to the defect is immediate in the contact-term language (both carry (Z alpha m R_E)^2) but is not stated as a relation | yes (K97) |
| 15 | J16 (iv): n-dependence (4/5)[psi(n) - ln n - 1/n] + profile constant ("post hoc") | KNOWN, for arbitrary n | K01 Eq. (26): D_n1(Nucl) contains +(4/3)(Z alpha)^2 [psi(n) - psi(1) - ln n + (n-1)/n - ...] (m R_E)^2 E_F, which is the corpus's form, (4/3)<r^2> = 4/5 R^2; K97 Eq. (8), same bracket printed with 2/3 (probably a misprint: its own log coefficient is 4/3, and K01 has 4/3; the corpus measures 4/3); KI02 Eq. (11) applies it at n = 2. Both K97 and K01 state that they used an effective non-relativistic theory: second-order perturbation theory with two delta-like potentials | yes (K97, K01, KI02) |
| 16 | J16: unregistered x^3 ln x term, n-independent, profile-dependent, constant unidentified | NOT FOUND | Searched K97, K01, KI02, EGS01 Sec. XXXIV, FP05, Vol06 App. A. Vol06 works to linear order in m R plus relativistic (Z alpha)^2 terms, and its (Z alpha)^3 m R ln(m R) terms are relativistic, not this non-relativistic x^3 ln x. EFT gives a candidate for the constant (2.3) | n/a |
| 17 | Paper 23's open note: does (iii) appear in EPJ D 23, 51 (2003)? | Answered: no, as far as the reproduction shows | Vol06 App. A (which states it reproduces [121] = Volotka-Shabaev-Plunien-Soff 2003) gives the non-relativistic correction only at linear order (Eq. 7.14, 7.23, 7.26) plus relativistic (Z alpha)^2 terms (Eq. 7.13, 7.25). The x^2 ln x term is in K97, not in Volotka et al. | via the thesis, not the EPJ D text |
| 18 | "Mechanism": short-distance moments consumed one per order, labels survive, logs from the composition | KNOWN as the EFT matching statement | HP11: NRQED Wilson coefficients fixed by matching; r_E defined through c_D (Eq. 4 and the "Bound state energies" paragraph); a structure-dependent contact d_2 "not determined by measured proton form factors". NRQED itself: Caswell-Lepage 1986 (cited by HP11 ref. [11] and EGS01 ref. [17]) | HP11 yes; Caswell-Lepage no |
| 19 | Transcendental bookkeeping (no pi at leading order; weight one from 4th order; the profile's class otherwise) | The content is in the known formulas; the atlas bookkeeping as such is NOT FOUND | PPY25 Eq. (30)-(32): ln, gamma, psi(n), and model logs ln(4/3), ln 2, ln 3. The classification is the corpus's framing | n/a |

---

## 2. Q2: is "one consumed moment per order plus a universal log" the EFT statement?

**What the EFT statement is.** The nucleus is much smaller than the Bohr radius (x = Z alpha m R much
less than 1). In NRQED with a structured nucleus [HP11; Caswell-Lepage], the nucleus enters the
electron's dynamics only through local operators: a contact delta(r) with coefficient proportional to
<r^2> (c_D), a nabla^2 delta(r) with <r^4>, the Friar contact with <r^3>_(2), and, for the hyperfine
factor, spin contacts with the Zemach and magnetic moments. Four consequences follow:

- (E1) Each coefficient is short-distance data that the long-distance theory cannot fix. This is the
  corpus's "consumed moments".
- (E2) The number of new numbers at a given order equals the number of independent operators at that
  dimension. It can be more than one.
- (E3) Logarithms come from UV divergences when lower-order operators are iterated (second-order
  perturbation theory with two delta potentials, exactly K97/K01's method). Their coefficients are fixed
  by products of lower-order coefficients, so a log is never new data.
- (E4) The constant beside each log is scheme- and profile-dependent. Only the log coefficient is
  universal. The n-dependence comes entirely from long-distance matrix elements (the reduced Coulomb
  Green function at E_n, which brings psi(n) - ln n), never from the coefficients.

**How the corpus statement compares.** On the mechanism it is the same statement, (E1) plus "labels
survive". It is *weaker* than EFT because it does not say (E3), that log coefficients are predicted by
lower orders. It is *stronger* than EFT, and here wrong, in two places: "one per order", against (E2),
and "a universal ln(2ZR) + 2 gamma_E", against (E4).

**What EFT predicts that the corpus statement does not, checked against the corpus's own numbers:**

1. *The x^5 log of the energy defect.* EFT says it is the cross term of the leading and Friar contacts,
   -2 delta_0 delta_3 x^5 ln x. With delta_0 = -2/5 and delta_3 = 32/63 this is 128/315, which is
   exactly the coefficient the corpus derived (finite_nucleus record, "Owed item": x^5 [(128/315) ln(2x)
   + ...]). The corpus has the number but does not state the rule.
2. *The hyperfine x^3 ln x coefficient (the "unregistered" term).* Iterating the contacts gives two
   pieces: (-2 delta_0)(-2 Z r_Z) from the Zemach operator crossed with the charge contact, plus
   -2 delta_3 from the charge's Friar contact, i.e. (4/5) x (Zemach coefficient) - 64/63. The corpus
   measured (4/5) x (Zemach coefficient) - 1.0147. The first piece is exactly the profile-dependent part
   it reports. The second is profile-independent, as measured, and 64/63 = 1.01587 lies within the stated
   1e-2 spread across column sets (1.2e-3 relative). **This is a candidate identification, not a result:**
   it needs the curve-fit audit and a fit with tighter columns (or an exact derivation) before it is used.
3. *An n-dependent x^4 energy defect.* (E4) predicts psi(n)-type n-dependence at the order where contacts
   are iterated, which for the energy is x^4. That is the same order at which J16 found psi(n) for the
   hyperfine factor (x^2 there corresponds to x^4 here). From PPY25 Eq. (30), using
   Delta E = -delta/n^3 - (3/2) delta^2/n^4 (Z = 1), the defect at x^4 is
   delta_4(n) = delta_4(inf) - delta_0^2 [psi(n) - ln n + 1/(2n)] - Z^4 <r^4>/(15 n^2).
   The -1/n inside PPY's bracket is exactly absorbed by the (3/2) delta_0^2/n^4 term, a non-trivial
   consistency check of this conversion. For the uniform sphere this gives approximately
   delta_4(1) - delta_4(inf) = -0.0162 x^4, and about -(8/525) x^4/n^2 at large n.
4. *A second consumed number at x^4.* <r^4> enters through nabla^2 delta with weight 1/n^2, so a
   threshold (n -> infinity) matching cannot separate it from the log-moment constant. Two numbers are
   consumed at fourth order, not one. Likewise, the hyperfine factor at x^2 consumes <r^2>_M (the
   1/(5n^2) piece, K01's R_M term) and a joint charge-magnetization constant (the corpus's profile
   constants 1.61973, 1.60095, 1.626), on top of the forced log.

---

## 3. Q3: does the corpus credit EFT for this mechanism?

**No.** A grep of papers/ (archive excluded) for NRQED, effective field theory, Wilson coefficient and
matching finds four hits, none of them for the finite-nucleus mechanism:

- Paper 32 (~line 5392): the "matched effective-field-theory chain (HQET -> NRQED -> pNRQED -> SCET;
  Fontes and Szafron)", for W2a (renormalization counterterms, LS-8a), not the finite nucleus. Its move
  ("not GeoVac-specific") is the template for J14/J16.
- Paper 34 (line 5581): "NRQED" only as a label for higher-order QED. Paper 51 (line 1717): a gravity
  table header. Paper 23 (line 1238): chiral EFT for nuclear forces. All unrelated.

Paper 23 sec:finite_nucleus credits Foldy, Friar 1979, Seaton, Pieper-Greiner, Zel'dovich-Popov, Zemach,
Bohr-Weisskopf, Friar-Payne 2005 and Paper 34 ("the physics is textbook; the contribution is the
accounting"). It does not credit (a) the EFT/NRQED mechanism; (b) Karshenboim for the psi(n) form or the
hyperfine log; (c) Friar 1979 / PPY25 for the x^4 log and gamma_E (the nuclear_profiles record, line 57,
says it is known; the paper does not). Friar-Payne 2005 covers the nuclear side (Zemach, Low,
polarizability), not the electron-side x^2 ln x. Seam J14/J16, Paper 57 sec:unity, the walls note and
memory/multi_focal_wall_pattern.md carry no EFT credit. Paper 57's Observation puts multi-focal
composition "outside the framework's native composition apparatus". But EFT matching is the standard
apparatus for composing two scales, and the price it charges, Wilson coefficients, is what the doors
measured.

---

## 4. Q4: verdict, what may be new, and re-scoping

**Known and to be credited.** The mechanism (E1-E4). The leading terms (items 2, 4, 12). The l = 1
defect (PPY25 Eq. 77). The x^4 log structure (PPY25 Eq. 30). R^(2 gamma) (Shabaev 1993). Z_cr. The
hyperfine log coefficient (K97 Sec. 2.5). The psi(n) - ln n - 1/n form (K97 Eq. 8; K01 Eq. 26), which
removes the "post hoc" status of J16 (iv): it has been published for arbitrary n since 1997/2001. The
magnetic-radius 1/n^2 term (K01 Eq. 26).

**Plausibly new (low novelty, accounting-level; none is a new mechanism):** (i) the hyperfine log stated
as -2 times the defect's leading coefficient (both are in the literature, proportional to
(Z alpha m R_E)^2, but the relation is not stated; it is one line in the contact picture); (ii) the
general-l product c_l prod(1 - s^2/n^2), not found beyond l = 1, though it is first-order perturbation
theory; (iii) the x^3 ln x term, and the candidate -64/63 if confirmed; (iv) the explicit finite-size
continuation with J13's factor; (v) the transcendental bookkeeping, the atlas's own framing.

**Overstatements found while checking (for the PI; I have not edited anything):**

1. **"Universal ln(2ZR) + 2 gamma_E" (Paper 23 "The price" and eq:defect_s text; seam J14; Paper 57 J14
   item; NP3 in nuclear_profiles.md).** Only the log coefficient, -delta_0^2, is universal. The constant
   beside it, gamma_E and ln 2 included, depends on the profile's log moment. PPY25 Eq. (30) has a single
   gamma for the exponential model, and EGS01 Eq. (188) averages <ln(2 m r Z alpha/n)> over the density.
   The sphere and the shell give 2 gamma because their log moments are rational. The Gaussian test in
   Paper 23 checked only the log coefficient (0.2%), so it does not support the gamma claim. NP3's
   "coefficient of gamma_E is -2 delta_0^2 for every profile" is, on this evidence, false in general.
2. **eq:defect_s "the same at every n" (Paper 23; seam J14 "the s-wave defect is the same for every n").**
   This is established through x^3 (n = 1 and 5) and at threshold (the 1.3e-18 check is the E = 0
   matching). At x^4, PPY25 Eq. (30) implies the n-dependence of Q2 point 3. That is a *falsifiable
   prediction* for the corpus's own matching code: at x = 1e-3 the n = 1 defect should differ from
   threshold by about -1.6e-14 for the sphere. Recommend running it before the claim is next touched. (My
   conversion from PPY's energy to a defect is my own algebra; it passed the 1/n consistency check.)
3. **"One consumed moment per order" / "then one joint moment per order" (Paper 23, Paper 57 HF-4 and
   J14 items, seam J16, the walls note).** At fourth order in the energy, <r^4> and a log-moment constant
   are both consumed. At second order in the hyperfine factor, <r^2>_M and a joint constant are both
   consumed. The accurate statement is EFT's: one new constant per independent operator at each order,
   with log coefficients fixed by lower orders.
4. **Paper 34 sec:proj_charge_density** says the Foldy/Friar correction "is rational in the charge
   form-factor moments at every order of the small-r expansion". That is contradicted by Paper 23 itself
   (logs and gamma_E from fourth order) and by PPY25 Eq. (30). This is an internal inconsistency.

**Recommended re-scoping (PI decision; nothing edited):**
- Wall statement (Paper 57 obs:multi_focal_wall, memory file, walls note): for J14 and HF-4, state the
  price as the EFT matching statement and cite Caswell-Lepage and Hill-Paz. The GeoVac-specific content
  is then only that the label skeleton has no native matching step (Paper 32 W2a is the template).
  J14/J16 support the wall's existence, not that it is GeoVac-specific.
- Seam J14/J16 prior-art columns: add EFT, PPY25 (Eq. 30, 77), K97/K01. Replace "post hoc" in J16 (iv)
  with "agrees with Karshenboim 1997/2001". Restrict "universal" to the log coefficient.
- Paper 23: the same credits; answer its open Volotka question; fix or test overstatements 1-3 (item 2
  is a cheap run). Paper 34: fix "rational at every order".

**Three most load-bearing citations:** PPY25 Eq. (30) and (77); K01 Eq. (26) with K97 Sec. 2.5 and
Eq. (8); HP11 (matching paragraphs, Eq. 2-4).
