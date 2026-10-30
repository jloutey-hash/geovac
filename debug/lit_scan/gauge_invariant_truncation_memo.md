# Gauge invariance on a spectral truncation: literature scan

Date 2026-10-29. Read-only scan. It was prompted by YC (`docs/predictions/2026-10-29_ccm_yukawa_coefficient_window.md`, GATE GO) and by `debug/graded_ac/graded_product_memo.md`. Both show that the compressed fluctuation A = a[D,b] + h.c. on the true CH window is not gauge invariant: with a = 1, the pure-gauge companion makes up most of k_A. Downloads: `debug/lit_scan/dl_gauge_truncation/`.

## Verdict

**KNOWN OBSTRUCTION** applies to the question as posed: local gauge invariance surviving the compression f ↦ PfP onto a fixed Dirac spectral window. The obstruction is elementary operator algebra, and it has been computed explicitly for the circle. **KNOWN FIXES exist only by changing the object.** Either truncate *after* fluctuating, or replace the operator system by its C*-envelope, which is a full matrix algebra with gauge group U(N). The fuzzy sphere is the one case where these two objects coincide. **OPEN:** no construction was found of a gauge theory on a *fixed* spectral window that has a nontrivial *local* gauge group, nor of a spectral-truncation analogue of link variables.

## 1. The obstruction

(a) **Compression and unitaries.** Let P be a projection and u a unitary. Then PuP is unitary on PH exactly when (1−P)uP = 0 and (1−P)u*P = 0, that is, when [u,P] = 0. This is the compression case of Choi's multiplicative-domain theorem. Choi (Illinois J. Math., "A Schwarz inequality for positive linear maps on C*-algebras") shows that a unital 2-positive map is multiplicative on a only if Φ(a*a) = Φ(a)*Φ(a). So a gauge function u acts on the window only if it commutes with the spectral projection of D. On the CH window the commutant of O^CH has dimension 1 (Paper 32 prop:odd(c), re-measured in the graded-product memo). The only gauge functions that act are therefore the constants.

(b) **The gauge group of an operator system, computed.** Dong (SIGMA 18 (2022) 060, arXiv:2111.13076) defines G(E) = {U ∈ U(C*(E)) : UEU* = E}. For the truncated circle (Toeplitz system Toep_n, whose C*-algebra is M_n), Dong proves Lemma 4.1: every U in G(Toep_n) is diagonal or anti-diagonal. Proposition 4.2 and Corollary 4.3 then give G(Toep_n) = U(1) × (U(1) ⋊ Z₂), the same for every n. What survives is the rotations and reflection of the circle plus a phase. The local gauge group C^∞(S¹, U(1)) is gone at every finite n. This fits Connes–van Suijlekom (arXiv:2004.14115, §3), who show that the truncation keeps the *isometries*, the unitaries with UDU* = D and UAU* = A, *because* such U commute with P. Inner (gauge) unitaries do not commute with D, and that argument does not cover them. Connes–van Suijlekom do not treat gauge fields (full-text grep).

(c) **What the YC miss is.** This is derived here, not taken from the literature. Check it before citing. Let χ be self-adjoint and U_ε = e^{iεχ}. Then U_ε D U_ε* = D + εB₁ + ε²B₂ + …, with B₁ = i[χ,D] and B₂ = −½[χ,[χ,D]]. Since Tr f(U D U*) = Tr f(D) for every f, the ε² coefficient of the *linear* pure-gauge perturbation is

  [ε²] Tr f(D + εB₁) = ½ Tr( f′(D) [χ,[χ,D]] ).

In the continuum, [χ,[χ,D]] = 0, because a function commutes with Clifford multiplication by its own differential. That is order zero plus order one, and it is the reason k^∞_pg = 0. On the window, D is diagonal and row i of the double commutator restricted to W is 2Σ_{j∈W}|χ_ij|²(λ_i−λ_j). In the continuum the full row sum vanishes, so

  k_pg(window) = − Σ_{i∈W, j∉W} |χ_ij|² f′(λ_i)(λ_i − λ_j).

This is an exact edge-crossing sum: a Toeplitz boundary term weighted by f′ at the window edge. With f = e^{−tx²}, f′(λ_edge) ∝ tλ e^{−τ}, where τ = tλ_edge². The term is therefore suppressed at fixed t as n_max grows, consistent with the fixed-t = 0.1 controls (ρ_pg 1.2e-2 → 1.4e-6). It is *not* suppressed at fixed τ, which is YC's scored regime.

Checked numerically (scratchpad, float64, finite-difference ε²). On the truncated antiperiodic circle with χ = compressed cos x and t = 0.1:

| R | FD | edge formula |
|---|---|---|
| 3 | −0.1183153 | −0.1183153 |
| 5 | −0.0430520 | −0.0430520 |
| 8 | −1.66193e-3 | −1.66193e-3 |
| 12 | −1.1404e-6 | −1.1402e-6 |

The double-commutator form ties the YC gauge breaking directly to the order-zero and order-one failure on the true window (Paper 32). CCvS: "the first order condition is equivalent to the linearity of the map from 1-forms to fluctuations" (arXiv:1304.7583).

## 2. Known fixes and what each would require on the CH window

**F1. Fluctuate, then truncate (sharp spectral action on D_A).** Chamseddine–Connes (hep-th/9606001) write the action as Trace χ(D/Λ), where D is the *fluctuated* operator (eqs. 1.23, 1.28). The inner automorphisms act as D ↦ uDu*, A ↦ u[D,u*] + uAu* (eq. 1.24). The action depends only on the spectrum, so it is gauge invariant for *any* χ, sharp or smooth. Smoothness is invoked only for the asymptotic expansion: "for any smooth cutoff function χ … one has" eq. (1.26). So sharp versus smooth does not matter for gauge invariance (task item 5). What matters is *which operator's* spectral projection defines the window. On the CH window this would mean replacing P_D by P_{D_A} = 1_{|D_A|≤Λ}: a window that moves with A. Exact invariance then needs the full D_A on L²(S³, S). Gauge invariance comes back, but the fixed finite-dimensional window is lost. Chamseddine–Connes eq. (1.27) proposes truncating to range χ(D/Λ) "while maintaining the commutation rule (1.20) for the algebra". Their programme is therefore a truncation that keeps an *algebra*, not a compressed operator system.

**F2. Promote to the C*-envelope and fluctuate nonlinearly.** Verhoeven (arXiv:2610.19791, Thm 4.10) proves that the C*-envelope of the Dirac spectral truncation S_J of SU(2), which is the CH window, is M_N(C) with N = dim P_J L²(S), and that 2 ≤ prop(S_J) ≤ 8. On M_N one has a genuine unitary group, and the CCvS perturbation semigroup applies: D ↦ Σ a_j D b_j with Σ a_j b_j = 1, including the quadratic A^(2) when first order fails (arXiv:1304.7583, Lemma 4 and §II). Pure gauge then costs exactly zero, because (u, u*) gives uDu*. Dong's Pert(E) is the operator-system version: a_i, b_i ∈ C*(E) with a_i E b_i* ⊂ E. **Cost on the window:** the gauge group becomes U(N), a non-local matrix-model symmetry. The fluctuations leave O^CH·[D,O^CH], and the local U(1) or SU(2) structure of S³ is not recovered. The gauge group G(O^CH) in Dong's sense has **not** been computed for the CH window. By analogy with Toep_n it is presumably the isometries times phases. That is unverified, and it is a cheap exact computation if wanted.

**F3. Fuzzy-sphere construction (task item 2).** Carow-Watamura and Watamura (CMP 212 (2000) 395, hep-th/9801195) take A_N ≅ M_{N+1}. Fields are A_N-bimodules, and "the local U(1) gauge transformation on the fuzzy sphere is identified with the left U(N+1) transformation of the field", with A′ = u d u* + u A u* (their eqs. 65–69). Steinacker (NPB 679 (2004) 66, hep-th/0307075) quantizes U(n) Yang–Mills on the fuzzy sphere as a matrix model. Gauge invariance is exact because the truncated function space *is* a matrix algebra: on S², ⊕_{l≤N}(2l+1) = (N+1)², so harmonics up to l = N fill End(V_{N/2}). This does **not** transfer to the CH window. There dim O^CH = 5, 30, 91 against N² = 16, 256, 1600, so the compressed functions are a proper subspace and the algebra they generate is the whole envelope (F2). For odd spheres, Ramgoolam (hep-th/0105006, abstract) states that higher fuzzy spheres are non-associative algebras obtained as projections of associative matrix algebras. So no direct fuzzy-S³ analogue keeps an associative function algebra.

**F4. Link variables (task item 3).** Marcolli–van Suijlekom, "Gauge networks in noncommutative geometry" (arXiv:1301.3480, J. Geom. Phys. 75 (2014) 71), build quivers of finite spectral triples with a lattice Dirac operator. The spectral action on a 4D lattice reduces to the Wilson action plus a lattice Higgs system. Pérez-Sánchez (arXiv:2508.17338) argues that the continuum limit is pure Yang–Mills with no Higgs. This is the NCG analogue of Wilson's links, and it keeps lattice gauge invariance exactly. It needs **sites**: a vertex algebra (⊕_v A_v) whose unitaries are the local gauge group, with gauge data carried by morphisms between vertex triples. The CH window has no sites, since the commutant of O^CH is one-dimensional, so the construction does not apply. **No spectral-truncation analogue of link variables was found.** This is OPEN.

## 3. Table

| source | what it constructs | keeps gauge invariance? | applies to a Dirac-spectral window? | verified? |
|---|---|---|---|---|
| Connes–van Suijlekom, arXiv:2004.14115 (CMP 383, 2021) | operator-system spectral triples (PAP, PH, PDP); Toeplitz and Fejér–Riesz; propagation number; state-space metric | no gauge treatment; keeps isometries only | yes, defines it | full text (abstract + §3 + grep) |
| Verhoeven, arXiv:2610.19791 (2026) | Peter–Weyl and Dirac truncations of SU(2); C*-env = M_N; 2 ≤ prop ≤ 8 | no gauge treatment | yes, it is the CH window | full text, Thm 4.10 |
| Dong, arXiv:2111.13076 (SIGMA 18, 060) | gauge group G(E) and perturbation semigroup of an operator system; G(Toep_n) = U(1)×(U(1)⋊Z₂) | only global/isometry part survives | definition yes; CH window not computed | full text |
| Choi, Illinois J. Math. (1974) | multiplicative domain of a 2-positive map | gives the obstruction ([u,P] = 0) | yes | title and author via Project Euclid; volume 18, pp. 565–574 SECONDARY |
| Chamseddine–Connes–van Suijlekom, arXiv:1304.7583 (JGP 73, 2013) | semigroup Pert(A); quadratic A^(2) without first order; D ↦ UDU* | yes, for an algebra | only via the envelope M_N (F2) | full text |
| Chamseddine–Connes, hep-th/9606001 | Tr χ(D_A/Λ); gauge action (1.24); truncation H_Λ = range χ(D/Λ) (1.27) | yes, for any χ (spectral); smooth χ only for the expansion | only if the window moves with A (F1) | full text; journal ref (CMP 186, 731) SECONDARY |
| Carow-Watamura & Watamura, hep-th/9801195 (CMP 212, 395) | U(1) gauge theory on fuzzy S²; gauge = left U(N+1) on bimodule | yes, exactly | no: needs truncated functions = matrix algebra | full text |
| Steinacker, hep-th/0307075 (NPB 679, 66) | U(n) YM on fuzzy sphere as random matrix model | yes | no (same reason) | abstract |
| Grosse–Klimčík–Prešnajder, hep-th/9602115 (CMP 180, 429) | truncated S⁴, scalar field, SO(5) kept | gauge not in abstract | no | abstract |
| Ramgoolam, hep-th/0105006 | higher/odd fuzzy spheres as non-associative projections of matrix algebras | n/a | indicates no associative fuzzy S³ | abstract; author and journal SECONDARY |
| D'Andrea–Lizzi–Martinetti, arXiv:1305.2605 (JGP 2014) | cutoff by projection P on D; Moyal plane, fuzzy sphere | gauge not in abstract | partly (cutoff yes, gauge no) | abstract |
| Marcolli–van Suijlekom, arXiv:1301.3480 (JGP 75, 71) | gauge networks; lattice Dirac; spectral action → Wilson + lattice Higgs | yes (lattice, links) | no: needs sites | abstract; edge placement of gauge data SECONDARY |
| van Nuland–van Suijlekom, arXiv:2104.09899 | spectral action expanded in Chern–Simons/YM cocycles, using its gauge invariance | yes, for an algebra | not addressed | abstract |

## 4. Consequences for the YC instrument

1. YC's reading ("a wall on the instrument, not a statement that the outer factor moves the 4") is consistent with the literature. Nothing found contradicts it, and §1(c) gives the mechanism in closed form: an edge-crossing boundary sum, equivalent to the order-zero and order-one defect.
2. A gauge-invariant redefinition that keeps the fixed window is not in the literature. Subtracting k_pg is not gauge invariance, since it removes one orbit direction, not the orbit.
3. The cheapest exact route is F2 restricted to unitaries: score the exponentiated companion Tr f(e^{iχ_w} D_w e^{−iχ_w}), which costs zero by construction, against the linear one. The difference *is* ½Tr(f′[χ_w,[χ_w,D_w]]), which isolates the defect. That is a diagnostic, not a fix.

Scope: commutative-factor compression only. Nothing here bears on WH1's open almost-commutative convergence or on Paper 32's H1 theorem.
