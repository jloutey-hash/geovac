# Sprint H1: Higgs from Inner Fluctuation on the GeoVac AC Extension

**Sprint:** Step 3 of the May 3-4 2026 PI three-step spectral-triple commitment.
**Status:** COMPLETE (positive-thin verdict).
**Files:** `geovac/almost_commutative.py` (~700 lines), `tests/test_almost_commutative.py` (38 tests, all passing), `debug/h1_sprint_analysis.py` (driver), `debug/data/h1_falsifier.json` (results), this memo.
**Date:** 2026-05-06.

---

## §0. Executive verdict

**POSITIVE-THIN** in the language of the scoping memo §0.

The minimal electroweak almost-commutative extension
$\mathcal{T}_{\mathrm{AC}} = \mathcal{T}_{\mathrm{GV}} \otimes \mathcal{T}_F$
with $\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H}$ is **constructively well-defined** at finite $n_{\max}$ on the truthful-CH GeoVac spectral triple. Connes' inner-fluctuation formula
$D \mapsto D + \omega + \epsilon' J \omega J^{-1}$
is well-defined; the gauge sector reproduces the Paper 25 / Paper 30 structure (as expected), and the Higgs sector $\Phi$ is **non-trivially populated** for any imposed non-zero Yukawa matrix $Y$ in $D_F$.

**However**, GeoVac structure does **not autonomously select** a non-trivial $Y$. The two candidate sources flagged in the scoping memo §3.3 are:
- **Candidate A** (offdiag CH chirality bridging): **RULED OUT** by Track 2's clean negative ($J D = +DJ$ residual ~2.0 on offdiag CH, see `debug/real_structure_finite_nmax_memo.md` §3D).
- **Candidate B** (Sturmian / DUCC bridging): not pursued in this sprint; would require committing to a non-trivial inter-shell coupling whose physical interpretation as a Yukawa is itself speculative.

So the verdict is what the scoping memo §3.4 named the **load-bearing case** and feared:
> If $D_F$ off-diagonal must be imposed by hand, then GeoVac is on the Marcolli–vS-without-Higgs side of the 2014/2024 distinction.

The construction works, but $Y$ is a free input. This is **exactly the structural reading** of Paper 32 §VIII.B G2 that the framing edits in this sprint sharpened to: gap G2 is **not** "inner fluctuations cannot define a Higgs" — they can — but rather "no GeoVac mechanism selects a non-trivial Yukawa structure."

---

## §1. Architecture, after Track 2 closure

The PI dispatch fixed the architecture as follows (post-Track-2):

1. **Use the truthful Camporesi–Higuchi Dirac $D_{\mathrm{GV}}$.** The offdiag CH Dirac is ruled out as the AC base by Track 2's verdict that $JD = +DJ$ fails on offdiag CH (residual ~2.0).
2. **Internal couplings live on the $M_n(\mathbb{C}) = \mathcal{A}_F$ factor itself.** $D = D_{\mathrm{GV}} \otimes \mathbb{1}_F + \gamma_{\mathrm{GV}} \otimes D_F$, where $\gamma_{\mathrm{GV}}$ is the chirality grading on $\mathcal{H}_{\mathrm{GV}}$ (defined as the SIGN of $D_{\mathrm{GV}}$ in our truthful convention) and $D_F$ is a Hermitian operator on $\mathcal{H}_F$.
3. **$\mathcal{H}_F$ is matter–antimatter doubled** following Connes–Marcolli 2008 Ch. 13: $\mathcal{H}_F = \mathcal{H}_F^{\mathrm{mat}} \oplus \mathcal{H}_F^{\mathrm{anti}} = \mathbb{C}^4 \oplus \mathbb{C}^4 = \mathbb{C}^8$. The algebra $\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H}$ acts on the matter sector only; the antimatter sector picks up the opposite-algebra action via $J$.
4. **$J = J_{\mathrm{GV}} \otimes J_F$**, with $J_F$ swapping matter $\leftrightarrow$ antimatter via complex conjugation:

  $J_F (\psi_{\mathrm{mat}}, \psi_{\mathrm{anti}}) = (\overline{\psi_{\mathrm{anti}}}, \overline{\psi_{\mathrm{mat}}}).$

  In matrix form $J_F = U_F K$ with $U_F = \sigma_x \otimes \mathbb{1}_4$.

This architecture is **canonical Connes–Chamseddine** (the matter–antimatter doubling is exactly the standard SM construction; see Connes–Marcolli 2008 §13.4); the only difference is that the GV factor here is the GeoVac S^3 truncation rather than the continuum Riemannian $S^3$.

### 1.1. Why an earlier (non-doubled) attempt failed

A naive first pass attempted $\mathcal{H}_F = \mathbb{C}^4$ undoubled, with $J_F = (i\sigma_2 K) \oplus (i\sigma_2 K)$ on $(L, R)$ blocks. Test failure: this $J_F$ does NOT satisfy $J_F D_F = +D_F J_F$ for general Yukawa $Y = \mathrm{diag}(y_\nu, y_e)$ — it requires $y_\nu = y_e$ (degenerate Yukawa) for the sign relation to hold. This is structurally because $i\sigma_2$ swaps the L-doublet entries with a sign, and a diagonal $Y$ is preserved by this swap only if its two entries agree.

The matter–antimatter doubling fixes this: on the doubled $\mathcal{H}_F$, $J_F = \sigma_x \otimes \mathbb{1}_4$ swaps SECTORS not entries within a sector, and so commutes with any $D_F$ of the form $\mathrm{block\_diag}(M, \overline{M})$. This is the standard Connes–Chamseddine convention, and it is the right convention for a *finite* electroweak triple over $\mathbb{C}$.

### 1.2. KO-dimension bookkeeping

| Factor | KO-dim | $J^2$ sign | $JD = \pm DJ$ |
|--------|:------:|:----------:|:-------------:|
| GV ($S^3$) | 3 | $-$ | $+$ |
| F ($\mathbb{C} \oplus \mathbb{H}$, doubled) | 6 | $+$ | $+$ |
| **Combined** | **3 + 6 = 9 ≡ 1 (mod 8)** | **$-$** | **$+$** |

So the combined real spectral triple sits at KO-dim 1, with $J^2 = -\mathbb{1}$ and $JD = +DJ$. Both signs verified to machine precision in `tests/test_almost_commutative.py::TestCombinedRealStructure`.

This shifts the KO-dim from the scoping memo §2.4 estimate (KO-dim 5 / 9). The discrepancy traces to the matter–antimatter doubling vs single-sector convention; the doubling is what makes $J_F D_F = +D_F J_F$ hold for arbitrary $Y$ rather than only for degenerate Yukawa, which is the structurally correct convention for the SM-flavored triple.

---

## §2. Inner fluctuation: gauge and Higgs decomposition

### 2.1. The omega calculation

For $\omega = \sum_i a_i [D, b_i]$ with $a_i, b_i \in \mathcal{A}_{\mathrm{GV}} \otimes \mathcal{A}_F$:

$[D, b] = [D_{\mathrm{GV}}, b^{\mathrm{GV}}] \otimes b^F + b^{\mathrm{GV}} \gamma_{\mathrm{GV}} \otimes [D_F, b^F]$

(using that $\gamma_{\mathrm{GV}}$ commutes with the scalar multiplier $b^{\mathrm{GV}}$, which is itself bosonic).

So $\omega$ decomposes as

$$\omega = \underbrace{\sum_i a_i^{\mathrm{GV}} [D_{\mathrm{GV}}, b_i^{\mathrm{GV}}] \otimes a_i^F b_i^F}_{\text{gauge piece}} + \underbrace{\sum_i a_i^{\mathrm{GV}} b_i^{\mathrm{GV}} \gamma_{\mathrm{GV}} \otimes a_i^F [D_F, b_i^F]}_{\text{Higgs piece}}.$$

The gauge piece reproduces Papers 25 and 30 (U(1) from the $\mathbb{C}$-summand, SU(2) from the $\mathbb{H}$-summand). The Higgs piece is non-zero exactly when $[D_F, b^F]$ has off-diagonal $\mathbb{C} \leftrightarrow \mathbb{H}$ matrix elements, which happens iff $D_F$ has a non-zero off-diagonal block — the Yukawa.

### 2.2. Matter / antimatter localization

A key structural property verified at machine precision:

> **Inner fluctuations preserve the matter / antimatter block decomposition.**

For any $a, b \in \mathcal{A}_{\mathrm{GV}} \otimes \mathcal{A}_F$, $\omega = a [D, b]$ has zero off-block between the matter and antimatter sectors of $\mathcal{H}_F$ (verified in `test_matter_antimatter_off_block_zero`). This is the well-known reason order-zero $[a, J b J^{-1}] = 0$ holds at the finite-sector level: $a$ is supported on matter only, $J b J^{-1}$ is supported on antimatter only (as it represents the opposite algebra), and they commute trivially.

### 2.3. Explicit decomposition (n_max = 2, sample)

| Scenario | $\|\omega\|$ | $\|\Phi\|$ Higgs | $\|\omega_{\mathrm{gauge}}\|$ | matter↔antimatter |
|----------|:-----------:|:----------------:|:-----------------------------:|:-----------------:|
| pure-Higgs ($a, b$ vary in $\mathcal{A}_F$, $a^{GV} = b^{GV} = I$) | 1.4422 | 1.4422 | 0.0000 | 0.0000 |
| pure-gauge ($b^F = I$, $b^{GV}$ varies) | 1.2732 | 0.0000 | 1.2732 | 0.0000 |
| combined | 1.2522 | 0.3143 | 1.2121 | 0.0000 |

(See `debug/data/h1_falsifier.json` for full data.)

The decomposition is clean: gauge and Higgs sectors are independent and complementary, and the matter↔antimatter off-block is identically zero.

---

## §3. Falsifier verdict

The scoping memo §5 statement (strong falsifier):

> Show that for *every* Hermitian $D_F$ on $\mathbb{C} \oplus \mathbb{H}$ derivable from GeoVac structure, the resulting inner fluctuation produces *only* gauge 1-forms.

| $D_F$ candidate | n_max | $\|\Phi\|_{\max}$ | $\|\omega_{\mathrm{gauge}}\|_{\max}$ | Falsifier |
|-----------------|:-----:|:------------------:|:-----------------------------------:|:---------:|
| $D_F = 0$ (zero Yukawa) | 1 | 0.000 | 0.000 | **HOLDS** |
| $D_F = 0$ | 2 | 0.000 | 0.389 | **HOLDS** |
| $D_F = 0$ | 3 | 0.000 | 0.580 | **HOLDS** |
| $y_e = 0.1$, $y_\nu = 0$ | 1 | 0.046 | 0.000 | FAILS |
| $y_e = 0.1$, $y_\nu = 0$ | 2 | 0.051 | 0.389 | FAILS |
| $y_e = 0.1$, $y_\nu = 0$ | 3 | 0.080 | 0.580 | FAILS |
| $y_e = 0.3$, $y_\nu = 0.2$ | 1 | 0.181 | 0.000 | FAILS |
| $y_e = 0.3$, $y_\nu = 0.2$ | 2 | 0.171 | 0.389 | FAILS |
| $y_e = 0.3$, $y_\nu = 0.2$ | 3 | 0.271 | 0.580 | FAILS |

**Reading.** The falsifier holds *iff* $D_F = 0$ (zero Yukawa). For any non-zero imposed Yukawa, the Higgs sector is non-trivially populated. The falsifier is **definitively NOT a structural property of GeoVac** — the construction admits a Higgs sector, it just requires the user to specify $D_F$.

This is **positive-thin** in the scoping memo §0 sense: positive (Higgs construction works) but thin (no GeoVac mechanism selects $D_F$).

### 3.1. Cross-validation: spectral action coherence

The Yukawa scan $\mathrm{Tr}(D^2)$ vs $y_e$ at $n_{\max} = 2$ gives the predicted polynomial:

$\mathrm{Tr}(D^2) = \mathrm{Tr}(D_0^2) + 4 \cdot \dim_{\mathrm{GV}} \cdot (|y_\nu|^2 + |y_e|^2)$

(verified in `test_Tr_D_squared_polynomial_in_yukawa`).

The cross-term $\mathrm{Tr}(D_{\mathrm{GV}} \gamma_{\mathrm{GV}} \otimes D_F) = \mathrm{Tr}(|D_{\mathrm{GV}}|) \cdot \mathrm{Tr}(D_F) = 0$ vanishes because $\mathrm{Tr}(D_F) = 0$ (Yukawa $D_F$ is off-diagonal in L↔R, hence traceless). This structural cancellation is the standard CC mechanism.

---

## §4. Why the strong falsifier doesn't fire as a clean negative

The strong falsifier in scoping memo §5 was designed to FIRE if either

(a) $\gamma$ commutes with every fiber-bridging $D_F$, or
(b) order-one forces all candidate $D_F$ to be block-diagonal in $\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H}$.

Neither (a) nor (b) holds:

- **(a) fails** because $\gamma_{\mathrm{GV}}$ only acts on the GV factor; it has no action on $\mathcal{A}_F$, so it neither commutes nor anti-commutes with $D_F$ in any non-trivial way. The factor of $\gamma_{\mathrm{GV}}$ in the Higgs piece of $\omega$ is just a chirality-sign flag on the GV index, not a $D_F$ constraint.
- **(b) fails** because order-one is automatically satisfied via matter–antimatter doubling: $a$ acts on matter, $J b J^{-1}$ on antimatter, and these are disjoint blocks. So order-one is trivial on this AC extension.

The construction **defines** the Higgs given any choice of $Y$. The question of *which* $Y$ is "right" is not addressable by inner fluctuations alone; in the SM that data is supplied by experiment.

---

## §5. What this means for Paper 32 §VIII.B G2

**G2 sharpens, but does NOT close.**

The G2 entry in Paper 32 §VIII.B (after the Phase-0 framing edit applied in this sprint) names the gap as:

> The natural extension is to $\mathcal{A}_{\mathrm{GV}} \otimes (\mathbb{C} \oplus \mathbb{H})$, in which inner fluctuations split into a gauge-1-form sector (Papers 25 / 30 recovered) and a candidate Higgs scalar sector. The off-diagonal $D_F$ ingredient that turns Marcolli–vS-without-Higgs into a Higgs construction is not yet selected from GeoVac structure.

This sprint confirms that diagnosis **and proves it constructively**:

- The construction goes through; inner fluctuations are well-defined; the Higgs sector is non-trivial for any non-zero Yukawa.
- GeoVac data does **not** select a non-zero Yukawa: zero is as natural a choice as any other from a purely GV-structural standpoint.

So Paper 32 G2 is **sharpened to "construction admits Higgs but does not autonomously generate it"** rather than closed. This is the positive-thin reading.

---

## §6. Cross-checks against scoping memo §6

### 6.1. Paper 18 (exchange constants) — neutral.

The Yukawa entries $y_\nu, y_e$ are free parameters in this construction; they are **calibration exchange constants** in Paper 18's sense (the same tier as $\pi$ in Paper 2). This is exactly what the scoping memo §6.1 anticipated: "if the GeoVac Higgs construction succeeds, the Higgs vev should appear as a Paper 18 calibration constant." It does. No new tier needed; no extension of Paper 32 §VIII case-exhaustion theorem needed (the Yukawa is not a $\pi$-bearing transcendental; it's a free real parameter).

### 6.2. Paper 35 (time as projection) — neutral.

The Higgs sector here is built without any temporal/spectral integration (we're at finite $n_{\max}$ throughout). So Paper 35's load-bearing principle is consistent: no $\pi$ enters because no continuous integration is performed. The Higgs vev would acquire $\pi$-content only if one performed a spectral-action heat-kernel computation in the continuum limit, which is exactly the temporal/spectral-window projection of Paper 35.

### 6.3. Paper 34 (projection taxonomy) — adds a 16th projection candidate.

Sprint H1 introduces a candidate **16th projection** on top of Paper 34's existing fifteen: **almost-commutative fiber promotion**. Variable axis: adds a fiber index $\sigma \in \{1, 2, ..., 8\}$ for the doubled electroweak fiber. Dimension axis: adds 0 spatial dimensions (the fiber is 0-dimensional). Transcendental axis: $\pi$-free at the bare-construction level (matter–antimatter doubling and Yukawa entries are free reals).

The 16th projection is documented here but NOT auto-applied to Paper 34 (per CLAUDE.md §13.5 paper-edit policy: this is a substantial structural addition; PI may apply if desired).

### 6.4. WH4 (four-way $S^3$ coincidence) — adds a fifth role.

The Higgs construction adds $S^3$ as the base manifold of the trivial $\mathbb{C}^2$ Higgs bundle. Five-way coincidence: Fock projection / Hopf base / CH spin carrier / SU(2) gauge manifold / Higgs base. Whether this is forced or accidental is unchanged by H1; the construction is *consistent with* but doesn't *prove* the four-way (now five-way) coincidence.

### 6.5. R3.5 chirality grading — does NOT close G3.

The R3.5 chirality grading $\chi = \pm 1$ on $\mathcal{H}_{\mathrm{GV}}$ in our construction is the SIGN of the truthful CH Dirac eigenvalue (so $\gamma_{\mathrm{GV}} = \mathrm{sign}(D_{\mathrm{GV}})$, commuting with $D_{\mathrm{GV}}$). This $\chi$ is independent of the SM weak-isospin chirality (the L/R distinction within the matter sector of $\mathcal{H}_F$). The two are NOT identified; weak-isospin chirality is a separate $\mathbb{Z}_2$-grading on $\mathcal{H}_F$, added externally.

So G3 remains open after H1: the natural GV chirality and the SM weak-isospin chirality are distinct gradings, and there is no GeoVac mechanism currently identifying them. Connes' SM uses the latter; GeoVac supplies the former.

---

## §7. Honest scope and limitations

**What H1 establishes:**

- The minimal electroweak almost-commutative extension exists at finite $n_{\max}$ on the truthful CH GeoVac triple.
- $J = J_{\mathrm{GV}} \otimes J_F$ satisfies the Connes axioms ($J^2 = -\mathbb{1}$, $JD = +DJ$, both verified at machine precision at $n_{\max} \le 3$).
- Inner fluctuations split into gauge (Papers 25/30) and Higgs sectors as expected from CC.
- For any non-zero imposed Yukawa $Y$, the Higgs sector is non-trivial.
- For $Y = 0$, no Higgs.

**What H1 does NOT establish:**

- A natural GeoVac mechanism selecting $Y \neq 0$.
- The Higgs vev as a derived quantity (it is a free input given $Y$ and $\Lambda$).
- A connection between the Higgs scalar potential coefficients and the GeoVac transcendental taxonomy (this would require a continuum spectral-action computation, deferred to Track 1 R2.5 L4 closure for the rigorous interpretation).
- Identification of the R3.5 chirality with weak-isospin chirality (G3 stays open).
- The full SM color factor $M_3(\mathbb{C})$, deferred per Paper 32 §VIII.B G4 (cross-manifold gap).

**What H1 explicitly closes:**

- The "can inner fluctuations be defined on the GeoVac AC extension at all?" question (yes).
- The "does the matter–antimatter doubling work cleanly at finite $n_{\max}$?" question (yes).
- The "is the candidate-A 'offdiag CH bridging' route open?" question (no — Track 2 closed it).

---

## §8. Files produced

- `geovac/almost_commutative.py` (~700 lines) — implementation with `ElectroweakFiniteTriple`, `AlmostCommutativeTriple`, `ElectroweakFiniteTriple.matter_dirac()`, `dirac_F()`, `real_structure_F()`, `algebra_action()`. Replaces Track 3 stub.
- `tests/test_almost_commutative.py` (38 tests, all passing) — covers the finite triple, algebra action, combined Dirac, real structure, inner fluctuations, falsifier, and spectral action sanity checks.
- `debug/h1_sprint_analysis.py` (driver) — runs the falsifier sweep at $n_{\max} \in \{1, 2, 3\}$ and the explicit decomposition.
- `debug/data/h1_falsifier.json` — raw data for memo §3 table.
- `debug/h1_ac_extension_memo.md` (this file).

---

## §9. Recommendation: Paper 32 §VIII.C addendum

Per the verdict (positive-thin), the right Paper 32 update is a **§VIII.C addendum** stating:

1. The construction works at finite $n_{\max}$.
2. The Higgs sector is well-defined.
3. GeoVac does not autonomously select a non-trivial Yukawa.
4. G2 is sharpened, not closed.

NOT a Paper 37 (the result is too thin to justify a self-contained paper; CLAUDE.md §1.5 rhetoric rule applies).

The §VIII.C addendum is applied alongside this memo; see `papers/group1_operator_algebras/paper_32_spectral_triple.tex` after this commit.

---

**End of memo.** Word count: approximately 3,200 words.
