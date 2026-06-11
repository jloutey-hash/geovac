# Sprint M-H1: Re-reading the H1 verdict through the dual modular propinquity lens

**Sprint:** Track M-H1 of the May 2026 dual-modular-propinquity sprint cluster (5 parallel tracks).
**Status:** COMPLETE.
**Verdict line:** **MIXED — SELF-ADJOINTNESS-REPHRASING with one structural sharpening of the falsifier (no new constraint on $Y$).**
**Files:** this memo. No production code or paper edits applied.
**Date:** 2026-05-23.

---

## §0. Executive summary

The PI hypothesis was that the dual modular Gromov–Hausdorff propinquity (Latrémolière 2018/2021, arXiv:1811.04534), by virtue of acting on metrized quantum vector bundles (Hilbert C*-bimodules carrying a D-norm), might impose a structural constraint on the Yukawa $Y$ that Sprint H1's algebra-level analysis missed. The Higgs cross-block of the inner-fluctuation $\omega = a\,dD\,b$ is, after all, a bimodule element between the $\mathbb{C}$-summand and the $\mathbb{H}$-summand of $\mathcal{A}_F$; if duality / Morita equivalence is to be respected, perhaps the cross-block is forced to vanish or to take a specific form.

The answer, on careful reading, is **no**: the dual modular propinquity does not produce a new constraint on $Y$ beyond what is already imposed by Connes' first-order condition, the Hermiticity of $D_F$, and the matter–antimatter doubling. The Higgs cross-block satisfies a "self-duality" relation $\Phi^\dagger \in \mathrm{Hom}(\mathcal{H}_\mathbb{H}, \mathcal{H}_\mathbb{C})$ as the right-action analog of the left-action element $\Phi \in \mathrm{Hom}(\mathcal{H}_\mathbb{C}, \mathcal{H}_\mathbb{H})$, but this is **strictly equivalent to $D_F^* = D_F$**, which is already imposed at the AC-extension level. The dual modular framework reads the same data as a duality requirement on the bimodule, which is more elegant vocabulary but does not pin $Y$.

There is **one structural sharpening** worth recording: the H1 falsifier, reformulated in bimodule language, asks whether every Hermitian bimodule endomorphism $\omega$ on $\mathcal{H}_{\mathrm{AC}}$ derivable from GeoVac structure has zero cross-block as an element of $\mathrm{Hom}_{\mathcal{A}_F\text{-bimod}}(\mathcal{H}_\mathbb{C}, \mathcal{H}_\mathbb{H})$. The bimodule reading makes one feature manifest that H1's matrix-level analysis treated as merely "matter-antimatter doubling makes order-zero automatic": the **left-action and right-action subspaces of the cross-block decouple completely**, which is the structural reason the cross-block doesn't pick up GeoVac-side constraints. The dual modular propinquity provides the right vocabulary to state this cleanly, but the conclusion is the same: $Y$ is free.

For G3 (electroweak chirality) and G4a (full Connes SM with $M_3(\mathbb{C})$), the bimodule reading is similarly **vocabulary-only** at the present sprint-scale. The G3 obstruction (independent $\mathbb{Z}_2$ gradings) survives the bimodule reading verbatim; the G4a obstruction (cross-manifold, blocked by Paper 24 §V Coulomb/HO asymmetry) is unaffected.

**Net:** H1's POSITIVE-THIN verdict is unchanged. The dual modular propinquity gives better mathematical vocabulary for stating what H1 already established. No new constraint on $Y$ emerges.

---

## §1. Bimodule reframing of inner fluctuations on the AC extension

### 1.1. Standard AC bimodule structure

The minimal AC extension built in Sprint H1 (Paper 32 §VIII.C, `geovac/almost_commutative.py`) has:

- **Algebra:** $\mathcal{A}_{\mathrm{AC}} = \mathcal{A}_{\mathrm{GV}} \otimes \mathcal{A}_F$ with $\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H}$.
- **Hilbert space:** $\mathcal{H}_{\mathrm{AC}} = \mathcal{H}_{\mathrm{GV}} \otimes \mathcal{H}_F$, where $\mathcal{H}_F = \mathcal{H}_F^{\mathrm{mat}} \oplus \mathcal{H}_F^{\mathrm{anti}} = \mathbb{C}^4 \oplus \mathbb{C}^4 = \mathbb{C}^8$ (matter–antimatter doubled per Connes–Marcolli 2008 Ch. 13).
- **Real structure:** $J = J_{\mathrm{GV}} \otimes J_F$ with $J_F = \sigma_x \otimes \mathbb{1}_4$ (the matter↔antimatter swap composed with complex conjugation), KO-dim 9 ≡ 1 (mod 8) for the combined triple ($J^2 = -\mathbb{1}$, $JD = +DJ$).
- **Dirac:** $D = D_{\mathrm{GV}} \otimes \mathbb{1}_F + \gamma_{\mathrm{GV}} \otimes D_F$ with $D_F = \mathrm{block\_diag}(M, \overline{M})$ and $M$ carrying the off-diagonal $L \leftrightarrow R$ Yukawa block.

The bimodule structure on $\mathcal{H}_{\mathrm{AC}}$ is:

- **Left action:** $\pi: \mathcal{A}_{\mathrm{AC}} \to B(\mathcal{H}_{\mathrm{AC}})$, supported on the matter sector of $\mathcal{H}_F$. Explicitly $\pi(a^{\mathrm{GV}} \otimes (\lambda, q)) = a^{\mathrm{GV}} \otimes \mathrm{diag}(M_{\mathrm{matter}}(\lambda, q), 0_4)$ where $M_{\mathrm{matter}}(\lambda, q) = q \oplus \mathrm{diag}(\lambda, \overline{\lambda})$ on the $\mathbb{C}^2_L \oplus \mathbb{C}^2_R$ decomposition of the matter sector. (See `almost_commutative.py` line 301ff.)

- **Right action:** $\pi^{\mathrm{op}}(a) := J\,\pi(a^*)\,J^{-1}$, supported on the antimatter sector. The matter–antimatter doubling realizes the opposite-algebra action $\mathcal{A}_{\mathrm{AC}}^{\mathrm{op}}$ as the antimatter-sector representation, which is exactly Connes' bimodule construction for finite spectral triples (Connes–Marcolli 2008 §13.4).

- **Bimodule:** $\mathcal{H}_{\mathrm{AC}}$ is a left-$\mathcal{A}_{\mathrm{AC}}$, right-$\mathcal{A}_{\mathrm{AC}}^{\mathrm{op}}$ module, with the two actions commuting because they are supported on disjoint sectors (this is the structural reason order-zero $[a, JbJ^{-1}] = 0$ holds automatically — H1 memo §2.2, Connes' order-zero condition).

### 1.2. Inner fluctuations as bimodule endomorphisms

Connes' inner fluctuation $\omega = \sum_i a_i [D, b_i]$ is a bimodule endomorphism of $\mathcal{H}_{\mathrm{AC}}$:

- The left action of $\mathcal{A}_{\mathrm{AC}}$ on $\omega$ is by multiplication on the left.
- The right action of $\mathcal{A}_{\mathrm{AC}}^{\mathrm{op}}$ (via $J$-conjugation) is by multiplication on the right.
- $\omega$ commutes with both actions modulo the gauge-fluctuation reading.

Crucially, $\omega$ itself is **not** in $\mathcal{A}_{\mathrm{AC}}$; it lives in $\Omega^1_D(\mathcal{A}_{\mathrm{AC}})$, the bimodule of Connes 1-forms generated by commutators $[D, b]$ for $b \in \mathcal{A}_{\mathrm{AC}}$. This is the standard setup of Connes–Chamseddine.

In the AC extension, the Connes 1-form decomposes as

$$\omega = \omega_{\mathrm{gauge}} + \omega_{\mathrm{Higgs}}$$

with the gauge piece living in $\Omega^1_{D_{\mathrm{GV}}}(\mathcal{A}_{\mathrm{GV}}) \otimes \mathcal{A}_F$ and the Higgs piece in $\mathcal{A}_{\mathrm{GV}} \cdot \gamma_{\mathrm{GV}} \otimes \Omega^1_{D_F}(\mathcal{A}_F)$ (H1 memo eq. (2.1)).

### 1.3. The metrized quantum vector bundle structure

Latrémolière's dual modular propinquity (arXiv:1811.04534) operates on **metrized quantum vector bundles**: Hilbert C*-bimodules over a Leibniz quantum compact metric space, equipped with a D-norm that generalizes the operator norm of a connection on a Riemannian manifold.

For the AC extension, the candidate metrized quantum vector bundle is the bimodule $(\mathcal{H}_{\mathrm{AC}}, \mathcal{A}_{\mathrm{AC}}, D, J)$ with the D-norm defined by

$$\|x\|_D := \|x\| + \|[D, x]\|.$$

(This is the operator-norm Lipschitz seminorm of Paper 46 §3, applied to the bimodule structure.)

The dual modular propinquity convergence requires that as one approximates a target bimodule by a sequence of finite-dimensional bimodules, the D-norms converge in a controlled way: specifically, the Latrémolière propinquity is bounded by a tunneling-pair construction (Berezin map + projection) that respects the bimodule structure.

**The structural question:** does this convergence requirement impose any non-trivial constraint on the Yukawa $Y$ in $D_F$?

---

## §2. Higgs cross-block in bimodule language

### 2.1. Cross-block as Hom-element

In the bimodule decomposition of $\mathcal{H}_F = \mathcal{H}_\mathbb{C} \oplus \mathcal{H}_\mathbb{H}$ (where $\mathcal{H}_\mathbb{C} = \mathbb{C}^2_R$ carries the $\mathbb{C}$-summand action via $\mathrm{diag}(\lambda, \overline{\lambda})$ and $\mathcal{H}_\mathbb{H} = \mathbb{C}^2_L$ carries the $\mathbb{H}$-summand action via $q$), the off-diagonal block of $D_F$ is:

$$\Phi := D_F^{\,RL} \in \mathrm{Hom}_\mathbb{C}(\mathcal{H}_\mathbb{H}, \mathcal{H}_\mathbb{C}) \quad \text{(an antimatter–antimatter sector restricted view)}.$$

Equivalently, $\Phi \in \mathrm{Hom}(\mathcal{H}_\mathbb{C}^2_L, \mathcal{H}_\mathbb{C}^2_R) \cong M_{2 \times 2}(\mathbb{C})$, which in the H1 implementation is $\Phi = Y = \mathrm{diag}(y_\nu, y_e)$ (the Yukawa matrix on the L-doublet).

The bimodule structure makes the Higgs an element of

$$\Phi \in \mathrm{Hom}\bigl(\mathcal{H}_\mathbb{H}^{\mathrm{left}} \to \mathcal{H}_\mathbb{C}^{\mathrm{left}}\bigr) \otimes \mathrm{Hom}\bigl(\mathcal{H}_\mathbb{H}^{\mathrm{right}} \to \mathcal{H}_\mathbb{C}^{\mathrm{right}}\bigr)$$

where the "left" superscript denotes left-$\mathcal{A}_F$-action support (matter sector) and "right" denotes right-$\mathcal{A}_F^{\mathrm{op}}$-action support (antimatter sector). The cross-block in the matter sector is $\Phi$; the cross-block in the antimatter sector is $\overline{\Phi}$ via $J_F$-conjugation; the matter↔antimatter cross-block is identically zero (H1 memo §2.2, machine-precision verified in `tests/test_almost_commutative.py::test_matter_antimatter_off_block_zero`).

### 2.2. Higgs as bimodule intertwiner

In bimodule terms, $\Phi$ is an intertwiner between two sub-bimodules of $\mathcal{H}_F$:

- The $\mathbb{C}$-sub-bimodule: $\mathcal{H}_\mathbb{C}^{\mathrm{mat}} \oplus \mathcal{H}_\mathbb{C}^{\mathrm{anti}}$, with left action $\lambda \in \mathbb{C}$ acting as $\mathrm{diag}(\lambda, \overline{\lambda})$ on the matter sector and as $\mathrm{diag}(\overline{\lambda}, \lambda)$ on the antimatter sector (via $J_F$-conjugation).
- The $\mathbb{H}$-sub-bimodule: $\mathcal{H}_\mathbb{H}^{\mathrm{mat}} \oplus \mathcal{H}_\mathbb{H}^{\mathrm{anti}}$, with left action $q \in \mathbb{H}$ acting as the quaternion-representation matrix on the matter sector and as its conjugate on antimatter.

$\Phi$ maps between these sub-bimodules. Connes' first-order condition $[[D, a], JbJ^{-1}] = 0$ implies $\Phi$ is a **bimodule intertwiner** — it commutes with the left action of $\mathcal{A}_F$ up to the right action of $\mathcal{A}_F^{\mathrm{op}}$.

The first-order condition is verified at the matter-antimatter level by the doubling: $a$ acts on matter only, $JbJ^{-1}$ acts on antimatter only, so $[[D, a], JbJ^{-1}] = 0$ is automatic for any $D_F$ supported within sectors (H1 memo §4(b)).

### 2.3. What the bimodule reading makes manifest

The matrix-level H1 analysis treated three facts as separate properties:

1. **Hermiticity of $D_F$:** $D_F^* = D_F$, hence $\Phi^* = \overline{\Phi}^T$ (matrix transpose conjugate).
2. **Matter-antimatter doubling:** $D_F = \mathrm{block\_diag}(M, \overline{M})$, hence the antimatter cross-block is $\overline{\Phi}^*$.
3. **Order-one (first-order) condition:** $[[D, a], JbJ^{-1}] = 0$, automatic by doubling.

The bimodule reading unifies these:

- **Property 1** is the requirement that $\Phi$, viewed as an element of $\mathrm{Hom}(\mathcal{H}_\mathbb{H}, \mathcal{H}_\mathbb{C})$ (left-action labeling), has a well-defined adjoint $\Phi^* \in \mathrm{Hom}(\mathcal{H}_\mathbb{C}, \mathcal{H}_\mathbb{H})$ compatible with the inner products on $\mathcal{H}_\mathbb{C}$ and $\mathcal{H}_\mathbb{H}$. This is the bimodule self-duality condition.
- **Property 2** is the requirement that the right-action picture (antimatter sector) is the complex conjugate of the left-action picture, which is the structural meaning of $J^2 = +I$ at KO-dim 6.
- **Property 3** is automatic for bimodules built by doubling: left action on one sector, right action on the other, commutation trivial.

The bimodule reading is more elegant but **does not add information**. Property 1 is what Connes calls "$D_F$ is self-adjoint as a bimodule operator"; Property 2 is the matter–antimatter doubling; Property 3 is order-zero. None of these constrains $|Y|$, $\arg Y$, or the structure of $Y$ beyond $\Phi^* = \Phi^\dagger$ as matrix.

---

## §3. Dual modular propinquity constraint on the cross-block

### 3.1. The D-norm structure

The dual modular propinquity convergence is governed by the D-norm on the bimodule. For the AC extension, the D-norm on the cross-block $\Phi$ is

$$\|\Phi\|_D = \|\Phi\| + \|[D, \Phi]\|_{\mathrm{commutator}},$$

where $[D, \cdot]$ is the graded commutator with $D = D_{\mathrm{GV}} \otimes \mathbb{1}_F + \gamma_{\mathrm{GV}} \otimes D_F$.

For $\Phi$ viewed as the constant matrix $Y$ on the GV side, $[D_{\mathrm{GV}} \otimes \mathbb{1}_F, \mathbb{1}_{\mathrm{GV}} \otimes \Phi] = 0$ identically (left action on GV side is by scalar multiplication), so the GV-side commutator contributes nothing. The fiber-side commutator $[\gamma_{\mathrm{GV}} \otimes D_F, \mathbb{1}_{\mathrm{GV}} \otimes \Phi] = \gamma_{\mathrm{GV}} \otimes [D_F, \Phi]$ encodes the Yukawa structure.

The dual modular propinquity convergence requires $\|\Phi\|_D < \infty$ on every finite-cutoff approximant and $\|\Phi\|_D \to \|\Phi_{\mathrm{limit}}\|_D$ as $n_{\max} \to \infty$.

### 3.2. Does this constrain $Y$?

**No.** The D-norm of $\Phi$ for a finite Yukawa $Y$ is bounded by $\|Y\| + \|[D_F, Y]\| \le \|Y\| + 2\|D_F\| \cdot \|Y\|$, which is finite for any finite $Y$ regardless of structure. The dual modular propinquity does not impose an upper or lower bound on $Y$.

What it **does** impose is a **duality requirement on the bimodule structure**: $\langle \cdot, \cdot \rangle_R = \overline{\langle \cdot, \cdot\rangle_L}$ for the right- and left-inner products on the bimodule. In the H1 construction, this is automatic by the matter–antimatter doubling: the right inner product is the complex conjugate of the left inner product via $J_F$.

Concretely, the duality requirement is

$$\langle \Phi x, y \rangle_{\mathcal{H}_\mathbb{C}} = \langle x, \Phi^* y \rangle_{\mathcal{H}_\mathbb{H}}$$

which is just self-adjointness of $D_F$ restricted to the cross-block, already imposed by Connes' axiom system. **It is exactly Hermiticity rephrased in bimodule language.**

### 3.3. Why one might have hoped for more

A naive reading might hope that the dual modular propinquity (which respects Morita equivalence) imposes that the cross-block be **Morita-trivial** in some sense — e.g., that $\Phi$ be expressible as $\Phi = u^* v$ for some bimodule generators $u, v$ derivable from the GeoVac side. If this were the case, the GeoVac data would in principle constrain $Y$ via the available $u, v$.

But this hope **does not survive the matter–antimatter doubling**. The cross-block $\Phi$ is supported in the matter sector only (the antimatter sector picks up $\overline{\Phi}$ via $J$). On the matter sector, the bimodule structure is just $\mathcal{H}_\mathbb{C} \oplus \mathcal{H}_\mathbb{H}$ with left action of $\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H}$. The cross-block $\Phi \in \mathrm{Hom}(\mathcal{H}_\mathbb{H}, \mathcal{H}_\mathbb{C})$ has no GeoVac-side label — it lives entirely in the fiber. The Morita equivalence question is between the $\mathbb{C}$-block and the $\mathbb{H}$-block in the fiber, both of which are finite-dimensional with no GeoVac structure.

So **there is no GeoVac-side bimodule data that could constrain $\Phi$**. The Yukawa lives in a sector of the bimodule that is structurally decoupled from the GeoVac index.

### 3.4. Connection to the inner-factor Mellin engine

This is consistent with the structural reading from Sprint 2026-05-07's inner-factor Mellin engine (CLAUDE.md §2): the inner factor $\mathcal{A}_F$ contributes via its own Yukawa Dirichlet ring $\mathbb{Q}[y_1^{-2s}, \ldots, y_n^{-2s}]$, categorically disjoint from the outer-factor (GeoVac) exchange constants. The bimodule reading confirms that this disjointness is at the bimodule-structural level, not just at the Mellin-engine level: the cross-block $\Phi$ has no GeoVac index, no Hopf-base measure $M_1$ content, no Seeley–DeWitt $M_2$ content, no vertex-parity $M_3$ content. It is pure inner-factor input data.

---

## §4. Comparison to H1's falsifier

### 4.1. H1's original falsifier (matrix form)

H1 §3 / scoping memo §5:

> Show that for every Hermitian $D_F$ on $\mathbb{C} \oplus \mathbb{H}$ derivable from GeoVac structure, the resulting inner fluctuation produces only gauge 1-forms.

This holds iff $D_F = 0$ (zero Yukawa) — H1 table at $n_{\max} \in \{1, 2, 3\}$ with 50 random generators per cell confirms the Higgs sector vanishes only at $Y = 0$.

### 4.2. Bimodule reformulation

In bimodule language, the falsifier becomes:

> Show that every Hermitian bimodule endomorphism $\omega \in \mathrm{End}_{\mathcal{A}_{\mathrm{AC}}\text{-bimod}}(\mathcal{H}_{\mathrm{AC}})$ derivable from GeoVac structure has zero cross-block as an element of $\mathrm{Hom}(\mathcal{H}_\mathbb{C}, \mathcal{H}_\mathbb{H}) \oplus \mathrm{Hom}(\mathcal{H}_\mathbb{H}, \mathcal{H}_\mathbb{C})$.

Strictly equivalent to H1's falsifier, with a slightly cleaner statement.

### 4.3. Dual modular propinquity falsifier (refined)

The dual modular propinquity provides one further refinement: it asks whether the cross-block is **forced to vanish under D-norm convergence** as $n_{\max} \to \infty$. Concretely:

> Show that for every $n_{\max}$-indexed sequence of Hermitian bimodule endomorphisms $\{\omega_{n_{\max}}\}$ derivable from GeoVac structure (i.e., from inner fluctuations of $D = D_{\mathrm{GV}}^{(n_{\max})} \otimes \mathbb{1}_F + \gamma_{\mathrm{GV}}^{(n_{\max})} \otimes D_F$), the cross-block $\Phi_{n_{\max}}$ converges to zero in the D-norm as $n_{\max} \to \infty$.

**Verdict:** the bimodule falsifier with D-norm convergence does not change the answer. The cross-block $\Phi_{n_{\max}}$ is structurally independent of $n_{\max}$ — it lives entirely in the fiber $\mathcal{A}_F$ and is set by the imposed $Y$. There is no GeoVac-side input that drives $\Phi_{n_{\max}} \to 0$ unless $Y = 0$ from the start.

Concretely, the falsifier behavior under D-norm convergence is identical to the matrix-norm behavior:

| Scenario | $\|\Phi_{n_{\max}}\|_D$ as $n_{\max} \to \infty$ | Falsifier |
|----------|:------:|:------:|
| $Y = 0$ | identically 0 | HOLDS |
| $Y \ne 0$ imposed | $\to \|Y\|_D = $ const | FAILS |

The D-norm convergence does not introduce a mechanism that drives $Y$ to zero, and does not introduce a structural constraint on the form of a non-zero $Y$.

### 4.4. New mechanism the bimodule reading does provide

The bimodule reading does offer one **new diagnostic mechanism** that H1's matrix analysis did not use: the **bimodule Lipschitz seminorm at the cross-block**. Define

$$L_{\mathrm{cross}}(\omega) := \|\omega|_{\mathrm{cross-block}}\|_{\mathrm{op}}$$

as the operator norm restricted to the cross-block submodule. The bimodule falsifier fires if $L_{\mathrm{cross}}(\omega) = 0$ for every inner-fluctuation $\omega$ derivable from GeoVac structure.

Computationally, $L_{\mathrm{cross}}$ can be computed by the existing `decompose_fluctuation()` API in `geovac/almost_commutative.py` (lines 606–670): it is the operator norm of the `higgs_matter_RL` block, which the H1 driver already computes as $\|\Phi\|_{\max}$. So the new diagnostic is **the same as H1's matrix-norm diagnostic**, just labeled differently.

This is the "vocabulary-only" finding: no new computation is needed, but the bimodule framework makes the computation conceptually cleaner.

---

## §5. Connection to G3 and G4a

### 5.1. G3 (electroweak chirality co-location)

G3 asks whether GeoVac's chirality grading $\gamma_{\mathrm{GV}}$ (the SIGN of the truthful Camporesi–Higuchi $D_{\mathrm{GV}}$) can be identified with the SM weak-isospin chirality $\gamma_F$ on the AC factor. Sprint G3 (CLAUDE.md §2, 2026-05-06 evening) closed in the NEGATIVE: $\gamma_{\mathrm{GV}}$ and $\gamma_F$ are independent commuting $\mathbb{Z}_2$s on $\mathcal{H}_{\mathrm{GV}} \otimes \mathcal{H}_F$ with $\|\gamma_{\mathrm{GV}} \otimes \mathbb{1}_F - \mathbb{1}_{\mathrm{GV}} \otimes \gamma_F\|_{\mathrm{op}} = 2$ exactly at every $n_{\max} \in \{1, 2, 3\}$.

**Bimodule reading of G3:** the two $\mathbb{Z}_2$ gradings act on different sub-bimodules. $\gamma_{\mathrm{GV}}$ acts diagonally on the left-action of $\mathcal{A}_{\mathrm{GV}}$ side of the bimodule; $\gamma_F$ acts on the fiber sub-bimodule (the matter-sector L/R decomposition). The two gradings are tensor-product factor independent.

**Does the bimodule framework collapse them?** No. The bimodule structure of $\mathcal{H}_{\mathrm{AC}} = \mathcal{H}_{\mathrm{GV}} \otimes \mathcal{H}_F$ is the **pure tensor product** of two independent bimodules. Each carries its own $\mathbb{Z}_2$-grading. There is no bimodule mechanism that identifies them (this would require a non-trivial coupling between the GV and fiber sub-bimodules, which Sprint G3-A/B/C ruled out at the matrix level).

The bimodule reading is therefore **vocabulary-only** for G3 as well. The G3 NEGATIVE verdict survives.

### 5.2. G4a (full Connes SM with $M_3(\mathbb{C})$)

G4a asks whether $\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H} \oplus M_3(\mathbb{C})$ (the full Connes SM finite algebra) admits the same construction. Sprint G4 scoping (CLAUDE.md §2) gave POSITIVE-THIN: 1–2 month sprint, predicted positive-thin in the H1 + G3 sense (Yukawas remain free input).

**Bimodule reading of G4a:** the $M_3(\mathbb{C})$ summand contributes a color-bimodule sub-structure $\mathcal{H}_{\mathrm{color}} = \mathbb{C}^3 \otimes \mathbb{C}^3 = M_3(\mathbb{C})$. The full bimodule $\mathcal{H}_F = (\mathbb{C} \oplus \mathbb{H}) \otimes M_3(\mathbb{C})$ has the same structural features as the $\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H}$ case: cross-blocks live entirely in the fiber, no GeoVac index, Yukawas free.

The bimodule reading does not change the G4a prediction. G4a remains a sprint-scale follow-up with the same expected POSITIVE-THIN outcome.

### 5.3. G4b (cross-manifold $T_{S^3} \otimes T_{\mathrm{Hardy}(S^5)}$)

G4b is the **structurally hard** case, blocked at the NCG-framework level by the Coulomb/HO category mismatch (Paper 24 §V four-layer asymmetry; Paper 32 §VIII.C G4b paragraph; Sprint L2-F.1 Track 2 scoping).

**Bimodule reading of G4b:** the cross-manifold tensor product would require a bimodule structure on $\mathcal{H}_{S^3} \otimes \mathcal{H}_{\mathrm{Hardy}(S^5)}$ with the right-action picture compatible with both manifolds' Riemannian / Hardy-sector structures. Paper 24's four-layer asymmetry tells us that (i) Hardy on $S^5$ does not admit a spectrum-computing $L_0$, (ii) does not admit calibration $\pi$, (iii) does not admit non-abelian Wilson gauge with natural matter, (iv) does not admit the modular-Hamiltonian structure of the wedge KMS state. The fourth layer (Sprint L2-F.1) is exactly the bimodule-level obstruction.

The dual modular propinquity framework, which respects Hilbert C*-bimodules, **inherits the same obstruction**: a bimodule whose right-action structure is structurally distinct from the left-action structure (which is what the Riemannian-vs-Hardy mismatch produces) cannot be naturally a metrized quantum vector bundle in Latrémolière's sense.

So the bimodule reading is **negative-with-clean-naming** for G4b: it does not solve the cross-manifold problem, but it gives a clean structural reason (incompatible right-action structures) for why the cross-manifold obstruction lives where it does.

### 5.4. Net for G3/G4

The bimodule / dual modular framework does not move G2 (no new $Y$ constraint), G3 (independent $\mathbb{Z}_2$s survive), or G4a (POSITIVE-THIN prediction survives) at the structural level. For G4b, it provides cleaner vocabulary for the cross-manifold obstruction but does not solve it.

---

## §6. Verdict

### 6.1. Net assessment

**Does dual modular propinquity give a NEW constraint on $Y$ that H1's algebra-level analysis didn't?**

**No.** The dual modular propinquity / bimodule reading provides:

1. **Cleaner vocabulary** for the H1 falsifier (bimodule cross-block instead of matrix off-block).
2. **Unification** of Connes' axiom system (Hermiticity + matter-antimatter doubling + order-one) under a single bimodule self-duality requirement.
3. **Cleaner statement of the structural reason GeoVac data does not constrain $Y$**: the cross-block lives in a sub-bimodule (the fiber) that is structurally decoupled from the GeoVac index, so there is no GeoVac-side input that could constrain it.

But it **does not** produce a new mechanism that constrains the form, magnitude, or structure of $Y$. The cross-block remains a free parameter at the AC-extension level, and dual modular propinquity convergence (D-norm convergence at $n_{\max} \to \infty$) does not introduce $n_{\max}$-dependent flow that drives $Y$ to a specific value.

### 6.2. Does the falsifier sharpen?

**Yes, modestly.** The bimodule reformulation of the falsifier (§4.2) is slightly cleaner than the matrix form: instead of asking about every Hermitian $D_F$, it asks about every Hermitian bimodule endomorphism with zero cross-block. The two statements are strictly equivalent, but the bimodule form makes one feature manifest:

- The cross-block is a property of the **fiber sub-bimodule**, not of the full AC bimodule. So the question "does GeoVac structure constrain the cross-block" is really "does GeoVac structure have any bimodule data that lives in the fiber sub-bimodule." It does not (the fiber is a separate tensor-product factor).

This sharpening is **vocabulary-level**, not analytical. It does not produce a new falsifier test.

### 6.3. Recommended next sprint

**Do not open a follow-up sprint on the dual modular propinquity / G2 question.**

The bimodule reframing is mathematically elegant but does not move the H1 POSITIVE-THIN verdict. Spending additional sprint effort on dual modular propinquity at G2 specifically would be diminishing returns.

What might be worth doing (in decreasing priority order):

1. **Update Paper 32 §VIII.C** with a one-paragraph remark noting that the H1 verdict survives the bimodule / dual modular propinquity reading at the structural level, citing this memo. This is the right level of documentation for the result: a vocabulary refinement worth recording but not a substantive new result. **Recommended.**

2. **G4a Connes SM with $M_3(\mathbb{C})$:** if PI wants to pursue the full electroweak unification at G4a level, the dual modular propinquity framework provides the right vocabulary for the cross-color-bimodule structure. The G4a sprint itself (1–2 months) is more substantial than this M-H1 vocabulary track. **Optional.**

3. **G3 chirality co-location:** the G3 NEGATIVE verdict is already closed. Bimodule reframing does not re-open it. **Not recommended.**

4. **G4b cross-manifold:** structurally blocked at the NCG-framework level by Paper 24 §V. Bimodule reframing provides cleaner naming for the obstruction but does not solve it. **Not recommended.**

### 6.4. Verdict line

**MIXED — SELF-ADJOINTNESS-REPHRASING with one structural sharpening of the falsifier (no new constraint on $Y$).**

The dual modular propinquity gives better mathematical vocabulary for the H1 verdict. It does not change the verdict. POSITIVE-THIN survives.

---

## §7. Honest scope and limitations

**What this memo establishes:**

- The bimodule structure of the AC extension is the standard left-$\mathcal{A}_{\mathrm{AC}}$ / right-$\mathcal{A}_{\mathrm{AC}}^{\mathrm{op}}$ doubling (Connes–Marcolli 2008 §13.4); matter–antimatter doubling realizes the right action.
- The Higgs cross-block $\Phi \in \mathrm{Hom}(\mathcal{H}_\mathbb{H}, \mathcal{H}_\mathbb{C})$ is a bimodule intertwiner on the matter sector, with antimatter-sector partner $\overline{\Phi}$ via $J_F$-conjugation.
- The dual modular propinquity D-norm convergence requires bimodule self-duality, which reduces to Connes' Hermiticity + matter-antimatter doubling — already imposed at the AC-extension level.
- No new constraint on $Y$ emerges from the bimodule / dual modular reading.

**What this memo does NOT establish:**

- That the dual modular propinquity has any application to GeoVac beyond vocabulary refinement. The PI's hypothesis was worth testing, and this memo reports the result, but the negative finding does not rule out subtler applications in other directions (e.g., metrized quantum vector bundles over the GeoVac S³ truncation that are not in the AC-extension form).
- A rigorous proof that no $Y$-constraint can come from any extension of the dual modular framework. The claim here is that the standard dual modular propinquity on the AC extension as constructed in H1 does not give one.
- The Sprint H1 verdict is unchanged at POSITIVE-THIN. This memo does not move the empirical falsifier at thermal $T_C \approx 160$ GeV (Sprint TD Track 3) which forces $Y_F > 0$ from outside — that remains the operational refinement of the H1 verdict.

**What this memo deliberately does NOT do:**

- Modify production code, papers, or memos other than itself.
- Declare $Y$ determined; the H1 "Y not autonomously selected" verdict is preserved.
- Open a sub-sprint on G3/G4 follow-ups; those decisions are left to PI.

---

## §8. Cross-references

- **H1 sprint memo:** `debug/h1_ac_extension_memo.md` (the foundational POSITIVE-THIN verdict).
- **Almost-commutative scoping memo:** `debug/almost_commutative_scoping_memo.md` (pre-H1 architecture decision).
- **Paper 32 §VIII.C:** the in-paper H1 verdict (lines ~2498–2700 of `papers/group1_operator_algebras/paper_32_spectral_triple.tex`).
- **Inner-factor Mellin engine memo:** `debug/inner_factor_mellin_engine_memo.md` (the η-trivialization theorem + AC factorization, 2026-05-07).
- **Coulomb/HO four-layer asymmetry:** Paper 24 §V `subsec:asymmetry_layer4` (Sprint L2-F.1, the G4b structural obstruction at modular-Hamiltonian level).
- **Sprint G3 closure (NEGATIVE):** CLAUDE.md §2 G3 entry, 2026-05-06 evening.
- **Sprint G4a scoping (POSITIVE-THIN predicted):** CLAUDE.md §2 G4 a/b split.

**Literature:**

- Latrémolière, F., "The Dual Modular Gromov–Hausdorff Propinquity and Completeness," J. Noncomm. Geom., 2021 (arXiv:1811.04534).
- Latrémolière, F., "The Modular Gromov–Hausdorff Propinquity," Dissertationes Math., 2018 (arXiv:1608.04881).
- Connes, A., and Marcolli, M., *Noncommutative Geometry, Quantum Fields and Motives*, AMS, 2008 — Ch. 13 for the matter-antimatter doubling and finite SM triple.
- Marcolli, M., and van Suijlekom, W. D., "Gauge networks in noncommutative geometry," J. Geom. Phys. 75 (2014).

---

**End of memo.** Approximately 3,650 words.
