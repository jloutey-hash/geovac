# Pythagorean-Orthogonality Extension: Scoping Memo

**Sprint:** Post-`h_local_residual_pslq` follow-up.
**Date:** 2026-05-17.
**Type:** Scoping (NO implementation).
**Verdict summary:** (a) SU(2) Wilson on S³ → **GO-WITH-PREREQS** (Wilson-coupled spinor matter sector + wedge inheritance). (b) SU(3) Bargmann on S⁵ → **NO-GO** (Paper 24 §V Coulomb/HO category mismatch resurfaces at the modular-Hamiltonian level, structurally analogous to G4b's blocker at the AC-extension level).

---

## 1. Motivation

The prior sprint (memo `debug/h_local_residual_pslq_memo.md`) found that on the BW-aligned hemispheric wedge $W \subset \mathcal{T}_{n_\max}$ of the truncated $S^3$ Camporesi–Higuchi spectral triple,

$$ \langle H_\text{local}, D_W^L \rangle_\text{HS} = \operatorname{Tr}(H_\text{local}^\dagger D_W^L) = 0 \quad \text{bit-exact} $$

at all tested $n_\max$ and $N_t$. The Pythagorean decomposition

$$ r^2 \;=\; \|H_\text{local}\|_F^2 + \|D_W^L\|_F^2 \;=\; \frac{S(n_\max)}{4\pi^2} + D(n_\max) $$

lives entirely in the master Mellin engine M1 ring (rational + $1/\pi^2$·rational), as predicted by the Paper 32 §VIII case-exhaustion theorem applied at $k=0$ (Hopf-base measure sector).

The structural question this opens: **is Pythagorean orthogonality $\langle H_\text{local}, D_W^L \rangle_\text{HS} = 0$ an S³-specific accident, or an architectural feature of compact spectral triples whose wedge construction is built from a Hopf-axis reflection?** The cleanest candidates for testing are the other two Wilson gauge settings unified in Paper 32 §VIII.B: SU(2) Wilson on the S³ Coulomb graph (Paper 30) and SU(3) Wilson on the Bargmann–Segal S⁵ graph (Sprint ST-SU3).

---

## 2. The four-question checklist on each setting

### 2.1 SU(2) Wilson on S³ = SU(2) (Paper 30)

**Q1. Natural matter Hilbert space H?**
The structurally correct answer is $\mathcal{H} = \mathcal{H}_\text{GV}^{(n_\max)} \otimes \mathbb{C}^{|E|}$ where $\mathcal{H}_\text{GV}^{(n_\max)}$ is the Camporesi–Higuchi spinor space at cutoff $n_\max$ (the full-Dirac basis indexed by $(n_\text{fock}, l, m_j, \chi)$, dim $= 2 \cdot \frac{2}{3}n_\max(n_\max+1)(n_\max+2)$) and $\mathbb{C}^{|E|}$ holds the gauge-link configuration. Wilson lattice gauge theory canonically couples to matter via $D_\text{cov}|_\text{edge}(i,j) = U_{i\to j} \cdot D_\text{GV}|_\text{edge}(i,j)$. CRITICAL: `geovac/su2_wilson_gauge.py` currently implements **pure gauge** (no matter sector wired). The Camporesi–Higuchi spinor space and the Wilson link variables both live on the same S³ Coulomb graph, but they do not currently talk to each other in production code.

**Q2. Does a hemispheric-wedge analog exist?**
YES, structurally inherited. The same `HemisphericWedge(axis="hopf")` from `geovac/modular_hamiltonian.py` lifts to the matter ⊗ gauge space by acting on the matter slot only: $P_W^\text{cov} = P_W \otimes I_{|E|}$. The Hopf-axis reflection $m_j \to -m_j$ is gauge-invariant (it acts on physical quantum numbers, not on gauge orbits), so it lands cleanly on the gauge-invariant Hilbert space after gauge fixing. **However**, the wedge construction must commute with the covariant Dirac $D_\text{cov}$ for the four-witness theorem to repeat. This commutation is non-trivial: gauge transformations preserve the wedge only if the gauge group action commutes with $R_\text{polar}$, which it generically does (SU(2) gauge transformations act on internal indices; $R_\text{polar}$ acts on $m_j$). Structurally clean.

**Q3. Candidate local modular Hamiltonian $K_\alpha^W$?**
The natural choice is the **same** $K_\alpha^W = \operatorname{diag}(\text{two\_m_j} > 0)$ acting on the matter slot only: $K_\alpha^{\text{cov}, W} = K_\alpha^W \otimes I_{|E|}$. Integer spectrum is preserved. The BW choice $H_\text{local} := K_\alpha^{\text{cov}, W}/\beta$ at $\beta = 2\pi$ gives a wedge KMS state $\rho_W^\text{cov} = e^{-K_\alpha^{\text{cov}, W}}/Z$ on the matter slot, with the gauge slot supporting a Haar measure (free / Wilson Boltzmann weight).

**Q4. Candidate Dirac D?**
The covariant Camporesi–Higuchi Dirac $D_\text{cov} = D_\text{GV} \otimes_{Wilson} U$. Concretely, where $D_\text{GV}$ has an off-diagonal matrix element between sites $i, j$ in the spinor multiplier representation, $D_\text{cov}$ has that same element multiplied by the link $U_{i \to j} \in$ SU(2). Equivalently: lift the spinor bundle from $S^3$ to its principal-SU(2)-bundle, then take the gauged Camporesi–Higuchi Dirac. This is the standard non-abelian-Wilson-Dirac construction. The wedge restriction $D_W^L|_\text{cov} = P_W^\text{cov} D_\text{cov} P_W^\text{cov}$.

**Q5. GO/NO-GO:** **GO-WITH-PREREQS.**

Prerequisites (each a 1–2 week sprint):
1. **Wire spinor matter into `geovac/su2_wilson_gauge.py`.** Currently the module is pure-gauge (link variables, plaquette enumeration, Wilson action, character expansion, Monte Carlo). Build `WilsonCovariantDirac(adjacency, links, ch_basis)` that returns the covariant Camporesi–Higuchi Dirac on $\mathcal{H}_\text{GV} \otimes (\text{gauge config})$.
2. **Verify Riemannian-limit recovery:** at trivial link configuration $U_e = I$ for all $e$, $D_\text{cov}$ must equal $D_\text{GV} \otimes I_{|E|}$ bit-exactly. Inherits the bit-exact Riemannian-limit-recovery discipline of Sprint L2-C.
3. **Check that the Hopf-polar reflection commutes with generic SU(2) gauge transformations** in the wedge sector. Expected: yes, gauge group and rotation group are independent.
4. **Choice of link configuration for the test:** the cleanest scoping probe is at the **trivial vacuum** $U_e = I$ for all $e$. In that limit, the covariant Dirac literally reduces to $D_\text{GV}$ on the matter slot, and the entire L2-E construction repeats verbatim. The Pythagorean orthogonality test then becomes: does it persist under Haar averaging over gauge configurations? This is the substantive question.

**Sketched construction at trivial vacuum (which is part of the Pythagorean check):** At $U_e = I$, $\langle H_\text{local}^{\text{cov}}, D_W^{\text{cov}, L}\rangle_\text{HS}$ factorizes as $\langle H_\text{local}, D_W^L\rangle_\text{HS}^{(S^3)} \times |\langle \text{vacuum gauge} | \text{vacuum gauge}\rangle|$, so the orthogonality at the trivial vacuum is **automatic** by the S³ result. The actual content is whether $\langle H_\text{local}^{\text{cov}}, D_W^{\text{cov}, L}\rangle_\text{HS}$ averaged over a non-trivial Wilson-Boltzmann distribution is still zero. Two structurally distinct possibilities:
- **Orthogonality survives Haar averaging:** the residual still lives in M1 ring with a Wilson-loop-expectation-valued prefactor. The Pythagorean structure is robust under gauging.
- **Orthogonality breaks at non-trivial $\beta_\text{Wilson}$:** the residual gains a $\langle W_\text{plaq} \rangle$-dependent off-diagonal contribution mixing with M3 (vertex-parity Hurwitz) sector through the plaquette structure of `geovac/su2_wilson_gauge.py`'s primitive non-backtracking walks. This would be the **substantive new finding**: gauging couples the wedge generator and the covariant Dirac.

**Verdict:** GO-WITH-PREREQS. Substantive (~3-4 week) sprint. The setup is structurally clean and inherits everything from L2-E except the gauge factor. The interesting science is whether Pythagorean orthogonality survives Haar averaging.

### 2.2 SU(3) Wilson on Bargmann–Segal S⁵ (Sprint ST-SU3)

**Q1. Natural matter Hilbert space H?**
The Bargmann–Segal graph has nodes $(N, l, m_l)$ with $N \geq 0$, $l \in \{N, N-2, \dots\}$, $m_l \in [-l, l]$, dim $= \dim(N,0)_{SU(3)} = (N+1)(N+2)/2$ per shell. This is a **scalar** (bosonic) Hilbert space: there is no spinor sector. The 3D HO has total spin $j = l$ (no spin-orbit, no Dirac), and the (N,0)-tower irreps of SU(3) are symmetric-traceless polynomial reps with no half-integer content.

In the Wilson construction of `geovac/su3_wilson_s5.py`, the link variables $U_e \in SU(3)$ live on edges between $(N, l, m_l) \leftrightarrow (N+1, l\pm1, m_l')$ — they are dipole-transition links, NOT a Dirac coupling.

**Q2. Does a hemispheric-wedge analog exist?**
**NO clean analog.** Three concrete obstructions:

(2a) **No half-integer spectrum to integer-spectrum convert.** Paper 42 / Sprint L2-E succeed because the spinor labels carry $\text{two\_m_j}$ which is an *odd integer* (so $e^{i \cdot 2\pi \cdot \text{two\_m_j}} = +1$ bit-exactly). The Bargmann graph has $m_l$ integer (bosonic), so the natural BW analog $K_\alpha^W = \operatorname{diag}(2 m_l > 0)$ would have *even* integer spectrum. The period closure becomes $e^{i\beta \cdot 2 m_l} = 1$ at $\beta = \pi$, not $\beta = 2\pi$. So the wedge KMS state at the BW canonical $\beta = 2\pi$ is no longer a "minimal" period — the BW thermal interpretation needs $\beta = \pi$, which is structurally different. Geometrically this corresponds to the fact that $S^5 / U(1)$ Hopf fibration over $\mathbb{CP}^2$ has a different equatorial structure than $S^3 / U(1)$ over $S^2$.

(2b) **The Bargmann graph has no Hopf-axis-distinguished hemispheric wedge.** $S^5$ does admit multiple Hopf foliations (Paper 24 §II), but none of them respects the (N, 0)-tower / Hardy-sector structure in the way the Hopf base $S^2 = S^3/S^1$ respects the spinor bundle on $S^3$. The natural "polar axis" on $S^5$ would be the $\mathbb{CP}^2$ symmetry axis, but the (N,0) reps of SU(3) decompose under $SU(2) \subset SU(3)$ in a non-trivial way (Clebsch–Gordan multiplicities that grow with N — exactly Sprint 5 Track S5's negative result on the (N,0)-tower SU(3)-adjacency action). Translating that into a "hemispheric wedge $P_W$ that respects the Hardy sector" leaves the wedge non-canonical.

(2c) **Coulomb/HO category mismatch resurfaces at the modular-Hamiltonian level.** Paper 24 §V argues that the Coulomb/HO asymmetry is multi-layered: spectrum-computing role of $L_0$ (yes S³, no S⁵), calibration $\pi$ (yes S³, no S⁵), Wilson gauge with natural matter coupling (yes S³, no S⁵ — Sprint ST-SU3 gauge-yes-matter-no), combinatorial vocabulary (universal). The Pythagorean orthogonality $\langle H_\text{local}, D_W^L\rangle_\text{HS} = 0$ is structurally a *modular-Hamiltonian* statement, which requires a thermal-time / KMS reading of the spectral triple. **The HO sector has no Dirac operator $D$ in the canonical NCG sense — the natural first-order operator on the (N,0) Hardy sector is the Euler operator $\hat{N} + 3/2$ (which is the Hamiltonian itself, not a Dirac).** So "$D_W^L$" has no analog: there is no second-order vs first-order distinction on the Hardy sector (Paper 24 §V Layer (i)/(iii)/(iv) asymmetry). The Pythagorean question, which is precisely the relation between the modular Hamiltonian (built from $K_\alpha$) and the wedge-restricted Dirac (built from $D$), evaporates because one of the two operands does not exist as an independent object.

**Q3. Candidate local modular Hamiltonian $K_\alpha^W$?**
The natural candidate is $K_\alpha^W = \operatorname{diag}(2 m_l > 0)$ on the Bargmann graph wedge, where the wedge is defined by an $m_l > 0$ projector. But this has only **even** integer spectrum (since $m_l$ is integer), so the BW closure period is $\pi$, not $2\pi$. The wedge KMS state at $\beta = 2\pi$ has fractional periodicity of $\sigma_t^\alpha$.

**Q4. Candidate Dirac D?**
None canonical. Three possible candidates, all problematic:
- **The Bargmann graph Laplacian L = D - A:** second-order (Laplacian, not Dirac). Not a Dirac.
- **The Kohn Laplacian on the Hardy space $H^2(S^5)$:** first-order in the holomorphic CR structure, but it acts on holomorphic functions in $\mathbb{C}^3$, NOT on a spinor bundle. The Kohn Laplacian is positive semidefinite (it's a Laplacian, not a Dirac), so its "square root" is not naturally a Dirac.
- **A shift operator $\hat{N} \pm 3/2$:** first-order in the Euler sense, but it has only integer eigenvalues (HO spectrum). It's the *Hamiltonian* itself, not a separate Dirac. So $\langle H_\text{local}, "D_W^L"\rangle_\text{HS}$ becomes $\langle K_\alpha^W/\beta, P_W \hat{N} P_W \rangle_\text{HS}$, which is a comparison of two diagonal operators (both in $(N, l, m_l)$). They live in the same diagonal subspace of operator space, so their HS inner product is **generically non-zero** (it's $\sum_v (2m_v > 0) \cdot (N_v + 3/2)$), and the Pythagorean structure collapses.

**Q5. GO/NO-GO:** **NO-GO at the modular-Hamiltonian level.**

Named structural obstructions:
1. **No spinor sector on Bargmann/Hardy.** The (N,0)-tower is bosonic; there is no spinor lift compatible with the holomorphic Hardy structure (would require breaking the Bargmann transform's $L^2$-isomorphism to $H^2(S^5)$). Spinors live on the full $S^5$ with the round-spinor bundle, not on the Hardy sector.
2. **No second-order vs first-order distinction.** Paper 24 §V Layer (i): the HO's spectrum-computing role is played by the *first-order* Euler operator on the Hardy sector, NOT by the second-order Laplacian. There is no "Dirac" sitting beside a "Laplacian" the way $D_\text{GV}$ sits beside $\Delta_\text{LB}$ on $S^3$. Without that distinction, "$H_\text{local}$ vs $D_W$" is not a well-defined comparison.
3. **No half-integer spectrum on the wedge.** $m_l$ is integer, so the natural BW closure period is $\pi$, not $2\pi$. This is a fundamental difference of architecture: the four-witness Wick-rotation theorem at $\beta = 2\pi$ is half-integer-spectrum-specific.
4. **Coulomb/HO category mismatch at the modular-Hamiltonian level.** This is the same blocker that killed G4b cross-manifold ($\mathcal{T}_{S^3} \otimes \mathcal{T}_\text{Hardy(S^5)}$) at the AC-extension level (Paper 32 §VIII.C / Paper 24 §V). It now resurfaces at the modular-Hamiltonian level: the Hardy-sector spectral triple does not have the structural ingredients (spinor bundle + half-integer wedge + second-order $\Delta$ alongside first-order $D$) that the Pythagorean orthogonality formula requires.

**Verdict:** NO-GO. The right reading is **dictionary content, not implementation target**: SU(3) Bargmann on $S^5$ joins the catalogue of structural-skeleton constructions where part of the L2-E architecture fails to lift, in a way that mirrors and refines the Paper 24 §V asymmetry. Adding this finding to Paper 24 §V is appropriate; opening a sprint is not.

---

## 3. Structural reading of the partition

The S³ Coulomb graph and the S⁵ Bargmann graph are both spectral triples in a generalized sense, but they sit in categorically different settings:

| Feature | $S^3$ Coulomb (Papers 25, 30) | $S^5$ Bargmann (Paper 24, ST-SU3) |
|:---|:---:|:---:|
| Manifold | Riemannian, parallelizable, $S^3 = SU(2)$ | Holomorphic Hardy sector of $S^5$ |
| Operator order | First-order Dirac $D_\text{GV}$ AND second-order $\Delta_\text{LB}$ on scalar bundle | First-order Euler $\hat{N}$ ONLY (no separate Dirac) |
| Spinor bundle | Yes, Camporesi–Higuchi, half-integer $m_j$ | No (Hardy sector is bosonic) |
| Hopf foliation supporting wedge | Yes, $S^2 = S^3/S^1$, preserves spinor bundle | $\mathbb{CP}^2 = S^5/U(1)$, breaks (N,0) tower |
| Wilson gauge | U(1), SU(2) (matter natural) | SU(3) (matter blocks at CG obstruction) |
| BW-aligned $K_\alpha^W$ spectrum | $2 m_j$ odd integer | $2 m_l$ even integer |
| Modular period at $\beta = 2\pi$ | $\sigma_{2\pi} = $ identity bit-exact | Would close at $\beta = \pi$ instead |
| Pythagorean $\langle H_\text{local}, D_W^L\rangle = 0$ | Yes (this sprint) | Vacuously, $D_W^L$ undefined |

The S³ side allows the L2-E construction because **first-order $D_\text{GV}$ and second-order $\Delta_\text{LB}$ coexist on the same parallelizable spinor bundle**, and the wedge respects the half-integer $m_j$ structure. The Bargmann S⁵ side has only a first-order operator (Euler = Hamiltonian), no spinor sector, and no half-integer wedge. The Pythagorean question requires precisely the S³ confluence.

This is **the same Coulomb/HO asymmetry pattern Paper 24 §V already identifies**, now extended to a fourth layer: **modular-Hamiltonian / wedge structure**. The asymmetry is consistent across all four layers documented in Paper 24:

1. Spectrum-computing role of $L_0$: S³ yes, S⁵ no.
2. Calibration $\pi$ from second-order Riemannian projection: S³ yes, S⁵ no.
3. Wilson lattice gauge with natural matter coupling: S³ yes, S⁵ no.
4. **Modular-Hamiltonian Pythagorean orthogonality on hemispheric wedge: S³ yes, S⁵ no.** (NEW from this scoping sprint.)

The four-layer asymmetry is now firmly an **architectural feature** of the framework: every modular / spectral-triple structural ingredient that is canonical on the second-order-Riemannian S³ has no clean Hardy-sector S⁵ analog.

---

## 4. SU(2) Wilson sub-sprint sketch (if GO is taken)

**Goal:** Test whether $\langle H_\text{local}^\text{cov}, D_W^{\text{cov}, L} \rangle_\text{HS} = 0$ persists when (a) Wilson links are non-trivial, (b) Haar-averaged at $\beta_\text{Wilson} \in \{0.5, 1, 2, 5\}$ matching Paper 30's tested range.

**Phase 1 (~1 week):** Build `WilsonCovariantDirac` in `geovac/su2_wilson_gauge.py`. API: take an SU(2) link configuration `links: Dict[OrientedEdge, np.ndarray]` and `n_max`, return the covariant Camporesi–Higuchi Dirac $D_\text{cov}$ on $\mathcal{H}_\text{GV} \otimes (\text{trivial gauge slot})$ at fixed link configuration. Verification: at $U_e = I_2$ for all $e$, $D_\text{cov} = D_\text{GV} \otimes I$ bit-exactly (load-bearing falsifier).

**Phase 2 (~1 week):** Wire in `LorentzianModularHamiltonian` of `geovac/modular_hamiltonian_lorentzian.py` to support a Wilson-link parameter. Verify the four-witness theorem at $U_e = I$ (must reproduce L2-E bit-exactly).

**Phase 3 (~1-2 weeks):** Compute $\langle H_\text{local}^\text{cov}, D_W^{\text{cov}, L} \rangle_\text{HS}$ at non-trivial link configurations. Two protocols:
- **Single-config probe:** pick a specific non-trivial $U_e$ pattern (e.g., $U_e = e^{i \theta_e \sigma_3 / 2}$ with $\theta_e$ on a chosen plaquette). Compute the inner product symbolically (sympy) and identify which sector it lands in (M1, M2, M3, or off-engine).
- **Haar-averaged probe:** use Paper 30's character-expansion machinery to compute $\langle\langle H_\text{local}^\text{cov}, D_W^{\text{cov}, L} \rangle_\text{HS}\rangle_\text{Wilson}$ at fixed $\beta_\text{Wilson}$.

**Phase 4 (~1 week):** Synthesis. If orthogonality survives, write up as Paper 32 §VIII.D Pythagorean-orthogonality corollary at finite cutoff *with* gauge coupling. If orthogonality breaks at finite $\beta_\text{Wilson}$, identify the off-orthogonal contribution and tag it to a Mellin-engine sector — the result becomes structurally **the way gauging couples the wedge generator to the covariant Dirac**, a new observation candidate.

**Total estimated effort:** 3–5 weeks, with `geovac/su2_wilson_gauge.py` extension as the load-bearing prerequisite. Aligned with Sprint L1 / L2-E timeline.

---

## 5. Recommended next move

**Two-sided answer:**

(A) **Accept SU(3) Bargmann on $S^5$ NO-GO as dictionary content.** Document the four-layer Coulomb/HO asymmetry extension in Paper 24 §V as a refinement (a new paragraph noting that the modular-Hamiltonian Pythagorean orthogonality is the fourth layer where S³ and S⁵ Hardy sectors diverge). Cross-reference from Paper 42 §10 / Paper 43 §10 as evidence that the H_local-vs-D_W distinction is structurally S³-specific. Cross-reference from Paper 32 §VIII.C as a sibling of the G4b cross-manifold blocker. **No implementation sprint.**

(B) **Open SU(2) Wilson on $S^3$ implementation sprint** as a *moderate* (3–5 week) follow-on. The setup is structurally clean, prerequisites are well-scoped, and the result is genuinely informative either way (positive = robustness of Pythagorean orthogonality under gauging; negative = explicit characterization of how gauging couples $K_\alpha$ and $D$). The sprint also produces side-product code (`WilsonCovariantDirac` couples Paper 30 SU(2) Wilson to Camporesi–Higuchi matter, which has been a structural gap since `geovac/su2_wilson_gauge.py` was first written as pure-gauge in Sprint 4 RH-Q).

**Sequencing recommendation:** Apply (A) first (~1 day, paper edits only, no code), then open (B) at PI discretion. Together they extend the Paper 24 §V asymmetry catalogue and clarify what the L2-E Pythagorean closure depends on — turning a single structural observation into a tested architectural feature of the framework.

---

## 6. Honest scope

- **This memo is scoping only.** No production code modified. No computational probe of the SU(2) Wilson case beyond a structural sketch at the trivial vacuum, where the orthogonality is automatic by the S³ result.
- **The NO-GO on SU(3) Bargmann is structural, not numerical.** It rests on the categorical absence of a spinor sector and a half-integer wedge on the Hardy sector. A different choice of matter Hilbert space (e.g., spinor bundle on the full $S^5$, NOT the Hardy sector, which the framework has not yet considered) could in principle revive the question — but that would break the Bargmann / Paper 24 graph identification and is out of scope for this sprint and likely for the framework as it stands.
- **The GO-WITH-PREREQS on SU(2) Wilson is contingent on the matter-coupling extension to `geovac/su2_wilson_gauge.py`.** The prerequisite work itself is independently valuable (closes a structural gap in the Wilson gauge module) and would land regardless of whether the Pythagorean probe returns positive or negative.
- **A formal proof of $\langle H_\text{local}, D_W^L\rangle_\text{HS} = 0$ at general $n_\max$ on S³** (currently empirical at $n_\max \leq 5$) is a parallel follow-on, independent of this scoping. Worth noting that the empirical certainty (5/5 panel cells bit-exact across $N_t$ values) is strong evidence the orthogonality is structural rather than coincidental.

End of memo.
