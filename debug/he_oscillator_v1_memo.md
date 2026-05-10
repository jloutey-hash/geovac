# Calc Track He-Oscillator: Helium 2¹P → 1¹S Oscillator Strength

**Sprint:** post-Sprint-Calc-rZG-extended, multi-track precision catalogue
extension (Track 1 of 5 §V.C decomposition targets, CLAUDE.md §1.8).
**Date:** 2026-05-09
**Status:** verification computation; Paper 34 §V.C.5 fill text drafted; **no production code or papers modified**.
**Verdict:** **Sturmian's continuum-closing property does NOT extend** to multi-electron transitions in a single-exponent hydrogenic basis. Both ω (transition energy) and |⟨z⟩|² (dipole matrix element) converge with n_max but to **wrong limits** — driven by the same excited-state non-variationality mechanism that broke the He 2¹S–2³S splitting (Sprint Calc-He). The framework reproduces f ≈ 0.45 vs Drake & Yan f = 0.27616, a +60% residual that does not improve at n_max = 7 (longest tested).

---

## 1. Reference value

| Source | Value | Notes |
|:-------|:------|:------|
| Drake & Yan 1992 (PRA 46, 2378) | f = 0.27616 | Length form, NR infinite-mass limit, 4-digit |
| Schiff & Pekeris 1964 (PR 134, A638) | f = 0.27617 | Variational, length-velocity gauge agreement |
| Theodosiou 1987 (At. Data Nucl. Data Tab. 36) | f = 0.276 | Compilation |
| Lewis 1989 experimental (PRL 63, 1397) | f ≈ 0.276 | Lifetime measurement |
| NIST ASD line list (2¹P¹ → 1¹S⁰, He I λ=584.33 Å) | A_{ij} = 1.798×10⁹ s⁻¹ | g_i f equivalent |

These are all infinite-mass NR limits; relativistic and recoil corrections are below 0.1% at He.

---

## 2. Focal-length formula

For a transition between an L=0 initial state and L=1 final state, the length-form oscillator strength is

```
f = (2/3) · ω · |⟨1¹S | r | 2¹P⟩|²
  = 2 · ω · |⟨1¹S | z | 2¹P_{m=0}⟩|²
```

with each piece traced to a Paper 34 §III projection:

| Piece | Source projection | What it is | Sturmian-closure-relevant? |
|:------|:------------------|:-----------|:---------------------------|
| τ₀ (atomic time) | Fock conformal (§III.1) | sets a₀, E_h jointly via p₀ = Z/n | foundational anchor |
| ω = E(2¹P) − E(1¹S) | Multi-electron CI on Fock graph | transition energy, framework-native | **YES — depends on Sturmian closure** |
| ⟨z⟩ multi-electron | Multi-electron CI + Wigner 3j (§III.8) + radial dipole | mixed angular + radial + correlation | **YES — depends on Sturmian closure** |
| (2/3) | Wigner 3j angular sum + photon polarization (§III.11) | rational angular average | always exact rational |

The framework-native pieces are ω and ⟨z⟩; both should converge cleanly with n_max if Sturmian closure extends to the multi-electron case. The single calibration constant in the assembly is α (which sits in atomic-unit conversion to SI); the structural transcendental class of f itself is rational + algebraic-extension.

---

## 3. GeoVac computation: convergence study

Using `geovac.casimir_ci.build_graph_native_fci`-style infrastructure with hydrogenic basis at exponent k_orb = Z = 2 (Coulomb Sturmian convention used in Sprint Calc-P). The 1¹S state lives in the (l₁=l₂=0, M_L=0) ss singlet sub-block; the 2¹P state lives in the (l₁=0, l₂=1, M_L=0) sp singlet sub-block. The transition dipole is computed via Slater-Condon rules:

```
⟨Ψ(1¹S) | z₁+z₂ | Ψ(2¹P)⟩ = Σ_{IJ} c_I^{1S} c_J^{2P} ⟨I,1S | z₁+z₂ | J,1P⟩
```

with each ⟨I|z|J⟩ reducing to single-particle dipoles ⟨a|z|c⟩ × delta on spectator electron. Single-particle dipoles use `geovac.dirac_matrix_elements.radial_matrix_element` (sympy-exact for ⟨n'l'|r|nl⟩) combined with the Gaunt-3j angular factor c₁(l_a, m_a, l_c, m_c) from `casimir_ci._gaunt_ck`.

### Path A: Coulomb Sturmian at k_orb = Z = 2

| n_max | configs | ω (Ha) | ⟨z⟩ (a.u.) | \|⟨z⟩\|² | f_length | err vs Drake | wall (s) |
|:-----:|:-------:|:------:|:----------:|:--------:|:--------:|:------------:|:--------:|
| 2 | 5 | 0.9044 | 1.0521 | 1.1069 | 2.0022 | +625% | 0.2 |
| 3 | 12 | 0.7045 | 0.6483 | 0.4203 | 0.5921 | +114% | 0.2 |
| 4 | 22 | 0.6930 | 0.6299 | 0.3968 | 0.5500 | +99.2% | 0.4 |
| 5 | 35 | 0.6891 | 0.5687 | 0.3234 | 0.4458 | +61.4% | 8.0 |
| 6 | 51 | 0.6893 | −0.5778 | 0.3339 | 0.4603 | +66.7% | 31.3 |
| 7 | 70 | 0.6893 | −0.5675 | 0.3221 | 0.4440 | +60.8% | 119.6 |
| **Drake reference** | — | **0.7799** | — | **0.1771** | **0.2762** | 0.00% | — |

Both ω and \|⟨z⟩\|² **plateau at wrong values** by n_max = 6:

* ω plateau: 0.6893 Ha (vs Drake 0.7799 Ha; **−11.6% low**)
* \|⟨z⟩\|² plateau: ~0.32 (vs needed ~0.18; **+82% high**)
* f plateau: ~0.45 (vs Drake 0.276; **+62% high**)

### Diagnostic: variational status of E(1¹S) and E(2¹P)

| n_max | E(1¹S) (Ha) | Δ vs Drake | E(2¹P) (Ha) | Δ vs Drake |
|:-----:|:-----------:|:----------:|:-----------:|:----------:|
| 2 | −2.886824 | +0.0169 ✓ | −1.982392 | +0.1415 ✓ |
| 3 | −2.889578 | +0.0141 ✓ | −2.185122 | **−0.0613** ✗ |
| 4 | −2.891487 | +0.0122 ✓ | −2.198485 | **−0.0746** ✗ |
| 5 | −2.892360 | +0.0114 ✓ | −2.203224 | **−0.0794** ✗ |
| 6 | −2.892855 | +0.0109 ✓ | −2.203568 | **−0.0797** ✗ |
| 7 | −2.893162 | +0.0106 ✓ | −2.203827 | **−0.0800** ✗ |
| Drake exact (NR ∞-mass) | −2.903724 | 0 | −2.123843 | 0 |

**The 2¹P state is variationally violated** by 0.08 Ha = 3.8% (E(2¹P)_GeoVac < E(2¹P)_Drake by ~80 mHa). The 1¹S state is variationally well-behaved (above Drake by ~10 mHa, the standard small-Z graph-validity-boundary residual). This is the signature of the same κ-induced over-binding mechanism documented in Sprint Calc-He 2¹S–2³S (CLAUDE.md §2): the graph κ = −1/16 adjacency adds a constant inter-shell coupling regardless of orbital extent, which over-binds Rydberg-like excited states.

### Path B: Z_eff conventions at n_max = 4

| k_orb | label | ω (Ha) | f_length | err vs Drake |
|:-----:|:-----:|:------:|:--------:|:------------:|
| 2.0 | bare Z | 0.6930 | 0.5500 | +99.2% |
| 1.6875 | variational 1s (27/16) | 0.8157 | 0.8308 | +201% |
| 1.69 | Slater rules 1s | 0.8147 | 0.8278 | +200% |
| 1.0 | screened 2p | 1.0942 | 3.0309 | +997% |

None of the natural Z_eff choices improve the result. k_orb = 27/16 lowers E(1¹S) (better) but leaves E(2¹P) at −2.23 Ha, *increasing* the over-binding gap relative to Drake. k_orb = 1 makes things dramatically worse (the basis is qualitatively wrong for the 1s electron). The single-exponent basis is **structurally inadequate** for representing the two distinct length scales (compact 1s, diffuse 2p) of He.

---

## 4. Comparison vs Drake/Theodosiou and residual attribution

The +62% residual at n_max = 7 cleanly decomposes:

| Component | Value | vs reference | Attribution |
|:----------|:------|:------------|:-----------|
| ω | 0.6893 Ha | 0.7799 Ha (−12%) | Excited-state E(2¹P) over-bound by 80 mHa (κ adjacency) |
| \|⟨z⟩\|² | 0.32 | 0.18 (+82%) | Same overlap construction sees the over-bound 2¹P wavefunction as too localized |
| f = 2ω·\|⟨z⟩\|² | 0.444 | 0.276 (+62%) | Compound effect; both factors wrong from same structural mechanism |

**Both effects share a single root cause**: the κ-induced inter-shell coupling in the graph one-body Hamiltonian h₁_off = κ · (−A) causes excited s and p Rydberg states to admix too strongly with low-n shells, giving over-bound and over-localized excited states. This is the same mechanism documented in CLAUDE.md §2 "Graph validity boundary" (Z_c ≈ 1.84; He at Z = 2 is just above) and in Sprint Calc-He (precision_catalogue_he_2s_singlet_triplet_memo.md, +11.86% on the 2¹S–2³S splitting).

**The structural reading**: Sturmian closure for one-electron polarizability (Sprint Calc-P, exact at N_basis = 2) works because the Dalgarno–Lewis function (r² + 2r) e⁻ʳ lives in span{S_{2,1}, S_{3,1}} — a finite-dimensional functional space the bare Sturmian basis spans by construction. For multi-electron transitions, the analogous "natural" wavefunction pair (1¹S and 2¹P CI states) does *not* live in any low-dimensional Sturmian sub-space, because the two states inhabit different effective one-body potentials (1¹S: doubly-screened nuclear attraction; 2¹P: single-electron-screened, much more diffuse). A single common Sturmian exponent cannot represent both.

---

## 5. Projection-chain breakdown

The He oscillator strength calculation invokes a **four-projection chain**:

```
f = (Fock conformal) ∘ (Wigner 3j angular) ∘ (vector-photon promotion)
       ∘ (multi-electron CI in Sturmian basis)
```

| Step | Projection | What it contributes | Status in this sprint |
|:-----|:-----------|:--------------------|:----------------------|
| 1 | **Fock conformal** (§III.1) | hydrogenic radial wavefunctions, atomic units, kinetic + h1 diagonal | **fails for excited 2¹P** (κ-induced over-binding, Z<Z_c+small artifact) |
| 2 | **Wigner 3j angular** (§III.8) | Δl = ±1 selection, c₁(l_a,m_a,l_c,m_c) Gaunt = 1/√3 for s↔p, exact rational | exact, no error |
| 3 | **Vector-photon promotion** (§III.11) | factor 2/3 polarization average, photon DOS | exact rational |
| 4 | **Multi-electron CI in single-exponent Sturmian basis** | configuration mixing, two-electron singlet sector | **fails** — single exponent inadequate for two length scales |

The cleanest reading is that the failure has **two compounded sources**:

(a) The **graph-validity-boundary at small Z** (CLAUDE.md §2) bites the *excited* state at He much harder than it bites the ground state. The 1¹S sits inside the Z > Z_c regime (graph perturbation small at Z=2 for tight 1s²); the 2¹P (with one diffuse 2p electron) sits closer to the boundary. **One-shot Sturmian closure fails because the "wrong" focal length on the 2p side propagates into both ω and ⟨z⟩.**

(b) The **single common exponent** cannot represent both compact 1s and diffuse 2p. This is structurally distinct from (a) but compounds with it: even if the graph were variationally well-behaved for 2¹P, the 1s² and (1s)(2p) configurations would benefit from different effective Z's, which is the very thing Sturmian convention forbids.

Either fix would be a sprint-scale extension:
* For (a): apply Sturmian closure within the **multi-focal architecture** (`geovac/cross_register_vne.py`-style multi-λ basis) — different λ_a per orbital block.
* For (b): use the **balanced coupled** architecture (`geovac/balanced_coupled.py`), which already supports multi-center / multi-exponent compositions for chemistry, generalized to internal (within-atom) multi-focal sectors.

---

## 6. Honest scope: what this sprint did and did not test

**Did test:**
* Whether Sturmian's continuum-closing property at single exponent extends from one-electron polarizability to two-electron transitions. **Answer: NO.** The structural mechanism that made Calc-P succeed (Dalgarno–Lewis function lives in a 2-dimensional Sturmian sub-space) does *not* generalize to two-electron transitions where the initial and final states want different effective Z.
* Whether the failure mode is angular (vector-photon promotion broken) or radial (Sturmian basis inadequate). **Answer: radial.** The angular Wigner 3j step is exact; both ω and ⟨z⟩ (radial-driven) plateau at wrong values.
* Whether the failure has the same signature as Sprint Calc-He 2¹S–2³S. **Answer: yes.** Excited-state variational violation by 3.8% on absolute energy, via the κ-induced inter-shell coupling.

**Did not test:**
* Multi-focal architecture (different λ per block within the atom). The "atomic split-Z" extension flagged in Sprint Calc-He memo §6 is the natural follow-on.
* Velocity-form vs length-form gauge consistency. (Without the multi-focal extension, both forms would have the same residual; the gauge-invariance check is meaningful only after the basis problem is fixed.)
* Explicitly-correlated trial functions (Hylleraas-type r₁₂). These could close the residual but would be a categorically different basis and outside the scope of "Sturmian closure".

**Two production-code follow-ups flagged for future tracks** (mirroring the Sprint Calc-He cleanup pattern):
* The CI matrix construction uses the same `_build_graph_h1` h1 + `two_electron_integral` machinery as Sprint Calc-He, which had two latent issues caught and fixed in the May-08 cleanup (`hypergeometric_slater.compute_rk_float` n≥6 catastrophic cancellation, fixed via threshold dispatch; `compute_he_spectrum` state mislabeling, fixed via subblock kwarg). Both fixes are now production; this sprint inherits the corrected pipeline. No new bugs surfaced.
* The single-particle dipole `radial_matrix_element(...)` returns sympy-exact symbolic results that we float-cast at the boundary. For n_max ≤ 7 this is safe; the cusp regime tested in cleanup (n ≥ 12) is well outside our range.

---

## 7. Pattern-finding: three-class tag

Per the §3 "Pattern-finding" rule, the failure decomposes:

| Class | Where the failure lives | Status |
|:------|:------------------------|:-------|
| Layer-2 input convention mismatch (`docs/curve_fit_audit_memo.md` class) | NO — angular machinery is exact, no Layer-2 itemization debate | Not the issue |
| GeoVac kernel approximation gap | YES — single-exponent basis is structurally too coarse for two-length-scale states | **Primary diagnosis** |
| Focal-length decomposition cataloguing (Paper 34 §V.C) | YES — this entry IS a §V.C autopsy of why one structurally simple operator (length-form dipole) fails when the underlying CI basis has the wrong focal-length structure | **Primary diagnosis** |

**The cleanest single-sentence takeaway**: Sturmian's continuum-closing is a *single-focal-length* property. For multi-electron transitions where initial and final states have structurally different effective focal lengths, single-exponent Sturmian basis does not close — the failure is in the CI sector (§III multi-electron projection), not in the angular (§III.8 Wigner 3j) or vector-photon (§III.11) projections. The natural extension is the multi-focal architecture (`cross_register_vne` for chemistry, internal split-Z for atomic excited states); this sprint identifies "internal multi-focal for atomic excited-state transitions" as a named follow-on track.

---

## 8. Proposed Paper 34 §V.C.5 fill text (DO NOT edit Paper 34 directly per sprint instructions)

```latex
\subsection{Helium $2{}^1\!P \to 1{}^1\!S$ oscillator strength}
\label{sec:catalog:he_oscillator}

The first multi-electron extension of Sprint Calc-L (hydrogen Lyman~$\alpha$,
\S\ref{sec:catalog:lyman_alpha}, $+0.055\%$ at depth-3) and Sprint Calc-P
(hydrogen polarizability, machine-exact at $N_{\rm basis}=2$ in the Sturmian
basis). The structural question: does the Sturmian continuum-closing property
extend to multi-electron transition matrix elements?

The framework-native length-form oscillator strength is computed via
$$
f = 2\omega\, |\langle 1{}^1\!S | z | 2{}^1\!P_{m=0}\rangle|^2
$$
in the graph-native CI sector at fixed $M_L$, with the $1{}^1\!S$ state in
the $(l_1=l_2=0)$ ss singlet sub-block, the $2{}^1\!P$ state in the
$(l_1=0, l_2=1, M_L=0)$ sp singlet sub-block, and the transition dipole
assembled by Slater-Condon expansion of the multi-electron CI vectors.

Convergence with the Sturmian basis (single exponent $k_{\rm orb} = Z = 2$):

\begin{tabular}{cccccc}
$n_{\max}$ & configs & $\omega$ (Ha) & $|\langle z\rangle|^2$ & $f_{\rm length}$ & err vs Drake \\
\hline
2 & 5 & 0.9044 & 1.107 & 2.002 & $+625\%$ \\
3 & 12 & 0.7045 & 0.420 & 0.592 & $+114\%$ \\
4 & 22 & 0.6930 & 0.397 & 0.550 & $+99\%$ \\
5 & 35 & 0.6891 & 0.323 & 0.446 & $+61\%$ \\
6 & 51 & 0.6893 & 0.334 & 0.460 & $+67\%$ \\
7 & 70 & 0.6893 & 0.322 & 0.444 & $+61\%$ \\
\hline
exact & --- & 0.7799 & 0.177 & 0.27616 & 0 \\
\end{tabular}

Both $\omega$ and $|\langle z\rangle|^2$ \emph{plateau} by $n_{\max} = 6$ at
\emph{wrong limits}: $\omega$ converges $-12\%$ low and $|\langle z\rangle|^2$
converges $+82\%$ high. The $2{}^1\!P$ state is variationally violated by
$\sim 0.08$ Ha (E(2$^1$P)$_{\rm GeoVac}$ $= -2.204$ vs Drake $-2.124$ Ha) —
the same $\kappa$-induced over-binding mechanism documented in
\S\ref{sec:catalog:he_2s_singlet_triplet} for the $2{}^1\!S - 2{}^3\!S$
splitting (small-$Z$ graph-validity-boundary artifact, $Z_c \approx 1.84$;
He at $Z = 2$ is just above). The angular (Wigner 3j) and vector-photon
(\S\ref{sec:projection:vector_photon}) projections are exact; the failure
is in the CI sector.

\textbf{Structural reading}: Sturmian's continuum-closing is a
\emph{single-focal-length} property. For multi-electron transitions where
initial and final states have structurally different effective focal lengths
(here, doubly-screened $1s$ and singly-screened $2p$), the single-exponent
Sturmian basis cannot represent both, and the framework's CI eigenvectors
overlap incorrectly even at large $n_{\max}$. The natural extension is the
\emph{multi-focal} architecture
(\verb|geovac/cross_register_vne.py| for chemistry; internal split-$Z$ for
atomic excited states is a named follow-on track per Sprint Calc-He
recommendation).

The $+61\%$ residual at $n_{\max} = 7$ does not improve at higher $n_{\max}$
within the single-exponent architecture and is therefore tagged error code
\texttt{B} (basis architecture) in Table~\ref{tab:catalog_off}. This entry
is a \emph{structural} negative result that complements the
\S\ref{sec:catalog:he_2s_singlet_triplet} positive observation: both rows
expose the small-$Z$ graph-validity-boundary at He as the dominant
limitation for excited-state precision, separable from the framework's
clean closures on hydrogenic and ground-state observables.
```

Suggested companion §V row (off-precision matches table) format:

| Match | Projection(s) | Vars | Dim | Trans. class | Match | Error code |
|:--|:--|:-:|:-:|:--|:-:|:-:|
| He $2{}^1\!P \to 1{}^1\!S$ oscillator strength | Fock $\circ$ Wigner 3j $\circ$ vector-photon $\circ$ multi-electron CI | $Z$, $\alpha$ | dimensionless | rational + algebraic | $+61\%$ vs Drake $f = 0.276$ at $n_{\max} = 7$ | B |

---

## 9. Files

* `debug/calc_track_he_oscillator_v1.py` — driver (initial Path A + Path B).
* `debug/he_osc_extended.py` — extended convergence n_max = 6, 7 + variational diagnostic.
* `debug/he_osc_sanity.py` — sanity check (Lyman α reproduction).
* `debug/data/he_oscillator_v1.json` — Path A + Path B results.
* `debug/data/he_oscillator_extended.json` — full n_max convergence + Drake comparisons.
* `debug/he_oscillator_v1.log`, `debug/he_osc_extended.log` — run logs.

No GeoVac production code or paper modified.

---

## 10. One-line summary

GeoVac multi-electron CI in single-exponent Sturmian basis at $k = Z = 2$ gives $f(2{}^1\!P \to 1{}^1\!S, \mathrm{He}) \approx 0.444$ at $n_{\max} = 7$ vs Drake $0.27616$ ($+61\%$, plateaued); failure mode is excited-state variational violation of $E(2{}^1\!P)$ by $0.08$ Ha (same $\kappa$-induced over-binding as Sprint Calc-He 2¹S-2³S), demonstrating that **Sturmian's continuum-closing property is single-focal-length and does NOT extend to multi-electron transitions** — the natural multi-focal extension (different $\lambda$ per orbital block) is the named follow-on track.
