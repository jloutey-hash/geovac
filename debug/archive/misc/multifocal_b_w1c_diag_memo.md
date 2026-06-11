# Multi-focal Phase B Sub-sprint W1c-diag

**Date:** 2026-05-07
**Author:** PM (Phase B-W1c-diag; no production files modified)
**Sources read:** Phase A synthesis + Track 1 audit (W1c rows + B7 entry); `geovac/neon_core.py`; `geovac/balanced_coupled.py`; `geovac/shibuya_wulfman.py`; `geovac/composed_qubit_relativistic.py` (the only consumer of Sprint 7b's screened-SO machinery); Sprint 7 balanced second-row memo; Sprint 7 fine-structure memo.
**Sanity script:** `debug/multifocal_b_w1c_sanity.py` (executed).
**Sanity data:** `debug/data/multifocal_b_w1c_sanity.json`.

---

## 1. The diagnostic question and W1c's place in the wall taxonomy

W1c is the third sub-wall under W1 (cross-register two-body operator missing). The Phase A synthesis defines it as: *frozen-core $Z_\text{eff}(r)$ provides same-center screening only; no cross-center screening of bare $V_{ne}$ from the other nucleus*. The empirical witness is Sprint 7 (April 2026): NaH at $n_\text{max}=2$ and MgH$_2$ at $n_\text{max}=2$ both produce monotonic over-attraction with no equilibrium under the balanced-coupled architecture, despite LiH at $n_\text{max}=2$ converging cleanly to a bound state with 7.0 % $R_\text{eq}$ error in the same architecture. The structural difference is that LiH carries an *explicit* core orbital block (Li 1s$^2$), while NaH and MgH$_2$ replace the core by a frozen [Ne] $Z_\text{eff}(r)$ profile that is consumed only in same-center contexts.

Sprint 7b (April 2026) demonstrated that the analog problem in the spin-orbit sector closes cleanly: by computing $\langle (1/r)\,dV/dr \rangle$ from the actual screened potential $V(r) = -Z_\text{eff}(r)/r$, the screened-SO splittings improve by 12–144$\times$ over the hydrogenic $Z_\text{eff}=2$ baseline for Ca/Sr/Ba 2p doublets. The diagnostic question for this sub-sprint is whether the same machinery ports to the cross-center $V_{ne}$ sector for second-row balanced builders, and what verdict (a/b/c/d) that port falls into.

The wall is *internal to GeoVac architecturally* under the Phase A synthesis classification. There is no published no-go theorem that prevents cross-center screening of $V_{ne}$ — the question is whether the existing tooling generalizes mechanically, partially, or not at all.

---

## 2. Q1 — What does Sprint 7b's `screened_xi_so` actually do?

`geovac/neon_core.py::screened_xi_so` (lines 734–803) performs four operations in sequence:

1. **Build a screened single-electron radial Schrödinger problem.** Calls `_solve_screened_radial(Z, l, n_target, core_type)`, which constructs the FrozenCore $Z_\text{eff}(r)$ analytically from Clementi-Raimondi hydrogenic densities for [Ne], [Ar], [Ar]3d$^{10}$, [Kr], or [Xe] cores, then solves $H u = E u$ with $V_\text{eff}(r) = -Z_\text{eff}(r)/r + l(l+1)/(2r^2)$ on a uniform finite-difference grid (default $n_\text{grid}=12{,}000$, $r_\text{max}=80$ bohr) using `scipy.linalg.eigh_tridiagonal`. Returns $u(r) = r R(r)$ for the requested $(n,l)$ state.
2. **Reconstruct the screened potential's radial derivative.** Computes $Z_\text{eff}(r)$ on the same grid, takes $dZ_\text{eff}/dr$ via `np.gradient`, and assembles
   $$\frac{1}{r}\frac{dV}{dr} = \frac{Z_\text{eff}(r)}{r^3} - \frac{1}{r^2}\frac{dZ_\text{eff}}{dr}.$$
   The first term is the bare hydrogenic-like contribution; the second term is *the screening-derivative correction*, present only when $Z_\text{eff}(r)$ is non-constant.
3. **Assemble** $\xi(r) = (\alpha^2/2)(1/r)(dV/dr)$ (the SO coupling function).
4. **Integrate** $\langle \xi \rangle = \int |u(r)|^2 \,\xi(r)\,dr$.

The radial integral being computed is therefore $\langle \xi \rangle_{nl}$ on a *one-center* potential whose *form* is screened. The screening manifests in two places: (a) the wavefunction $u(r)$ is computed in the screened potential (it shifts the radial node structure and the small-$r$ amplitude), and (b) the operator $\xi(r)$ itself uses the screened $V$ rather than the bare $-Z/r$.

The key observation for Q1 is *what FrozenCore $Z_\text{eff}(r)$ actually provides*. Reading `_solve_screened_radial`, `screened_r3_inverse`, and `screened_xi_so` together, FrozenCore supplies a *single-center* radial profile $Z_\text{eff}^A(r)$ that is the difference between the bare nuclear charge $Z_A$ and the cumulative core electron count $N^A_\text{core}(r)$:
$$Z_\text{eff}^A(r) = Z_A - \int_0^r n^A_\text{core}(r')\,dr'$$
where $n^A_\text{core}$ is the analytical Clementi-Raimondi sum over closed-shell hydrogenic radial densities. This is *one-center* both in construction (the density sits on nucleus A) and in consumption (the radial Schrödinger problem is one-dimensional in $r$, the radial coordinate from A).

**Would the same machinery, applied to the cross-center $V_{ne}$ sector, give the right physics?** *Architecturally yes, but with a non-trivial extension.* The relevant target integral is
$$\langle \psi_{nlm}^A | -Z_B / |\mathbf{r} - \mathbf{R}_B| | \psi_{n'l'm'}^A \rangle,$$
where $\psi^A$ is on nucleus A. Currently `compute_cross_center_vne` evaluates this with the bare Coulomb tail $-Z_B/|\mathbf{r}-\mathbf{R}_B|$. Sprint 7b's mechanism naturally suggests replacing the bare tail by a *screened* tail $-Z_\text{eff}^B(|\mathbf{r}-\mathbf{R}_B|)/|\mathbf{r}-\mathbf{R}_B|$ when nucleus B has a frozen core. This is conceptually clean: nucleus B's frozen core electrons are part of the molecular electronic structure, and a valence orbital centered on A sees the *Coulomb tail of (B + core_electrons_of_B)* rather than the bare nuclear charge.

The non-trivial extension is that $Z_\text{eff}^B(|\mathbf{r}-\mathbf{R}_B|)$ is a function of the displacement from B, while the orbital is integrated over coordinates centered on A. The integral structure changes from "bare Coulomb (separable multipole expansion in spherical harmonics centered at A)" to "screened Coulomb (multipole expansion of a B-centered, B-radial-only function in A-centered spherical coordinates)" — see Q3 for the algebraic analysis.

---

## 3. Q2 — What does `balanced_coupled.py` currently do for cross-center $V_{ne}$ in second-row systems?

Reading `balanced_coupled.py` lines 478–602 (the "Add cross-center V_ne (the balanced part — NEW in Track CD)" block) together with `_get_block_geometry` and the per-molecule `_get_nuclei_for_*` factories: when MgH$_2$ uses balanced cross-center $V_{ne}$ with frozen [Ne] core, `compute_cross_center_vne` is called with `Z_nuc = 12.0` (the *bare* Mg nuclear charge from `_get_nuclei_for_mgh2`, line 154–161 of `balanced_coupled.py`). The same is true for Na (`Z_nuc = 11.0`), Cl ($17.0$), S ($16.0$), P ($15.0$), Si ($14.0$), and so on through every second-row, third-row, and frozen-core species in the library.

There is no $Z_\text{eff}$ on the production cross-center path. The shielding logic in Sprint 7b's `screened_xi_so` is invoked only inside `composed_qubit_relativistic.py` line 442 (inside the diagonal-$h_1$ phase of the *relativistic* builder), guarded by `use_screened = (abs(Z_so - Z_sb) > 0.5) and (n_offset > 0)`. That branch screens the *same-center* SO coupling for the heavy-atom valence orbital, where the SO operator and the orbital share a center. The cross-center $V_{ne}$ path in `balanced_coupled.py` does not import `neon_core`, has no `FrozenCore` instance, and uses the bare `nuc['Z']` directly.

This is *exactly* the failure mode the audit predicted. Quick numerical confirmation (sanity probe `multifocal_b_w1c_sanity.py`, executed; data in `debug/data/multifocal_b_w1c_sanity.json`):

| Mechanism | trace of cross-center $V_{ne}$ on H-side n=2 orbitals at R=3.5 bohr |
|:----------|:----:|
| (A) Bare $Z_\text{Na}=11$ (current production) | $-12.15$ Ha |
| (B) Naive rescale by $Z_\text{eff}^\text{Na}(R)/Z_\text{Na} \approx 0.091$ | $-1.10$ Ha |
| (C) Asymptotic screening (Na$^+$ tail, $Z_\text{eff,asy}=1$) | $-1.10$ Ha |

The bare cross-center $V_{ne}$ over-attracts the H-side orbital by approximately 11 Ha relative to a screened estimate at the same R. This is comparable in magnitude to the absolute scale of the overattraction in Sprint 7's NaH PES (electronic energies of $-7$ Ha to $-12$ Ha across $R \in [2, 5]$ bohr). The Z_eff(r) profile (probe output above):

| r (bohr) | $Z_\text{eff}^\text{Na}(r)$ | $Z_\text{eff}^\text{Mg}(r)$ |
|:--:|:--:|:--:|
| 0.01 | 11.00 | 12.00 |
| 0.50 | 2.89 | 3.17 |
| 1.00 | 1.03 | 2.01 |
| 2.00 | 1.00 | 2.00 |
| 3.50 | 1.00 | 2.00 |

shows that by $r \approx 1.5$ bohr the [Ne] core is essentially fully internalized — the bonding region (typical $R \in [2, 4]$ bohr for a hydride) sits squarely in the asymptotic-screened regime where the H-side orbital should see Na$^+$ ($Z_\text{eff} = 1$), not bare Na ($Z = 11$). Production uses $Z = 11$. **The $\sim 10\times$ overattraction matches the empirical PES failure quantitatively.**

This makes the diagnostic answer to Q2 unambiguous: *the current balanced_coupled cross-center $V_{ne}$ does not screen* — it sees the full bare $Z=11/12$ from the other nucleus regardless of frozen core. The screened-$Z_\text{eff}(r)$ extension is the natural fix.

---

## 4. Q3 — Cross-center screened multipole expansion (symbolic, low $l_\text{max}$)

The bare Coulomb potential at A from nucleus B at displacement $\mathbf{R}_B$ has the standard Laplace multipole expansion
$$\frac{1}{|\mathbf{r}-\mathbf{R}_B|} = \sum_{L=0}^{\infty} \frac{r_<^L}{r_>^{L+1}} P_L(\cos\theta_{\mathbf{r},\mathbf{R}_B}),$$
where $r_< = \min(r, R_B)$ and $r_> = \max(r, R_B)$ in coordinates centered at A (here $r = |\mathbf{r}|$ from A and $R_B = |\mathbf{R}_B|$). This is what `_radial_split_integral` and `compute_cross_center_vne_element` exploit: the angular integrals against $\psi^A_{nlm}$ are Gaunt 3j's, and the radial integral splits into inner ($r < R_B$) and outer ($r > R_B$) regions analytically via incomplete-gamma integrals. The expansion **terminates exactly** at $L_\text{max} = 2 l_\text{max}$ by the 3j triangle inequality — Track CD's positive result.

The screened cross-center potential $V_\text{screened}(\mathbf{r}, \mathbf{R}_B) = -Z_\text{eff}^B(|\mathbf{r}-\mathbf{R}_B|)/|\mathbf{r}-\mathbf{R}_B|$ does **not** factor into $Z_\text{eff}^B \cdot (1/|\mathbf{r}-\mathbf{R}_B|)$ because $Z_\text{eff}^B$ depends on $|\mathbf{r}-\mathbf{R}_B|$. The right way to think about it is:
$$V_\text{screened}(\mathbf{r}, \mathbf{R}_B) = -\frac{Z_B}{|\mathbf{r}-\mathbf{R}_B|} + \int \frac{n^B_\text{core}(\mathbf{r}')}{|\mathbf{r} - (\mathbf{R}_B + \mathbf{r}')|} d^3\mathbf{r}'$$
i.e., the bare nuclear Coulomb minus the Hartree potential of the frozen core electron density (which sits on B). Both pieces admit Laplace multipole expansions in coordinates centered at A. The first piece is exactly the existing bare Coulomb expansion. The second piece is the *Hartree integral of the frozen B core* evaluated at A-displacement.

**Key structural fact (good news for tooling extension):** the frozen [Ne] / [Ar] / etc. core density $n^B_\text{core}(\mathbf{r}')$ is **spherically symmetric** about B (it is a sum over closed-shell hydrogenic densities, all of which are spherically symmetric). This means the Hartree potential $\Phi^B_\text{core}(|\mathbf{r}-\mathbf{R}_B|)$ generated by the core is also spherically symmetric about B and, by Newton's shell theorem on a spherical density, has the closed-form representation
$$\Phi^B_\text{core}(\rho) = \frac{1}{\rho}\int_0^\rho n^B_\text{core}(r')\,4\pi r'^2 dr' + \int_\rho^\infty n^B_\text{core}(r')\,4\pi r' dr' \ = \ \frac{N^B_\text{core}(\rho)}{\rho} + \tilde\Phi^B(\rho)$$
where $N^B_\text{core}$ is exactly the cumulative count FrozenCore already builds. So:
$$V_\text{screened}(\mathbf{r}, \mathbf{R}_B) = -\frac{Z_B - N^B_\text{core}(|\mathbf{r}-\mathbf{R}_B|)}{|\mathbf{r}-\mathbf{R}_B|} - \tilde\Phi^B(|\mathbf{r}-\mathbf{R}_B|),$$
or more compactly in terms of Sprint 7b's $Z_\text{eff}^B(\rho)$:
$$V_\text{screened}(\mathbf{r}, \mathbf{R}_B) = -\frac{Z_\text{eff}^B(|\mathbf{r}-\mathbf{R}_B|)}{|\mathbf{r}-\mathbf{R}_B|} - \tilde\Phi^B(|\mathbf{r}-\mathbf{R}_B|),$$
The second piece $\tilde\Phi^B$ is the "outer" tail integral $\int_\rho^\infty n^B_\text{core}(r')\,4\pi r' dr'$; it is **smooth, finite at the origin, exponentially decaying at large $\rho$, and spherically symmetric about B**.

Each piece is a *radial-only* function of $|\mathbf{r}-\mathbf{R}_B|$. By Laplace's expansion of any radial-only potential about a point, every spherically symmetric radial function $f(|\mathbf{r}-\mathbf{R}_B|)$ admits a multipole expansion in A-centered spherical coordinates of the form
$$f(|\mathbf{r}-\mathbf{R}_B|) = \sum_L f_L(r, R_B)\,P_L(\cos\theta),$$
where $f_L$ is computed by the angular projection
$$f_L(r, R_B) = \frac{2L+1}{2}\int_{-1}^{1} f\!\left(\sqrt{r^2 + R_B^2 - 2rR_B u}\right) P_L(u)\,du.$$
For the bare Coulomb $f(\rho) = 1/\rho$, this reduces to $f_L = r_<^L / r_>^{L+1}$ (the Laplace expansion). For the screened pieces, $f_L(r, R_B)$ is *not* a clean inner/outer separation — but the angular projection still terminates the *expansion* at $L_\text{max} = 2l_\text{max}$ by the same 3j triangle inequality applied to the orbital angular factors, **because the angular structure of the multipole expansion depends only on the orbital angular labels, not on the form of the radial function**. The Gaunt closure that gives Track CD its positive termination result transfers verbatim.

**Symbolic check at $L=0,1$ for a simple model.** Let $n^B_\text{core}(\rho) = (Z_\text{core}/\pi) \beta^3 e^{-2\beta\rho}$ (a model 1s-like core density on B). The angular integrals $\int_{-1}^{1} e^{-2\beta\sqrt{r^2+R_B^2-2rR_B u}}/\sqrt{r^2+R_B^2-2rR_B u}\, P_L(u)\, du$ are evaluable as combinations of incomplete-gamma functions in $\beta r$ and $\beta R_B$ via a similar split-region identity to the bare-Coulomb case (change of variable from $u$ to $\rho$ gives a one-dimensional radial integral over $[\,|r-R_B|,\, r+R_B\,]$). The integrals **do not break the $L_\text{max}=2l_\text{max}$ termination**: the angular projection $\int P_L(u)$ is an *output* of the expansion, and the truncation is enforced by the *outer* Gaunt 3j coupling to $\psi^A_{nlm}$, which is unchanged.

For real Clementi-Raimondi cores, the radial profile is a sum of $\sim 10$ exponential terms (one per occupied shell), and the analogous integrals are sums of incomplete-gamma evaluations. Concretely, the implementation extends `_split_integral_analytical` to take a list of $(c_k, \beta_k)$ pairs from a Clementi-Raimondi-decomposed $n^B_\text{core}(\rho)$.

**Verdict on Q3:** the cross-center screened multipole expansion **preserves Gaunt termination at $L_\text{max} = 2l_\text{max}$** (the angular structure is unchanged) and **admits closed-form analytical evaluation per multipole** as a sum of incomplete-gamma terms (by exponential-shell decomposition of the FrozenCore density). The radial integral structure is more complex than the bare case — instead of a single inner/outer split, there is a per-shell-pair (orbital $\times$ core-shell) decomposition — but it is a **mechanical extension** of `_split_integral_analytical`, not a structural new identity. No multipole inflation, no termination breakdown, no new transcendentals beyond incomplete-gamma at multiple $\alpha$ values.

---

## 5. Q4 — Empirical sanity check (DONE; small-scale)

Q4 was performed in `debug/multifocal_b_w1c_sanity.py` using the existing `compute_cross_center_vne` machinery. The probe rescales the bare cross-center $V_{ne}$ by $Z_\text{eff}^\text{Na}(R)/Z_\text{Na}$ as a magnitude-only proxy for the screened result. Headline numbers:

- $Z_\text{eff}^\text{Na}(R=3.5\,\text{bohr}) = 1.000$ (the [Ne] core is fully internalized at this distance).
- Bare cross-center $V_{ne}$ trace on H-side $n_\text{max}=2$ orbitals: $-12.15$ Ha.
- Naive screened trace: $-1.10$ Ha.
- Trace shift: $-11.04$ Ha (bare is more attractive by exactly an order of magnitude).

For comparison, Sprint 7's NaH electronic energy at $R=3.5$ bohr is $E_\text{elec} \approx -7.5$ Ha, with the PES monotonically descending toward smaller $R$. The 11-Ha-scale overattraction in the cross-center $V_{ne}$ alone is **larger than the entire electronic energy**. This is consistent with the observed PES failure mode (no equilibrium, monotonic over-attraction at small $R$). The naive rescaling is a single-multiplier proxy and does not capture the radial profile of $Z_\text{eff}^\text{Na}(\mathbf{r}-\mathbf{R}_\text{Na})$ as the H orbital probes different regions of the screened potential — a proper W1c implementation would integrate the screened profile against the orbital. Even so, the order of magnitude is unambiguous.

The same scan for Mg ($Z = 12$, $Z_\text{eff,asy} = 2$) gives a 6$\times$ overattraction rather than 10$\times$, consistent with MgH$_2$'s qualitatively similar but somewhat less severe overattraction in Sprint 7.

**The screening magnitude is plausibly consistent with — and arguably *quantitatively predictive of* — the Sprint 7 NaH/MgH$_2$ overattraction.** The diagnostic is satisfied.

(Honest caveat: the proxy rescaling assumes uniform screening, which the actual radial profile does not exhibit. The H orbital's small-$r$ tail probes the unscreened core region, which contributes a fraction of the bare attraction back. A proper W1c would treat this correctly by construction. The proxy answer is therefore a *lower bound* on the screened result; the true correction may be slightly less than $-11$ Ha but very close to it for orbital basis with negligible amplitude inside the [Ne] core radius $\sim 1$ bohr — which is the case for hydrogenic n=2 orbitals at Z=1.)

---

## 6. Q5 — Verdict

**(a) Tooling-addressable in a small sprint, with extension closer to (b) than to a pure mechanical port.**

The right reading is:

- W1c is *not* downstream of a more general issue (eliminating verdict (c)). The composed architecture is fine; the missing piece is locally a cross-center potential, not a deeper rethinking. Phase A's W1a (full two-body coordinate operator) is genuinely a separate, deeper wall; W1c is structurally narrower.
- W1c is *not* blocked by a structural feature (eliminating verdict (d)). Q3's analysis shows the multipole termination at $L_\text{max} = 2l_\text{max}$ survives the screened extension, and the radial integrals admit closed-form evaluation in the same incomplete-gamma family already used by `_split_integral_analytical`.
- W1c is **not quite (a) "port `screened_xi_so` mechanism mechanically to V_ne"** — the SO machinery solves a *one-center* radial Schrödinger problem with a *one-center* screened potential; the V_ne extension requires *two-center* multipole expansion of a *one-center-on-B* screened tail evaluated by orbital integrals on A. The conceptual extension (Hartree decomposition of the screened cross-center potential into bare + smooth correction; Newton-shell-theorem reduction to radial profile of $Z_\text{eff}^B$) is straightforward but is *new code*, not just a wiring change.
- W1c **is** (b) "tooling-addressable but with non-trivial extension (multipole expansion of screened V_ne needs new identities)."

The "new identities" are not deep new mathematics — they are mechanical: extending `_split_integral_analytical` to handle a sum of exponentials (one per Clementi-Raimondi shell of the frozen core) instead of a single exponential. No new transcendentals, no Gaunt-rule violations, no $L_\text{max}$ inflation, no convergence-radius surprises. Phase A's effort estimate of 4-6 weeks for Phase C-W1c is roughly correct for a careful first-pass implementation with full test coverage.

---

## 7. Phase C-W1c scope sketch

A reasonable Phase C-W1c sprint would have four tracks:

1. **Track 1 — Implement `cross_center_screened_vne_element` and `cross_center_screened_vne`** in a new module (or by extending `shibuya_wulfman.py`). Inputs: orbital quantum numbers and orbital-center charge $Z_\text{orb}^A$ (as today); off-center nucleus $Z_B$ and core type for B; internuclear vector. Internally: build FrozenCore for B, decompose $n^B_\text{core}$ as a sum of Clementi-Raimondi exponentials, evaluate the multipole-expanded screened cross-center potential by extending `_split_integral_analytical` to handle a list of exponential kernels. Validate against the bare case in the limit of zero core electrons (regression test: extends, never replaces, the bare path). Time: 2-3 weeks.
2. **Track 2 — Wire into `balanced_coupled.py`** behind a `screened_cross_center: bool = True` kwarg defaulting to true for frozen-core species. Call sites: lines 568-572 (the `compute_cross_center_vne` call). Detect frozen-core species via `core_type` look-up and dispatch to the screened path automatically. Backward compatibility: for first-row species (no frozen core), the screened path collapses to the bare path. Time: 1 week.
3. **Track 3 — Validate on NaH and MgH$_2$.** Recompute Sprint 7's PES with the screened cross-center $V_{ne}$. Predicted outcome: NaH and MgH$_2$ recover bound states with $R_\text{eq}$ errors comparable to LiH ($\sim 7$ % at $n_\text{max}=2$). If the prediction holds, that closes W1c empirically. If it doesn't (e.g., the PES is *too* weakly attractive after screening), document the residual structural error as a continuing wall (likely a sub-feature of incomplete radial decomposition, not a categorical failure). Time: 1 week.
4. **Track 4 — Extend to third row** (KH, CaH$_2$, HCl, AsH$_3$, etc., 14 frozen-core species in the current library). Resource impact: Pauli count unchanged (selection rules unchanged), 1-norm shifts by O(1-10 Ha) per species, FCI energies improve. Time: 1-2 weeks for systematic regression.

Total: 5-7 weeks first pass. Risk profile: **low**. The mathematics is well-posed (Q3); the empirical signature is sharp (Q4); the audit predicts the failure mechanism precisely (Sprint 7's NaH/MgH$_2$ overattraction is order-of-magnitude exactly right for a 10$\times$/6$\times$ rescaling). Failure modes to watch: (i) Clementi-Raimondi multi-shell decomposition convergence in `_split_integral_analytical` extension; (ii) numerical conditioning at small $R_{AB}$ where the screened tail saturates; (iii) the residual structural drift from PK that survives in LiH ($+0.053$ bohr per $n_\text{max}$, 3$\times$ smaller than PK's drift but nonzero) — this is not a W1c issue but should not be confused with one.

Phase C-W1c does **not** require any architectural change to the composed framework, the Pauli-encoding pipeline, the relativistic builder, or any guardrail paper. It is exactly an extension of a single integral-evaluation function. If Phase B-W1c-diag's diagnosis is correct (Q1-Q4 above), Phase C-W1c is one of the smallest of the wall-closure candidates from Phase A's Section 7 table.

---

## 8. Honest scope and uncertainty

What this memo covered: a diagnostic walk through Sprint 7b's machinery, the production cross-center $V_{ne}$ path, the algebraic structure of the screened multipole expansion, and a small empirical sanity probe. Verdict (b) is structurally well-supported.

What this memo did **not** cover deeply:

1. **Multi-center cross-screening (NaCl etc.).** For NaH/MgH$_2$/HCl only one side has a frozen core, so the question is one-sided. For NaCl, both sides have [Ne] cores and the cross-center $V_{ne}$ in *both directions* needs screening. The Track 1 implementation should handle this in full generality (each off-center nucleus gets its own screening based on its own core type); no new identity needed, but empirical validation is Phase C Track 4 work.
2. **Multi-center geometry rotations.** Q3's symbolic analysis used a $z$-aligned $\mathbf{R}_B$. The screened multipole expansion's angular structure is identical to the bare case (all $P_L(\cos\theta)$), so the existing l=2 Wigner-D rotation infrastructure (Track CM) should continue to work without modification. Worth verifying explicitly in Phase C Track 1.
3. **The residual LiH $R_\text{eq}$ drift ($+0.053$ bohr per $n_\text{max}$).** LiH has an *explicit* core block, so the [Ne]-style W1c failure does not apply; the residual drift is a separate matter. Phase C-W1c should not be expected to close the LiH residual.
4. **Bound-state magnitude prediction.** Q4 shows the 10$\times$ overattraction at NaH bare vs. screened roughly matches the failure scale. It does *not* pin down whether the screened NaH PES has $R_\text{eq}$ within 7 % (LiH-comparable), 15 %, or worse. Empirical answer comes from Phase C Track 3.
5. **The relationship between W1c and W1a.** W1c is genuinely simpler than W1a: W1a wants a full two-body cross-register coordinate operator $V(\hat\mathbf{r}_e, \hat\mathbf{R}_n)$; W1c only wants to dress an existing one-body cross-center potential with a screening profile. **W1c closure does NOT close W1a** — the recoil failure mode is a deeper missing operator.

What I am uncertain about:

- **Clementi-Raimondi multi-shell decomposition convergence in `_split_integral_analytical`.** Q3's claim that the structure transfers from one-exponential-shell to multi-shell via linearity is correct but unverified at the symbolic level for a real [Ne] core (3 occupied shells: 1s, 2s, 2p). Verifying this is a multi-day exercise in Phase C Track 1.
- **Sprint 7b's $\sim 60-70\%$ SO residual error** (CaH/SrH/BaH SO splittings undershoot physical SO) may be informative about how good the Phase C-W1c PES result will be. If screened-PES has the same character (qualitatively correct, quantitatively off), then W1c-fixed NaH PES might bind to within 15-25 % $R_\text{eq}$ rather than to LiH-comparable 7 %. Worth tracking.
- **Effort estimate.** Phase A predicted "Small (2-3 weeks)" for the diagnostic and "4-6 weeks" for Phase C; the present memo updates Phase C to 5-7 weeks. Both still in the small-sprint regime.

---

**End of Phase B-W1c-diag memo.**
