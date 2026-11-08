# Sprint HF Track 1 — Bohr-Fermi $A_{\mathrm{hf}}$ for hydrogen 21 cm

**Date:** 2026-05-07
**Goal:** Derive the hyperfine constant $A_{\mathrm{hf}}$ for the ground state of atomic hydrogen at the Bohr-Fermi level (point-particle Dirac proton, Pauli electron at $g_e = 2$, no recoil, no electron anomalous moment, no nuclear structure, no multi-loop QED) from GeoVac's native machinery, and compare to the experimental 21 cm transition $\nu_{\mathrm{hf}} = 1\,420\,405\,751.768$ Hz.
**Status:** Closed-positive; one honest caveat about continuum embedding flagged.

---

## 1. Why this sprint exists (refresh from briefing)

GeoVac is exact at single-focal-length resolution: selection rules, transcendental signatures, scaling laws, and structural prefactors are derived correctly without fits. Every accumulated residual (Lamb shift $-0.534\%$, alpha-conjecture $8.8\times 10^{-8}$, composed-PK errors $5$–$26\%$, G4b cross-manifold obstruction) sits at a *multi-focal seam* where the framework's composition machinery has limits. Hydrogen 21 cm is the cleanest possible multi-focal observable: it requires a continuum hydrogenic wavefunction (electronic focal length), a Dirac/Pauli spin coupling (relativistic focal length), and a nuclear magnetic moment (proton focal length), and the experimental value is known to 12 digits. Track HF-1 calibrates the Bohr-Fermi piece and locates the residual as a known physics target for HF-2 onward.

## 2. The Bohr-Fermi formula

The Fermi-contact Hamiltonian for an electron-proton bound state in the non-relativistic limit, in atomic units, is

$$ H_F \;=\; \tfrac{8\pi}{3} \, g_e \, g_p \, \mu_B \, \mu_N \; \mathbf{S}_e \!\cdot\! \mathbf{I}_p \; \delta^{3}(\mathbf{r}). $$

Taking the expectation in the 1s state and folding $\langle \delta^3(\mathbf{r}) \rangle = |\psi_{1s}(0)|^2$ into a coefficient of $\mathbf{I}\!\cdot\!\mathbf{S}$, one defines the hyperfine constant

$$ A_{\mathrm{hf}} \;=\; \tfrac{8\pi}{3} \, g_e \, g_p \, \mu_B \, \mu_N \; |\psi_{1s}(0)|^{2}. $$

In atomic units with $c = 1/\alpha$, $\mu_B = \alpha/2$ and $\mu_N = (m_e/m_p)\,\alpha/2$, so

$$ A_{\mathrm{hf}} \;=\; \tfrac{2}{3} \, g_e \, g_p \, \alpha^{2} \, \tfrac{m_e}{m_p} \, Z^{3} \quad \text{Hartree}. $$

For hydrogen ($Z = 1$) at $g_e = 2$ this gives the canonical Bohr-Fermi result

$$ A_{\mathrm{hf}}^{\mathrm{BF}} \;=\; \tfrac{4}{3} \, g_p \, \alpha^{2} \, \tfrac{m_e}{m_p} \quad \text{Hartree}. $$

(The PI briefing wrote the prefactor as $(4/3)\,\alpha^{2}\,g_e\,(g_p/2)\,(m_e/m_p)$. Read literally this differs from the canonical result by a factor of two; in context "$(g_p/2)$" is shorthand for $\mu_p/\mu_N \approx 2.793$, the proton magnetic moment in nuclear magnetons. The Fermi-contact *operator* carries the full $g_e g_p$, and the briefing's expected outcome 1418.84 MHz is computed from the canonical $g_e g_p$ form. We use the canonical form throughout.)

## 3. Per-factor decomposition: what came from where

| Factor | Value | Source | Category |
|--------|-------|--------|----------|
| $\frac{8\pi}{3}$ | exact | Standard Pauli/Dirac NR-limit Fermi contact (Bethe-Salpeter Eq. 22.18) | B (standard physics) |
| $g_e$ | $= 2$ | Tree-level Dirac equation; on $S^3$ the Camporesi–Higuchi Dirac operator gives the same value at tree level, but this is the gamma-matrix algebra, not S^3 specific | B (standard physics) |
| $g_p$ | $5.5856946893$ | Measured proton g-factor (CODATA 2018); QCD-internal, NOT computable in GeoVac | C (external calibration) |
| $\mu_B = \alpha/2$ | atomic-unit identity | Standard | B + C ($\alpha$ external) |
| $\mu_N = (m_e/m_p)\alpha/2$ | atomic-unit identity | Standard | B + C ($m_e/m_p$ external) |
| $\alpha$ | $7.2973525693\times 10^{-3}$ | CODATA 2018 (Paper 2 in Observations is the framework's open conjecture about deriving $\alpha$; we do NOT use that conjecture here) | C (external calibration) |
| $m_e/m_p$ | $1/1836.15267$ | CODATA 2018 | C (external calibration) |
| $\lvert\psi_{1s}(0)\rvert^{2}$ | $Z^{3}/\pi$ | **A: From GeoVac.** Symbolic evaluation of `geovac.dirac_matrix_elements._hydrogenic_radial_wavefunction(n=1, l=0, r, Z)` at $r=0$, multiplied by the standard $|Y_{0,0}|^{2} = 1/(4\pi)$. Reproduces $R_{1,0}(0) = 2 Z^{3/2}$ symbolically | A (with caveat — see §4) |

Three categories used: **A** = derivable from existing GeoVac machinery; **B** = standard physics result that GeoVac does not contradict (the Dirac equation works the same on $S^3$ as in flat space at tree level); **C** = external CODATA calibration not derivable from current GeoVac.

The only piece that genuinely runs through GeoVac code is $|\psi_{1s}(0)|^{2}$. Everything else is either a CODATA input or a piece of standard atomic-unit / Pauli-Dirac physics. That is honest and exactly what one would expect for a single-focal-length hydrogenic observable.

## 4. The $|\psi_{1s}(0)|^{2}$ caveat: continuum embedding, not graph-native

The hydrogenic radial wavefunction $R_{n,l}(r)$ in `geovac.dirac_matrix_elements._hydrogenic_radial_wavefunction` is the *continuum* solution to the Schrödinger equation in atomic units. The Fock projection makes the $S^3$ graph Laplacian's spectrum exactly $n^2 - 1$ (Paper 7, 18 symbolic proofs); but the spatial *profile* $\psi_{n,l}(\mathbf{r})$ in three-dimensional Euclidean space is the framework's continuum embedding, not a graph quantity. CLAUDE.md §4 says the Fock projection re-parameterizes $\mathbb{R}^{3}$ onto $S^{3}$ via stereographic projection at $p_0 = \sqrt{-2E_n}$; the wavefunction at the origin in $\mathbb{R}^{3}$ corresponds (after this projection) to the south-pole of $S^{3}$, and its amplitude in the embedded graph would be $R_{1,0}(0) \cdot Y_{0,0}$ — the same number we computed. So the answer is consistent with the framework, but it is a *continuum-embedding* answer.

The framework does not yet have a *graph-native* notion of $|\psi(0)|^{2}$ that does not reach back to the continuum 1s wavefunction. A pure-graph definition would presumably be the amplitude on the $n=1$ node summed over angular degrees of freedom, normalized appropriately — but that quantity is dimensionless on the bare graph (no metric), so converting it to a $1/a_0^{3}$ density requires the same Fock-projection conformal factor $\Omega(p) = (1+p^{2})/(2p_{0})$ that re-introduces the continuum. This is precisely the "Layer 1 (bare graph) + Layer 2 (named projections)" architecture of Paper 34: the mass scale and metric come from a Layer-2 projection, and the Fermi-contact density falls in that bucket. It is an honest consequence of the framework's structural-skeleton scope (see CLAUDE.md memory file `geovac_structural_skeleton_scope_pattern.md`), not a deficiency to hide.

## 5. Numerical results

With CODATA 2018 inputs $\alpha = 7.2973525693\times 10^{-3}$, $g_p = 5.5856946893$, $m_e/m_p = 1/1836.15267343$, and $1$ Ha $= 6.5797\times 10^{15}$ Hz:

| Quantity | Hartree | Hz | MHz |
|----------|---------|----|----|
| Strict Bohr-Fermi ($g_e=2$, no recoil) | $2.15992\times 10^{-7}$ | $1.42116\times 10^{9}$ | **1421.16** |
| Reduced-mass-corrected (recoil included) | $2.15639\times 10^{-7}$ | $1.41884\times 10^{9}$ | 1418.84 |
| Experimental (CODATA / NIST) | $2.15894\times 10^{-7}$ | $1.42041\times 10^{9}$ | 1420.4058 |

**Strict Bohr-Fermi residual:** $A_{\mathrm{hf}}^{\mathrm{BF, strict}} - \nu_{\mathrm{hf}}^{\mathrm{exp}} = +0.754$ MHz $= +531$ ppm.

**Reduced-mass-corrected residual:** $A_{\mathrm{hf}}^{\mathrm{BF, recoil}} - \nu_{\mathrm{hf}}^{\mathrm{exp}} = -1.566$ MHz $= -1102$ ppm.

The PI briefing's quoted value 1418.84 MHz / +1.57 MHz / +1100 ppm corresponds to the **reduced-mass-corrected** Bohr-Fermi number. Without recoil, BF over-shoots by 531 ppm. The signs of the strict and recoil residuals are opposite because two known corrections to BF — recoil $(\sim -1635$ ppm$)$ and the electron anomalous moment $a_e \sim +1160$ ppm — partially cancel. The headline residual the PI flagged $(+1100$ ppm) is dominated by $a_e$, with the recoil already absorbed.

The two cleanly separable residual contributions, using known physics for context only:

| Correction | Magnitude | Sign | Source |
|------------|-----------|------|--------|
| $a_e \approx \alpha/(2\pi)$ | $\sim 1160$ ppm | $+$ (raises $A_{\mathrm{hf}}$) | Electron anomalous magnetic moment (Schwinger 1948) |
| $(1 + m_e/m_p)^{-3} - 1$ | $\sim -1635$ ppm | $-$ | Reduced-mass / recoil |
| Two-loop QED + Lamb-like + nuclear-structure | $\sim 33$ ppm | $\pm$ | Sub-residual physics (out of HF-2 scope) |

The strict-BF $+531$ ppm and the recoil-included $-1102$ ppm bracket the truth. The framework reproduces both reference points correctly given its inputs; the residual is exactly where standard QED + nuclear physics says it should be.

## 6. Three things to test next (input to Track HF-2)

1. **Derive $a_e$ from GeoVac's vertex correction machinery and reduce the residual to $\sim 33$ ppm.** The framework already has `geovac.qed_anomalous_moment` and the graph-native $F_2(\kappa)$ pipeline (Paper 28 §curvature_coefficients); HF-2 should plug the GeoVac-derived $a_e^{\mathrm{GV}}$ into the Bohr-Fermi formula in place of $g_e = 2$, giving $g_e \to 2(1 + a_e)$, and report the residual against experiment. Falsifier: if the GeoVac $a_e$ at $n_{\max} = 3, 4, 5, 6$ either undershoots or systematically misidentifies the $\alpha/(2\pi)$ Schwinger coefficient, HF-2 returns honest-negative.
2. **Recoil from the Fock projection's reduced-mass conformal factor.** The reduced-mass factor $(1 + m_e/m_p)^{-3}$ for $|\psi(0)|^{2}$ should be derivable from making the Fock projection two-body: the conformal factor at $p_0 = \sqrt{-2 E_{\mathrm{red}}}$ uses the reduced mass instead of $m_e$, and $E_{\mathrm{red}} = E_{m_e} \cdot (1 + m_e/m_p)^{-1}$. Track the conformal factor through to $|\psi(0)|^{2}$; the result should be $(1 + m_e/m_p)^{-3}$ exactly. If GeoVac can produce that factor algebraically, that's a meaningful structural result (the framework correctly handles two-body recoil in its native projection); if it cannot, that's a clean statement of where the framework is missing a piece.
3. **Bound the QCD content of $g_p$ from GeoVac's gauge-truncation scope.** CLAUDE.md memory `bertrand_sm_gauge_truncation` says GeoVac's gauge content is upper-bounded by $U(1) \times SU(2) \times SU(3)$. The proton $g$-factor lives in $SU(3)$ confinement; the framework cannot derive $g_p$, but it can in principle predict the *parametric form* of the deviation from the Dirac-point value $g_p = 2$. A scoping memo on whether the GeoVac $SU(3)$ Wilson construction (Paper 30 + Sprint ST-SU3) reproduces the right structural form $g_p = 2 \cdot \text{(QCD strong-form factor)}$ would be a valuable structural-skeleton result, even if the numerical value remains a calibration input.

## 7. Honest blockers / limitations

- **No graph-native $|\psi(0)|^{2}$.** As discussed in §4, the wavefunction-at-origin amplitude is a continuum-embedding quantity. This is not a bug; it is an explicit instance of GeoVac's structural-skeleton scope (Layer 1 + Layer 2 projections, Paper 34). It does mean that any physical observable involving a contact density inherits the continuum embedding's metric.
- **No autonomous derivation of $\alpha, g_p, m_e/m_p$.** All three are CODATA inputs. Paper 2 in Observations is the framework's open conjecture for $\alpha$; $g_p$ is QCD; $m_e/m_p$ is presumably bound by deeper structure not currently in the framework.
- **Reduced-mass factor not derived.** The strict-BF result without recoil is $+531$ ppm above experiment; the recoil-corrected is $-1102$ ppm below. The recoil correction is a known one-line standard-physics result, but it is not currently produced by the framework's machinery (HF-2 idea 2 above is a way to address this).
- **Convention warning logged.** The briefing's $(g_p/2)$ in the prefactor reads literally as a factor-of-2 error; we used the canonical $g_p$ form. The PI may want to revise the briefing language for HF-2.

## 8. Files

- `debug/sprint_hf_track1.py` — computation (this sprint).
- `debug/data/sprint_hf_track1.json` — structured outputs.
- `debug/sprint_hf_track1_memo.md` — this memo.

No production GeoVac code was modified.
