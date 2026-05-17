# Sprint L1-C — Four-witness physical-setup specification on $\mathcal{T}_{n_{\max}}$

**Date:** 2026-05-16
**Sprint:** L1-C (third of three parallel scoping tracks for Sprint L1 modular-Hamiltonian closure)
**Sibling tracks:** L1-A (modular-Hamiltonian architecture); L1-B (operator-system infrastructure / reusability)
**Builds on:** Sprint TD Track 1 (`geovac/thermal_tensor_triple.py`); Sprint TD Track 4 (`debug/sprint_td_track4_memo.md`); Track D (`debug/bisognano_wichmann_track_d_memo.md`); Unruh pendant (`debug/unruh_pendant_memo.md`); Sprint L0 audit (`debug/lorentzian_l0_audit_memo.md`); Paper 32 §VIII Rem `rem:bisognano_wichmann_reading` + §VIII.E; Paper 34 §III.27; Paper 38 (WH1 PROVEN, GH-convergence on $S^3$).
**Status:** Scoping memo. **No production code or papers modified.**
**Scope:** Per-witness physical-setup specification — what each face's GeoVac mapping must reproduce at the operator-system level, with explicit ambiguities named.

---

## §1. Four-witness theorem recap

The Hawking–Sewell–Bisognano-Wichmann–Unruh theorem (codified in
`debug/unruh_pendant_memo.md` §4) is one Wick-rotation statement with four
physical instantiations: a state restricted to a region $\mathcal{R}$ bounded
by a bifurcate Killing horizon (Schwarzschild exterior, Rindler wedge, etc.)
with surface gravity $\kappa_g$ is KMS at inverse temperature
$\beta = 2\pi / \kappa_g$ with respect to the Killing-time evolution. The
single $2\pi$ across all four faces is $\mathrm{Vol}(S^1)$ — the regularity
period of the Euclidean rotation of the observer's causal-domain coordinates.
On the GeoVac side, this $2\pi$ is the M1 Hopf-base measure signature of the
master Mellin engine (Paper 32 §VIII `rem:master_mellin_domain`), instantiated
on a temporal $S^1$. Sprint L1's goal is to lift the structural correspondence
to literal operator-system identification: compute the modular Hamiltonian
$K_{n_{\max}}$ on the truncated operator system $\mathcal{O}_{n_{\max}}$ for
a wedge-restricted Camporesi–Higuchi state and verify
$\sigma_{i \cdot 2\pi}(O) = O$ at finite $n_{\max}$.

---

## §2. Per-witness specification

### §2.1 Hartle–Hawking

| Field | Specification |
|:------|:--------------|
| **Lorentzian setup** | Eternal Schwarzschild exterior region $\mathcal{R} = \{r > 2M\}$ of mass-$M$ black hole; static external observer at fixed Schwarzschild $r$. State: Hartle–Hawking thermal vacuum on Kruskal-extended spacetime, restricted to $\mathcal{R}$. |
| **Surface gravity $\kappa_g$** | $\kappa_g = 1/(4M)$ (closed form). |
| **Inverse temperature $\beta$** | $\beta_{\mathrm{HH}} = 2\pi / \kappa_g = 8\pi M$. |
| **Continuum $K$** | $K_{\mathrm{HH}} = \beta_{\mathrm{HH}} \cdot H_{\xi}$ where $H_{\xi}$ generates the Killing-time evolution along the static timelike Killing vector $\xi^\mu = (\partial_t)^\mu$ in the exterior. Modular flow is the surface-gravity-rescaled Killing flow (Sewell 1982). |
| **Wedge / restriction** | Static exterior $r > 2M$, bounded by the bifurcate horizon $S^2 \times \{r = 2M\}$. |
| **KMS state $\omega_\beta$** | $\omega_{\mathrm{HH}}(O) = Z^{-1} \mathrm{Tr}(e^{-\beta H_\xi} O)$ on the wedge algebra. |
| **GeoVac mapping (a) $P_{\text{wedge}}$** | **Candidate H-α (recommended):** the wedge is an "exterior" of a hypothetical horizon $S^2 \subset S^3$ at some fixed polar angle $\theta_0$; $P_{\text{wedge}}$ projects $\mathcal{H}_{n_{\max}}$ onto modes with support on $\theta > \theta_0$. *Ambiguity*: $\theta_0$ is not forced by GeoVac data — it parameterises the "horizon size." Natural choice: $\theta_0 = \pi/2$, i.e. the horizon is the equator of $S^3$, giving equal-volume wedges. **Candidate H-β:** drop horizon embedding entirely and define $P_{\text{wedge}}$ via the spectral cut of a Killing-like generator on $T_{S^3} \otimes T_{S^1_\beta}$ (Sprint TD Track 1's tensor construction). This avoids the embedded-$S^2$ ambiguity but at the cost of leaving "wedge" purely algebraic. Recommendation: try H-α first because it is geometrically interpretable; fall back to H-β if the horizon embedding produces no clean modular structure. |
| **GeoVac mapping (b) $K_{n_{\max}}$** | $K_{n_{\max}} = \beta_{\mathrm{HH}} \cdot H_{\xi, n_{\max}}$ where $H_{\xi, n_{\max}}$ is the truncated Killing generator. *Ambiguity*: $S^3 = \mathrm{SU}(2)$ has six Killing vectors (so(4) algebra); choose the one along the polar axis defining the horizon $\theta = \theta_0$ in H-α (so $H_\xi = -i \partial_\phi$ for $\phi$ the azimuthal angle around the polar axis). |
| **GeoVac mapping (c) Truncated KMS state** | $\omega_{n_{\max}}(O) = Z_{n_{\max}}^{-1} \mathrm{Tr}_{P_{\text{wedge}} \mathcal{H}_{n_{\max}}}(e^{-\beta_{\mathrm{HH}} H_{\xi, n_{\max}}} O)$. |
| **GeoVac mapping (d) Predicted $\sigma_{i\cdot \beta}$** | $\sigma_{i \cdot 8\pi M}(O) = e^{2\pi i \cdot H_{\xi, n_{\max}} / \kappa_g} \cdot O \cdot e^{-2\pi i \cdot H_{\xi, n_{\max}} / \kappa_g}$. After unit normalisation $\kappa_g \to 1$ (see §3), predicted equality $\sigma_{i \cdot 2\pi}(O) = O$ to machine precision (bit-exact in exact-arithmetic mode; $\lesssim 10^{-12}$ in floating point) at $n_{\max} \in \{2, 3, 4\}$. |

### §2.2 Sewell

| Field | Specification |
|:------|:--------------|
| **Lorentzian setup** | Globally hyperbolic Lorentzian manifold with a bifurcate Killing horizon (the framework Sewell 1982 generalised from Schwarzschild). For GeoVac concreteness: same Schwarzschild exterior as HH. State: any KMS state with respect to Killing-time evolution at $\beta = 2\pi / \kappa_g$. |
| **Surface gravity $\kappa_g$** | $\kappa_g$ = Killing surface gravity of the bifurcate horizon. For Schwarzschild $= 1/(4M)$. |
| **Inverse temperature $\beta$** | $\beta_{\mathrm{Sew}} = 2\pi/\kappa_g$ (general; Sewell's content is identifying this $\beta$ as a KMS-Tomita-Takesaki theorem, not deriving a new $\beta$). |
| **Continuum $K$** | $K_{\mathrm{Sew}} = 2\pi \cdot M_{\xi}$ where $M_\xi$ is the boost-class generator built from $\xi$. Sewell's theorem: the modular automorphism is the surface-gravity-rescaled Killing flow. |
| **Wedge / restriction** | Static exterior region of the bifurcate horizon (Schwarzschild exterior $r > 2M$ for the concrete case). |
| **KMS state $\omega_\beta$** | Same as HH (the same state object; Sewell's contribution is the *KMS-modular* interpretation of HH's path-integral derivation). |
| **GeoVac mapping** | **Sewell's mapping is HH's mapping at the operator-system level.** Sewell does not add new variables or operators beyond HH; it adds a *modular-theoretic interpretation* of the same construction. The operator-system test for Sewell is therefore the same operator-system test as for HH; the *content* that distinguishes Sewell is the identification of $K_{n_{\max}}$ with the modular Hamiltonian of $P_{\text{wedge}}$, which is established via Tomita's theorem once $K_{n_{\max}}$ has been verified to satisfy $\sigma_{i \cdot \beta}(O) = O$. Sewell is therefore a **derived corollary** at the operator-system level, not an independent setup. |
| **Predicted $\sigma_{i\cdot \beta}$** | Identical to HH. Sewell's verdict is automatic once HH closes. |

### §2.3 Bisognano–Wichmann

| Field | Specification |
|:------|:--------------|
| **Lorentzian setup** | Minkowski spacetime; right Rindler wedge $W_R = \{x > |t|\}$; restriction of the Minkowski vacuum $\Omega$ to the wedge-localised algebra $\mathcal{A}(W_R)$. State: Minkowski vacuum restricted to $W_R$. |
| **Surface gravity $\kappa_g$** | In natural Lorentz units, $\kappa_g = 1$ (the boost rapidity rate). The rapidity rate of the boost subgroup mapping $W_R$ to itself. For an accelerated observer at proper acceleration $a$ this becomes $\kappa_g = a$; the BW theorem is stated in rapidity-units convention. |
| **Inverse temperature $\beta$** | $\beta_{\mathrm{BW}} = 2\pi$ (in rapidity units; for a specific observer at proper acceleration $a$, $\beta = 2\pi/a$ — see Unruh below). |
| **Continuum $K$** | $K_{\mathrm{BW}} = -2\pi \cdot M_{\text{boost}}$ where $M_{\text{boost}}$ is the boost generator preserving $W_R$, i.e. the generator of the one-parameter Lorentz boost subgroup that maps $W_R$ to itself. *Bisognano–Wichmann theorem*: the Tomita modular operator $\Delta$ of $(\mathcal{A}(W_R), \Omega)$ satisfies $\Delta = e^{-K_{\mathrm{BW}}}$. |
| **Wedge / restriction** | Right Rindler wedge of Minkowski. |
| **KMS state $\omega_\beta$** | Minkowski vacuum $\Omega$ restricted to $\mathcal{A}(W_R)$, which is KMS at $\beta = 2\pi$ w.r.t. boost evolution. |
| **GeoVac mapping (a) $P_{\text{wedge}}$** | **Recommended:** half-$S^3$ algebra (one hemisphere under the polar-axis decomposition $S^3 = D^3_+ \cup_{\text{equator}} D^3_-$). $P_{\text{wedge}}$ is the projection onto modes with support on $D^3_+$ (a particular hemisphere). This is the "Riemannian-side analog of a Rindler wedge" already named in Track D §3.1 step 1. *Ambiguity*: which polar axis (which Killing direction) — the answer is gauge-equivalent under SO(4), but a specific choice must be fixed before any computation. Natural choice: align with the polar axis defining the Hopf base $S^2$ in Paper 25 (so the wedge respects the Hopf fibration). |
| **GeoVac mapping (b) $K_{n_{\max}}$** | $K_{n_{\max}} = -2\pi \cdot M_{\text{boost}, n_{\max}}$ where $M_{\text{boost}, n_{\max}}$ is a "boost-class" generator on $\mathcal{H}_{n_{\max}}$. *This is the deep ambiguity*: GeoVac has no actual Lorentz boost on a compact $S^3$, only an analog. Candidates: (BW-α) the rotation generator $J_z$ around the polar axis defining the wedge — this is a *spatial* rotation, not a boost, but it plays the analogous role of "the Killing generator preserving the wedge"; (BW-β) the dilation generator $\rho \partial_\rho$ in a chosen stereographic chart on $S^3$ — this is closer in form to a boost in the conformal compactification reading; (BW-γ) the modular Hamiltonian of $(\mathcal{O}_{n_{\max}, +}, \omega)$ extracted directly from Tomita's $S = J \Delta^{1/2}$ construction without identifying any geometric generator first. Recommendation: BW-γ as the *definitional* choice (let Tomita's theorem define $K_{n_{\max}}$ on the operator system); use BW-α as a sanity check (compare the extracted $K_{n_{\max}}$ to $-2\pi J_z$ and see whether they agree). BW-β is unlikely to add insight given the compact $S^3$ has no natural dilation. |
| **GeoVac mapping (c) Truncated KMS state** | $\omega_{n_{\max}}(O) = \langle \Omega_{n_{\max}} | O | \Omega_{n_{\max}} \rangle$ for $O \in P_+ \mathcal{O}_{n_{\max}} P_+$, where $\Omega_{n_{\max}}$ is the Camporesi–Higuchi ground state on $\mathcal{H}_{n_{\max}}$. |
| **GeoVac mapping (d) Predicted $\sigma_{i\cdot \beta}$** | $\sigma_{i \cdot 2\pi}(O) = e^{2\pi K_{n_{\max}} / (2\pi)} O e^{-2\pi K_{n_{\max}} / (2\pi)} = e^{K_{n_{\max}}} O e^{-K_{n_{\max}}}$. Under Tomita's theorem this equals $O$ identically when $K = -\log \Delta$ for the correctly identified $\Delta = S^* S$. Predicted exact at machine precision (bit-exact in exact arithmetic; $\lesssim 10^{-12}$ floating point) at $n_{\max} \in \{2, 3, 4\}$. *This is the principal Sprint L1 test*. |

### §2.4 Unruh

| Field | Specification |
|:------|:--------------|
| **Lorentzian setup** | Same Rindler wedge as BW, but with a specific observer at fixed proper acceleration $a$. State: Minkowski vacuum restricted to the accelerated observer's causal patch (i.e., $\mathcal{A}(W_R)$ as in BW). |
| **Surface gravity $\kappa_g$** | $\kappa_g = a$ (proper acceleration). |
| **Inverse temperature $\beta$** | $\beta_U = 2\pi / a$. |
| **Continuum $K$** | $K_U = -(2\pi/a) \cdot M_{\text{boost}} = -\beta_U \cdot M_{\text{boost}}$. Same modular Hamiltonian as BW, rescaled by $1/a$ to express in proper-time units of the accelerated observer. |
| **Wedge / restriction** | Same right Rindler wedge as BW. |
| **KMS state $\omega_\beta$** | Same wedge-restricted Minkowski vacuum; thermal at $T_U = a/(2\pi)$ in the accelerated observer's proper-time frame. |
| **GeoVac mapping** | **Identical to BW at the operator-system level**, but parameterised by $a$ instead of by a unit rapidity. $P_{\text{wedge}}$ = same half-$S^3$ projection; $K_{n_{\max}}$ = $-\beta_U \cdot M_{\text{boost}, n_{\max}} = -(2\pi/a) \cdot M_{\text{boost}, n_{\max}}$. The factor of $1/a$ is an overall rescaling that cancels in $\sigma_{i \cdot \beta_U}$: $e^{i \cdot i \beta_U \cdot K_{n_{\max}}/(-\beta_U)} = e^{-K_{n_{\max}}/\beta_U \cdot \beta_U \cdot i \cdot i / (\ldots)}$ — once unit-normalised, falls onto the BW test exactly. |
| **Predicted $\sigma_{i\cdot \beta}$** | Identical to BW after rapidity rescaling. Unruh is the *parameterised* form of the BW test, with $a$ as the variable observer-specific surface gravity. |

---

## §3. Cross-witness unification

The four-witness theorem says these are four faces of one statement; the
operator-system tests should be four parameterised instances of one
underlying construction. Cross-witness unification table:

| Witness | $\kappa_g$ | $\beta$ | $P_{\text{wedge}}$ | $K_{n_{\max}}$ | Test |
|:--------|:-----------|:--------|:-------------------|:---------------|:-----|
| HH | $1/(4M)$ | $8\pi M$ | exterior of horizon $S^2$ at $\theta = \theta_0$ | $\beta_{\mathrm{HH}} \cdot H_\xi$ | $\sigma_{i \cdot 8\pi M}(O) = O$ |
| Sew | $\kappa_{\text{Killing}}$ | $2\pi/\kappa$ | same as HH | $2\pi \cdot M_\xi$ | corollary of HH |
| BW | $1$ (rapidity) | $2\pi$ | hemisphere of $S^3$ | $-2\pi \cdot M_{\text{boost}}$ | $\sigma_{i \cdot 2\pi}(O) = O$ |
| Unruh | $a$ | $2\pi/a$ | same as BW | $-(2\pi/a) \cdot M_{\text{boost}}$ | rescaling of BW |

**Unifying construction.** There is a single underlying object
$\mathfrak{M}(\mathcal{R}) = -\beta \cdot \xi_{\text{Killing}}(\mathcal{R})$
parameterised by the choice of wedge $\mathcal{R}$ and the natural inverse
temperature $\beta = 2\pi / \kappa_g(\mathcal{R})$. The four witnesses are
four choices of $(\mathcal{R}, \xi)$:

- HH/Sew: $\mathcal{R}$ = Schwarzschild exterior, $\xi$ = static-time Killing field
- BW: $\mathcal{R}$ = Rindler wedge, $\xi$ = boost generator (rapidity)
- Unruh: $\mathcal{R}$ = Rindler wedge, $\xi$ = boost generator scaled by $a$ (parameterised BW)

On the GeoVac side, the unification is sharper: under recommended choices
H-α (hemisphere bounded by Hopf-base equator) and BW-γ (Tomita-defined
$K_{n_{\max}}$), HH and BW collapse to the **same operator-system test**
applied to different wedges. The horizon-embedded Schwarzschild test (HH)
reduces to a hemispheric wedge test in the limit where the embedded $S^2$
becomes the Hopf base; BW is exactly this hemispheric wedge test at unit
rapidity; Unruh is the BW test at rapidity $a$; Sewell is the abstract
KMS-modular reading of either.

**Unit normalisation.** The GeoVac graph is dimensionless. The
surface-gravity scalars $\kappa_g = 1/(4M)$, $\kappa_g = a$ have dimensions
of inverse length in Lorentzian theory and need to be reinterpreted. Three
unit-normalisation candidates:

- **U-1 (natural-units):** set $\kappa_g = 1$ throughout (the BW
  rapidity-units convention). All four witnesses then test the same
  $\sigma_{i \cdot 2\pi}(O) = O$ identity. Recommended for Sprint L1: it
  makes the test independent of any specific physical scale and emphasises
  the structural-period-closure nature of the claim.
- **U-Z (Coulomb-anchored):** set $\kappa_g = Z$ to preserve Coulomb-anchored
  GeoVac conventions. Test then reads $\sigma_{i \cdot 2\pi/Z}(O) = O$,
  which is a $Z$-dependent rescaling of U-1's test. Not recommended as
  primary: introduces a Coulomb-specific scale into a structural test that
  should hold independent of GeoVac's physics anchoring.
- **U-obs (observable-matched):** choose $\kappa_g$ to reproduce a specific
  physical observable (e.g., $\kappa_g$ such that $T_H$ matches a chosen
  Schwarzschild mass). Useful for catalogue placement but not for the
  structural period-closure test.

**Recommendation:** Sprint L1 uses U-1 (natural units). The witnesses
collapse to the same period-closure test $\sigma_{i \cdot 2\pi}(O) = O$
modulo wedge choice; the four-witness theorem becomes a single operator-system
identity verifiable at finite $n_{\max}$.

---

## §4. Predicted falsifiers per witness

| Witness | What can fail | Diagnostic |
|:--------|:--------------|:-----------|
| **HH** | (i) The horizon embedding $S^2 \subset S^3$ is unnatural — no GeoVac data selects $\theta_0$; (ii) Killing-time Hamiltonian $H_\xi$ on $\mathcal{H}_{n_{\max}}$ does not commute with the wedge-restriction projector, breaking the modular-flow construction; (iii) the truncation $P_{n_{\max}}$ does not respect the hemispheric structure (a $C^\infty$ function supported on $\theta > \theta_0$ has Fock content extending across all shells). | (i) is the deep ambiguity — recommendation H-α (equatorial horizon $\theta_0 = \pi/2$) makes the construction maximally symmetric but is a *choice*, not a derivation. (ii) Test: compute $[H_\xi, P_{\text{wedge}}]$ on $\mathcal{H}_{n_{\max}}$ at small $n_{\max}$; verify it is small in operator norm relative to $\|H_\xi\|$. (iii) Test: project $f(\theta) = \chi_{\theta > \pi/2}$ onto the Fock basis at $n_{\max}$ and measure the residual outside the wedge. If non-zero at $O(1)$, the wedge is leaky at finite cutoff. |
| **Sew** | Sewell is automatic given HH. The only failure mode is HH's failures. | If HH passes but the operator system $\mathcal{O}_{n_{\max}}$ is not Tomita-cyclic (which would prevent the modular-Hamiltonian identification from going through), Sewell's reading fails even though HH's spectrum-level test passes. Test: verify cyclicity of the truncated KMS state on $\mathcal{O}_{n_{\max, +}}$. |
| **BW** | (i) GeoVac has no actual Lorentz boost on a compact $S^3$, so $M_{\text{boost}}$ is genuinely ambiguous (BW-α/β/γ). (ii) The hemispheric algebra may not be cyclic/separating for the truncated Camporesi–Higuchi ground state at finite $n_{\max}$; Tomita's $S$ operator may not be densely defined. (iii) The period closure $\sigma_{i \cdot 2\pi}(O) = O$ may fail not by an order-1 amount but by an $O(1/n_{\max})$ amount — a "soft" failure indicating asymptotic convergence to the BW period in the GH limit (Paper 38) but not literal closure at finite cutoff. | (i) Compare BW-α (spatial rotation $J_z$), BW-γ (Tomita-defined). If they agree at machine precision, the modular generator is unique. If they disagree by $O(1)$, the BW reading is geometrically ambiguous and the recommendation is to take BW-γ as the definitional choice. (ii) Test: compute $S = $ closure of $O \mapsto O^*$ on $\mathcal{O}_{n_{\max}, +} \Omega_{n_{\max}}$; verify that its polar decomposition exists at $n_{\max} = 2, 3, 4$. (iii) Test: measure $\|\sigma_{i \cdot 2\pi}(O) - O\|_{\text{op}}$ as a function of $n_{\max}$; if it scales as $\gamma_{n_{\max}} \to 0$ from L2 of Paper 38, the literal-identification claim is recovered in the GH limit but not at finite $n_{\max}$. This is the *most likely failure mode* and may be the actual outcome. |
| **Unruh** | Same as BW (Unruh is parameterised BW). Specific failure mode: the $a$-rescaling may not commute with $P_{n_{\max}}$ if the truncation introduces an implicit length scale. | Compute the period closure for at least two values of $a$ (e.g., $a = 1$ and $a = 2$) at the same $n_{\max}$. If the residual scales with $a$, the truncation is introducing a hidden length scale and the period-closure test is breaking the rescaling invariance. |

---

## §5. Sprint L0 / M3-trivialization connection

Sprint L0 predicted (`debug/lorentzian_l0_audit_memo.md` §3 & §4) that the
M3 sub-mechanism of the master Mellin engine (vertex-parity Hurwitz /
Dirichlet-$L$ content from Paper 28's $\zeta(s, 3/4)$ / $\zeta(s, 5/4)$
reductions producing Catalan's $G$ and $\beta(4)$) becomes structurally
empty at signature $(3, 1)$ because the BBB Table 2 entry flips the
$\{J, \gamma_5\}$ anticommutation relation.

**Does this affect the L1 witness specs?** No, at L1's scope.

The modular Hamiltonian construction of Sprint L1 operates on the
*Riemannian* truncated spectral triple $\mathcal{T}_{n_{\max}}$ at
signature $(3, 0)$. The witness specs above are entirely Riemannian: the
"wedge" is a Riemannian hemisphere or polar-cap region of $S^3$; $K_{n_{\max}}$
is built from Riemannian Killing generators or spatial rotations; the
truncated KMS state lives on the Riemannian Camporesi–Higuchi spectrum.
The bridge to the Lorentzian witnesses is the published Wick-rotation
chain (HH 1976 → Sewell 1982 → BW 1976 → Unruh 1976) applied to the
metric-functional output, not a re-derivation in Krein signature.

The Sprint L0 prediction that M3 trivializes at $(3, 1)$ is a Sprint L2
target — it requires actually constructing the Krein-space spectral
triple at $(3, 1)$ per BBB 2018. Sprint L1 does not enter this signature.
The four-witness M1 mechanism ($k = 0$, Hopf-base measure / $\mathrm{Vol}(S^1)$)
is uniformly load-bearing across all four witnesses; M3 trivialization
does not bite because the witnesses' transcendental signature is M1, not
M3.

**Consistency check.** If Sprint L1 succeeds, the lifted operator-system
identification at $(3, 0)$ would say: the framework's M1 mechanism on
$\mathcal{T}_{n_{\max}}$ IS the modular-flow / boost-orbit period of all
four faces simultaneously, verified internally without recourse to the
published Wick-rotation chain. The M3 trivialization at $(3, 1)$ would
then be a *separate* Sprint L2 consistency check on the Krein-lifted
framework: does the Lorentzian extension preserve the M1-only character
of the four-witness theorem, or does it introduce new content via the
flipped $\{J, \gamma_5\}$? The expected answer is "preserves M1 only"
— but verifying it is L2 work, not L1.

---

## §6. Verdict and recommendations for L1 implementation

**Verdict: clean spec; one deep ambiguity per witness, all named.**

The four witnesses collapse to **one operator-system test** at the
recommended unit-normalisation U-1 (natural rapidity units, $\kappa_g = 1$):
verify $\sigma_{i \cdot 2\pi}(O) = O$ for $K_{n_{\max}}$ extracted from
Tomita's $S = J \Delta^{1/2}$ construction on $(\mathcal{O}_{n_{\max}, +}, \omega_{n_{\max}})$
at $n_{\max} \in \{2, 3, 4\}$, where $\mathcal{O}_{n_{\max}, +}$ is the
hemispheric sub-algebra of the truncated operator system.

**Recommended L1 implementation order:**

1. **BW first.** Hemispheric wedge on $S^3$ is the cleanest geometric
   object; BW's mapping is the most parsimonious. If $\sigma_{i \cdot 2\pi}(O) = O$
   passes for BW, all four witnesses lift simultaneously via the cross-witness
   unification of §3.

2. **HH/Sew as parameterised BW.** Reinterpret HH's "exterior of horizon
   $S^2$" as the hemispheric wedge with a specific Killing-axis choice.
   Sewell follows automatically from HH's modular-Hamiltonian extraction.

3. **Unruh as scaled BW.** Trivial parameterisation; the BW test at unit
   rapidity, evaluated at acceleration $a$, gives Unruh directly.

**Recommended falsifier prioritisation:**

- **Primary test:** BW with BW-γ (Tomita-defined $K_{n_{\max}}$). Either
  closes literal identification across all four witnesses or returns an
  $n_{\max}$-decay residual whose scaling identifies the cutoff-recovery
  mode.
- **Sanity check:** compare BW-γ to BW-α (the spatial rotation $J_z$
  around the polar axis). Agreement at machine precision elevates the
  identification from "Tomita-defined" to "geometrically natural."
- **Stress test:** verify Unruh's $a$-rescaling invariance at two values
  of $a$ to confirm the truncation does not introduce a hidden length
  scale.

**Honest scope:** the most likely outcome is the "soft" failure mode
(BW iii in §4): $\|\sigma_{i \cdot 2\pi}(O) - O\|_{\text{op}}$ scales as
the L2 mass-concentration moment $\gamma_{n_{\max}} \to 0$ from Paper 38,
giving literal identification only in the GH limit, not at finite
$n_{\max}$. This would still be a substantial result — it would land
Sprint L1 as a positive-with-cutoff-rate, matching Paper 38 WH1's
qualitative-rate maturity, and would unify the BW reading with the
GH-convergence machinery.

**Three open ambiguities flagged for L1-A architecture track:**

1. **Horizon-axis choice for HH.** Whether to fix $\theta_0 = \pi/2$
   (equatorial, Hopf-base-aligned) or treat $\theta_0$ as a parameter
   over which to scan.
2. **Boost-generator choice for BW.** BW-γ (Tomita-defined) is the
   recommended definitional choice, but BW-α ($J_z$) and BW-β
   (stereographic dilation) deserve at least sanity-check comparison at
   $n_{\max} = 2$ before committing.
3. **Hemispheric truncation correction.** The Fock basis does not
   respect hemispheric support exactly; the projection $P_{\text{wedge}}$
   on $\mathcal{H}_{n_{\max}}$ leaks across the equator at finite
   $n_{\max}$. The leakage scales with $n_{\max}$ in a way that needs
   to be quantified before the modular construction can be trusted.

**Files referenced:** `debug/unruh_pendant_memo.md`,
`debug/bisognano_wichmann_track_d_memo.md`,
`debug/lorentzian_l0_audit_memo.md`,
`papers/synthesis/paper_32_spectral_triple.tex` (§VIII
`rem:bisognano_wichmann_reading`, §VIII.E),
`papers/observations/paper_34_projection_taxonomy.tex` (§III.27),
`geovac/thermal_tensor_triple.py` (Sprint TD Track 1 apparatus),
`geovac/operator_system.py` (Connes-vS truncated operator system),
`geovac/connes_distance.py` (existing SDP framework, extensible to
modular flow).

**Literature anchors (re-verified this sprint):**

1. Hartle, J.B. and Hawking, S.W. "Path-integral derivation of black-hole
   radiance." *Phys. Rev. D* **13**, 2188–2203 (1976). [APS link confirmed.]
2. Bisognano, J.J. and Wichmann, E.H. "On the duality condition for
   quantum fields." *J. Math. Phys.* **17**, 303–321 (1976). [ADS / NASA
   confirmed; nLab Bisognano-Wichmann theorem entry cross-references.]
3. Sewell, G.L. "Quantum fields on manifolds: PCT and gravitationally
   induced thermal states." *Annals of Physics* **141**, 201–224 (1982).
   [Verified in Track D bibliography.]
4. Unruh, W.G. "Notes on black-hole evaporation." *Phys. Rev. D* **14**,
   870–892 (1976). [Verified in Unruh-pendant bibliography.]

No phantom citations.
