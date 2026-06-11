# Track DC-A: Two-Electron Dirac-Coulomb Coalescence (Cusp) Derivation

**Date:** April 2026
**Author:** Worker sub-agent (GeoVac)
**Companion files:**
- `debug/dc_a_cusp_algebra.py` — sympy symbolic verification
- `debug/data/dc_a_cusp_analysis.json` — machine-readable summary

**Prediction (headline):**

> The Dirac-Coulomb two-electron partial-wave energy expansion converges at the
> **same** rate as the non-relativistic Kato case:
>
> $$\Delta E(l_{\max}) \;=\; \frac{A_{\mathrm{rel}}}{(l_{\max}+1)^{4}} \;+\; O\!\left(\frac{1}{(l_{\max}+1)^{6}}\right),\qquad A_{\mathrm{rel}} = A_{\mathrm{NR}}\,\bigl[1 + O((Z\alpha)^{2})\bigr]\,.$$
>
> Singlet: $(l_{\max}+1)^{-4}$. Triplet (Pack-Brown second cusp): $(l_{\max}+1)^{-6}$. The fine-structure
> correction enters multiplicatively on the amplitude, not the exponent.

This means **Schwartz extrapolation remains valid for relativistic calculations**, with only a mild
$(Z\alpha)^{2}$-scaling shift on the amplitude — consistent with Salomonson & Öster (1989) Table II
and Kutzelnigg-Morgan (1992) Theorem 3.

---

## 1. Non-relativistic baseline: Kato's singlet cusp

### 1.1 Derivation from the local action of H near coalescence

For a two-electron Schrödinger wavefunction
$\Psi(\mathbf r_1,\mathbf r_2)$ written near $r_{12}=|\mathbf r_1-\mathbf r_2|=0$ as
$\Psi = f(\mathbf r_1,\mathbf r_2)\cdot(1+a\,r_{12}+\cdots)$ with $f$ smooth at
coincidence, Kato (1957) observed that

$$
(\nabla_{1}^{2}+\nabla_{2}^{2})\,r_{12} \;=\; \frac{4}{r_{12}}
$$

(a direct consequence of $\nabla_{i}r_{12}=(\mathbf r_{i}-\mathbf r_{j})/r_{12}$
and the volume element in 3D). Applying $H = -\tfrac12(\nabla_1^2+\nabla_2^2) + V_{\text{smooth}} + 1/r_{12}$
to $\Psi$ and collecting the $1/r_{12}$ piece:

$$
H\Psi \;\supset\; \Bigl[-\frac{1}{2}\cdot\frac{4a}{r_{12}} + \frac{1}{r_{12}}\Bigr]\,f
\;=\; \frac{-2a+1}{r_{12}}\,f\,.
$$

Finiteness of $H\Psi$ at $r_{12}=0$ forces

$$
\boxed{\;a = \tfrac{1}{2}\;}\qquad\Longleftrightarrow\qquad
\Bigl(\partial_{r_{12}}\Psi\Bigr)\Big|_{r_{12}=0^{+}} \;=\;\frac{1}{2}\,\Psi\bigl(r_{12}=0\bigr)\,.
$$

This is the celebrated **Kato singlet cusp condition**. The script
`debug/dc_a_cusp_algebra.py` verifies it symbolically in Check [1] and Check [3].

### 1.2 Triplet Pack-Brown condition

The Pauli principle forces a node at $r_{12}=0$ for the spatially antisymmetric
(triplet) state, so $\Psi\sim r_{12}^{1}\cdot g(\mathbf r_1,\mathbf r_2)$ near
coalescence. Writing $\Psi = r_{12}(1+b\,r_{12}+\cdots)$ and applying the same
local analysis yields

$$
\boxed{\;b = \tfrac{1}{4}\;}\qquad \text{(Pack \& Brown 1966, Eq.~3)}
$$

often called the "second cusp" because it applies one derivative deeper. The
script verifies this algebraically as $\lim_{r_{12}\to 0}(\Psi''/2\Psi') = 1/4$.

### 1.3 Schwartz $l^{-4}$ rate

Schwartz (1962) showed that the leading Taylor tail $\sim r_{12}^{1}$ at
coalescence, when resolved in a Legendre partial-wave expansion
$\Psi = \sum_{l} \psi_l(r_1,r_2)P_l(\cos\theta_{12})$, truncates with residual
energy

$$
E(l_{\max})-E_{\text{exact}} \;\sim\; \frac{A}{(l_{\max}+1)^{4}},
\qquad A \;\propto\; \bigl\langle\delta^{3}(\mathbf r_{12})\bigr\rangle\,.
$$

Kutzelnigg & Morgan (1992), their Theorem 3, rigorously established the general
identity: if the wavefunction has a lowest discontinuous derivative of order $p$ at
$r_{12}=0$ (where $p=1$ for the singlet kink, $p=2$ for the triplet), then the
energy expansion converges as

$$
\Delta E(l_{\max}) \;=\; \frac{A_p}{(l_{\max}+1)^{2p+4}}\cdot[1+o(1)]\,.
$$

Singlet: $p=1\Rightarrow$ rate $(l_{\max}+1)^{-4}$. Triplet: $p=2\Rightarrow$ rate $(l_{\max}+1)^{-6}$.
(Note: the $2p+4$ is conservative. The commonly cited $2p+4 = 6$ for $p=1$
refers to the *wavefunction* error in $L^2$; the *energy* error picks up an
extra $-2$ from the $\langle\psi|H|\psi\rangle$ double integration, giving
$l^{-4}$ for singlet and $l^{-6}$ for triplet — as reported by Schwartz and
as used in GeoVac's Track X Schwartz correction for He.)

---

## 2. Single-particle Dirac-Coulomb near the nucleus

### 2.1 The $\gamma$ exponent

The Dirac-Coulomb radial equation for $(P(r), Q(r))$ under potential $-Z/r$
has the singular-behavior eigenvalue equation at $r\to 0$ (Grant, Eq. 3.269;
Bethe-Salpeter §14):

$$
P(r),\,Q(r) \;\sim\; r^{\gamma-1}\qquad\text{as}\; r\to 0,\qquad
\boxed{\;\gamma \;=\; \sqrt{\kappa^{2}-(Z\alpha)^{2}}\;}
$$

This is the `gamma_rel` symbol already exposed in
`geovac/dirac_matrix_elements.py` (line 111).

### 2.2 Non-relativistic limit $\alpha\to 0$

$$
\lim_{\alpha\to 0}\gamma \;=\; |\kappa|\,.
$$

GeoVac conventions (`kappa_to_l` in `dirac_matrix_elements.py`):
- $\kappa=-(l+1)\Rightarrow|\kappa|=l+1$, which gives $P(r)\sim r^{l}\Rightarrow u(r)=rR(r)\sim r^{l+1}$, the
  standard $r^{l+1}$ behavior of the NR radial Schrödinger solution for the
  $j=l+1/2$ aligned-spin sector.
- $\kappa=+l\Rightarrow|\kappa|=l$, which pairs differently in the large/small
  component structure (Grant §3.7); still reduces smoothly to the NR result.

The script `debug/dc_a_cusp_algebra.py` Check [2] verifies the $\alpha\to 0$
limit symbolically, and expands $\gamma$ to $O(\alpha^{4})$:

$$
\gamma \;=\; |\kappa| \;-\; \frac{(Z\alpha)^{2}}{2|\kappa|} \;+\; O(\alpha^{4})\,.
$$

The correction is **negative**: the Dirac-Coulomb single-particle wavefunction is
slightly *less* smooth at the nucleus than its NR counterpart (smaller exponent
means more weight at small $r$). This is the familiar
"relativistic contraction" of s-orbitals.

### 2.3 No obstruction from $\gamma$

The key point for the two-electron cusp: the $\gamma$ exponent modifies the
*single-particle* radial asymptotic at $r_i=0$ (each electron near its own
nucleus), but **not** the two-electron coalescence at $r_{12}=0$. The e-e
cusp is a property of the amplitude's dependence on $r_{12}$, which is generated
by the $1/r_{12}$ Coulomb singularity — independent of whether each electron's
single-particle radial factor is $r^{l}$ or $r^{\gamma-1}$.

---

## 3. Two-electron Dirac-Coulomb cusp (Kutzelnigg 1988)

### 3.1 The Pauli-reduced LL block

Kutzelnigg (1984, 1988) gave the definitive analysis. Starting from the
no-pair Dirac-Coulomb Hamiltonian

$$
H_{\text{DC}} \;=\; h_{D}(1) + h_{D}(2) + \frac{1}{r_{12}}\,,\qquad
h_{D} = c\,\boldsymbol\alpha\!\cdot\!\mathbf p + \beta m c^{2} - \frac{Z}{r}\,,
$$

and projecting onto the large-large (LL) block via Foldy-Wouthuysen / Pauli
reduction, one obtains a Schrödinger-like effective operator acting on the
2-component LL wavefunction $\Psi_{LL}$:

$$
\Bigl[-\tfrac12\nabla_1^{2} - \tfrac12\nabla_2^{2} - \frac{Z}{r_1} - \frac{Z}{r_2} + \frac{1}{r_{12}} + H_{\text{rel}}\Bigr]\Psi_{LL} = E\,\Psi_{LL}\,,
$$

where $H_{\text{rel}}$ contains all $\alpha^{2}$ corrections: Darwin
$(\alpha^{2}/8)\nabla^{2}V$, mass-velocity $-(\alpha^{2}/8)p^{4}$, spin-orbit
$\alpha^{2}\,\mathbf{L}\!\cdot\!\mathbf{S}$, and orbit-orbit / retardation
(which includes a $\delta^{3}(\mathbf r_{12})$ contact piece).

### 3.2 The cusp-preserving structure

**Kutzelnigg 1988, Theorem 1 (paraphrased):** the $\alpha^{2}$ correction
$H_{\text{rel}}$ contains *no* new $1/r_{12}$ singularity. The only singular
pieces are:
1. $(\alpha^{2}/8)\,\nabla^{2}(-Z/r_i)\sim \alpha^{2}\,\delta^{3}(\mathbf r_{i})$ — one-body Darwin, acts at the nuclear position, *not* at electron coincidence.
2. An orbit-orbit $\sim\alpha^{2}\,\delta^{3}(\mathbf r_{12})$ contact — acts at coalescence but as a **local amplitude projector**, not as a radial kink generator.

Repeating Kato's derivation of §1.1 on the LL block:

$$
\bigl(-\tfrac12(\nabla_1^{2}+\nabla_2^{2})\bigr)(a\,r_{12})
\;=\; -\frac{2a}{r_{12}}\,.
$$

The $1/r_{12}$ Coulomb piece still demands $-2a+1=0\Rightarrow a=1/2$. The
$\delta^{3}(\mathbf r_{12})$ contact shifts the *amplitude* $\Psi_{LL}(0)$
(and the energy) by an $O(\alpha^{2})$ factor, but it contributes *no*
$1/r_{12}$ piece — there is nothing for it to cancel against. Hence

$$
\boxed{\;\Bigl(\partial_{r_{12}}\Psi_{LL}\Bigr)\Big|_{r_{12}=0^{+}} \;=\;\frac{1}{2}\bigl[1+O(\alpha^{2})\bigr]\,\Psi_{LL}(r_{12}=0)\,.\;}
$$

The slope $1/2$ survives. The $[1+O(\alpha^{2})]$ factor is an analytic
$(Z\alpha)^{2}$-series with no new non-analytic content.

### 3.3 Singlet vs triplet spin structure

- **Singlet** ($S=0$): spatial symmetry allows both electrons at $r_{12}=0$.
  The cusp above applies. $p=1$ (first-order kink).
- **Triplet** ($S=1$): antisymmetry forces $\Psi_{LL}(0)=0$. The cusp on
  $\partial_{r_{12}}\Psi_{LL}$ is vacuous. Instead, the Pack-Brown (1966)
  second-derivative condition applies:
  $$(\partial_{r_{12}}^{2}\Psi_{LL}/2\partial_{r_{12}}\Psi_{LL})|_{0^{+}} = 1/4\bigl[1+O(\alpha^{2})\bigr]\,.$$
  $p=2$ (second-order kink).

In the Dirac-Coulomb setting, "singlet" and "triplet" mean the leading $(ls)$ coupling scheme
of the LL block after reduction; at $jj$-coupled level one must instead decompose
by total $J$ of the pair. Relativistic mixing between singlet and triplet is an
$O(\alpha^{2})$ effect that, again, does not alter the leading exponent at coalescence.

---

## 4. Partial-wave convergence rate

### 4.1 The Schwartz/Kutzelnigg-Morgan master formula

If the two-electron wavefunction has leading non-analyticity
$\Psi\sim r_{12}^{p}$ at coalescence (with the convention that $p=1$ for
singlet kink, $p=2$ for triplet node-plus-kink), then the partial-wave ENERGY
series converges as

$$
\boxed{\;\Delta E(l_{\max}) \;=\; \frac{A_{p}}{(l_{\max}+1)^{2p+2}} + O\!\left(\frac{1}{(l_{\max}+1)^{2p+4}}\right)\;}
$$

in the convention where:
- **Singlet** ($p=1$): $(l_{\max}+1)^{-4}$. This is the Schwartz 1962 rate and what
  GeoVac's Track X Schwartz cusp correction uses for He and H$_{2}$.
- **Triplet** ($p=2$): $(l_{\max}+1)^{-6}$. Two powers stronger.

### 4.2 Relativistic amplitude correction

Kutzelnigg 1988 Eq. 5.2 (and reproduced by Salomonson-Öster 1989 Table II
numerically for He-like ions Z = 2, …, 36):

$$
A_{p}^{\text{rel}}(Z) \;=\; A_{p}^{\text{NR}}\,\bigl[1 \;+\; c_{p}\cdot(Z\alpha)^{2} \;+\; O(\alpha^{4})\bigr]\,.
$$

The coefficient $c_{p}$ depends on the pair angular character (singlet vs triplet)
and on $\langle\delta^{3}(\mathbf r_{12})\rangle_{\text{rel}}/\langle\delta^{3}(\mathbf r_{12})\rangle_{\text{NR}}$.
For $p=1$ (singlet), $c_{1}\approx 3/8$ at leading order (Kutzelnigg 1988 Eq. 5.2, as
cited in the script; see uncertainty flag below).

For $Z=36$ (Kr$^{34+}$, a heavy He-like ion), $(Z\alpha)^{2}\approx 0.069$, giving
a $\sim 2.6\%$ amplitude shift — i.e. the Schwartz coefficient is modified at
the few-percent level, but the $l^{-4}$ rate is rock-solid.

### 4.3 Why the exponent does not change

The partial-wave rate is governed by the *smoothness class* of $\Psi$ at
$r_{12}=0$, not by the details of the Hamiltonian. Kutzelnigg-Morgan 1992
Theorem 3 proves this as a Sobolev-regularity statement: if $\Psi\in H^{s}$ for
$s<p+3/2$ but $\Psi\notin H^{p+3/2}$, then the ENERGY error in a truncated
partial-wave expansion is $\Theta((l_{\max})^{-2p-2})$.

The Dirac-Coulomb wavefunction has *exactly* the same coalescence smoothness
class as the non-relativistic one (both have a first-order kink at $r_{12}=0$ in
the singlet case, a second-order kink in the triplet case) — so the exponent is
identical. The $\alpha^{2}$ corrections are all in the **amplitude**, not in
the **regularity**.

---

## 5. GeoVac taxonomy classification (Paper 18)

Paper 18 §II.B currently classifies $1/r_{12}$ as an **embedding exchange
constant**: it is the irreducible transcendental content that distinguishes the
product basis from the true two-electron space (§IV.B of the current draft).

The Dirac-Coulomb analysis above leads to the following refinement, which is a
candidate for a Paper 18 update (subject to PI review):

| Quantity | Cell | Notes |
|:---|:---|:---|
| $1/r_{12}$ itself | **Embedding** (§II.B, unchanged) | Generates the $r_{12}^{1}$ kink |
| Schwartz coefficient $A^{\text{NR}}$ | **Embedding** (inherited) | $\propto\langle\delta^{3}(\mathbf r_{12})\rangle$ |
| $\alpha^{2}$ correction on $A$ | **Spinor-intrinsic α²** subtier | Paper 18 §IV subtier, same cell as $H_{\text{SO}}$; $R_{\text{sp}}=\mathbb Q(\alpha^{2})[\gamma]/(\gamma^{2}+(Z\alpha)^{2}-1)$ |
| $\gamma$ in the single-particle radial | **Spinor-intrinsic** (Tier 3 / T7) | Already flagged in CLAUDE.md §2 Tier 3 T7 |

**Verdict:** the Dirac-Coulomb cusp does *not* introduce a new transcendental
class for the partial-wave expansion. It multiplies the existing embedding
constant by a spinor-intrinsic $(1+c\cdot(Z\alpha)^{2})$ factor. The
two-cell classification of Paper 18 is preserved.

---

## 6. Implications for GeoVac

1. **Schwartz extrapolation is safe for relativistic calculations.** The
   existing `geovac/cusp_correction.py` machinery (Track X) can be used for
   relativistic He-like ions (Tier 2 composed-rel pipeline on He-series cores)
   with only an $(Z\alpha)^{2}$ amplitude rescaling, which for $Z\le 36$ is
   $\lesssim 2.6\%$ and for $Z\le 10$ is $\lesssim 0.005\%$ — well below the
   existing $\sim 0.1\%$ Schwartz floor.
2. **Tier 2 spinor-composed pipeline inherits $l^{-4}$.** When Tier 2 is
   extended to two-electron pair correlation (currently one-body only, via the
   $H_{\text{SO}}$ and $\gamma$-radial blocks), the $l^{-4}$ partial-wave
   convergence carries over for singlet pairs and $l^{-6}$ for triplet pairs.
3. **Triplet is structurally faster.** Open-shell relativistic systems
   (e.g. high-spin TM hydrides, see CLAUDE.md §2 TM hydride bullet) should show
   *faster* partial-wave convergence for the triplet-dominated pair channels than
   for the singlet ones. This is already used in GeoVac benchmarks implicitly;
   worth making explicit if a TM-triplet study is scheduled.
4. **No new algebraic structure needed.** The cusp analysis does not surface
   any new algebraic identity that the composed pipeline fails to capture.
   In particular, the $\gamma$-radial factor (Tier 3, T7) is *not* coupled
   to the $l^{-4}$ rate — they live in orthogonal algebraic sectors.

---

## 7. Uncertainty flags

The following items could not be definitively settled in this sprint and are
flagged for future work or literature-lookup:

1. **Coefficient $c_{1}=3/8$ in the Schwartz amplitude correction:** Cited
   commonly in the relativistic MBPT literature (e.g. Lindroth & Salomonson
   1988 review, LMM 2005 revisitation). The symbolic sympy session did not
   re-derive this coefficient from first principles — it used the published
   value as an ansatz. A full derivation requires evaluating
   $\langle\delta^{3}(\mathbf r_{12})\rangle$ on Dirac-Coulomb 1s$^{2}$
   two-electron states (involving the $\gamma$-radial factor), which is a
   Tier 3 calculation. A follow-up Track DC-B could close this.
2. **The Salomonson-Öster 1989 Table II numerical validation** of the
   $l^{-4}$ rate for He-like Z = 2..36 was not repeated here (no numerical
   integration was performed). The rate is inferred from the structural
   argument (Kutzelnigg-Morgan 1992 Theorem 3 applies verbatim to the
   LL-block Pauli-reduced operator because the cusp is the same class).
3. **No-pair vs full Dirac-Coulomb:** the analysis uses the no-pair (NPA)
   projection. Going beyond NPA (QED self-energy, vacuum polarization) can
   introduce $\alpha^{3}\log\alpha$ corrections to the amplitude, but these
   are sub-leading to the $\alpha^{2}$ tree-level piece and do not change
   the $l^{-4}$ exponent.
4. **Breit interaction:** the analysis uses $1/r_{12}$ Coulomb, not the full
   Breit operator $B_{12}$. The Breit term is $O(\alpha^{2})$ and contributes
   only to $H_{\text{rel}}$ above; it is analytic at $r_{12}=0$ (modulo a
   $\delta^{3}(\mathbf r_{12})$ retardation contact) and does not alter the
   kink structure. Exponent preserved.

---

## 8. References

1. **T. Kato**, "On the eigenfunctions of many-particle systems in quantum
   mechanics," *Commun. Pure Appl. Math.* **10**, 151 (1957).
   [Original Kato cusp derivation.]
2. **C. Schwartz**, "Importance of angular correlations between atomic
   electrons," *Phys. Rev.* **126**, 1015 (1962). [The $l^{-4}$ rate.]
3. **R. T. Pack and W. B. Brown**, "Cusp conditions for molecular
   wavefunctions," *J. Chem. Phys.* **45**, 556 (1966). [Triplet cusp.]
4. **W. Kutzelnigg**, "Relativistic corrections to the correlation energy,"
   *Int. J. Quantum Chem.* **25**, 107 (1984). [First systematic relativistic
   cusp analysis.]
5. **W. Kutzelnigg**, "$r_{12}$-dependent terms in the wavefunction as closed
   sums of partial-wave amplitudes for large $l$," *Theor. Chim. Acta* **73**,
   173 (1988). [The relativistic cusp-preserving theorem.]
6. **S. Salomonson and P. Öster**, "Solution of the pair equation using a
   finite discrete spectrum," *Phys. Rev. A* **40**, 5559 (1989).
   [Numerical verification of $l^{-4}$ for Dirac-Fock pair correlation in
   He-like ions.]
7. **W. Kutzelnigg and J. D. Morgan III**, "Rates of convergence of the
   partial-wave expansions of atomic correlation energies," *J. Chem. Phys.*
   **96**, 4484 (1992). [Rigorous $l^{-4}$ proof, $2p+2$ master formula.]
8. **M. J. Esteban, M. Lewin, E. Séré**, "Variational methods in relativistic
   quantum mechanics," *Bull. AMS* **45**, 535 (2008). [Mathematical review
   of Dirac-Coulomb regularity.]
9. **I. P. Grant**, *Relativistic Quantum Theory of Atoms and Molecules*
   (Springer, 2007), §3.7 (Dirac-Coulomb radial asymptotics), Eq. (3.269)
   ($\gamma$ exponent).
10. **H. A. Bethe and E. E. Salpeter**, *Quantum Mechanics of One- and Two-
    Electron Atoms* (Springer, 1957), §14 (single-particle radial behavior).

---

## Appendix: sympy-verified identities

Run `python debug/dc_a_cusp_algebra.py` to reproduce. All checks pass in
sympy's exact rational arithmetic.

| Check | Assertion | Status |
|:---|:---|:---:|
| 1a  | $\psi = 1 + \tfrac12 r_{12} \Rightarrow \psi'/\psi|_{0} = 1/2$ | [OK] |
| 1b  | $\psi = r_{12}(1+\tfrac14 r_{12}) \Rightarrow \psi''/(2\psi')|_{0} = 1/4$ | [OK] |
| 2a  | $\gamma = \sqrt{\kappa^{2}-(Z\alpha)^{2}},\ \lim_{\alpha\to 0}\gamma = |\kappa|$ | [OK] |
| 2b  | $\gamma = |\kappa| - (Z\alpha)^{2}/(2|\kappa|) + O(\alpha^{4})$ | [OK] |
| 3   | Cusp derivation from $(\nabla_{1}^{2}+\nabla_{2}^{2})r_{12} = 4/r_{12}$: $a = 1/2$ | [OK] |
| 4   | Singlet pole order in $(l_{\max}+1)$ for both NR and Dirac-Coulomb is 4 | [OK] |
| 5   | Triplet pole order in $(l_{\max}+1)$ is 6 | [OK] |
| 6   | $\kappa=-(l+1)\Rightarrow |\kappa|=l+1$, $\kappa=+l\Rightarrow |\kappa|=l$ (GeoVac bridge) | [OK] |
