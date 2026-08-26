# Coulomb-Sturmian & Generalized-Sturmian Secular Equations — Equations Sheet

**Purpose:** exact classical formulation of J. Avery's generalized Coulomb-Sturmian method, for
implementation + validation of a quantum algorithm. Every equation is tagged with a **verified
primary source** (fetched and read, not inferred from an abstract). Extraction date: 2026-08-18.

> Notation warning — three symbols for one object. The literature uses different letters for the
> single common Sturmian scale: **k** (Herbst–Avery–Dreuw; Avery theses §10 / book ch. 6),
> **p₀** (Aquilanti/Coletti group), and **pκ** (Avery, many-electron atomic case). They are the
> same quantity: the isoenergetic scale with **E = −½ k² = −½ p₀² = −½ pκ²**. A cross-walk is in §6.

---

## 0. Sources (each item tagged with the tag used below)

| Tag | Reference | Access | Used for |
|-----|-----------|--------|----------|
| **[HAD]** | M. F. Herbst, J. E. Avery, A. Dreuw, *"Quantum chemistry with Coulomb Sturmians: construction and convergence of Coulomb Sturmian basis sets at Hartree–Fock level."* arXiv:1811.05777 (2018); publ. Mol. Phys. | **OPEN** — read via ar5iv HTML | Item 1 (one-electron CS def, normalization, k↔Z/n, E) |
| **[CAL]** | D. Calderini, S. Cavalli, C. Coletti, G. Grossi, V. Aquilanti, *"Hydrogenoid orbitals revisited: from Slater orbitals to Coulomb Sturmians."* **J. Chem. Sci. 124(1), 187–192 (2012).** | **OPEN** — full PDF read (ias.ac.in) | Item 1 (Sturmian ODE eq 6; potential-weighted orthonormality eq 10) |
| **[BK6]** | J. E. Avery & J. S. Avery, *"The Generalized Sturmian Method"*, World Scientific book **chapter 6** (typeset 2010-07-28), reproduced **as an appendix** in [PhD] below. | **OPEN** — read via [PhD] appendix | Item 2 (secular eq 6.5–6.35), Item 5 (Tables 6.2/6.3, He −2.90250) |
| **[PhD]** | J. E. Avery, PhD thesis *"New Computational Methods in the Quantum Theory of Nano-Structures"*, Univ. Copenhagen (DIKU), 2011. | **OPEN** — `hjemmesider.diku.dk/~avery/thesis.pdf` | Item 1 (§10.1), Item 3 (Shibuya–Wulfman §10.4–10.5) |
| **[MSc]** | J. E. Avery, MSc thesis *"The Generalized Sturmian Method: Development, Implementation and Applications in Atomic Physics"*, Univ. Copenhagen, 2008. | **OPEN** — `hjemmesider.diku.dk/~avery/speciale.pdf` | Item 4 (interelectron integrals §3.2.2), Item 5 (§4.1 benchmarks) |
| [SW65] | T. Shibuya & C. E. Wulfman, *"Molecular orbitals in momentum space"*, 1965. | origin ref only (not re-read) | Item 3 provenance |
| [AA07] | J. E. Avery & J. S. Avery, *"Generalized Sturmians and Atomic Spectra"*, World Scientific (2006/2007). | **PAYWALLED** (Google-Books snippets only) | its content reached indirectly via [BK6]/[PhD]/[MSc] |
| [AA04] | J. Avery & J. Avery, *"Generalized Sturmian solutions for many-particle Schrödinger equations"*, J. Phys. Chem. A **108**, 8848 (2004). | **PAYWALLED** (abstract only) | canonical secular-eq paper; content = [BK6] |

`[MSc]` cites the manual as `[Av08]` and the monograph as `[A06]/[AA07]`; equation numbers of the
form (3.2.x)/(6.x)/(10.x.y) below are the theses'/book-chapter's own numbering.

---

## 1. One-electron Coulomb Sturmian — isoenergetic relation, normalization, weighted orthonormality

**Radial function** (position space) — [PhD eq 10.1.2] / [HAD eq 4–5] / [CAL eq 7]:

    R_nl(r) = N_nl (2kr)^l e^{-kr} · F(l+1−n | 2l+2 | 2kr)               [confluent hypergeometric]
            = N_nl (2kr)^l e^{-kr} · L^{2l+1}_{n−l−1}(2kr)               [equivalent Laguerre form, HAD eq 4]

**Potential-weighted normalization constant** — [PhD eq 10.1.3] / [HAD eq 5]:

    N_nl = ( 2 k^{3/2} / (2l+1)! ) · sqrt[ (l+n)! / ( n (n−l−1)! ) ]

**Governing (Sturm–Liouville / "conjugate eigenvalue") equation** — [PhD eq 10.1.5] / [BK6 eq 6.5]:

    [ −½∇²  −  n k / r  +  ½ k² ] φ_nlm(x) = 0

i.e. the hydrogenlike equation with the replacement **Z/n → k**, held at fixed energy.

**Isoenergetic / "k = Z/n" quantization** — [CAL, text after eq 4; BK6 text after eq 6.9]:
every basis member shares one k, hence one energy

    ε = −½ k²          (all n,l,m in the set)                            [PhD eq 10.1.7; BK6 6.5]

The Sturmian *becomes the physical hydrogenlike orbital* exactly when **k → Z/n**; equivalently the
"weight/charge eigenvalue" is β = n·p₀ = n·k, varied to satisfy the boundary condition
([CAL eq 6], where E₀ = −½p₀² is the *fixed* energy and β_nl = n p₀ is the eigenvalue of the charge).
[HAD] state the same as **β_n = k n / Z** with **E = −k²/2** ([HAD eq 2, 6]).

**POTENTIAL-WEIGHTED ORTHONORMALITY** (the defining property; weight = the bare-Coulomb 1/r) —
[PhD eq 10.1.6] = [BK6 eq 6.7] = [CAL eq 10]:

    ∫ d³x  φ*_{n'l'm'}(x) · (1/r) · φ_nlm(x)  =  (k / n) · δ_{n'n} δ_{l'l} δ_{m'm}

Equivalent unit-normalized ("Shibuya–Wulfman overlap") form — [BK6 eq 6.8]:

    ∫ d³x  φ*_{μ'}(x) · ( (−∇² + k²) / (2k²) ) · φ_μ(x)  =  δ_{μ'μ}

Table of the first R_nl is [BK6 Table 6.1] / [PhD Table 9.7]: e.g. R_10=2k^{3/2}e^{−kr},
R_20=2k^{3/2}(1−kr)e^{−kr}, R_21=(2/√3)k^{3/2} kr e^{−kr}.

---

## 2. Many-electron GENERALIZED STURMIAN secular equation (the full statement)

Source: **[BK6] §6.1.2–6.1.5** (Avery & Avery book ch. 6), cross-checked against [MSc]. This is the
"Goscinskian" atomic case — V₀ = bare-nucleus attraction.

**N-electron Hamiltonian pieces** — [BK6 eq 6.22–6.24, 6.26]:

    Ĥ = −½ Σ_j ∇²_j + V(x),      V(x) = V₀(x) + V′(x)
    V₀(x) = − Σ_{j=1}^N Z / r_j                 (nuclear attraction, the "weighted"/zeroth-order potential)
    V′(x) = Σ_{j>i} Σ_i 1 / r_ij                 (interelectron repulsion)
    x ≡ (x₁,…,x_N) including spin.

**Generalized-Sturmian basis** = isoenergetic solutions of the *weighted* wave equation
(Goscinski 1968) — [BK6 eq 6.9 / 6.25]:

    [ −½ Σ_j ∇²_j + β_ν V₀(x) − E_κ ] Φ_ν(x) = 0,     E_κ = −½ pκ²,   pκ ≡ √(−2E_κ)

with the weight β_ν chosen so **all** Φ_ν share the one energy E_κ (isoenergetic).

**Goscinskian configurations** = Slater determinants of hydrogenlike spin-orbitals with a common
weighted charge Q_ν — [BK6 eq 6.27–6.29]:

    Φ_ν = |φ_{μ1} φ_{μ2} … φ_{μN}|,     φ_{nlm,ms}(x_j) = R_nl(r_j) Y_lm(θ_j,φ_j) χ_{ms}
    weighted charge:   Q_ν = Z β_ν = pκ / R_ν
    configuration "root":   R_ν = sqrt( 1/n₁² + 1/n₂² + … + 1/n_N² )
    energy:            E_κ = −½ pκ² = −½ Q_ν² R_ν² = −( Q_ν²/2n₁² + … + Q_ν²/2n_N² )      [BK6 6.30]

**Matrix definitions** — [BK6 eq 6.13]:

    T⁰_{ν'ν} ≡ −(1/pκ) ∫ dτ  Φ*_ν'(x) V₀(x) Φ_ν(x)        (nuclear-attraction / "kinetic-derived" matrix)
    T′_{ν'ν} ≡ −(1/pκ) ∫ dτ  Φ*_ν'(x) V′(x) Φ_ν(x)        (interelectron-repulsion matrix)

**Nuclear-attraction matrix is DIAGONAL and p-independent** — [BK6 eq 6.34] (= [MSc] ⟨V₀⟩=−ZR·p):

    T⁰_{ν'ν} = Z R_ν · δ_{ν'ν}

**THE GENERALIZED STURMIAN SECULAR EQUATION** — [BK6 eq 6.19], atomic form [BK6 eq 6.35]:

    Σ_ν  [  T⁰_{ν'ν}  +  T′_{ν'ν}  −  pκ δ_{ν'ν}  ]  B_ν  =  0
    ⇔   Σ_ν  [  Z R_ν δ_{ν'ν}  +  T′_{ν'ν}  −  pκ δ_{ν'ν}  ]  B_ν  =  0        (Goscinskian, atoms)

- **Eigenvalue = the scaling parameter pκ** (NOT an energy). Energies recovered by
  **E_κ = −pκ² / 2**  [BK6 eq 6.21]. (So yes: eigenvalue = pκ = √(−2E), and E = −pκ²/2.)
- **Eigenvector B_ν** = configuration-mixing (CI) coefficients; Φ_κ = Σ_ν Φ_ν B_ν [BK6 6.20].
- **T′ is the "energy-independent interelectron-repulsion matrix"**: its elements are **pure numbers
  that depend only on the electron number N** — independent of pκ, of E, and of Z — so a single T′
  serves an entire isoelectronic series and all its states [BK6 eq 6.35 + boxed remark].
  (Mechanism: the two-electron integrals scale ∝ Q_ν = pκ/R_ν; the −1/pκ prefactor cancels the pκ.)

**The six structural features** (verbatim gist, [BK6] after eq 6.21):
(1) kinetic-energy term has vanished; (2) the V₀ (nuclear) matrix is diagonal; (3) roots are not
energies but scaling parameters pκ ∝ √(binding energy); (4) only the *shapes* of Φ_ν are known
before solving; (5) the solution yields a near-optimal basis + the states + energies at once;
(6) "the Hamiltonian formalism is nowhere to be seen." The Slater exponents Q_ν=pκ/R_ν that drop
out are automatically near-optimal (no preliminary HF, no exponent guessing).

> General (non-atomic) form: for an arbitrary V₀, eq (6.19) `[T⁰+T′−pκ 1]B=0` still holds; only the
> convenient facts "T⁰ diagonal = ZR_ν" and "T′ = pure numbers" are special to the bare-nucleus
> (Goscinskian) choice.

---

## 3. The Shibuya–Wulfman matrix (one-electron, many-center)

Source: **[PhD] §10.4–10.5**. Role: it is the **overlap-type ("kinetic-derived") matrix** of the
*isoenergetic one-electron* secular problem when Coulomb Sturmians are placed on several centers
(the molecular / many-center analogue of §2's T⁰). Originates with [SW65].

**One-electron molecular orbital**, expanded in many-center CS — [PhD eq 10.4.1–10.4.6]:

    [ −½∇² + v(x) − ε_j ] ψ(x) = 0,     v(x) = − Σ_a Z_a / |x − x_a|,   ε_j = −½ k²
    ψ(x) = Σ_μ φ_μ(x) C_μ,    φ_μ(x) := φ_{nlm}(x − x_a),    μ := (n, l, m, a)   [a = center index]

**Shibuya–Wulfman integral / matrix** S — [PhD eq 10.4.8]:

    S_{μ',μ}  :=  (1/k²) ∫ d³x  φ*_{μ'}(x) [ −½∇² + ½ k² ] φ_μ(x)

**Wulfman (potential) integral** W — [PhD eq 10.4.9]:

    W_{μ',μ}  :=  −(1/k) ∫ d³x  φ*_{μ'}(x) v(x) φ_μ(x)

**Secular equation** (roots are k, not energies) — [PhD eq 10.4.10]:

    Σ_μ  [ W_{μ',μ}  −  k · S_{μ',μ} ]  C_μ  =  0,     then  ε_j = −½ k²      [PhD text after 10.4.10]

**Momentum-space (Fock-projection) evaluation** — [PhD §10.5]. With the Fock map of momentum p to
the unit 4-sphere u (eq 10.5.11) and the transformed CS  t̄_μ(p)=M(p) Y_{n−1,l,m}(u), M(p)=4k^{5/2}/(k²+p²)²
(eq 10.1.10–10.1.11, 10.5.6):

    Sturmian overlap:      m_{μ',μ} = ∫ d³p e^{ip·R} M(p)² Y*_{μ'}(u) Y_μ(u)               [PhD 10.5.12]
    Shibuya–Wulfman:       S_{μ',μ} = ∫ d³p e^{ip·R} ( 2k/(k²+p²) )³ Y*_{μ'}(u) Y_μ(u)     [PhD 10.5.13]
    with R = x_{a'} − x_a  (inter-center displacement).

Both reduce to hyperangular integrals over 4-D hyperspherical harmonics
Y_{λ,l,m}(u) = N_{λ,l} C^{1+l}_{λ−l}(u₄) Y_{l,m}(u₁,u₂,u₃) [PhD 10.5.6, Gegenbauer C]. A closed radial
reduction S_{μ',μ} → f_{n,l}(kR) Y_{lm}(R̂) with a 3-term R_{n,l} recurrence is [PhD eq 10.5.16].

---

## 4. Interelectron-repulsion matrix in the CS basis for a 2-electron atom (He)

Source: **[MSc] §3.1–3.2.2** (Slater–Condon reduction + explicit Gaunt/3j + radial Slater integral).
T′ of §2 is assembled from these two-electron integrals.

**Determinant-level (Slater–Condon for a 2-electron operator)** — [MSc eq 3.1.11–3.1.13]:

    ⟨Φ|Ŵ|Φ⟩ = Σ_{i<j, k<l} (−1)^{i+j+k+l} C_{ij;kl} |S_{ij;kl}|,
    C_{ij;kl} = ⟨φ_μi φ_μj | w₁₂ | φ_μk φ_μl − φ_μl φ_μk⟩,      w₁₂ = 1/|x₁ − x₂|  (Coulomb − exchange)

For **He 1s²** with orthonormal orbitals this collapses (rule 3.1.13) to the single Coulomb–exchange
pair term C_{ij;ij}.

**Radial/angular separation of each two-electron integral** — [MSc eq 3.2.18–3.2.22]:

    J = Σ_l a_l I_l                                                                    [3.2.20]
    with the Laplace/Legendre expansion  1/r₁₂ = Σ_l ( r_<^l / r_>^{l+1} ) P_l(x̂₁·x̂₂)  [3.2.19]

**Angular factor a_l = product of Gaunt integrals = products of Wigner 3j** — [MSc eq 3.2.11, 3.2.23–3.2.26]:

    ∫ dΩ Y_{l1 m1} Y_{l2 m2} Y_{l m}
        = sqrt[ (2l1+1)(2l2+1)(2l+1) / 4π ] · ( l1 l2 l ; 0 0 0 ) · ( l1 l2 l ; m1 m2 m )
    a_l ∝ [∫ Y_{l,−mi−mk} Y_{li mi} Y_{lk mk}] · [∫ Y_{l,−mj−ml} Y_{lj mj} Y_{ll ml}]     [3.2.25]

**Selection rules** (kill the infinite l-sum → finite) — [MSc eq 3.2.12–3.2.13, 3.2.26]:

    m_i + m_j = m_k + m_l,     (−1)^{l_i + l_k} = (−1)^{l_j + l_l},     l ∈ {|l1−l2|,…,l1+l2} (step 2)

**Radial factor I_l = the Slater integral R^l** (exponential-type orbitals) — [MSc eq 3.2.22, 3.2.9]:

    I_l = ∫₀^∞ dr₁ r₁² p₁(r₁) e^{−α₁ r₁} ∫₀^∞ dr₂ r₂² p₂(r₂) e^{−α₂ r₂} ( r_<^l / r_>^{l+1} )

(evaluated in closed form via the polynomial-times-exponential identity ∫₀^∞ r^j e^{−αr} dr = Γ(j+1)/α^{j+1},
[MSc eq 3.2.8–3.2.9]; p = product polynomial, α = α₁+α₂.)

**Coupled (LS) configurations → 6j.** When the determinants are angular-momentum-coupled, the products
of 3j's recombine into **Wigner 6j** coefficients (the hyperspherical/"λ-basis" recoupling), per
[CAL §4] ("connected … through a generalized 6j coefficient") and [PhD Tables 10.1–10.2] (integrals =
pure functions of the scaled separation). **Pure-number property**: because I_l ∝ Q = pκ/R and T′ carries
the −1/pκ prefactor, T′_{ν'ν} is pκ/E/Z-independent (§2, [BK6 6.35]).

**Concrete He leading element (standard hydrogenic value, structure per above):** the 1s² diagonal
Coulomb integral for hydrogenlike charge Q is ⟨1s²|1/r₁₂|1s²⟩ = (5/8) Q; with Q = pκ/R and R=√2 for 1s²,
the corresponding **T′ entry is the pure number (5/8)/√2 ≈ 0.4419**, independent of Z — illustrating
[BK6 6.35]. (The (5/8)Q value is the well-known hydrogenic result; only the *pure-number reduction* is
Avery's claim I verified.)

---

## 5. Concrete PUBLISHED numbers to validate against

All from **[BK6]** (Avery & Avery, book ch. 6), which cites exact nonrel. references
Nakashima–Nakatsuji (2008) and NIST/Ralchenko (2008).

**(a) He ground state 1s² ¹S** — [BK6 §6.2, verbatim]:
> "Values calculated using **102 isoenergetic configurations based on Coulomb Sturmians** are given by
> us in [AA07], Table F.1. For the **helium ground state we obtain −2.90250 Hartrees**."
  - Compare exact nonrelativistic He: **−2.90372 Ha** (Nakashima–Nakatsuji). Error ≈ 1.2×10⁻³ Ha.
  - Note [BK6]: the *Goscinskian* basis is deliberately poor for the He **ground** state (best for
    excited states); the CS-isoenergetic basis (102 configs) is the accurate one → use (a) as the
    ground-state validation target.

**(b) He excited ¹S series, 40 Goscinskians** — [BK6 Table 6.2] (Hartrees):
  | state | GS method | Nakatsuji (exact nonrel) |
  |-------|-----------|--------------------------|
  | 1s2s ¹S | −2.1429 | −2.1460 |
  | 1s3s ¹S | −2.0603 (≈, table) | −2.0611 |
  Whole 2-electron isoelectronic table (He…N⁵⁺) computed "in a few milliseconds."

**(c) He excited ³S, 36 Goscinskians** — [BK6 Table 6.3]:
  | state | GS method | Nakatsuji |
  |-------|-----------|-----------|
  | 1s2s ³S | −2.1736 | −2.1752 |

**(d) Two-electron (He-like) isoelectronic series, 98-function Goscinskian basis `N2medium1S`**
— [MSc §4.1, Fig 4.2]: ground-state relative error vs NIST **< 0.5 %**, dropping below the
relativistic-correction magnitude by Z≈8; large-Z few-electron approximation reaches 10⁻⁴–10⁻⁶
relative error [MSc §4.1 text]. (Curves, not a single tabulated Ha value.)

**Recommended primary validation point:** **He 1s² ground state = −2.90250 Ha with 102 Coulomb-Sturmian
isoenergetic configurations** [BK6 §6.2 / AA07 Table F.1]; secondary: 1s2s ¹S = −2.1429 Ha (40 Goscinskians)
[BK6 Table 6.2].

---

## 6. Notation cross-walk (implementers read this)

| Concept | [HAD]/[PhD]/[BK6] | [CAL] (Aquilanti) | [BK6] many-electron |
|---|---|---|---|
| common Sturmian scale | k | p₀ | pκ |
| energy | E = −½k² | E₀ = −½p₀² | E_κ = −½pκ² |
| charge/weight eigenvalue | β_n = kn/Z | β_nl = n p₀ | Q_ν=Zβ_ν=pκ/R_ν |
| becomes physical orbital when | k = Z/n | p₀ = λ = Z/n | Q_ν = Z (single config) |
| secular eigenvalue | k (1-electron, §3) | — | pκ (N-electron, §2) |
| weighted orthonormality weight | 1/r (⇒ k/n) | 1/r (⇒ p₀/n) | V₀ (⇒ −pκ² diag) |

**One-line implementation summary.** Fix E<0 ⇒ scale s ∈ {k, pκ}=√(−2E). Build the two matrices
(nuclear T⁰: diagonal ZR_ν for atoms, or Wulfman W for molecules; interelectron T′: pure-number
3j/6j·Slater-integral array). Solve the **linear (generalized) eigenproblem for the scale s**, not for
E: `[T⁰ + T′ − s·1]B = 0` (atoms) or `[W − s·S]C = 0` (molecules). Read back E = −s²/2. The kinetic
operator never appears; T′ (the only nontrivial matrix) is energy- and Z-independent.

---

## 7. Verification status (honest)

| Item | Status | Basis |
|---|---|---|
| 1. one-electron CS (def, N_nl, k↔Z/n, E=−½k², **weighted orthonormality ∫φ*(1/r)φ=(k/n)δ**) | **PRIMARY-VERIFIED (×3)** | [HAD] eq 1–6 (ar5iv), [CAL] eq 3–10 (open PDF), [PhD] eq 10.1.2–10.1.7 |
| 2. N-electron generalized-Sturmian secular eq `[ZR_ν+T′−pκ]B=0`, eigenvalue pκ, E=−pκ²/2, T′ pure numbers | **PRIMARY-VERIFIED** | [BK6] eq 6.9–6.35 (Avery book ch.6, appended in [PhD]); cross-checked [MSc] |
| 3. Shibuya–Wulfman matrix `S_{μ'μ}=(1/k²)∫φ*(−½∇²+½k²)φ`, `[W−kS]C=0`, momentum-space form | **PRIMARY-VERIFIED** | [PhD] eq 10.4.1–10.4.10, 10.5.1–10.5.16 |
| 4. He interelectron matrix = Slater–Condon → Gaunt/3j (·6j coupled) × radial Slater integral I_l | **PRIMARY-VERIFIED** (structure) | [MSc] eq 3.1.11–3.2.26; (5/8)Q leading value is standard hydrogenic, only its pure-number reduction is [BK6 6.35] |
| 5. He ground state **−2.90250 Ha / 102 CS configs**; 1s2s ¹S −2.1429/40 Gosc.; series <0.5% | **PRIMARY-VERIFIED** | [BK6] §6.2, Tables 6.2/6.3; [MSc] §4.1 |

**Paywall-blocked (abstract-only), NOT used as equation source:** the monograph [AA07] "Generalized
Sturmians and Atomic Spectra" (World Scientific 2006) and the JPCA paper [AA04] — but their content is
the very [BK6]/[PhD]/[MSc] material read above, so no gap remains for items 1–5. The only number that
lives *only* behind a paywall is the full 102-config Table F.1 of [AA07]; its **He entry (−2.90250)** is
quoted openly in [BK6] and is captured here. [SW65] (original Shibuya–Wulfman paper) was not re-read;
its definition is taken from [PhD] §10.4, which attributes it.

**Un-derived / flagged:** exact Gaunt→6j recombination coefficients for specific coupled He configs are
stated structurally (per [CAL]/[PhD]); the explicit 6j tables live in [AA07] Ch.3–4 / Appendix A
(paywalled) — reproduce them from the Slater–Condon+3j primitives in [MSc §3.2] rather than trusting a
secondary transcription.
