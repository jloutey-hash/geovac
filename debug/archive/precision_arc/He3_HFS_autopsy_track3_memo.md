# Sprint memo: Helium-3 2³S₁ hyperfine — first multi-electron HFS autopsy

**Date:** 2026-05-18
**Track:** Multi-track Roothaan-autopsy sprint, Track 3
**Sibling tracks:** D HFS (Track 5), H 21cm (v1), μH Lamb (v1), He 2³P (v1), Ps 1S-2S (Track 4)
**Status:** CLOSED — five-component decomposition + two structural findings

---

## 1. Task and reference

Decompose the helium-3 2³S₁ metastable HFS splitting at I=1/2, J=1
into a five-component framework-native autopsy.

**Reference (experimental):** $\nu_\text{HFS}(\text{He-3},\, 2{}^{3}\!S_1) = 6739.701177(16)$ MHz
(Schuessler–Fortson–Dehmelt PR 187, 5 (1969); Prior–Wang PRA 16, 2071 (1977);
Rosner–Pipkin PRA 1, 571 (1970)). The state is metastable
(t_lifetime ~ 8000 s); the splitting is between F=3/2 and F=1/2 sub-levels.

**Architectural novelty:** This is the **first multi-electron HFS** test in
the Paper 34 §V.C catalogue. The previous five autopsies (H 21cm, μH Lamb,
He 2³P, D HFS, Ps 1S-2S) are all single-electron-active or atomic-fine-
structure observables; the He-3 2³S₁ HFS forces the framework to compose
a multi-electron contact density at the nucleus.

---

## 2. Five-component framework-native architecture

| # | Component | Projection chain | Output (MHz) | Residual (ppm) |
|---|-----------|------------------|--------------|----------------|
| 1 | Bohr–Fermi 2³S₁ (multi-electron) | §III.1 ∘ §III.7 ∘ §III.8 ∘ §III.22 ∘ §III.20 | 6650.31 | −13260 |
| 2 | + Schwinger $a_e$ | §III.6 (Parker-Toms) | 6658.03 | −12118 |
| 3 | + Reduced-mass recoil ($m_p \to m_{He3}$) | §III.14 | 6654.40 | −12657 |
| 4 | + §III.18 Zemach ($r_Z = 1.965$ fm, $Z=2$ scaling) | §III.18 | 6653.41 | −12804 |
| 5 | (§III.20 PK orthogonality diagnostic) | §III.20 | (10.5% spread across paths) | DOMINANT |

**Headline:** framework-native cumulative residual is **−12804 ppm = −1.28%**
at the Slater-screened contact-density prescription. The dominant Layer-2
content is §III.20 multi-electron screening (10.5% spread) — *not* multi-loop
QED (LS-8a wall) as for H/D/μH HFS.

---

## 3. Architectural content (the substantive new findings)

### 3.1 First-time content: multi-electron contact density at the nucleus

The 2³S₁ state has electronic configuration (1s, 2s) with triplet spin
coupling (S=1, M_S ∈ {-1, 0, 1}). The spatial wavefunction is the
antisymmetric Slater determinant:

$$ \Psi_T(r_1, r_2) = \frac{1}{\sqrt{2}}\,[\varphi_{1s}(r_1)\varphi_{2s}(r_2) - \varphi_{2s}(r_1)\varphi_{1s}(r_2)] $$

The Fermi-contact Hamiltonian acts via $\sum_i \delta^3(r_i)\,\mathbf{s}_i\!\cdot\!\mathbf{I}$.
Evaluating in the |J=S=1, M_S=1⟩ state:

$$ \Big\langle \sum_i \delta^3(r_i)\,s_{i,z} \Big\rangle_{M_S=1} = \frac{|\varphi_{1s}(0)|^2 + |\varphi_{2s}(0)|^2}{2} $$

The factor 1/2 is a **structural multi-electron projection** that does
NOT appear in single-electron HFS calculations. Mechanism: the
antisymmetric Slater determinant has $\langle\delta^3(r_i)\rangle =
(|\varphi_{1s}|^2 + |\varphi_{2s}|^2)/2$ for each electron, and the
spin operator $\langle s_{i,z}\rangle = 1/2$ for both electrons in the
M_S=1 state. Summing two electrons gives total/2 (each contributes
density/2 × 1/2 × 2 electrons = density/2).

This is a new architectural ingredient relative to the single-electron
HFS template used in Sprint MH (μH HFS, Mu HFS, H/D HFS, etc.). The
generalization to multi-electron HFS in atomic precision catalogues is
**a §III.22 bipolar-harmonic content** (the |J=S, M_S⟩ basis projection
of the per-electron spin-density operator).

### 3.2 Convention finding: BF formula mass slot is m_p, not m_N

A surprising load-bearing finding surfaced during diagnostic. The
`bohr_fermi_a_constant` formula's mass slot uses m_e/m_p **always**
(because the nuclear magneton is referenced to the proton mass μ_N =
e\hbar/(2m_p), independent of the actual nucleus). Track 5 (D HFS)
used m_d in this slot with a doubled g-factor (g_d_atomic = 2·μ_d/I)
which happens to give the correct splitting by a *factor cancellation*
specific to I=1 (the m_d/m_p factor cancels the 2/2 g-factor doubling).

For I=1/2 nuclei (He-3, H, Mu) the cancellation does NOT hold. The
standard convention (g_N = μ/I, m_p in mass slot) is required. The
He-3+ 1s hydrogen-like ion sanity check confirms this: standard
convention matches experiment at **+0.05%** (predicted A = -8661.25
MHz vs experimental -8665.65 MHz), while Track 5's convention is off
by 3/2.

**This is a §V.D candidate convention exposure:** D HFS Track 5
convention is system-dependent and cannot be transferred directly to
He-3 (or any I=1/2 multi-electron HFS). The standard convention
generalizes correctly across all I-values.

### 3.3 §III.18 magnetization-density operator at I=1/2 He-3

Reproduces Eides analytic leading-order Zemach to **machine precision**
(0.0 ppm absolute, 0.0% of LO shift) at r_Z(He-3) = 1.965 fm. Profile
(Gaussian vs exponential) independence preserved at machine precision.
The §III.18 module's Z=1 basis output is multiplied by Z=2 (He nucleus
charge) to give the He-3 leading-order Zemach shift.

This is the second I=1/2-nucleus operator-level §III.18 test
(after H 21cm v1) and the first multi-electron one. The operator-level
machinery extends without modification.

### 3.4 §III.20 Phillips-Kleinman / core-valence orthogonality detectability

**THE HEADLINE FINDING.** The 2³S₁ contact density depends sensitively
on the screening prescription for the outer 2s electron:

| Screening prescription | Z_eff(1s) | Z_eff(2s) | ν_HFS (MHz) | Residual (ppm) |
|------------------------|-----------|-----------|-------------|----------------|
| Hydrogenic (unscreened) | 2.00 | 2.00 | 7311.34 | +84817 |
| Slater (Z-σ_1s = 2-0.85) | 2.00 | 1.15 | 6653.41 | −12804 |
| Full screen (Z_eff,2s = 1) | 2.00 | 1.00 | 6600.52 | −20651 |

**Spread across paths: 105468 ppm = 10.5%.** This is the FIRST
precision-catalogue entry where the framework's multi-electron screening
prescription is the dominant Layer-2 systematic. For H/D HFS the
dominant wall is LS-8a multi-loop QED (~10-300 ppm); for He-3 2³S₁
HFS, the §III.20 PK convention dominates (10.5% spread, well above
both LS-8a and the experimental precision).

**§III.20 PK orthogonality is structurally detectable at He-3 2³S₁
HFS precision.** The framework's composed-architecture treatment with
an actual §III.20 PK ab initio calculation at finite n_max would
provide a unique answer; the Slater Z_eff=1.15 prescription is a
sprint-scope placeholder.

---

## 4. Architectural decisions

### 4.1 Contact-density prescription used

We adopted **Slater's-rule-screened** Z_eff for the 2s electron
(Z_eff(2s) = Z - σ_{1s} = 2 - 0.85 = 1.15). This is the most defensible
cleanly-specified prescription at sprint scope. Alternative paths
documented:

- **Hydrogenic unscreened** (Z_eff = 2 for both): physically incorrect
  for the outer 2s (overshoot residual by +85000 ppm); included as
  diagnostic to demonstrate §III.20 detectability.
- **Full screen** (Z_eff(2s) = 1): aggressive screening; modest
  improvement vs Slater but no clear physical justification.
- **Graph-native CI** at n_max=2: gives identical contact density to
  hydrogenic-Slater because the 2³S₁ state is the unique (1s,2s)
  triplet configuration in the n_max=2 graph-native ss-sector
  (verified: triplet ss-block at n_max=2 is 1×1, eigenvalue
  E(2³S) = -2.124143 Ha at 2.35% vs Drake).
- **Hylleraas r12** (geovac/hylleraas_r12.py): the named follow-on
  closure path; explicitly NOT used in this autopsy because the
  single-alpha Hylleraas basis is known non-variational for the
  2³S₁ state (+209% at omega=5 per §V.C.4 closure path). The
  Hylleraas-Eckart double-alpha extension is the named structural
  follow-on (~2-4 week sprint).

### 4.2 Five-component framework-native subtotal at Slater screening

ν_HFS = 6653.41 MHz, residual −12804 ppm = **−1.28%**.

For comparison, the existing precision-catalogue entries with
multi-loop-QED-bounded residuals:
- H 21 cm: +18 ppm (Eides Tab 7.3 bounded)
- D 1S HFS: +286 ppm (Pachucki–Yerokhin Layer-2 bounded)
- μH 1S HFS: +2 ppm (Sprint MH B closure)
- Mu 1S HFS: +199 ppm (cleanest LS-8a empirical anchor)

The He-3 2³S₁ residual of −1.28% is 50–500× larger than these,
**because the dominant systematic is §III.20 screening, not LS-8a
multi-loop QED**. This is the first precision-catalogue entry where
the framework's multi-electron Layer-2 prescription is load-bearing.

---

## 5. I·J multiplet structure at I=1/2, J=1

Operator-level I·J Hamiltonian construction reproduces the F=3/2 vs
F=1/2 splitting structure exactly:

- He-3 (I=1/2, J=1) eigenvalues of I·J: {−1.0, +0.5}, splitting 3/2
- D (I=1, J=1/2) eigenvalues of I·S: {−1.0, +0.5}, splitting 3/2
- Multiplicity ratio He-3/D = 1.000000000000 (bit-exact)

**I↔J swap symmetry verified at machine precision.** This is a structural
consequence of Wigner 3j symmetry: the F-multiplet decomposition of
$\mathbf{I}\otimes\mathbf{J}$ at total angular momentum F is symmetric
under I↔J swap when the total multiplet structure (F=3/2, F=1/2) is
the same. Pauli encoding identical structure: I=1/2 nuclear (1 qubit) +
J=1 electronic (2 qubits for the (2J+1)=3-dim Hilbert) for He-3, vs
I=1 nuclear (2 qubits) + J=1/2 (1 qubit) for D — same total Hilbert
dimension, same Pauli decomposition.

---

## 6. Layer-2 attribution

Cumulative residual −12804 ppm (Slater path) attributes as:

| Wall | Magnitude (ppm) | Tier |
|------|-----------------|------|
| §III.20 multi-electron screening (1s screening of 2s) | ±50000 spread | §III.20 PK |
| Multi-electron QED Pachucki-Drachman | few hundred ppm | LS-8a multi-particle |
| He-3 nuclear polarizability (α-particle-like binding) | few ppm | W3 inner-factor |
| Multi-loop QED $\alpha^2(Z\alpha)$ | tens of ppm | LS-8a renormalization |

**The dominant wall is §III.20 multi-electron screening**, NOT LS-8a
multi-loop QED. This restructures the catalogue's understanding of
where the structural skeleton's residuals sit for multi-electron HFS:
single-electron HFS (H, D, μH, Mu) has LS-8a as the dominant wall;
multi-electron HFS (He-3 2³S₁, future Li-7, Cs-133) has §III.20 PK as
the dominant wall.

---

## 7. Honest scope and named follow-ons

### 7.1 Limitations of this autopsy

- The contact density extraction at Slater-screened Z_eff is a
  sprint-scope placeholder. Ab initio §III.20 PK at finite n_max
  with the actual He 2s radial wavefunction would give the unique
  framework-preferred value. **Named follow-on:** wire `screened_psi_origin_squared`
  (`geovac/neon_core.py`) to compute |φ_2s(0)|² for He without a
  frozen core (Z=2, no [Ne] or other registered core); test as
  §III.20 ab initio path.

- Multi-loop QED corrections (Pachucki-Drachman) are tens of ppm
  and structurally outside the framework's one-loop scope (LS-8a wall).

- He-3 nuclear polarizability is small (~few ppm, alpha-particle-like
  binding stronger than deuteron's n+p) and sits in W3 inner-factor.

- Hylleraas-Eckart double-alpha extension is the named structural
  closure to access the proper many-body contact density (currently
  blocked at single-alpha Hylleraas by the §V.C.4 negative result).

### 7.2 §V.D candidate convention exposures from this autopsy

1. **BF formula mass slot convention**: D HFS Track 5 uses m_d (cancels
   for I=1 only); standard convention uses m_p. **§V.D.7 candidate:**
   document this as a class-(i) inter-track / inter-system itemization
   convention.

2. **§III.20 PK screening prescription (hydrogenic-vs-Slater-vs-ab-initio)**:
   first instance in the catalogue where the framework's own screening
   convention is the dominant Layer-2 systematic. **§V.D.8 candidate:**
   document as a framework-internal class-(ii) convention (similar to
   Friar RMS-vs-first-moment §V.D.4) — different prescriptions yield
   structurally distinct Layer-2 attributions.

### 7.3 Coverage matrix update

Adding He-3 2³S₁ to the multi-focal precision catalogue (now 9 systems):

| System | Mass hierarchy | Nuclear spin | Multi-electron | Residual |
|--------|---------------|--------------|----------------|----------|
| H 21cm | $m_e \ll m_p$ | I=1/2 | single-e | +18 ppm |
| D HFS | $m_e \ll m_d$ | I=1 | single-e | +286 ppm |
| μH HFS | $m_\mu \sim m_e \times 207$ | I=1/2 | single-e | +2 ppm |
| Mu HFS | $m_e \ll m_\mu$ | I=1/2 | single-e | +199 ppm |
| Ps HFS | equal mass | I=1/2 | single-e | +4900 ppm |
| He 2³P | $m_e \ll m_\alpha$ | I=0 (4He) | multi-e | −0.014/−0.20% |
| Mu 1S-2S | $m_e \ll m_\mu$ | I=1/2 | single-e | −0.11 ppm |
| μH Lamb | overlap | I=1/2 | single-e | +0.013% |
| **He-3 2³S₁ HFS** | $m_e \ll m_{He3}$ | I=1/2 | **multi-e + HFS** | **−1.28%** |

He-3 2³S₁ HFS adds the **multi-electron + HFS** intersection point.
The framework-native residual is 50–500× larger than single-electron
HFS systems, attributable cleanly to §III.20 PK rather than LS-8a.

---

## 8. Architecture used (per CONSTRAINT clause)

- Framework architecture for multi-electron density: **hydrogenic-Slater**
  (Z_eff = 2 - 0.85 for the 2s electron), the most defensible cleanly-
  specified prescription at sprint scope.
- §III.18 magnetization-density module: used at the operator level
  (geovac/magnetization_density.py) for the Zemach component — no
  modification.
- Hylleraas r12 module: NOT used (single-alpha known non-variational
  for 2³S₁ per §V.C.4 closure path).

### 8.1 Is §III.20 PK orthogonality structurally detectable?

**YES.** The 10.5% spread across screening prescriptions is well above
both the ppm-scale experimental precision (16 Hz / 6.7 GHz = 2.4 ppb)
and the existing precision-catalogue framework residuals (tens to
hundreds of ppm for single-electron HFS). §III.20 PK orthogonality is
structurally detectable and is the **first §V.D candidate exposure**
to surface from a multi-electron HFS autopsy.

---

## 9. Files

- `debug/He3_HFS_autopsy_track3.py` — five-component decomposition driver
- `debug/He3_HFS_autopsy_track3_memo.md` — this memo
- `debug/data/He3_HFS_autopsy_track3.json` — structured results

## 10. Verdict

- Framework-native residual: **−1.28%** at Slater-screened contact density.
- §III.18 operator-level at I=1/2 He-3: **machine-precision reproduction**
  of Eides analytic leading-order Zemach.
- §III.20 PK orthogonality: **structurally detectable at 10.5% spread**.
  **First precision-catalogue entry where multi-electron screening
  prescription is the dominant Layer-2 wall.**
- New architectural finding: **multi-electron triplet contact density
  carries a structural factor 1/2 vs single-electron, from |J=S=1,M_S=1⟩
  spin projection**. This is §III.22 bipolar-harmonic content for HFS
  observables.
- New convention exposure (§V.D candidate): **BF formula mass slot
  is m_p, not m_N**; Track 5 (D HFS) convention is system-dependent.
- Named follow-on: §III.20 ab initio screening for the 2s via
  `screened_psi_origin_squared` (geovac/neon_core.py, no-frozen-core
  Z=2 case) — would give the framework-preferred unique value within
  the spread.
