# Track SM-E: Consistency check — is Δ = 1/40 electromagnetic/Coulomb-specific?

**Verdict: PASS.** Δ is absent from both Paper 24 (HO/Bargmann–Segal) and Paper 23
(nuclear shell), and the existing universal/Coulomb-specific partition in
Papers 22 and 24 is structurally consistent with Δ being electromagnetic.

## 1. Paper 24 (Bargmann–Segal lattice for 3D HO on S⁵)

Theorem 1 (Bargmann graph rationality): every diagonal and every adjacency weight
of the HO graph at any N_max is an **exact rational number**. The verification
script at N_max=5 reports:
- `all_diagonal_rational: True`
- `all_adjacency_rational: True`
- `irrationals_encountered: []`

Theorem 3 (Structural origin of calibration π, dual of Fock rigidity):
the HO's first-order complex-analytic projection (linear spectrum
E_N = ℏω(N + 3/2), linear in N) **requires no calibration exchange constant**.
The Coulomb/HO asymmetry table (Sec. on asymmetry) explicitly lists:

| Feature              | Coulomb (S³)     | HO (S⁵ Hardy)    |
|----------------------|------------------|------------------|
| Projection constant  | p₀ = Z/n         | ℏω (universal)   |
| Transcendentals      | π (calibration)  | **none**         |

Paper 24 Sec. 7.1 explicitly states that the π in K = π(B + F − Δ)
"is structurally Coulomb-specific: it comes from the spectral zeta values
of S³ … there is no analog in any other dynamical-symmetry-classified
central potential." There is nowhere in the HO framework that could host a
Δ-like rational correction of transcendental origin — the entire HO graph
is already bit-exactly rational, so there is nothing to "correct away."

Grep for `1/40`, `Delta`, `vacuum polarization`, `fine structure`, `calibration
constant` in paper_24 returns only: (a) HO graph selection-rule symbols
ΔN, Δl, Δm_l (unrelated to Δ = 1/40), (b) explicit references to Paper 2's
combination B + F − Δ in Sec. 7.1 where it is classified as Coulomb-specific,
(c) framing statements that calibration constants are tied to second-order
Riemannian operators.

## 2. Paper 23 (deuteron, He-4, composed nuclear-electronic)

Grep for `\pi`, `transcendental`, `calibration`, `fine-structure`, `1/40`,
`Delta` on paper_23_nuclear_shell.tex returns only:
- Line 642: narrative mention of "the fine-structure constant" in a
  forward-looking "Implications" paragraph — no formula, no calculation.
- Line 799: bibliography entry for Paper 2.
- A single δV(r) symbol (perturbation notation in the Fock-rigidity discussion)
  and the standard ΔE notation for finite-nucleus correction (1.01 × 10⁻¹⁰ Ha),
  which is real physics and ~10⁸ times smaller than the Δ = 1/40 effect
  being investigated.

The deuteron and He-4 Hamiltonians are built from:
1. HO single-particle energies ℏω(2 n_r + l + 3/2) — diagonal, exact in ℏω;
2. Minnesota NN potential (Gaussian parameters V_R, V_s, V_t and exponents κ_R,
   κ_s, κ_t) — fitted constants, not transcendentals of geometric origin;
3. Moshinsky–Talmi brackets and two-species Jordan–Wigner encoding — both
   rational/combinatorial.

There is no place in this construction where a "calibration π" would enter,
and no one-loop-like correction sits between the geometric operator and the
physical spectrum. The nuclear Hamiltonian does not require (and has no room
for) any Δ-like term.

## 3. Paper 22 (potential-independent angular sparsity)

Theorem 2 gives a potential-independent upper bound on ERI density depending
only on l_max. Sec. "Universal vs Coulomb-Specific Structure" explicitly
partitions:

- **Universal** (Coulomb, HO, Woods-Saxon, square well, Yukawa):
  angular sparsity, Gaunt selection rules, composed block architecture.
- **Coulomb-specific**: the universal constant K = -1/16, the S³ Fock
  projection, the π-free graph structure of S³, and the spectral-zeta
  machinery on S³.

Quote: "The universal–specific partition is sharp. Angular sparsity, Gaunt
selection rules, and composed block architectures are universal across
spherical fermion systems … \[the S³ Fock machinery\] is Coulomb-specific
and breaks under any modification of the radial potential."

Δ = 1/40 lives in Paper 2's Hopf-bundle / S³ construction, which sits
firmly on the Coulomb-specific side of this partition. This is fully
consistent with the SM-running hypothesis: Δ should appear only in the
electromagnetic / Coulomb sector and should have no analog in the
HO / nuclear sector.

## 4. α_s running (sign check)

The SM-running hypothesis requires β > 0 for QED (α runs **upward** with
energy, from ~1/137 at m_e to ~1/127 at M_Z), and β < 0 for QCD (asymptotic
freedom, α_s runs **downward** with energy). Any hypothetical α_s analog of
Δ would therefore have the **opposite sign** (Δ_s > 0 as a positive additive
correction, not a subtractive one). No such α_s analog appears anywhere in
Papers 22–24. This is consistent — the nuclear / strong-interaction sector
is built on the HO geometry in GeoVac, and that geometry is already rational
and needs no Δ.

## Surprises

- None. The check is clean and the existing paper framework already
  anticipated this partition.
- Paper 24 Sec. 7.1 already classifies the π in K = π(B + F − Δ) as
  Coulomb-specific via Theorem 3. Track SM-E confirms that the *subtractive*
  piece Δ shares this classification: if it has any physical origin (SM
  running), that origin must also be Coulomb/electromagnetic-specific,
  consistent with Δ being absent from the HO and nuclear constructions.

## One-line verdict

PASS — Δ = 1/40 is absent from Papers 22–24 nuclear/HO constructions; the
universal/Coulomb-specific partition already established in Paper 22 and
sharpened by Paper 24's Theorem 3 is fully consistent with Δ being an
electromagnetic (Coulomb/S³) object with no HO/nuclear analog.
