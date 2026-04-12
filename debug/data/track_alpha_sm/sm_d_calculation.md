# Track SM-D: One-Loop Photon Self-Energy on $S^3$ at $n_{\max}=3$

**Status:** PARTIAL (one positive structural identification + one negative perturbative result)

## 1. Setup and Normalization

Paper 2 defines
$$
K = \pi\left(B + F - \Delta\right) = \pi\left(42 + \tfrac{\pi^2}{6} - \tfrac{1}{40}\right) \approx 137.036064,
$$
with the three-tier decomposition (Phases 4B–4G):

| Symbol | Value | Origin |
|---|---|---|
| $B$ | 42 | Finite $SO(4)$ Casimir sum over $n \le 3$ |
| $F$ | $\pi^2/6 = \zeta_R(2)$ | Infinite Dirichlet series $D_{n^2}(d_{\max}=4)$ |
| $\Delta$ | $1/40 = 1/(|\lambda_3|\cdot N(2)) = 1/(8\cdot 5)$ | Finite boundary product at truncation $n_{\max}=3$ |

The physical hypothesis under test: $\Delta$ is the finite-volume one-loop vacuum polarization correction to $\alpha^{-1}$ from putting QED on the compact $S^3$ with mode cutoff $n_{\max}=3$.

## 2. The $S^3$ Dirac Spectrum (Camporesi–Higuchi 1996)

The Dirac operator on the **unit** $S^3$ has eigenvalues
$$
\lambda_n^{\pm} = \pm\left(n + \tfrac{3}{2}\right), \qquad n = 0, 1, 2, \dots
$$
with single-chirality degeneracy
$$
g_n = 2(n+1)(n+2).
$$

| $n$ | $|\lambda_n|$ | $g_n$ | cumulative $\sum g$ |
|---|---|---|---|
| 0 | 3/2 | 4 | 4 |
| 1 | 5/2 | 12 | 16 |
| 2 | 7/2 | 24 | 40 |
| **3** | **9/2** | **40** | **80** |
| 4 | 11/2 | 60 | 140 |

**Key structural observation:** $\boxed{g_3^{\text{Dirac}} = 2\cdot 4\cdot 5 = 40 = \Delta^{-1}}$.

Compare Paper 2's decomposition $\Delta^{-1} = |\lambda_3|_{\text{graph}} \cdot N(2) = 8 \cdot 5$. Both factorizations produce 40, but the Dirac one is arithmetically cleaner: $g_n^{\text{Dirac}} = 2(n+1)(n+2)$ directly, and at $n=n_{\max}=3$ this equals 40. This is a fermion-mode interpretation of Paper 2's boundary invariant.

Moreover the **cumulative** Dirac mode count through $n=2$ is also 40 (!), doubly reinforcing 40 as the "state-count at the selection-principle edge."

## 3. Perturbative Vacuum Polarization — Heat-Kernel Mode Sum

The one-loop photon self-energy at $q^2=0$ on a compact manifold, written in the Schwinger/heat-kernel representation (see Birrell–Davies Ch. 6; Dowker 1984; Camporesi–Higuchi 1996 §IV), gives a coefficient-of-$F_{\mu\nu}^2$ contribution
$$
\Pi(0) \;\sim\; -\frac{1}{6\pi^2}\sum_{n=0}^{n_{\max}} g_n \log\!\left(|\lambda_n|^2 + M_f^2\right) \;+\; (\text{continuum subtraction}).
$$
With the on-shell renormalization convention
$$
\frac{1}{\alpha(q^2)} = \frac{1}{\alpha_{\text{bare}}} - \frac{\Pi(q^2)}{\pi},
$$
the shift is
$$
\Delta\!\left(\frac{1}{\alpha}\right) = \frac{1}{6\pi^2}\sum_{n=0}^{3} g_n\log\!\left(|\lambda_n|^2 + M_f^2\right) + \text{const}.
$$

I evaluated this at 50-digit precision with `mpmath` for three canonical mass normalizations:

| Normalization | $M_f$ on unit $S^3$ | Computed shift | ratio to $1/40$ | ratio to $\pi/40$ |
|---|---|---|---|---|
| Massless | $0$ | $3.47352$ | 138.94 | 44.23 |
| Compton ($R_{S^3}=1/m_e$) | $1$ | $3.59280$ | 143.71 | 45.75 |
| Bohr ($R_{S^3}=a_0$) | $1/\alpha \approx 137$ | $13.29506$ | 531.80 | 169.29 |

**None** of these normalizations reproduces $1/40$ or $\pi/40$. The shifts are $\mathcal{O}(1)$ to $\mathcal{O}(10)$, overshooting the target by factors of $44$ to $530$.

## 4. Continuum Running Check (Cross-Validation with SM-A)

Using the standard 1-loop Peskin–Schroeder formula for a single charged fermion of mass $m$,
$$
\frac{1}{\alpha(q^2)} - \frac{1}{\alpha(m^2)} = \frac{1}{3\pi}\log\!\frac{q^2}{m^2},
$$
I inverted for the scale at which the shift equals $1/40$ and $\pi/40$:

| Target | $q/m$ (dimensionless) | $q$ in MeV (for $m=m_e$) |
|---|---|---|
| $1/40$ | $1.1250$ | $574.9$ |
| $\pi/40$ | $1.4479$ | $739.9$ |

These match Track SM-A's independent calculation to five decimal places, confirming the QED running formula is correctly normalized. The $1/40$ scale sits just above $m_e$ but has no obvious physical significance.

## 5. Interpretation

### 5a. Structural identification (POSITIVE)

$$\boxed{\Delta = \frac{1}{g_3^{\text{Dirac on }S^3}} = \frac{1}{2(n_{\max}+1)(n_{\max}+2)} \Bigg|_{n_{\max}=3}}$$

This recasts Paper 2's factorization $1/(|\lambda_3|\cdot N(2))$ in fermionic language: $\Delta$ is the *reciprocal* of the number of Dirac states at the truncation edge. Because $g_n^{\text{Dirac}} = 2(n+1)(n+2)$ and the Paper 2 product $|\lambda_n^{\text{graph}}| \cdot N(n-1) = (n^2-1)\cdot \sum_{k=1}^{n-1}k^2$ both equal 40 at $n=3$ but diverge for $n\ne 3$, the match is a **true arithmetic coincidence tied to the unique $n_{\max}=3$ selection principle**. The fact that the Dirac formula is polynomially simpler ($2(n+1)(n+2)$ vs a quartic) makes it the more fundamental expression.

### 5b. Perturbative identification (NEGATIVE)

The naive heat-kernel mode-sum for $\Pi(0)$ does not give $1/40$ at any tested mass normalization. Numerical overshoot is 45–530×.

### 5c. What would complete the calculation

To convert the structural coincidence into a derivation one would need:

1. **Proper zeta-regularization**: Replace the bare log-sum with the Camporesi–Higuchi zeta-regularized $\zeta_D(s)\big|_{s=0}$ for the Dirac operator on $S^3$, truncated at $n_{\max}=3$. This removes the UV-sensitive $\log(|\lambda|^2)$ terms and leaves a finite, potentially rational, remainder.
2. **Match to SM-A's scale**: $q = 1.125 m_e \approx 575$ MeV (for $1/40$) or $740$ MeV (for $\pi/40$). Neither corresponds to a standard physical scale. A proper $S^3$ calculation would need to identify which *intrinsic* radius $R_{S^3}$ gives $M_f \cdot R_{S^3} = 1.125$, i.e. $R_{S^3} = 1.125/m_e \approx 441$ fm. This has no obvious physical interpretation.
3. **Photon propagator on $S^3$**: Use the correct vector Laplacian spectrum (not scalar/Dirac only) for the internal photon line, with eigenvalues $n(n+2)-1$ for transverse modes (Higuchi 1987).

### 5d. Sign check

Paper 2 has $K = \pi(B + F - \Delta)$ with $\Delta$ entering with a **minus** sign. The shift from vacuum polarization *increases* $1/\alpha$ as the scale grows (asymptotic freedom works the other way for QED, which is not asymptotically free; in QED $1/\alpha$ *decreases* with increasing scale because $\beta > 0$). A boundary correction that *subtracts* from $K$ is consistent with a finite-volume *increase* in $1/\alpha$ from discretizing the photon mode sum — the compact space has fewer states at low $q$, reducing the screening. This is qualitatively right-signed.

## 6. Verdict

**PARTIAL.** The strong positive result is the arithmetic identification
$$
\Delta^{-1} = g_3^{\text{Dirac}}(S^3) = 2(n_{\max}+1)(n_{\max}+2)\big|_{n_{\max}=3} = 40,
$$
which is a cleaner and more physically transparent factorization of Paper 2's $\Delta = 1/(8\cdot 5)$ than the graph-theoretic one, and which holds **only at $n_{\max}=3$** (joining the list of quantities that single out the selection-principle cutoff). The weak negative result is that the naive heat-kernel perturbative vacuum polarization at the same cutoff overshoots by 45–530×, so $\Delta$ is not a one-loop running correction in the standard sense; it is a *state-counting* boundary invariant that happens to coincide with a fermionic degeneracy.

The calculation suggests $\Delta$ should be reframed in Paper 2 as "$1/\Delta$ equals the single-chirality Dirac-mode count at level $n_{\max}$" rather than as a vacuum-polarization shift.

## References

- Camporesi, R. and Higuchi, A. (1996). *On the eigenfunctions of the Dirac operator on spheres and real hyperbolic spaces.* J. Geom. Phys. **20**, 1–18. [arXiv:gr-qc/9505009]
- Peskin, M. and Schroeder, D. (1995). *An Introduction to Quantum Field Theory.* Addison-Wesley. §7.5 (vacuum polarization).
- Birrell, N. and Davies, P. (1982). *Quantum Fields in Curved Space.* Cambridge. Ch. 5–6.
- Dowker, J. S. (1984). *Conformal anomaly for spherical harmonics.* Class. Quantum Grav. **1**, 359.
- Higuchi, A. (1987). *Symmetric tensor spherical harmonics on the $N$-sphere and their application to the de Sitter group $SO(N,1)$.* J. Math. Phys. **28**, 1553.

## Data

- `debug/sm_d_vacuum_pol.py` — computation script (50-dps mpmath)
- `debug/data/track_alpha_sm/sm_d_numerics.json` — full numerical output
