# EP-3: Derivation of S_n = k log A_n with A_n ∝ n^4

## 1. Paper 5 passage (verbatim, lines 487-498 of `papers/conjectures/Paper_5_Geometric_Vacuum.tex`)

> \subsection{Entropy and the second law}
> The holographic entropy of shell $n$ is (Paper 3, 4):
> \begin{equation} S_n = k \ln A_n + \text{const}, \end{equation}
> where $A_n \propto n^4$ is the "surface area" (number of plaquettes). Since $A_n$ increases with $n$, entropy is monotonic:
> \begin{equation} \frac{dS}{dn} = \frac{k}{A} \frac{dA}{dn} = \frac{4k}{n} > 0. \end{equation}

Also line 702: `S_n ~ ln(n^4) = 4 ln n` (labeled "2D CFT regime").

The claim is asserted, not derived, in Paper 5. Paper 0 contains only the one-body shell count $N_k = 2(2k-1)$ (per-shell) and cumulative $2n^2$ — no n^4.

## 2. Microstate counts on S^3 Coulomb packing

| Quantity | Formula | n=1 | 2 | 3 | 4 | 5 | 10 | asymptotic |
|---|---|---|---|---|---|---|---|---|
| Orbital deg g_n^orb | n² | 1 | 4 | 9 | 16 | 25 | 100 | n² |
| With spin g_n | 2n² | 2 | 8 | 18 | 32 | 50 | 200 | 2n² |
| Cumulative N(n) | n(n+1)(2n+1)/3 | 2 | 10 | 28 | 60 | 110 | 770 | (2/3)n³ |
| Pair within shell g_n(g_n−1)/2 | n²(2n²−1) | 1 | 28 | 153 | 496 | 1225 | 19900 | 2 n⁴ |
| Pair product g_n² | 4n⁴ | 4 | 64 | 324 | 1024 | 2500 | 40000 | 4 n⁴ |
| Adjacent-shell pairs g_n·g_{n−1} | 4n²(n−1)² | 0 | 16 | 144 | 576 | 1600 | 32400 | 4 n⁴ |

## 3. Shannon entropies

**Single-shell uniform entropy** (one-body):
S_n = log g_n = log(2n²) = **2 log n + log 2**

This is the entropy of a single electron randomly placed in shell n. It does NOT give 4 log n. Factor is **off by 2**.

**Cumulative entropy** (up to shell n):
S_cum(n) = log N(n) ~ **3 log n + log(2/3)** asymptotically.

Also not 4 log n.

## 4. Pair-counting entropy (two-body)

**Unordered pairs within shell n** (two-electron antisymmetrized):
S_pair(n) = log[g_n(g_n−1)/2] = log[n²(2n²−1)] = **4 log n + log 2 + O(1/n²)**

**Product/ordered pairs** g_n²:
S_pair(n) = log(4n⁴) = **4 log n + 2 log 2**

**Adjacent-shell pairs** g_n · g_{n−1}:
log ~ **4 log n + 2 log 2** + O(1/n)

All three two-body counts give `4 log n`. The factor of 4 is the SAME in all pair interpretations: it comes from log(n²) + log(n²) = 2·2 log n, i.e. each electron contributes 2 log n from its own g = 2n², and there are two of them.

## 5. Derivation sketch

On the S^3 packing lattice:
- Shell n carries g_n = 2n² states (Paper 0 cumulative, equivalent to sum_{l=0..n-1} 2(2l+1)).
- A **plaquette** is a two-site elementary face — i.e., a (site, site) pair. The natural counting on shell n is either g_n² (directed plaquettes) or g_n(g_n−1)/2 (undirected).
- Either way, A_n = #plaquettes_n ∝ n⁴. log A_n = 4 log n + const. ✓

The area-law scaling A_n ∝ n⁴ is therefore **intrinsically two-body**. It cannot be recovered from single-particle (one-body) microstate counting, which gives at most 2 log n per shell or 3 log n cumulatively.

Hopf bundle / antipodal pairing does NOT produce a second factor of 2 (antipodal identification halves, not doubles, the count). Spin is already included in g_n. The only clean source of the doubling is **two electrons = two factors of g_n in the microstate count**.

## 6. Verdict

**Paper 5's claim S_n = k log A_n with A_n ∝ n^4 is derivable from microstate counting — but ONLY as a two-body statement.** The single-shell single-particle Shannon entropy is log(2n²) = 2 log n + log 2, which is a factor of 2 short.

The factor of 4 emerges naturally and uniquely from pair-counting on the (n, l, m, spin) lattice:
    A_n = #{plaquettes on shell n} ≈ g_n² ∝ n⁴.

This **strengthens the EP-1 thesis**: the area law is a two-body phenomenon, consistent with the view that entropy in GeoVac lives in the electron-electron potential V_ee (off-diagonal pair structure), not in the one-body graph Laplacian h₁.

## 7. Implications for Paper 27

The holographic entropy scaling is a statement about **pair states**, not single particles. Paper 27 should present A_n = g_n² (plaquette = ordered pair of sites on shell n) as the operational definition, note that A_n ∝ n⁴ then follows immediately from the Paper-0 shell degeneracy g_n = 2n², and frame the area law as a consequence of two-body microstate counting — consistent with the EP sprint thesis that entropy lives off the diagonal of V_ee.
