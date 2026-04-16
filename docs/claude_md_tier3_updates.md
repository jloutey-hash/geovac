# CLAUDE.md Tier 3 Updates (v2.11.0 -> v2.12.0)

## Edit 1: §1 version bump
Change: `v2.11.0 (April 15, 2026)` -> `v2.12.0 (April 15, 2026)`

## Edit 2: §2 Tier 3 summary
Append after any existing Tier 2 summary in §2 (Current Development Frontier):

"Tier 3 (v2.12.0): gamma radial corrections (T7, exact for n_r=0), Darwin+mass-velocity completing alpha^4 ladder (T8, Dirac formula exact), squared Dirac D^2 spectral zeta theorem (T9: zeta_{D^2}(s) = polynomial in pi^2 at every integer s; 4th cell of Paper 18 grid degenerate with calibration; operator order is the transcendental discriminant). Fine-structure splitting accuracy unchanged (66-211%, multi-electron SS/SOO needed)."

## Edit 3: §3 new failed approach row
Add to the table in Section 3:

| Darwin+MV for He/Li/Be 2p-doublet improvement | 1 | NEGATIVE. Both 2p states share l=1, so Darwin=0 for l>=1 and MV cancels in the splitting. Residual 66-211% errors trace to multi-electron SS/SOO. Tier 3 Track T8. |

## Edit 4: §10 new benchmark rows
Add to the Validation Benchmarks table:

| Fine-structure Dirac formula (n<=4) | exact symbolic | Dirac formula verified for all (n,l,j) through n=4 |
| Dirac accidental degeneracy | 6/6 pairs confirmed | E_FS depends on (n,j) only |
| gamma radial <1/r> NR limit | Z/n^2 | Hellmann-Feynman all-state formula |
| gamma radial <1/r^2> n_r=0 NR limit | Z^2/(n^3(l+1/2)) | Pochhammer ratio |
| gamma radial <1/r^3> n_r=0 NR limit | Z^3/(n^3 l(l+1/2)(l+1)) | Pochhammer ratio |
| zeta_{D^2}(2) | pi^2 - pi^4/12 | Squared Dirac spectral zeta (T9) |
| zeta_{D^2}(s) pi^{even} only | theorem | No odd-zeta content at any integer s |

## Edit 5: §11 new topic mappings
Add to the Topic-to-Paper Lookup table:

| Lichnerowicz formula on S^3 | 18 | §IV | Core |
| Squared Dirac spectral zeta | 18 | §IV (T9 result) | Core |
| Darwin term | 14 | §V | Core |
| Mass-velocity term | 14 | §V | Core |
| Operator-order transcendental discriminant | 18, 24 | §IV, §V | Core |
| Dirac fine-structure formula | 14 | §V | Core |
| Kramers-Pasternak recursion (deferred) | 18 | §IV | Core |

## Edit 6: §12 gamma radial status
In the Level 5 (Composed Geometries) algebraic registry section, add or update:

| gamma radial <r^s> (s=-1, all states) | algebraic | Hellmann-Feynman from exact Dirac energy (T7, v2.12.0) |
| gamma radial <r^s> (s=-2,-3, n_r=0) | algebraic | Pochhammer ratios from single-term wavefunction (T7, v2.12.0) |
| gamma radial <r^s> (s=-2,-3, n_r>=1) | numerical-required | Kramers-Pasternak three-term recursion needed (T7, deferred) |
