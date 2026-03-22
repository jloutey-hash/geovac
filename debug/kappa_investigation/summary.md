# kappa = -1/16 Investigation: Results Summary
# Date: 2026-03-22

## Key Finding: The Two "8"s are NOT the Same at n_max = 3

CRITICAL: At n_max = 3, lambda_max(graph) = 5.0, NOT 8.

The graph Laplacian max eigenvalue only reaches 8 in the ASYMPTOTIC limit
(n_max -> infinity). The "8" in kappa = -1/16 is from:
  kappa = E_1 / lambda_max(graph, n_max->inf) = -0.5 / 8 = -1/16

The S3 eigenvalue |lambda_3| = 3^2 - 1 = 8 is a DIFFERENT 8.

## Numerical Data: lambda_max(graph) vs n_max

| n_max | N_states | lambda_max(D-A) | lambda_max(D-A+W) | n^2-1 | E_ground  | err%   |
|:-----:|:--------:|:---------------:|:-----------------:|:-----:|:---------:|:------:|
|    3  |     14   |        5.000000 |          4.821853 |     8 | -0.301366 | 39.727 |
|    5  |     55   |        6.618034 |          6.551403 |    24 | -0.409463 | 18.107 |
|   10  |    385   |        7.611436 |          7.592472 |    99 | -0.474529 |  5.094 |
|   15  |   1240   |        7.821269 |          7.812430 |   224 | -0.488277 |  2.345 |
|   20  |   2870   |        7.898179 |          7.892703 |   399 | -0.493294 |  1.341 |
|   25  |   5525   |        7.934293 |          7.930740 |   624 | -0.495671 |  0.866 |
|   30  |   9455   |        7.954094 |          7.951632 |   899 | -0.496977 |  0.605 |

## lambda_max(graph) vs |lambda_{n_max}(S3)| Comparison

| n_max | lambda_max(graph) | n^2-1 | match? | ratio |
|:-----:|:-----------------:|:-----:|:------:|:-----:|
|    2  |            3.0000 |     3 |    YES | 1.000 |
|    3  |            5.0000 |     8 |     no | 0.625 |
|    4  |            6.0000 |    15 |     no | 0.400 |
|    5  |            6.6180 |    24 |     no | 0.276 |
|    6  |            7.0322 |    35 |     no | 0.201 |
|    8  |            7.4200 |    63 |     no | 0.118 |
|   10  |            7.6114 |    99 |     no | 0.077 |
|   15  |            7.8213 |   224 |     no | 0.035 |
|   20  |            7.8982 |   399 |     no | 0.020 |

## Analysis

### What the data shows:

1. **lambda_max(D-A) -> 8 asymptotically.** This is a 2D grid property.
   The lattice is bipartite with max degree 4. For bipartite graphs,
   lambda_max <= 2 * d_max = 8. The bound is approached but never reached
   at finite n_max.

2. **At n_max = 3: lambda_max(graph) = 5.0, NOT 8.** The graph is too small
   at n_max = 3 to reach the asymptotic limit. Only 14 states.

3. **The only exact match is at n_max = 2:** lambda_max(graph) = 3.0 = |lambda_2| = 3.
   For all other n_max, the two quantities diverge.

4. **kappa = -1/16 works because of the ASYMPTOTIC lambda_max = 8**, not
   because of the S3 eigenvalue at the Hopf cutoff.

### The Two "8"s:

- **Graph 8:** 2 * d_max = 2 * 4 = 8 (bipartite 2D grid with 4 neighbors)
- **S3 8:** |lambda_3| = 3^2 - 1 = 8 (Laplace-Beltrami on unit S3 at n=3)

Are they structurally related? The graph is a discrete approximation to S3.
But the graph lambda_max = 8 comes from the GRID TOPOLOGY (degree bound),
while the S3 lambda_3 = -8 comes from the CASIMIR OPERATOR (representation
theory of SO(4)). These are genuinely different mathematical structures.

### Number-theoretic coincidence:

The equation 2 * d_max = n_max^2 - 1 with d_max = 4 gives:
  8 = n_max^2 - 1  =>  n_max = 3

And d_max = 4 because the lattice has 2 transition types (radial, angular)
each allowing +/- steps = 4 total neighbors.

So: 2 * (2 * dim_transition) = n_max^2 - 1 => 4 * dim_transition = n_max^2 - 1
With dim_transition = 2: n_max = 3 is the unique solution.

This IS suggestive: the number of transition dimensions in the discrete lattice
(2: radial and angular) determines d_max, and this d_max picks out n_max = 3
as the S3 eigenvalue that matches the grid bound. Since the lattice is built
from quantum numbers that live on S3, the 2 transition dimensions are related
to the S3 geometry.

### Is kappa derivable from the Hopf bundle?

**Not directly as initially hypothesized.** The formula:
  kappa = -1/(2 * |lambda_{n_max}(S3)|) with n_max = 3 from Hopf

gives the right answer (-1/16) but for the WRONG reason at finite n_max.
At n_max = 3, the graph has lambda_max = 5.0, giving E_ground = -0.301
(39.7% error). The formula works because |lambda_3(S3)| = 8 happens to
equal the asymptotic graph lambda_max, not because the graph at n_max = 3
faithfully represents S3.

**However**, there is a deeper structural argument:
- kappa = E_1 / lambda_max(graph, n_max -> inf)
- lambda_max(graph, n_max -> inf) = 2 * d_max = 2 * 4 = 8
- d_max = 4 because 2 transition types * 2 directions
- The Hopf selection n_max = 3 satisfies n^2 - 1 = 2 * d_max
- So kappa = -1/(2(n_Hopf^2 - 1)) = -1/16 is a VALID derivation,
  provided one accepts that n_Hopf^2 - 1 = lambda_max(graph) is structural.

### Self-consistent kappa test:

kappa_self(n_max) = E_exact / lambda_max(n_max):
  n_max =  3: kappa_self = -0.10369458, ratio to -1/16 = 1.6591
  n_max =  5: kappa_self = -0.07631953, ratio to -1/16 = 1.2211
  n_max = 10: kappa_self = -0.06585471, ratio to -1/16 = 1.0537
  n_max = 15: kappa_self = -0.06400058, ratio to -1/16 = 1.0240
  n_max = 20: kappa_self = -0.06334965, ratio to -1/16 = 1.0136
  n_max = 25: kappa_self = -0.06304582, ratio to -1/16 = 1.0087
  n_max = 30: kappa_self = -0.06288018, ratio to -1/16 = 1.0061

### Alternative kappa formulas:

| n_max | kappa = -1/(2(n^2-1)) | E_ground (at that n_max) | err%   |
|:-----:|:---------------------:|:------------------------:|:------:|
|   2   |  -1/6 = -0.1667      |  -0.4583                 |  8.3%  |
|   3   |  -1/16 = -0.0625     |  -0.3014                 | 39.7%  |
|   4   |  -1/30 = -0.0333     |  -0.1958                 | 60.8%  |
|   5   |  -1/48 = -0.0208     |  -0.1365                 | 72.7%  |

Note: kappa = -1/6 at n_max = 2 actually gives BETTER accuracy at n_max = 2
than kappa = -1/16 does at n_max = 3. But kappa = -1/16 is the correct
asymptotic value for all n_max -> infinity.
