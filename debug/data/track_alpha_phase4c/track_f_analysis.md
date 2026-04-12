# Track alpha-F Analysis: Delta = 1/(B - N_init)?

## Setup
Tests whether Paper 2's Delta = 1/(|lambda_{n_max}| * N(n_max - 1))
coincides algebraically with 1/(B(n_max) - N_init), where N_init = 2
(Paper 0 Axiom 1), and whether any coincidence is polynomial in n_max.

Closed forms:
- |lambda_n| = n^2 - 1
- N(m) = m(m+1)(2m+1)/6
- B(m) = sum_{n=1}^{m} sum_{l=0}^{n-1} (2l+1) l(l+1)
- B(m) closed form: m**5/10 + m**4/4 - m**2/4 - m/10
- B(3) = 42 (verified against Paper 2 Eq. box line 226)

## Subtask 1: Algebraic equivalence

Compare LHS(m) = |lambda_m| * N(m-1) against RHS(m) = B(m) - N_init.

| m | \|lambda_m\| | N(m-1) | LHS | B(m) | RHS = B(m) - 2 | LHS - RHS | agree |
|---|-----|--------|-----|------|------|------------|-------|
| 1 | 0 | 0 | 0 | 0 | -2 | 2 | no |
| 2 | 3 | 1 | 3 | 6 | 4 | -1 | no |
| 3 | 8 | 5 | 40 | 42 | 40 | 0 | YES |
| 4 | 15 | 14 | 210 | 162 | 160 | 50 | no |
| 5 | 24 | 30 | 720 | 462 | 460 | 260 | no |
| 6 | 35 | 55 | 1925 | 1092 | 1090 | 835 | no |
| 7 | 48 | 91 | 4368 | 2268 | 2266 | 2102 | no |
| 8 | 63 | 140 | 8820 | 4284 | 4282 | 4538 | no |

### Symbolic polynomial check

- LHS(m) = m**5/3 - m**4/2 - m**3/6 + m**2/2 - m/6
- RHS(m) = m**5/10 + m**4/4 - m**2/4 - m/10 - 2
- LHS(m) - RHS(m) = 7*m**5/30 - 3*m**4/4 - m**3/6 + 3*m**2/4 - m/15 + 2

**Polynomial identity in m: DOES NOT HOLD**

**Verdict for Subtask 1: COINCIDENCE AT m=3 ONLY (fails at m = [1, 2, 4, 5, 6, 7, 8])**

## Subtask 2: Alternative decompositions of Delta = 1/40

### Manual candidates
| Formula | Value | Equals 1/40 |
|---------|-------|-------------|
| 1/(B - N_init) = 1/40 | 1/40 | YES |
| kappa^2 * B | 21/128 | no |
| |kappa|/B | 1/672 | no |
| B * |kappa|^2 | 21/128 | no |
| 1/(|lambda_3| * N(2)) [Paper 2] | 1/40 | YES |
| N_init / (B * N(2)) | 1/105 | no |
| 1/(d_max^2 * N(2)) | 1/80 | no |
| 1/(8*5) | 1/40 | YES |
| 1/(2*20) | 1/40 | YES |
| 1/(4*10) | 1/40 | YES |
| 1/(d_max * 2 * N(2)) | 1/40 | YES |
| 1/(|lambda_3| * (B - N_init)/|lambda_3|) | 1/40 | YES |
| |kappa| * (|lambda_3|/N(2)) | 1/10 | no |
| 2*|kappa| * (|lambda_3|/N(2)) | 1/5 | no |

### Systematic single-quantity search (a/(b*Q) = 1/40)
Found 15 hits (including trivial reorderings):

- 8/(1*N(2)) (Q = 5)
- 5/(1*|lambda_3|) (Q = 8)
- 10/(2*|lambda_3|) (Q = 8)
- 1/(1*B-N_init) (Q = 40)
- 2/(2*B-N_init) (Q = 40)
- 3/(3*B-N_init) (Q = 40)
- 4/(4*B-N_init) (Q = 40)
- 5/(5*B-N_init) (Q = 40)
- 6/(6*B-N_init) (Q = 40)
- 7/(7*B-N_init) (Q = 40)
- 8/(8*B-N_init) (Q = 40)
- 9/(9*B-N_init) (Q = 40)
- 10/(10*B-N_init) (Q = 40)
- 10/(1*d_max) (Q = 4)
- 10/(11*B+N_init) (Q = 44)

### Two-quantity products with value 40
- 2/(d_max*N(2)) = 1/10
- 1/(|lambda_3|*N(2)) = 1/40
- 1/(N(2)*|lambda_3|) = 1/40
- 2/(N(2)*d_max) = 1/10


## Verdict

**Subtask 1:** NOT a polynomial identity; see m-table for failing rows
**Subtask 2:** 7 manual decompositions equal 1/40.
The forms 1/(|lambda_3| * N(2)) and 1/(B - N_init) both equal 1/40, but for distinct reasons unless the m-polynomial identity holds.

**Classification of Delta:**

At m = 3 both quantities equal 40, but they disagree elsewhere.
This means the Phase 4B observation Delta = 1/(B - N_init) is a single-point
coincidence, not a structural identity. Delta remains an independent ingredient
in Paper 2's K formula.

Failing rows (m != 3): [1, 2, 4, 5, 6, 7, 8]

**Recommendation: CLOSE this line of investigation.** Delta is not reducible to
B and N_init except at the specific cutoff m = 3. Phase 4B's observation was
real at the data point but does not generalize.
