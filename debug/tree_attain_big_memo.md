# tree_attain_big: Lemma 1 of prop:star_attain on larger trees, modulo two primes

Date 2026-11-19. Scripts `debug/tree_attain_big_modp.py`, `debug/tree_attain_big_fire.py`; data
`debug/data/tree_attain_big_*.json`, logs `debug/data/tree_attain_big_modp*.log`. No paper, test, doc or
tracked file edited.

## Verdict

**Every target CONFIRMED, by the least-order operator (not the lemma's product construction), at two primes.**
For each tree the least-theta-order a in A = C[z]<theta> with a H pure was computed; its symbol at (1, -t) splits
completely over the integers with every root of absolute value <= n-2, it is nonzero at t = -(n-1), and a H is
pure and nonzero, by exact polynomial arithmetic mod p. Every requested tree finished in under 90 s per prime; the
lemma's product construction was never needed.

## Method (and why it is fast)

Same module model as the reviewer's `debug/tree_attain_review_modp.py` and `tests/_tree_attain.py` (rules
re-typed, not imported). Two changes:

1. **Orbit reduction.** H is invariant under the automorphisms of the rooted tree and theta is equivariant, so
   A H sits in the invariant part. Coordinates are orbit sums (orbit = decorated rooted tree up to isomorphism, via
   an AHU-type canonical string, plus the D-power i); theta acts on them by an integer matrix (alpha + beta z),
   integrality asserted. a H is pure iff every non-pure orbit coefficient vanishes (orbits have disjoint
   supports). Star of eight: 7290 keys -> 120 coordinates.
2. **Relation over C(z), then clear denominators.** The least theta-order N is the first j at which
   theta^0 H .. theta^j H become dependent modulo the pure span (rank at two random points). The relation
   is one-dimensional over C(z); its primitive polynomial representative is the least-z-degree a (every
   polynomial solution is a C[z]-multiple of it), i.e. the same object as the reviewer's. It is obtained by
   solving at m random points with c_N = 1, rational reconstruction of each c_j(z), lcm of denominators, division
   by the content gcd. m grows by 1.6x until the **exact** check passes: a H computed as polynomials mod p on every
   coordinate, pure and non-pure. The reconstruction is therefore never trusted; the exact check is the gate.
(The reviewer's script solved one large nullspace for a directly; that cost stopped it on the branch-rooted fork.)

## Controls (reproduce earlier exact / mod-p results)

| tree, root | N | z-deg | weight | roots (multiplicity) | earlier record |
|:--|:-:|:-:|:-:|:--|:--|
| star 3 | 7 | 5 | 13 | 0:3, +-1:2 | memo exact: -4 t^3 (t-1)^2 (t+1)^2 |
| path 4, end | 12 | 14 | 32 | 0:4, +-1:2, +-2:2 | memo exact, reviewer |
| star 4 | 16 | 25 | 56 | 0:6, +-1:3, +-2:2 | memo exact (also without orbit reduction) |
| path 4, inner | 23 | 49 | 105 | 0:7, +-1:4, +-2:4 | memo exact c t^7 (t-+2)^4 (t-+1)^4 |
| path 5, centre | 34 | 134 | 276 | 0:8, +-1:6, +-2:5, +-3:2 | reviewer (also without orbit reduction) |
| fork 5, root 0 | 36 | 162 | 330 | 0:6, +-1:7, +-2:4, +-3:4 | reviewer |
| path 6, end | 30 | 122 | 250 | 0:6, +-1:4, +-2:4, +-3:2, +-4:2 | reviewer (orbit and full) |

All match in every column. The branch-rooted fork also matches with the orbit reduction switched off.

## Results (p = 2147483629 and p' = 2147483587; identical at both primes in every column)

Fork of five = edges 0-1, 1-2, 1-3, 3-4. Spider 6 = 0-1, 1-2, 1-3, 3-4, 3-5. Broom 6 = path 0-1-2-3 with leaves
4, 5 on 3. Path 7 centre = root at the middle vertex.

| tree, root | n | coords (non-pure) | N | z-deg | weight | symbol roots (multiplicity) | at -(n-1) | a H | time | verdict |
|:--|:-:|:-:|:-:|:-:|:-:|:--|:-:|:-:|:-:|:-:|
| fork 5, branch vertex | 5 | 74 (69) | 53 | 335 | 681 | 0:11, +-1:10, +-2:7, +-3:4 | nonzero | pure, 5 nonzero | 1 s | CONFIRMED |
| star 6 | 6 | 56 (50) | 50 | 329 | 670 | 0:12, +-1:8, +-2:6, +-3:3, +-4:2 | nonzero | pure, 6 | 1 s | CONFIRMED |
| spider 6, root 0 | 6 | 88 (82) | 74 | 808 | 1628 | 0:12, +-1:10, +-2:10, +-3:7, +-4:4 | nonzero | pure, 6 | 9 s | CONFIRMED |
| broom 6, root 0 | 6 | 40 (34) | 34 | 162 | 330 | 0:6, +-1:5, +-2:4, +-3:3, +-4:2 | nonzero | pure, 6 | 0.4 s | CONFIRMED |
| path 7, end | 7 | 49 (42) | 42 | 263 | 532 | 0:6, +-1:6, +-2:4, +-3:4, +-4:2, +-5:2 | nonzero | pure, 7 | 0.6 s | CONFIRMED |
| path 7, centre | 7 | 100 (93) | 93 | 1326 | 2667 | 0:15, +-1:12, +-2:12, +-3:8, +-4:5, +-5:2 | nonzero | pure, 7 | 18 s | CONFIRMED |
| star 7 | 7 | 84 (77) | 77 | 859 | 1733 | 0:15, +-1:12, +-2:8, +-3:6, +-4:3, +-5:2 | nonzero | pure, 7 | 10 s | CONFIRMED |
| path 8, end | 8 | 64 (56) | 56 | 500 | 1008 | 0:8, +-1:6, +-2:6, +-3:4, +-4:4, +-5:2, +-6:2 | nonzero | pure, 8 | 4 s | CONFIRMED |
| star 8 | 8 | 120 (112) | 112 | 1946 | 3912 | 0:20, +-1:15, +-2:12, +-3:8, +-4:6, +-5:3, +-6:2 | nonzero | pure, 8 | 47 s | CONFIRMED |
| star 8, rooted at a leaf | 8 | 114 (106) | 106 | 1793 | 3600 | 0:14, +-1:15, +-2:12, +-3:8, +-4:6, +-5:3, +-6:2 | nonzero | pure, 8 | 40 s | CONFIRMED (p only) |

| spider 8 at its centre (legs 2, 2, 3: 0-1-2, 0-3-4, 0-5-6-7) | 8 | 408 (400) | 264 | 11691 | 23416 | 0:34, +-1:30, +-2:30, +-3:24, +-4:18, +-5:9, +-6:4 | nonzero | pure, 8 | 34 min | CONFIRMED (p only) |
| caterpillar 8 (0-1, 1-2, 1-3, 3-4, 3-5, 5-6, 5-7), root 0 | 8 | 370 (362) | 318 | - | - | - | - | - | stopped | NOT DONE (extra) |

The last three rows go beyond the request (branched n = 8 trees) and were run at one prime only. Where the time
went on the spider of eight: the z-degree (11,691) is far above the ratio seen elsewhere (d/N about 4 to 17;
here 44), so the point count climbed 2112 -> 35424 over six rounds (about 1100 s in failed rounds), and the last
round spent 650 s in 35,424 small dense solves (264 x 264) and 250 s in rational reconstruction. The caterpillar
(N = 318, no symmetry to exploit beyond one swap) would need the same or more and was stopped to stay inside the
cap; it is not a failure, only not computed. The lemma's product construction was not needed for any requested
tree.

N = least theta-order; weight = max(2k + j) over terms z^k theta^j. In every case the symbol has degree N in t
and the listed multiplicities sum to N (no other factor mod p); it is even in t. "pure, k" = no non-pure
coordinate of a H survives and k pure coordinates K0^(n-1-i) D^i are nonzero.

## Fire test (`debug/tree_attain_big_fire.py`, data `tree_attain_big_fire.json`)

On star 4, path 5 (centre) and the fork rooted at its branch vertex:
- W1, drop one coefficient of a: a H is no longer pure (2, 2, 3 non-pure coordinates left). Rejected.
- W3, left-multiply a by theta^2 - (n-1)^2 z (symbol gains t^2 - (n-1)^2): a H **stays pure** (the pure span is
  theta-stable), so only the symbol test can reject it; it does: roots +-(n-1) found, value at -(n-1) zero.
- The true a passes both tests on the same code path.

## What this does and does not establish

- It confirms Lemma 1's conclusion on every requested tree (n = 7, 8 and the branch-rooted fork included) by a
  least-order operator. Mod-prime evidence, not a proof; agreement at two primes in every column makes a
  reduction accident unlikely, not excluded.
- Free-generator module model; the true module is a quotient, which cannot break purity. a H nonzero as a
  function is not claimed (the proof gets it by contradiction).
- Nothing here touches Lemma 2, the irreducibility of F_n (inherited for n <= 8), or l > 0.
- Observation, not used: the largest |root| is exactly n-2 in every tree; its multiplicity is 2 for stars,
  paths (end or centre) and the broom, 4 for the fork of five (both rootings checked) and the spider of six.

## Files

- `debug/tree_attain_big_modp.py` (driver; `TREE_BIG_PRIME` env var selects the prime, `--full` disables the
  orbit reduction), `debug/tree_attain_big_fire.py`
- `debug/data/tree_attain_big_<tree>.json`, `..._p2147483587.json`, `..._full.json`
- logs: `debug/data/tree_attain_big_modp_run1.log` (first batch, earlier Vandermonde interpolation),
  `tree_attain_big_modp_star8.log`, `tree_attain_big_modp_p2.log`, `tree_attain_big_modp_extra8.log`
- `debug/tree_attain_big_memo.md`
