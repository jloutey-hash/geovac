# Banana dictionary for Paper 59's tree ladder (2026-10-17)

**Question.** Is Paper 59's tree ladder F_n(t) = int_0^inf K0(x)^(n-1) K0(t x) dx the known
banana / Bessel-moment ladder in disguise, at the level of differential operators? What in it is
prior art?

**Verdict: GO, but not as an equivalence.** F_n is not gauge-equivalent to any banana operator.
It is related to the equal-mass banana by one exact, standard operation: a Hadamard product with
the Legendre family (A) = binom(2m,m)^2. For n = 2, 3, 4, 5 the resulting operators are
verbatim entries of the published tables. Every step is checked in exact rational arithmetic.

## 1. The exact map

Notation. theta = t d/dt. The operators are those of `debug/data/linking_ceiling_B_ode.json`,
re-derived independently (`banana_dict_ops.py`: F == paper for n = 2..5).
Banana: Vanhove's I_N(s) = 2^(N-1) int x I0(sqrt(s) x) K0(x)^N dx (arXiv:1401.6438, eq. 9.1,
Table 1). This is the equal-mass (N-1)-loop banana with N propagators in D = 2.

**Statement.** Put y = 1/(16 t^2) and F_n = t^(-1) g(y). Then

  L_{F_n}(y) = (A) * B_{n-1},

where * is the Hadamard product of Almkvist and Zudilin (math/0402386, Sec. 7). Here
(A) = theta^2 - 4y(2theta+1)^2 is the operator of sum binom(2m,m)^2 y^m, the Legendre K. And
B_N is the operator of sum_k W_N(2k) y^k. W_N(2k) is the sum of squared multinomials
(Broadhurst 2009, eqs. 13, 60, 61). It is the large-s holomorphic period of the N-propagator
banana.

Why it holds (two lines). The Mellin transform in t of F_n is B(z) A_{n-1}(1-z), with
B(z) = 2^(z-2) Gamma(z/2)^2. Here A_m is the Mellin transform of K0^m. For the weight-x moment
G_n(t) = t int x K0(x)^(n-1) K0(t x) dx the transform is B(z+1) A_{n-1}(1-z). The ratio is
r(z) = Gamma(z/2)^2 / (2 Gamma((z+1)/2)^2). It obeys r(z+2)/r(z) = z^2/(z+1)^2, which is the
term ratio of binom(2m,m)^2. So, as an identity of functions,

  F_n(t) = (2/pi) int_0^1 K(sqrt(1-u^2)) G_n(t/u) du/u .

[Transcendentals: K is the complete elliptic integral, the Legendre period. The 2/pi is the
Gamma(1/2)^2 = pi of the Mellin ratio, an angular-measure pi (Paper 18 calibration class).]

A second form of the same fact needs no Mellin transform. Since int_0^inf J0(p x) dp = 1/x,

  F_n(t) = int_0^inf dp  int_0^inf x J0(p x) K0(t x) K0(x)^(n-1) dx .

The inner integral is the n-propagator 2-D banana at Euclidean momentum p, with masses
(t, 1, ..., 1). So **F_n is the n-propagator banana integrated over a line of external momentum.**

### Checks run

| Check | Script | Result |
|:--|:--|:--|
| Mellin re-derivation of L_{F_n} equals Paper 59's, n = 2..5 | `banana_dict_ops.py` | exact equality |
| Hadamard-derived recurrence = c(z) x recurrence of L_{F_n} (rational left factor only), n = 2..5 | `banana_dict_divide.py` | remainder exactly 0 |
| Recurrence of L_{F_n} for the t = inf series is proportional to (A) * Broadhurst's W_{n-1} recurrence; the series binom(2m,m)^2 W_{n-1}(2m) satisfies it, m < 40 | `banana_dict_hadamard.py` | exact, residual 0 |
| Vanhove's banana operator (N = 2..5 propagators) is a right factor of the operator of int x K0(x)^N K0(t x) dx, s = t^2 | `banana_dict_vanhove.py` | exact |
| Vanhove's banana operator kills sum W_N(2k) s^(-k-1), N = 2..5; fire-tested (one coefficient perturbed: fails) | `banana_dict_walks_banana.py` | exact |
| L_{F_n} in y equals the AZ / AESZ operators verbatim (below) | `banana_dict_aesz.py` | exact equality, n = 2..5 |
| Kernel identity, n = 2, t = 0.5, 1.3 | `banana_dict_kernel.py` | agree to 2e-19 (n = 3 nested quadrature timed out at 600 s; not checked numerically) |
| Line-integral identity, n = 2 (closed-form bubble), t = 0.5, 1.3 | `banana_dict_line.py` | agree to working precision (25 digits) |

Data: `debug/data/banana_dict_{ops,divide,hadamard,vanhove,walks_banana,aesz,kernel,line}.json`.

## 2. The dictionary, rung by rung

| n | F_n in y = 1/(16t^2) | = (A) * banana with n-1 propagators | literature entry |
|:-:|:--|:--|:--|
| 2 | theta^2 - 4y(2theta+1)^2 | (A) * tadpole (W_1 = 1) | Legendre K, classical |
| 3 | theta^3 - 8y(2theta+1)^3 | (A) * bubble (W_2 = binom(2m,m)) | sum binom(2m,m)^3, Clausen 3F2(1/2,1/2,1/2;1,1;64y) |
| 4 | theta^4 - 4y(2theta+1)^2(10theta^2+10theta+3) + 144y^2(2theta+1)^2(2theta+3)^2 | (A) * sunrise, AZ sequence (c) = W_3 | **AESZ #58**, listed as (A)*(c) |
| 5 | theta^5 - 8y(2theta+1)^3(5theta^2+5theta+2) + 2^10 y^2(theta+1)(2theta+1)^2(2theta+3)^2 | (A) * three-loop banana, AZ sequence (alpha) = W_4 | **the fifth-order operator D' whose pullback is AESZ #39**, listed as (A)*(alpha) and (c)*(e) |

Notes on the table. AZ sequence (c) = sum binom(m,k)^2 binom(2k,k) equals W_3(2m) (Broadhurst
2009, eq. 13). AZ (alpha) = sum binom(m,k)^2 binom(2k,k) binom(2m-2k,m-k) equals W_4(2m): group
the four multinomial indices in pairs. Almkvist (arXiv:0808.1480) independently identifies (alpha)
as Verrill's four-term sum. The fifth-order D' is printed verbatim under AESZ #39. Its
superseeker row is N0 = 3, |N1| = 96, |N3| = 26208. Row #58 is N0 = 3, |N1| = 48, |N3| = 11056,
which matches Paper 59's instanton numbers 16 and 11056/3 under N0 = 3.

## 3. The off-by-one, resolved exactly

In D = 2 the banana is a weight-x moment: int x I0 K0^N. That is Vanhove eq. 9.1, and Almkvist's
main theorem uses the odd moments c_{m,2n+1}. Paper 59's F_n is a weight-x^0 moment. The two
weights differ by the Mellin factor r(z). That factor is a Hadamard product with the Legendre
family. A Hadamard product of a CY-d1 operator with a CY-d2 operator is of CY type d1 + d2 + 1, and
Legendre has d = 1. So it raises the Calabi-Yau dimension by 2 relative to the banana with n-1
propagators. Relative to the banana with all n Bessel factors as propagators, it raises it by 1:
that is the momentum-line integral. With your counting (only the K0(x) factors as propagators):

- 3 factors: the sunrise is elliptic; F_3 is K3. F_3 is the sunrise integrated over p, and (A) * bubble.
- 4 factors: the 3-loop banana is K3; F_4 is CY3 (AESZ #58 = (A) * sunrise).
- 5 factors: the 4-loop banana is CY3 (AESZ #34, the 5-step walk); F_5 is (A) * three-loop banana
  (pulls back to AESZ #39).

So the rung of F_n matches the banana with n+1 propagators, but the family does not. The
singular points show it at once. F_n's singular values of t^2 are {1}, {4}, {1,9}, {4,16} for
n = 2..5. These are the thresholds of the banana with n-1 propagators, not of the banana with
n+1 propagators (e.g. {1,9,25} for 5 propagators). Gauge equivalence in a fixed variable keeps
true singular points, so no gauge map exists. A rational pullback was not searched. The Hadamard
statement makes that search moot.

## 4. What is prior art in Paper 59's ladder

- **Tree identity (eq. tree).** A one-line Fourier argument. I found no source stating it for
  molecular densities. The integral it produces, int prod K0(a_v x) dx, is the standard
  Bessel-moment family (BBBG 2008; not re-read here).
- **n = 2, 3 operators.** Classical: Legendre K and Clausen's 3F2.
- **n = 4 operator.** Fully prior art: AESZ #58, with its Hadamard factorization (A)*(c), its
  holomorphic coefficients, and its superseeker integers. Paper 59 already credits this.
- **The Hadamard description for all n.** The construction (A) * L and the n = 4, 5 tables are
  Almkvist and Zudilin's (math/0402386, Sec. 7). The identification of W_N with the banana is
  standard (Broadhurst 2009 for W_N and lattice Green functions; Vanhove 2014 and Almkvist 2008
  for banana operators). **Not found in any source I read:** the statement that the weight-x^0
  moment with one rescaled argument equals (A) * banana for every n. Its proof is the two-line
  Mellin computation above, so it is a corollary of classical facts, not a new theorem.
- **n = 5.** Prior art, and it **answers Paper 59's open question.** The paper says whether F_5
  "is new relative to n = 4 is open". Its operator is verbatim the AESZ fifth-order D' for #39,
  with N0 = 3, |N1| = 96, |N3| = 26208. That differs from #58 (|N1| = 48), so F_5 is tied to a
  different Sp4 operator, #39, not to the F_4 rung.
- **Irreducibility and quadric counts (0, 1, 0, 1).** These follow from the listed operators being
  CY operators (n = 4) or exterior-square-type (n = 5). Paper 59's own computation is a consistent
  re-derivation.
- **Ring (eq. star_triangle).** The 3-D uniqueness relation is classical. Preti's STR paper
  (arXiv:1811.04935, eq. 2.7) gives the general-D formula. At D = 3 with unit indices it gives
  exactly Paper 59's coefficient: pi^(3/2) Gamma(1/2)^3 / Gamma(1)^3 = pi^3, with exponent 1 on
  each |x_ij|. [Transcendental: pi^3, three Gamma(1/2) factors plus pi^(D/2), angular-measure pi.]
  Preti attributes the first 3-D use to D'Eramo, Parisi and Peliti, Lett. Nuovo Cim. 2, 878
  (1971), and the method to Kazakov, Phys. Lett. 133B, 406 (1983) and TMF 58, 223 (1984). Paper 59
  cites none of these. The fibre reading (Struve moment, L_G, the factorization into elementary
  and Legendre pieces) was not found in the sources read. Under the dictionary, its elliptic piece
  is the n = 2 rung, (A) itself.

## 5. Citations

VERIFIED (read in the primary source this session):
- Vanhove, arXiv:1401.6438: eq. 9.1 (banana = 2^(N-1) int x I0(sqrt t x) K0^N); Table 1
  operators N = 2..6; recursion of Borwein-Salvy. Used as data in `banana_dict_vanhove.py`.
- Broadhurst, "Bessel moments, random walks and Calabi-Yau equations" (preprint, 2 Nov 2009,
  carmamaths.org): W_n definition, eq. 13 (squared multinomials), eq. 14 and P_{n,j}, eq. 60-62
  (Green functions G_n; G_3 tied to three-particle decays), D = 4 diamond = AESZ #34. It works
  with J0^n moments and W_n(2k). It does **not** discuss binom(2m,m)^2 Hadamard products or
  K0(t x).
- Almkvist, van Enckevort, van Straten, Zudilin, arXiv:math/0507430: entry #58 (operator and
  A_n), entry #39 (the pullback of D', printed), superseeker rows for #58 ((A)*(c)) and #39
  ((c)*(e), (A)*(alpha)).
- Almkvist, Zudilin, arXiv:math/0402386, Sec. 7: definitions of (A), (c), (alpha); the
  Hadamard formulas for orders 4 and 5; the tables (c) x (A) = 58 and (alpha) x (A) = 39.
- Almkvist, arXiv:0808.1480: main theorem (odd moments c_{m,2n+1} and Verrill's A_n^(m) share a
  PF equation); m = 4 gives (alpha), m = 5 gives #34.
- Fresan, Sabbah, Yu, arXiv:2006.02702, Sec. 1: M_k is symplectic for k even and orthogonal for
  k odd, k = number of Bessel factors. This matches Paper 59's parity sentence.
- Bonisch et al., arXiv:2008.10574 (abstract only): banana linked to CY manifolds of dimension
  loop order minus one.
- Preti, arXiv:1811.04935: eq. 2.7 and its references [14], [17], [18].

NOT VERIFIED (not read here): D'Eramo-Parisi-Peliti 1971 and Kazakov 1983/84 (seen only as
Preti's citations); BBBG 2008 and its Eq. 33; Bloch-Kerr-Vanhove 2015; Broadhurst 2016; Zhou
2018; Gegelia-van Straten 2024; CYDB entry 2.9.

## 6. Defects found in Paper 59

1. **sec:three_electron, "the repulsion integral pairs with the two-loop sunrise and the chain
   with the three-loop banana."** True of the rung, false of the family. F_2 = (A) * tadpole is
   the momentum-line integral of the one-loop bubble. F_3 = (A) * bubble is the line integral of
   the sunrise. The chain's K3 is binom(2m,m)^3 (Clausen). The three-loop banana's K3 is W_4
   (Verrill, AZ (alpha)). Suggest: "sits on the same rung as".
2. **sec:four_density, "whether it is new relative to n = 4 is open."** Stale against the
   literature. F_5's operator is AESZ's printed D' for #39, distinct from #58.
3. **sec:four_density, "though the operators differ."** True but understated. The operators are
   related exactly: L_{F_n} = (A) * banana_{n-1}. The Hadamard construction and the n = 4, 5 tables
   should be credited to Almkvist-Zudilin (math/0402386). The paper credits the Hadamard form to
   Broadhurst's walk preprint, which treats W_n but not this product.
4. **sec:density_ring, "the uniqueness relation is classical."** It carries no citation. Add
   D'Eramo-Parisi-Peliti (1971) and Kazakov (1983/84) after checking them in the primary sources.
5. **Possible, not audited: sec:literature / abstract wording "on the sunrise elliptic curve".**
   At the D = 0 fibre the two-density period is the Legendre family, i.e. the line integral of the
   bubble, not the sunrise. Whether the shifted-argument object of eq. besselmoment is on a
   sunrise curve was not checked.

Nothing outside debug/ was modified. Nothing was committed.
