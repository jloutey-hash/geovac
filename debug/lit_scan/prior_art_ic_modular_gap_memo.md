# Prior-art check: the IC instrument's conjugate-root gap (2026-11-01, PI-asked)

**Question.** The planted control 4(c - w^2) y'' - 4 w y' + y = 0, with solutions (w - s)^(1/2) and (w + s)^(1/2)
where s^2 = c, is triangular, yet `debug/intclass_lib.py` calls it STUCK when c is not a square mod p. Is the
equation, or the gap, already known?

## The equation

Put w = sqrt(c) x. It becomes (1 - x^2) y'' - x y' + (1/4) y = 0, which is Chebyshev's equation of order n = 1/2.
With x = cos t it is y'' + n^2 y = 0, so the solutions are cos(t/2) = sqrt((1+x)/2) and sin(t/2) = sqrt((1-x)/2).
Classical. Nothing new here. Checked by direct substitution, and the trigonometric form is standard
(Encyclopedia of Mathematics, "Chebyshev equation"; ProofWiki, general solution).

## The gap

Known, and handled in the literature, for the case where the base field is not algebraically closed.

- Johansson, Kauers and Mezzarobba, arXiv:1301.2486, Sec. 1 and 2.1 (read via ar5iv). The classical method
  computes local solutions at each singular point and tests combinations. It assumes the constant field C is
  algebraically closed, and the singular points are the roots of the leading coefficient in C. Its Sec. 1 credits
  the exact method to van Hoeij, which needs number fields whose degree can grow. It cites van der Put and Singer,
  Ch. 4. Our search instead uses one exponent per F_p-irreducible factor and never extends the constants. That is
  the step this literature shows is not enough.
- Cluzeau and van Hoeij, J. Symb. Comput. 38 (2004) 1043-1076, "A modular algorithm for computing the
  exponential solutions of a linear differential operator" (DOI 10.1016/j.jsc.2003.11.004). It is the modular
  method itself. Its abstract says it avoids "unnecessarily large algebraic extensions of the constants". The
  author PDF has no text layer, so its body was NOT read. What follows rests on the companion paper.
- Cluzeau and van Hoeij, "Computing Hypergeometric Solutions of Linear Recurrence Equations" (author preprint,
  2005-06-15), Sec. 8, "Computing a hard solution". Read in full text. Citing the differential case
  ([CH04, Sec. 6.1, Lemma 6.3]), it says solutions may have to be searched for over algebraic extensions. Its
  Lemma 3 says the conjugates of such a solution are also solutions. Of minimal algebraic degree, distinct
  conjugates are of distinct type and linearly independent. The final algorithm returns solutions "up to
  conjugation". Its last example has two solutions conjugate over Q and needs the "hard case". It also notes
  that Mathematica's Petkovsek implementation finds no solution there, because it "is limited to the case where
  a0 and an have no irreducible factors of degree > 2". That is the same class of miss.

## Reading

The gap and its fix (allow the exponent choice to vary over the conjugates of a degree >= 2 factor, which means
working over the extension F_p[w]/f) are standard computer algebra. What is ours: the planted control, the
telltale test (a degree >= 2 factor with a non-integer exponent class, or exponents outside F_p), and its result
on the scored residuals. No telltale on any clean STUCK residual: N1 (6 of 6, p2), W2 (6 of 6, p2).
W1 (6 of 6, both primes). Downloads: `debug/lit_scan/dl_ic_modular/`, `dl_ic_modular2/` (untrusted, excluded from git).
