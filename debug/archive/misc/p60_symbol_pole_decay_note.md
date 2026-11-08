# The other pole: an envelope law for the Paper 60 cross-block symbol

*2026-09-11, PM main-session probe (conversation, not a sprint). Driver:
`debug/p60_symbol_pole_decay.py`. Not yet captured into any paper -- held
pending the Toeplitz literature scan (`debug/lit_scan/toeplitz_finite_section_memo.md`).*

## Context

Paper 60 `eq:sigma_law` governs the **IR pole** of the two-center metric symbol
a(chi) = j0(kR cot(chi/2)): sup a = 1 attained at chi = pi (which is p = 0 under
the Fock map p = k cot(chi/2)), quadratically, giving
1 - sigma_max = (kR)^2/24 * pi^2/n^2 and cond(S) ~ n^2.

The **UV pole** chi -> 0 (p -> infinity) was unexamined. There a(chi) is a chirp:
amplitude ~ chi/(2kR), phase ~ 2kR/chi.

## Result

Stationary phase on the chirp (stationary point chi* = sqrt(2kR/j)) predicts,
with no fitted constants:

    |c_j| = (2 pi)^(-1/2) * 2^(-3/4) * (kR)^(-1/4) * j^(-5/4)
            * |sin(2 sqrt(2 kR j) + pi/4)|  + o(j^(-5/4))

where c_j = (1/pi) int_0^pi cos(j chi) a(chi) dchi, and the metric's off-diagonal
entries are <n|a|m> = c_(n-m) - c_(n+m). So **the off-diagonal envelope exponent
is exactly -5/4**, modulated by a phase that also comes out of the prediction.

Verified by three independent quadrature routes (mpmath quadosc in the
s = kR cot(chi/2) variable; a split fast-phase/slow-phase route; a deterministic
fixed-order Gauss-Legendre route on phase-resolved panels), j = 64..65536,
kR in {1, 2, 5}. Agreement 1-3% except within a few percent of the modulating
sine's zeros, where the leading term is small by construction. The first two
routes fail intermittently at large j (they disagree with each other and return
values that *grow* with j, impossible for a continuous symbol); the GL route is
stable to ~1e-11 under doubling the truncation radius and is the one to trust.
The sign pattern of c_j is reproduced by the predicted phase -- a parameter-free
check, not a fit.

## Consequence: the banded-Loewdin lever is closed, negative

S^(-1/2) is the multiplication operator with symbol (1 +/- a)^(-1/2), so its
locality is set by that symbol's smoothness. Both sectors fail, at different poles:

- **Gerade** (1 + a)^(-1/2) is bounded and smooth at the IR pole -- Paper 60 has
  cond(I + C) -> 2.555041..., flat in basis and in R. But at the UV pole it
  inherits the chirp linearly, (1 + a)^(-1/2) = 1 - a/2 + O(a^2), so its
  coefficients are exactly -1/2 times the raw symbol's. **Measured ratio
  -0.4988, -0.5018, -0.5002 at j = 4096, 16384, 65536** (the j = 1024 point sits
  on a zero of the modulation and is uninformative). Same j^(-5/4) decay.
  l1 band-truncation error therefore falls only as b^(-1/4): bandwidth ~4e5 for
  1e-2, ~4e9 for 1e-3, ~4e21 for 1e-6. Not a lever.
- **Ungerade** (1 - a)^(-1/2) is worse and not merely slowly-decaying. Near
  chi = pi, cot(chi/2) -> (pi - chi)/2 and 1 - j0(x) ~ x^2/6, so
  1 - a ~ (kR)^2 (pi - chi)^2 / 24 and the symbol ~ 2 sqrt(6) / (kR (pi - chi)).
  Then int |symbol| dchi diverges logarithmically: the symbol is not in L^1, so
  its Fourier coefficients do not exist. This is the conditioning blow-up seen
  from the symbol side. (Elementary; asserted analytically, not measured -- the
  numeric probe in the driver cuts the wrong endpoint and does not test it.)

## The sharp statement

**The gerade sector is perfectly conditioned and still not local.** Conditioning
fails at the IR pole; locality fails at the UV pole; they are different poles, so
an operation that fixes one buys nothing for the other. That is a mechanism for
the observed pattern that every repair *relocates* the non-orthogonality cost
rather than removing it (cost conservation; dual-basis theorem, Artacho-del Bosch
PRA 43, 5770, 1991) -- there is no single operation with both poles in its reach.

## Status

MEASURED + SYMBOLIC. Load-bearing claim is the -5/4 envelope and the -1/2
inheritance ratio. Not in any paper yet. Open: whether the -5/4 law is a
rediscovery (stationary phase on a chirp is classical; the question is whether
this symbol class is named in the Toeplitz finite-section literature).
