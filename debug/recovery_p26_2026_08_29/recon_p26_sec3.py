r"""RECONSTRUCTION step 2: Paper 26 SS III, rebuilt in one pass.

Source of truth: debug/recovery_p26_2026_08_29/paper_26_PRE_TRUNCATION_text.txt
All SS III values RE-VERIFIED post-order-fix (identical -- counts, densities
and the theta-profile are order-independent):
  107/625 = 17.1%   canonical 41/225 = 18.2%   smallest nonzero 0.0172
  theta-profile (seed 42): 0.427 / 0.702 / 0.830 / 1.000
  spans: theta=1e-8 0.38-0.46 ; theta=1e-4 0.76-0.83
  n_max=4: 57,700/810,000 = 7.12% ; canonical 15,293/216,225 = 7.07%
  density improves ~ M^-0.49 over M = 5, 14, 30
"""
import io

P = "papers/group6_precision_observations/paper_26_entanglement.tex"
s = io.open(P, encoding="utf-8").read()

old = r"""At $n_{\max} = 2$ (five spatial orbitals: $1s, 2s, 2p_{-1}, 2p_0,
2p_1$), the graph basis produces 265 nonzero ERIs out of
$5^4 = 625$ possible four-index combinations, corresponding to an
ERI density of $42.4\%$.  Rotating the basis by any angle
$\theta > 10^{-6}$ fills the ERI tensor to $\sim 620/625$ nonzero
entries ($99.2\%$ density).

The transition is effectively a step function at the identity:
\begin{equation}
  D(\theta) =
  \begin{cases}
    0.424, & \theta = 0 \\
    0.992, & \theta > 10^{-6}.
  \end{cases}
  \label{eq:step}
\end{equation}

This result is not an artifact of the small basis.  At
$n_{\max} = 4$ (30 spatial orbitals), the graph basis produces
79,465 nonzero ERIs, and any rotation beyond $\theta \sim 10^{-6}$
fills the tensor to near-complete density.

\subsection{$Z$-independence of ERI count}

The nonzero ERI count is identical across all $Z$ at fixed
$n_{\max}$: 79,465 nonzero ERIs at $n_{\max} = 4$ for
$Z = 2, 3, 4, 5, 10$.  This confirms that the sparsity pattern
is purely geometric---determined entirely by angular momentum
conservation (Gaunt selection rules)---and independent of the
radial wavefunctions.  Changing $Z$ rescales the radial functions
$R_{nl}(r) \to Z^{3/2} R_{nl}(Zr)$ but does not alter which
angular integrals vanish."""

new = r"""At $n_{\max} = 2$ (five spatial orbitals: $1s, 2s, 2p_{-1}, 2p_0,
2p_1$), the graph basis produces 107 nonzero ERIs out of
$5^4 = 625$ possible four-index combinations, corresponding to an
ERI density of $17.1\%$.  Rotating the basis by any nonzero angle
fills the ERI tensor: the density rises to $42.7\%$ at
$\theta = 10^{-8}$, $70.2\%$ at $10^{-6}$, $83.0\%$ at $10^{-4}$,
and saturates at $100\%$ ($625/625$) for $\theta \gtrsim 10^{-2}$.

\emph{Corrected 2026-08-29.}  This subsection previously reported
$265/625 = 42.4\%$.  That count came from an evaluator that omitted
the Coulomb selection rule $m_a + m_b = m_c + m_d$, so $158$ of the
$265$ entries ($59.6\%$) were physically zero.  The correction
\emph{strengthens} the sparsity claim and reverses the basis-size
trend reported below; the qualitative content --- that the sparse
point is isolated and any angular-mixing departure is dense --- is
unchanged.

The discontinuity is at the identity itself; the approach to
saturation is rapid but graded, not a step at a finite $\theta$:
\begin{equation}
  D(\theta) =
  \begin{cases}
    0.171,  & \theta = 0 \\
    0.427,  & \theta = 10^{-8} \\
    0.702,  & \theta = 10^{-6} \\
    0.830,  & \theta = 10^{-4} \\
    1.000,  & \theta \gtrsim 10^{-2}.
  \end{cases}
  \label{eq:step}
\end{equation}
(Corrected 2026-08-28: this equation previously reported a
two-branch form with $D = 0.992$ for $\theta > 10^{-6}$.  Neither
value reproduces --- saturation is complete at $625/625$, not
$620/625$.  The qualitative content is unchanged and if anything
stronger.)

\textbf{Generator dependence.}  The two endpoints --- $0.171$ at the
identity and $1.000$ saturation for $\theta \gtrsim 10^{-2}$ --- are
properties of the basis and are reproduced by every random
antisymmetric generator tested.  The three intermediate densities
are \emph{not}:\ across five seeds of the same construction the
$\theta = 10^{-8}$ value spans $0.38$--$0.46$ and the
$\theta = 10^{-4}$ value spans $0.76$--$0.83$ (the quoted spans are
sampling-dependent, not generator-variation limits).  The values
tabulated above are for one representative generator (seed 42) at
the $10^{-10}$ nonzero-entry cutoff used throughout; they should be
read as a profile shape, not as basis constants.  The two endpoints
are cutoff-robust:\ the identity-point zeros are \emph{exactly}
$0.0$ and the smallest nonzero entry is $0.0172$, so the $107/625$
count is independent of any cutoff below that value; saturation is
complete ($625/625$) at every $(\theta, \text{seed})$ tested.

This result is not an artifact of the small basis.  At
$n_{\max} = 4$ (30 spatial orbitals), the graph basis produces
$57{,}700$ nonzero entries out of $30^4 = 810{,}000$, a full-tensor
density of $7.12\%$.  The density therefore \emph{improves} with
basis size, from $17.1\%$ at $n_{\max} = 2$ to $7.12\%$ at
$n_{\max} = 4$.  Counting instead over canonical unique tuples
($p \le q$, $r \le s$) gives $15{,}293$ of $216{,}225$, a density of
$7.07\%$ --- and the same improvement is visible within that
convention alone ($18.2\%$ at $n_{\max} = 2$), so the trend is not
an artifact of mixing the two counts.  Any rotation again fills the
tensor to near-complete density.

\emph{Corrected 2026-08-29.}  The retired text reported $318{,}720$
($39.3\%$) and $79{,}465$ ($36.8\%$) and concluded that the density
was ``essentially flat \dots\ does not degrade with basis size.''
Those counts carried the same missing selection rule as the
$n_{\max} = 2$ figure, and at $n_{\max} = 4$ it inflated them by a
factor of $5.5$.  With the rule imposed the conclusion is not merely
restored but reversed in the framework's favour:\ the density falls
with basis size (empirically $\sim\!M^{-0.49}$ over $M = 5, 14, 30$).

\subsection{$Z$-independence of ERI count}

The nonzero ERI count is identical across all $Z$ at fixed
$n_{\max}$: $15{,}293$ canonical-unique nonzero ERIs at
$n_{\max} = 4$ for $Z = 2, 3, 10$ (and $107$ full-tensor at
$n_{\max} = 2$ for $Z = 2, 6, 10$).  This confirms that the sparsity
pattern is purely geometric---determined entirely by angular momentum
conservation (Gaunt selection rules)---and independent of the
radial wavefunctions.  Changing $Z$ rescales the radial functions
$R_{nl}(r) \to Z^{3/2} R_{nl}(Zr)$ but does not alter which
angular integrals vanish."""

assert s.count(old) == 1, s.count(old)
io.open(P, "w", encoding="utf-8").write(s.replace(old, new))
print("  ok  SS III reconstructed (sparsity transition + Z-independence)")
