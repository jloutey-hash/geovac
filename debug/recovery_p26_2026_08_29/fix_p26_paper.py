r"""Paper 26 Sec III: correct the ERI sparsity numbers.

RAW STRINGS ONLY -- LaTeX-bearing edit (project hard rule: never through a
bash heredoc, which halves backslashes in this harness).

Corrected (measured 2026-08-29):
  n_max=2 full-tensor      107/625      = 17.1%   (was 265/625 = 42.4%)
  n_max=2 canonical-unique  41/225      = 18.2%
  theta profile (seed 42)  0.4272 / 0.7024 / 0.8304 / 1.000
  spans  theta=1e-8 0.38-0.46 ; theta=1e-4 0.76-0.83
  smallest nonzero |ERI|   0.017163     (was 0.0122)
  n_max=4 full-tensor    57,700/810,000 = 7.12%   (was 318,720 = 39.3%)
  n_max=4 canonical      15,293/216,225 = 7.07%   (was 79,465 = 36.8%)
"""
import io

P = "papers/group6_precision_observations/paper_26_entanglement.tex"
s = io.open(P, encoding="utf-8").read()
n = 0


def rep(old, new, label):
    global s, n
    assert s.count(old) == 1, f"{label}: {s.count(old)} matches"
    s = s.replace(old, new)
    n += 1
    print(f"  ok  {label}")


rep(r"$42.4\%\!\to\!100\%$ basis-intrinsic sparsity transition (with a",
    r"$17.1\%\!\to\!100\%$ basis-intrinsic sparsity transition (with a",
    "intro transition figure")

rep(r"""2p_1$), the graph basis produces 265 nonzero ERIs out of
$5^4 = 625$ possible four-index combinations, corresponding to an
ERI density of $42.4\%$.  Rotating the basis by any nonzero angle
fills the ERI tensor: the density rises to $72.2\%$ at
$\theta = 10^{-8}$, $91.4\%$ at $10^{-6}$, $97.8\%$ at $10^{-4}$, and
saturates at $100\%$ ($625/625$) for $\theta \gtrsim 10^{-2}$.""",
    r"""2p_1$), the graph basis produces 107 nonzero ERIs out of
$5^4 = 625$ possible four-index combinations, corresponding to an
ERI density of $17.1\%$.  Rotating the basis by any nonzero angle
fills the ERI tensor: the density rises to $42.7\%$ at
$\theta = 10^{-8}$, $70.2\%$ at $10^{-6}$, $83.0\%$ at $10^{-4}$, and
saturates at $100\%$ ($625/625$) for $\theta \gtrsim 10^{-2}$.

\emph{Corrected 2026-08-29.}  This subsection previously reported
$265/625 = 42.4\%$.  That count came from an evaluator that omitted the
Coulomb selection rule $m_a + m_b = m_c + m_d$, so $158$ of the $265$
entries ($59.6\%$) were physically zero.  The correction \emph{strengthens}
the sparsity claim and reverses the basis-size trend reported below; the
qualitative content --- that the sparse point is isolated and any
angular-mixing departure is dense --- is unchanged.""",
    "n_max=2 count + profile")

rep(r"""    0.424,  & \theta = 0 \\
    0.722,  & \theta = 10^{-8} \\
    0.914,  & \theta = 10^{-6} \\
    0.978,  & \theta = 10^{-4} \\""",
    r"""    0.171,  & \theta = 0 \\
    0.427,  & \theta = 10^{-8} \\
    0.702,  & \theta = 10^{-6} \\
    0.830,  & \theta = 10^{-4} \\""",
    "eq:step values")

rep(r"""\textbf{Generator dependence.}  The two endpoints --- $0.424$ at the
identity and $1.000$ saturation""",
    r"""\textbf{Generator dependence.}  The two endpoints --- $0.171$ at the
identity and $1.000$ saturation""",
    "endpoint value")

rep(r"""$\theta = 10^{-8}$ value spans $0.72$--$0.84$ and the $\theta = 10^{-4}$
value spans $0.95$--$0.98$ (a 20-seed sweep widens these to
$0.65$--$0.84$ and $0.94$--$0.98$; the quoted spans are
sampling-dependent, not generator-variation limits).""",
    r"""$\theta = 10^{-8}$ value spans $0.38$--$0.46$ and the $\theta = 10^{-4}$
value spans $0.76$--$0.83$ (the quoted spans are
sampling-dependent, not generator-variation limits).""",
    "seed spans")

rep(r"""nonzero entry is $0.0122$, so the $265/625$ count is independent of any
cutoff below that value;""",
    r"""nonzero entry is $0.0172$, so the $107/625$ count is independent of any
cutoff below that value;""",
    "cutoff robustness")

rep(r"""This result is not an artifact of the small basis.  At
$n_{\max} = 4$ (30 spatial orbitals), the graph basis produces
$318{,}720$ nonzero entries out of $30^4 = 810{,}000$, a full-tensor
density of $39.3\%$ --- essentially flat against the $42.4\%$ at
$n_{\max} = 2$, so the sparsity does not degrade with basis size.
Counting instead over canonical unique tuples ($p \le q$, $r \le s$)
gives $79{,}465$ of $216{,}225$, a density of $36.8\%$ in that
convention.  Both conventions are quoted here because the $79{,}465$
figure is canonical-unique while the $42.4\%$ above is full-tensor;
comparing them directly without the totals invites a spurious
impression that the density improves with basis size.  Any rotation
again fills the tensor to near-complete density.""",
    r"""This result is not an artifact of the small basis.  At
$n_{\max} = 4$ (30 spatial orbitals), the graph basis produces
$57{,}700$ nonzero entries out of $30^4 = 810{,}000$, a full-tensor
density of $7.12\%$.  The density therefore \emph{improves} with basis
size, from $17.1\%$ at $n_{\max} = 2$ to $7.12\%$ at $n_{\max} = 4$.
Counting instead over canonical unique tuples ($p \le q$, $r \le s$)
gives $15{,}293$ of $216{,}225$, a density of $7.07\%$ --- and the same
improvement is visible within that convention alone ($18.2\%$ at
$n_{\max} = 2$), so the trend is not an artifact of mixing the two
counts.  Any rotation again fills the tensor to near-complete density.

\emph{Corrected 2026-08-29.}  The retired text reported $318{,}720$
($39.3\%$) and $79{,}465$ ($36.8\%$) and concluded that the density was
``essentially flat \dots\ does not degrade with basis size.''  Those
counts carried the same missing selection rule as the $n_{\max} = 2$
figure, and at $n_{\max} = 4$ it inflated them by a factor of $5.5$.
With the rule imposed the conclusion is not merely restored but
reversed in the framework's favour: the density falls with basis size
(empirically $\sim\!M^{-0.49}$ over $M = 5, 14, 30$).""",
    "n_max=4 counts + trend reversal")

rep(r"""$n_{\max}$: $79{,}465$ canonical-unique nonzero ERIs at $n_{\max} = 4$ for
$Z = 2, 3, 4, 5, 10$.""",
    r"""$n_{\max}$: $15{,}293$ canonical-unique nonzero ERIs at $n_{\max} = 4$ for
$Z = 2, 3, 10$ (and $107$ full-tensor at $n_{\max} = 2$ for
$Z = 2, 6, 10$).""",
    "Z-independence counts")

io.open(P, "w", encoding="utf-8").write(s)
print(f"\n{n} edits applied to {P}")
