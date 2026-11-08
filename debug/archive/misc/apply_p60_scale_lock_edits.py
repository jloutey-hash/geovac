"""Apply the Paper-60 scale-lock corrections: registry entries + the He passage.

Scope (PI-approved 2026-09-08): fix the -2.90250 citation and register the
numbers it now cites.  Per CLAUDE.md Sec.15 rule 2 the surrounding sentence is
re-read and corrected too -- it carried the mechanism claim ("genuine basis
incompleteness") that the span-vs-posing diagnostic refutes.

NOT in scope here: the Sec.4 floor paragraph ("cannot reach chemical accuracy at
any K", the [OPEN] high-l question).  Held pending the excited-state ladder.
"""
import io
import re

REG = "debug/qa/numeric_registry.py"
PAP = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"

# ---------------------------------------------------------------- registry
NEW_ENTRIES = '''    "p60_span_deficit_spdf": dict(
        value=1.28, convention="constant: mHa above the exact He ground state "
                               "reached by a VARIATIONAL CI over the identical "
                               "Goscinskian span (l_max=3, K=130) with the global "
                               "scale lambda optimized -- the span's own deficit, "
                               "with the isoenergetic scale-lock removed",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_scale_scan.py. Companion "
                   "locked-scale value at the same K is 7.46 mHa. s-sector "
                   "counterpart (against the known exact s-limit -2.879029 Ha) is "
                   "0.15 mHa at K=136 vs 4.43 locked. Pipeline unit-tested at K=1, "
                   "where both postings coincide and return -2.8476562 = "
                   "-(2-5/16)^2 exactly.",
        aliases={1.640: "l_max=3, K=100", 2.232: "l_max=3, K=74"}),
    "p60_posing_cost_ground": dict(
        value=4.21, convention="constant: mHa, E_iso - min_lambda E_var over the "
                               "SAME span, He ground state, s-only, K=105. "
                               "Reference-free: needs no known limit",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_posing_cost_by_state.py, "
                   "lambda by GRID scan (Brent found a local minimum at nmax=4 and "
                   "reported a NEGATIVE cost, which the variational bound forbids). "
                   "Grows with K: 3.52 (K=36), 4.13 (K=78), 4.21 (K=105).",
        aliases={3.524: "K=36", 4.125: "K=78"}),
    "p60_posing_cost_exc": dict(
        value=0.98, convention="constant: mHa, same quantity as "
                               "p60_posing_cost_ground but for He 2^1S (the second "
                               "root of M), s-only, K=105",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_posing_cost_by_state.py. "
                   "4.3x SMALLER than the ground state and the ratio widens with K "
                   "(0.681/3.524 at K=36 -> 0.983/4.212 at K=105). This is the "
                   "measured form of Avery's split-shell mechanism: a Goscinskian "
                   "1s^2 configuration pins both electrons to one exponent, an "
                   "excited configuration gets two free from n_a != n_b.",
        aliases={0.681: "K=36", 0.885: "K=78"}),
'''

src = io.open(REG, encoding="utf-8").read()
anchor = '    "p60_cond_S_converged": dict('
assert anchor in src, "registry anchor not found"
assert "p60_posing_cost_ground" not in src, "entries already present"
src = src.replace(anchor, NEW_ENTRIES + anchor, 1)
io.open(REG, "w", encoding="utf-8").write(src)
print("registry: 3 entries added before p60_cond_S_converged")

# ---------------------------------------------------------------- paper
OLD = r"""The residual $\sim\!7$~mHa is genuine basis incompleteness, not a missing metric:
the Goscinskian basis is deliberately poor for the helium ground state---Avery and
Avery reach $-2.90250$ with $102$ optimized Coulomb--Sturmian
configurations~\cite{avery2006}, still $1.2$~mHa short---and every point stays
above the exact value (no variational overshoot).  The machinery is correct,
including the higher-$\ell$ angular correlation (Gaunt / Wigner-$3j$ multipole
coupling)."""

NEW = r"""The residual is \emph{not} a deficiency of the configuration span.  A
variational CI over the \emph{identical} span, with the global scale optimized,
reaches $\gvq{p60_span_deficit_spdf}{1.28}$~mHa of the exact energy at $K=130$,
where the isoenergetic solution of that same span stalls at $7.46$~mHa;  in the
$s$-sector, where the exact $s$-limit $-2.879029$~Ha is known independently, the
same comparison is $0.15$~mHa against $4.43$~mHa at $K=136$.  The residual is
therefore the price of the \emph{posing}, not of the basis:  the metric-free form
requires the basis scale to be locked to the eigenvalue, $\lambda=p_\kappa$, and
that lock is not the variational optimum.  \textbf{[INTERNAL THEOREM]} The same
identity makes the bound automatic rather than fortunate---the isoenergetic root
is the lowest root of the fixed-scale generalized problem $H(p_\kappa)C=E\,SC$,
so every point lies above the exact value by construction, not by observation.

\textbf{[ESTABLISHED, from Avery]} The mechanism is Avery's, not ours.  The
helium ground state needs \emph{in-out} radial correlation, which variational
treatments supply through split-shell ($1s\,1s'$) functions carrying two
independent exponents;  a Goscinskian $1s^{2}$ configuration places both
electrons at $n=1$ and therefore at the identical exponent $Q_\nu/n_j$---for
\emph{any} weighting potential, since a separable $V_0$ still assigns one
$\beta_\nu$ per configuration.  Excited configurations obtain two scales free
from $n_a\neq n_b$, and the state-dependence this predicts is measurable:  the
posing cost is $\gvq{p60_posing_cost_ground}{4.21}$~mHa for the ground state
against $\gvq{p60_posing_cost_exc}{0.98}$~mHa for $2\,^{1}S$ at $K=105$, and the
ratio widens with $K$.  The floor is a ground-state pathology rather than a
property of the method.

\textbf{[OBSERVATION]} A comparison to the $-2.90250$~Ha reached with $102$
optimized Coulomb--Sturmian configurations~\cite{avery2006} is \emph{withdrawn}
as a characterization of the present method:  that figure ($1.2$~mHa short)
matches our scale-optimized ladder ($1.64$~mHa short at $K=100$) rather than our
locked-scale one ($7.70$~mHa at the same $K$), and which posing produced it could
not be confirmed against the primary source.  The machinery is correct, including
the higher-$\ell$ angular correlation (Gaunt / Wigner-$3j$ multipole coupling)."""

tex = io.open(PAP, encoding="utf-8").read()
assert OLD in tex, "paper passage not found verbatim"
tex = tex.replace(OLD, NEW, 1)
io.open(PAP, "w", encoding="utf-8").write(tex)
print("paper: He-ladder passage replaced (%d -> %d chars)" % (len(OLD), len(NEW)))
