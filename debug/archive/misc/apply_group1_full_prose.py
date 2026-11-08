r"""group1 FULL-cert remediation, part 1: paper + synthesis prose.

All findings verified against primary text before writing. Every fix moves a
claim toward its true tier (the honest-scope direction); none asserts anything
new.

  P42 M1  L481 heading "Propinquity convergence to the continuum" contradicts
          its own body ("state-space GH") and mis-credits the withdrawn Paper 38
          L5 propinquity leg. -> "State-space GH convergence to the continuum".
  P43 M2  L331-337 "The Lorentzian content ... literally satisfied at (3,1) ...
          not merely a structural correspondence" -- un-swept pre-descope voice;
          the abstract (L150-154) and Paper 42 Sec.10 settle it signature-blind.
          Add the same inline caveat Paper 42 Sec.10 carries.
  P43 eps caveat: the numerical-verification paragraph states residual
          ~ O(sqrt(dim) eps_machine) ~ 1e-16; the nmax=1 table rows read 1e-32,
          which invited a "fabricated?" flag. The code reviewer independently
          derived the mechanism: a degenerate single-eigenvalue wedge gives a
          SCALAR modular unitary, so the O(eps) term cancels and the residual is
          O(eps^2). One-line caveat added.
  P53 M1  Theorem 5.6 (thm:plane_propinquity) asserts convergence in the
          strictly stronger Latremoliere pointed/proper propinquity via
          max(reach, height) -> 0, while its own governing remark
          (rem:height_constant) says the height leg is L_s>1, Lambda-independent,
          and "has not been done here". Scope the theorem to what is established
          (state-space GH) and mark the Latremoliere form conditional -- matching
          the paper's own remark.
  Synth M1 L1474 "The inner arrow is propinquity convergence (Paper 45/46 main
          theorem applied)" -- present tense for a descoped arrow; the Status
          note 20 lines above says it is descoped. Add the inline caveat.
  Synth M2 L1516 "analytically AND empirically established" -- contradicts Sec.6's
          own "automatic formula identity, not an operator-level falsifier".
          Reword to the spatial-rate-formula reading.
  Retired-label NITs (state-space GH, not "propinquity", for GeoVac's OWN
          Paper 38/39/40 results): P42 open-question loci, P44 bucket, synthesis
          L1565.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P42 = "papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex"
P43 = "papers/group1_operator_algebras/paper_43_lorentzian_extension.tex"
P44 = "papers/group1_operator_algebras/paper_44_lorentzian_operator_system.tex"
P53 = "papers/group1_operator_algebras/paper_53_disk_propinquity.tex"
SY = "papers/synthesis/group1_operator_algebras_synthesis.tex"
EDITS = []


def edit(path, old, new, label):
    EDITS.append((path, old, new, label))


# ---- P42 M1: the section heading -----------------------------------------
edit(P42,
     r"""\subsection{Propinquity convergence to the continuum}""",
     r"""\subsection{State-space GH convergence to the continuum}""",
     "P42 M1: heading -> state-space GH (matches its own body)")

# ---- P43 M2: the pre-descope "Lorentzian content" framing -----------------
edit(P43,
     r"""\emph{Lorentzian} content of Theorem~\ref{thm:main_intro} is that this
identification persists through the Krein-space lift and the
BBB Connes axiom audit at $(4, 6)$ — the four-witness theorem, which
is a Lorentzian statement in the Wightman axiomatic framework, has its
operator-system analog literally satisfied at signature $(3, 1)$ at
finite cutoff, not merely as a structural correspondence via the
Wick-rotation chain.""",
     r"""\emph{Lorentzian} content of Theorem~\ref{thm:main_intro} is that this
identification persists through the Krein-space lift and the
BBB Connes axiom audit at $(4, 6)$ — the four-witness theorem, which
is a Lorentzian statement in the Wightman axiomatic framework, has its
operator-system analog satisfied on a signature-$(3, 1)$-labeled Krein
substrate at finite cutoff.  \emph{This is signature-blind:}\ as
Paper~42~\S 10 records for the same construction, the finite-cutoff
closure carries no metric statement that distinguishes the signature,
and the identification is a Krein-substrate carrier choice rather than
distinguishing Lorentzian content.  The construction is genuine;\ what
it does not yet supply is a Lorentzian metric statement (open, see
Paper~45).""",
     "P43 M2: add the signature-blind caveat (matches abstract + Paper 42 Sec.10)")

# ---- P43 eps caveat -------------------------------------------------------
edit(P43,
     r"""\varepsilon_{\mathrm{machine}})$. We have computed the maximum""",
     r"""\varepsilon_{\mathrm{machine}})$ for the non-degenerate wedge
($n_{\max}\ge 2$).  At $n_{\max}=1$ the wedge generator has a single
degenerate eigenvalue, so the modular unitary $e^{i2\pi K_\alpha}$ is a
scalar phase and its conjugation cancels the $O(\varepsilon)$ term,
leaving an $O(\varepsilon^{2})\sim10^{-32}$ residual;\ this is why the
$n_{\max}=1$ rows below sit far below the $10^{-16}$ floor.  We have
computed the maximum""",
     "P43 eps: caveat for the sub-epsilon nmax=1 residual (degenerate scalar unitary)")

# ---- P53 M1: scope the theorem to match its own governing remark ----------
edit(P53,
     r"""\begin{theorem}[Plane propinquity convergence]
\label{thm:plane_propinquity}
The band-limited truncated triples $\Tcal_{\Lambda}=(P_{\Lambda}C_{0}(\R^{2}
_{\alpha})P_{\Lambda},\,L,\,0)$ converge to the plane triple
$\Tcal_{\infty}=(C_{0}(\R^{2}_{\alpha}),\,L,\,0)$ in the Latr\'emoli\`ere
pointed/proper propinquity, with
\[
\Lprop(\Tcal_{\Lambda},\Tcal_{\infty})\;\le\;C\,\gamma_{\Lambda}\;\to\;0,
\]
where $\gamma_{\Lambda}$ is the Bochner--Riesz reconstruction rate of
Theorem~\ref{thm:plane} (measured $\gamma_{\Lambda}\sim\Lambda^{-0.88}$ for
a Gaussian observable).
\end{theorem}""",
     r"""\begin{theorem}[Plane convergence]
\label{thm:plane_propinquity}
The band-limited truncated triples $\Tcal_{\Lambda}=(P_{\Lambda}C_{0}(\R^{2}
_{\alpha})P_{\Lambda},\,L,\,0)$ converge to the plane triple
$\Tcal_{\infty}=(C_{0}(\R^{2}_{\alpha}),\,L,\,0)$ in the van~Suijlekom
state-space Gromov--Hausdorff distance at the qualitative rate
$\gamma_{\Lambda}$ of Theorem~\ref{thm:plane}
(measured $\gamma_{\Lambda}\sim\Lambda^{-0.88}$ for a Gaussian observable).
Convergence in the strictly stronger Latr\'emoli\`ere pointed/proper
propinquity, $\Lprop \le C\gamma_{\Lambda}\to 0$, holds
\emph{conditionally} on the extent rebuild of
Remark~\ref{rem:height_constant}:\ the height leg carries the
$\Lambda$-independent Lebesgue constant $\Lebesgue_{s}>1$, which does not
vanish, so the pointed/proper statement is not established here (see the
remark).
\end{theorem}""",
     "P53 M1: theorem scoped to state-space GH; Latremoliere form marked conditional")

# ---- Synthesis M1: the inner-arrow present-tense descoped clause ----------
edit(SY,
     r"""The inner arrow is propinquity convergence at
coupled cells (Paper~45/46 main theorem applied along the scaling);
the outer arrow is norm-resolvent convergence of the Lorentzian Dirac""",
     r"""The inner arrow was proposed as propinquity convergence at
coupled cells (Paper~45/46 main theorem applied along the scaling), but
that arrow is \emph{descoped} by Paper~45's degeneracy theorem;\
the outer arrow is norm-resolvent convergence of the Lorentzian Dirac""",
     "Synth M1: inner arrow marked descoped inline")

# ---- Synthesis M2: "established" -> spatial-rate-formula reading -----------
edit(SY,
     r"""upgrades the L3c-$\alpha$ theorem from ``analytical corollary of
Paper~45'' to ``analytically AND empirically established.''""",
     r"""confirms the $T$- and $N_t$-independence of the surviving
\emph{spatial} rate formula.  Under the degeneracy analysis this is an
automatic formula identity (Sec.~\ref{sec:degeneracy_reading}), not an
operator-level convergence, so it does not upgrade the descoped
L3c-$\alpha$ panel to an established Lorentzian result.""",
     "Synth M2: reworded to the spatial-rate-formula reading (matches Sec.6)")

# ---- Retired-label NITs (GeoVac's own results are state-space GH) ---------
edit(P44,
     r"""\item[(a)] \textbf{Riemannian propinquity convergence on spectral""",
     r"""\item[(a)] \textbf{Riemannian state-space GH convergence on spectral""",
     "P44 NIT: bucket label -> state-space GH")

edit(SY,
     r"""k-fold extension of the Paper~39 two-factor tensor-product propinquity""",
     r"""k-fold extension of the Paper~39 two-factor tensor-product state-space GH bound""",
     "Synth NIT: L1565 propinquity -> state-space GH bound")

by_path = {}
for path, old, new, label in EDITS:
    by_path.setdefault(path, []).append((old, new, label))

applied, failed = [], []
for path, items in by_path.items():
    with io.open(path, encoding="utf-8") as fh:
        t = fh.read()
    for old, new, label in items:
        if old in t:
            t = t.replace(old, new, 1)
            applied.append(label)
        else:
            failed.append("%s   [%s]" % (label, path))
    with io.open(path, "w", encoding="utf-8") as fh:
        fh.write(t)

print("applied %d of %d" % (len(applied), len(EDITS)))
for a in applied:
    print("  +", a)
if failed:
    print("")
    print("UNMATCHED (%d):" % len(failed))
    for f in failed:
        print("  -", f)
    sys.exit(1)
