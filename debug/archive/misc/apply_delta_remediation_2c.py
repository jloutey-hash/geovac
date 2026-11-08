"""DELTA remediation 2c: the Paper 12 edits still outstanding after 2a/2b.

2a wrote per-path only when every edit for that path matched, so P13 and P15
landed while P12 did not; 2b then applied F1 and F11.  This applies the rest,
against the file's ACTUAL current wrapping (read, not reconstructed).

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P12 = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


# ---- F2a abstract
edit(
    r"""no six-dimensional quadrature and no fitted parameters; the sole
residual numerical integration is a single one-dimensional adaptive
quadrature for the log-singular second-kind Legendre moment $B_l$.""",
    r"""no six-dimensional quadrature and no fitted parameters; in the
$\sigma$ sector the sole residual numerical integration is a single
one-dimensional adaptive quadrature for the log-singular second-kind
Legendre moment $B_l$.  (The $|m| \le 1$ extension reported below
keeps the same kernel but evaluates its ordered-$\xi$ integral by a
two-dimensional spectral quadrature rather than by those recurrences,
so the quadrature-free property is a $\sigma$-sector statement
throughout this paper.)""",
    "F2a: abstract quadrature claim scoped to sigma")

# ---- F2b abstract
edit(
    r"""them in the same basis with the same algebraic $V_{ee}$ gives
$99.09\%$ of $D_e$ at $|m| \le 1$,""",
    r"""them in the same basis with the same Neumann kernel gives
$99.09\%$ of $D_e$ at $|m| \le 1$,""",
    "F2b: abstract 'same algebraic V_ee' -> 'same Neumann kernel'")

# ---- F2c Sec VII.B
edit(
    r"""recurrence relations on quantum number labels, the sole residual
numerical step being a single one-dimensional adaptive quadrature
for the log-singular second-kind moment $B_l$.""",
    r"""recurrence relations on quantum number labels, the sole residual
numerical step in the $\sigma$ sector being a single one-dimensional
adaptive quadrature for the log-singular second-kind moment $B_l$
(Sec.~\ref{sec:azimuthal} states what the $|m| \ge 1$ extension
costs).""",
    "F2c: Sec. VII.B scoped")

# ---- F2d conclusion
edit(
    r"""six-dimensional quadrature for electron-electron repulsion, leaving
only a single one-dimensional quadrature for the log-singular $B_l$
moment.""",
    r"""six-dimensional quadrature for electron-electron repulsion, leaving
only a single one-dimensional quadrature for the log-singular $B_l$
moment in the $\sigma$ sector.""",
    "F2d: conclusion scoped")

# ---- F2e conclusion item 6
edit(
    r"""channels at $|m| \le 1$ in the same basis with the same
    algebraic $V_{ee}$ (Sec.~\ref{sec:azimuthal}), a gain of""",
    r"""channels at $|m| \le 1$ in the same basis with the same
    Neumann kernel (Sec.~\ref{sec:azimuthal}), a gain of""",
    "F2e: conclusion item 6 wording")

# ---- F8 the italicised 'central thesis'
edit(
    r"""86.8\% with $p \leq 2$ basis functions.  This confirms
the central thesis: \emph{the cusp does not need $r_{12}^p$
basis functions---it needs exact $V_{ee}$ computed from
quantum number algebra.}""",
    r"""86.8\% with $p \leq 2$ basis functions.  What that comparison
establishes is narrower than it looks:\ \emph{an exact $V_{ee}$ beats
$r_{12}$ basis functions carrying grid error}.  It says nothing about
the cusp in either direction.  (Earlier versions drew the stronger
conclusion that ``the cusp does not need $r_{12}^p$ basis
functions---it needs exact $V_{ee}$'';\ that goes with the rest of the
withdrawn cusp reading, and
Sec.~\ref{sec:azimuthal_diagnosis} gives the cusp's actual cost here.)""",
    "F8: 'central thesis' cusp misattribution withdrawn")

# ---- F9a Eq. basis_mu
edit(
    r"""  u_{j l \mu}(\xi,\eta,\varphi) =
    \xi^{\,j}\,\eta^{\,l}\,
    (\xi^2-1)^{\mu/2}(1-\eta^2)^{\mu/2}\,e^{-\alpha\xi},
  \label{eq:basis_mu}""",
    r"""  u_{j l \mu}(\xi,\eta) =
    \xi^{\,j}\,\eta^{\,l}\,
    (\xi^2-1)^{\mu/2}(1-\eta^2)^{\mu/2}\,e^{-\alpha\xi},
  \qquad
  \Phi_\mu = u\,u'\,\cos\mu(\varphi_1-\varphi_2),
  \label{eq:basis_mu}""",
    "F9a: basis_mu is a one-electron factor; the phi pairing made explicit")

# ---- F9b 'recovers exactly'
edit(
    r"""Setting $\mu = 0$ recovers
Eq.~\eqref{eq:basis_function} exactly.""",
    r"""At $\mu = 0$ the factor $u$ reduces to the one-electron factor of
Eq.~\eqref{eq:basis_function} and $\Phi_0$ to its symmetrised
product, so the $\mu = 0$ sector is exactly the original basis.""",
    "F9b: 'recovers exactly' type mismatch corrected")

# ---- F10a kept dimensions
edit(
    r"""$(3,2)$ &  46 & 92.37 &  92 & 99.00 \\
$(3,3)$ &  72 & 92.42 & 144 & 99.09 \\""",
    r"""$(3,2)$ &  46 & 92.37 &  92 (76) & 99.00 \\
$(3,3)$ &  72 (65) & 92.42 & 144 (115) & 99.09 \\""",
    "F10a: kept dimensions added to the table")

# ---- F10b threshold
edit(
    r"""scanned, canonical orthogonalization throughout.
  $D_e^{\rm exact} = 0.174475$~Ha.}""",
    r"""scanned, canonical orthogonalization throughout --- eigenvectors of
  $S$ with eigenvalue below $10^{-11}\lambda_{\max}$ are discarded, and
  the surviving dimension is given in parentheses where it differs
  from $N$.  $D_e^{\rm exact} = 0.174475$~Ha.}""",
    "F10b: orthogonalisation threshold stated")

with io.open(P12, encoding="utf-8") as fh:
    text = fh.read()

applied, failed = [], []
for old, new, label in EDITS:
    if old in text:
        text = text.replace(old, new, 1)
        applied.append(label)
    else:
        failed.append(label)

# write whatever matched, then report -- partial application is visible and
# re-runnable, which is better than the all-or-nothing that split 2a.
with io.open(P12, "w", encoding="utf-8") as fh:
    fh.write(text)

print("applied %d edits" % len(applied))
for a in applied:
    print("  +", a)
if failed:
    print("STILL UNMATCHED (%d):" % len(failed))
    for f in failed:
        print("  -", f)
    sys.exit(1)
