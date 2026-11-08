"""STEP 5: the two-way direction -- where the backing proves MORE than the prose.

The QA gate returns a two-way verdict by design, and this run produced two
upgrades on the results the whole week rests on.  Both were found by the code
reviewer re-deriving the claims independently rather than checking the tests
pass.

U1.  eq:W_diagonal is stated as "verified entrywise to 5e-11" -- a MEASUREMENT
sitting inside an [INTERNAL THEOREM] block.  It is an exact identity, provable in
three cases, and the middle one is the substantive one:  hermiticity of (T - E)
applied to the two Sturmian equations gives (Q_mu - Q_nu) W_{mu,nu} = 0, so the
off-diagonal vanishes whenever the weighted charges differ.  The reviewer
confirmed it symbolically (sympy, closed-form hydrogenic integrals, two charges,
no grid:  exact 0).

The part worth stating plainly:  this is the SAME generalized-eigenvalue
orthogonality the paper already invokes molecularly ("automatic for eigenvectors
of a common Sturmian problem with distinct beta").  The atomic case is that
statement's own special case and was never written down -- the paper measured
numerically what it had already argued structurally eighty lines later.

The degenerate branch was checked rather than assumed:  R_mu = R_nu at the same
l first collides at n_max = 35 (none at n_max <= 17, the largest computed), and
even there W = 0, because equal R forces equal charge and the two configurations'
n-multisets are then necessarily disjoint.  The identity is unconditional.

U2.  The variational bound is asserted via an intermediate the paper never
proves and does not need ("E_iso is the LOWEST root").  Given W diagonal,

    H(lambda) + (lambda^2/2) S = lambda * (lambda * 1 - M)    exactly,

and S is positive definite, so by Sylvester's law of inertia the number of pencil
roots strictly below -lambda^2/2 equals #{k : lambda_k(M) > lambda}.  At
lambda = lambda_max(M) that count is 0 -- the bound -- and at lambda = lambda_k
it is exactly k, which ALSO proves the root-by-root correspondence that
test_c4 currently only asserts numerically and that every excited-state number
in Sec.4 depends on.  Three lines replace a sampled claim, and prove more.

CODE-M3 is folded in here because it is the same equation:  the sentence
declares Q_nu = lambda/R_nu and then writes W = R_nu delta, but at that charge
W_{nu,nu} = Q_nu R_nu^2 = lambda R_nu.  The next clause says "at unit scale" for
T, so the intent is clear;  the labelled equation just never carried the
qualifier.  Verified: W/R_nu = 1.000000 / 2.000000 / 2.397696 at lambda = 1 / 2 /
p_kappa.
"""
import io

PAP = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
s = io.open(PAP, encoding="utf-8").read()

# ---- U1 + CODE-M3: state the scale convention, and prove the identity.
OLD = ("Because each $\\Phi_\\nu$ is hydrogenic at $Q_\\nu=\\lambda/R_\\nu$, the\n"
       "one-body Coulomb metric is exactly diagonal,\n"
       "\\begin{equation}\n"
       "  W_{\\mu\\nu}\\equiv\\Big\\langle\\Phi_\\mu\\Big|\\sum_j r_j^{-1}\\Big|\\Phi_\\nu\\Big\\rangle\n"
       "    = R_\\nu\\,\\delta_{\\mu\\nu},\n"
       "  \\label{eq:W_diagonal}\n"
       "\\end{equation}\n"
       "verified entrywise to $5\\times10^{-11}$.")
NEW = ("Because each $\\Phi_\\nu$ is hydrogenic at $Q_\\nu=\\lambda/R_\\nu$, the\n"
       "one-body Coulomb metric is exactly diagonal;\\ written at unit scale\n"
       "($\\lambda=1$, the convention used for $T$ below, so that at general $\\lambda$\n"
       "the diagonal reads $\\lambda R_\\nu$),\n"
       "\\begin{equation}\n"
       "  W_{\\mu\\nu}\\equiv\\Big\\langle\\Phi_\\mu\\Big|\\sum_j r_j^{-1}\\Big|\\Phi_\\nu\\Big\\rangle\n"
       "    = R_\\nu\\,\\delta_{\\mu\\nu}.\n"
       "  \\label{eq:W_diagonal}\n"
       "\\end{equation}\n"
       "\\textbf{[SYMBOLIC]} This is exact, not numerical.  The diagonal is\n"
       "$W_{\\nu\\nu}=Q_\\nu\\sum_j n_j^{-2}=Q_\\nu R_\\nu^{2}=R_\\nu$ from\n"
       "$\\langle 1/r\\rangle_{n}=Q/n^{2}$.  Off the diagonal there are three cases.  If\n"
       "$l_\\mu\\neq l_\\nu$ the integral vanishes by angular orthogonality.  If the\n"
       "$l$ agree and $R_\\mu\\neq R_\\nu$, apply $(T-E)$ to the two Sturmian equations\n"
       "and use its hermiticity:\\ the cross terms give\n"
       "$(Q_\\mu-Q_\\nu)\\,W_{\\mu\\nu}=0$, so $W_{\\mu\\nu}=0$ whenever the weighted\n"
       "charges differ.  The remaining case, $R_\\mu=R_\\nu$ at equal $l$, does not\n"
       "occur below $n_{\\max}=35$ (the largest basis computed here is\n"
       "$n_{\\max}=17$), and when it does occur $W$ still vanishes:\\ equal $R$ forces\n"
       "equal charge, and two distinct configurations sharing a charge have disjoint\n"
       "$n$-multisets, so the radial factors are orthogonal.  The identity is\n"
       "therefore unconditional on this family, and the entrywise agreement to\n"
       "$5\\times10^{-11}$ is a check on the implementation rather than the evidence\n"
       "for the claim.  \\emph{This is the same potential-weighted orthogonality that\n"
       "Sec.~\\ref{sec:molecular} invokes for the molecular problem --- automatic for\n"
       "eigenvectors of a common Sturmian problem with distinct $\\beta$ --- read back\n"
       "in its atomic special case.}")
assert OLD in s, "U1 locus not found"
s = s.replace(OLD, NEW, 1)

# ---- U2: replace the unproven intermediate with the inertia argument.
OLD2 = ("The variational bound is automatic rather than\n"
        "fortunate --- $E_{\\rm iso}$ is the lowest root of $H(p_\\kappa)C=E\\,SC$, so no\n"
        "point can fall below the exact value.")
NEW2 = ("The variational bound is automatic rather than\n"
        "fortunate, and the argument gives more than the bound.  Because $W$ is\n"
        "diagonal, $H(\\lambda)+\\tfrac12\\lambda^{2}S=\\lambda\\,(\\lambda\\mathbb{1}-M)$\n"
        "\\emph{exactly};\\ $S$ is positive definite, so by Sylvester's law of inertia\n"
        "the number of roots of the pencil lying strictly below $-\\tfrac12\\lambda^{2}$\n"
        "equals $\\#\\{k:\\lambda_k(M)>\\lambda\\}$.  Taking\n"
        "$\\lambda=\\lambda_{\\max}(M)$ that count is zero, which is the bound;\\ taking\n"
        "$\\lambda=\\lambda_k$ it is exactly $k$, which is the root-by-root\n"
        "correspondence the excited-state results of this section rely on.  Neither\n"
        "step needs $E_{\\rm iso}$ to be identified as the lowest root in advance.")
assert OLD2 in s, "U2 locus not found"
s = s.replace(OLD2, NEW2, 1)

io.open(PAP, "w", encoding="utf-8").write(s)
print("U1: eq:W_diagonal MEASURED -> [SYMBOLIC], three-case proof, scale convention stated")
print("U2: variational bound -> Sylvester inertia; root-by-root correspondence now proved")
