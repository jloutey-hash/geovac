"""STEP 4: the prose findings -- four classes, swept claim-wide rather than locus
by locus (the habit that cost this corpus twelve avoidable passes in the v5.4.4 arc).

CLASS 1 (claims M2 = synthesis F1).  "reverses to superlinear".  Two dimensions
found this independently from opposite ends -- synthesis at the citer, claims at
the owner's abstract.  The body measures TWO objects and the bridging sentence
merges them:  under full-shell growth the pure-number block reaches K^1.07, but
the TOTAL exponent only rises 0.867 -> 0.911.  0.911 < 1:  no measured full-shell
point is superlinear.  The conclusion (L989) already says only "reverses under
full-shell growth" and is the model.

CLASS 2 (claims M1).  "ill-conditioned" survives at six loci -- abstract, three
in the body, one in the molecular section, one in the synthesis -- after Sec.2
withdraws exactly that characterization ("ordinary, not ill-conditioned").  The
sharpest instance calls cond 5.8 ill-conditioned while Sec.2 calls cond 56.3
ordinary:  same paper, factor of ten, opposite verdicts.  The surviving
mechanism is not conditioning at all -- it is that Loewdin's S^{-1/2} is DENSE,
and that the higher basis members grow progressively linearly dependent so the
coefficient distribution spreads.  That is what the prose should say, because it
is what the paper measured.

CLASS 3 (claims M4).  The abstract calls -2.873 Ha "the s-sector limit".  It is
the locked-posing s-sector LADDER VALUE;  the exact s-limit is -2.879029 and the
body says so twice.  The 6.0 mHa between them is precisely what the scale-lock
argument turns on, so the mislabel erases the paper's own evidence.

CLASS 4 (claims M7, the UNDERCLAIM).  The Acknowledgments say the paper "adds
only its reading as a quantum-encoding target and the resource measurements" --
handing four of its own theorems to Avery, two of which the DoD's C8.13 names
explicitly:  "the identification of that matrix AS the V0-weighted overlap is
ours, and the variational bound for the fixed-scale metric-free problem is not
in his canon.  Crediting either to Avery = MATERIAL."  The sentence was correct
when softened on 2026-08-18 and went stale when this week's results landed --
the strengthening-mirror class Sec.9 names as NOT covered by the retraction
machinery, so nothing could have caught it but a reader.
"""
import io

PAP = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
SYN = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"
REG = "docs/claims_register.md"


def patch(path, pairs, label):
    s = io.open(path, encoding="utf-8").read()
    for old, new in pairs:
        assert old in s, "%s: not found: %.65s" % (path, old.replace("\n", " "))
        s = s.replace(old, new, 1)
    io.open(path, "w", encoding="utf-8").write(s)
    print("  %-46s %s (%d)" % (label, path.split("/")[-1], len(pairs)))


print("CLASS 1 -- 'reverses to superlinear' (owner + citer + register):")
patch(PAP, [
 # abstract
 ("holds when $l_{\\max}$ is fixed and higher-$n$ configurations are added, and\n"
  "reverses to superlinear when the basis is grown as full hydrogenic shells.",
  "holds when $l_{\\max}$ is fixed and higher-$n$ configurations are added, and\n"
  "weakens when the basis is grown as full hydrogenic shells --- there the total\n"
  "exponent rises toward $1$ ($0.867\\to0.911$) while remaining below it, and only\n"
  "the pure-number block $T'$ itself goes superlinear."),
 # the bridging sentence that manufactures the merge
 ("construction gives $\\|T'\\|_1\\sim K^{1.07}$ and a total exponent that\n"
  "\\emph{rises} with $K$ ($0.867\\to0.911$ over $K=56$--$220$).  Superlinear\n"
  "behaviour and a rising exponent are therefore real --- for a different\n"
  "growth rule.",
  "construction gives $\\|T'\\|_1\\sim K^{1.07}$ --- superlinear --- while the\n"
  "\\emph{total} exponent rises but stays below $1$ ($0.867\\to0.911$ over\n"
  "$K=56$--$220$).  A superlinear block and a rising total are therefore both\n"
  "real for a different growth rule;\\ but no measured full-shell point makes the\n"
  "total itself superlinear, and we do not claim one."),
], "abstract + bridging sentence")

patch(SYN, [
 ("not to the construction}:\\ it holds at fixed $l_{\\max}$ and reverses to\n"
  "$K^{1.07}$ under full-shell growth, and the cheap rule carries an accuracy",
  "not to the construction}:\\ it holds at fixed $l_{\\max}$ and weakens under\n"
  "full-shell growth, where the pure-number block reaches $K^{1.07}$ while the\n"
  "total exponent rises only to $0.911$ --- still sublinear.  The cheap rule\n"
  "carries an accuracy"),
], "citer")

patch(REG, [
 ("reverses to K¹·⁰⁷ under full-shell growth",
  "weakens under full-shell growth (the pure-number block reaches K¹·⁰⁷ while the "
  "total exponent rises only to 0.911 — still sublinear)"),
], "claims register row 27")

print("\nCLASS 2 -- the 'ill-conditioned' cluster (6 loci):")
patch(PAP, [
 ("hydrogenic basis), driven by the growing ill-conditioning of the overlap.",
  "hydrogenic basis), driven by the dense $S^{-1/2}$ that L\\\"owdin\n"
  "orthogonalization introduces as the basis members grow progressively linearly\n"
  "dependent."),
 ("($3.0,5.8,13.9,32.2$ for $N=2,3,5,8$), and the dense $S^{-1/2}$ built from an\n"
  "ill-conditioned Gram inflates $\\lambda$.",
  "($3.0,5.8,13.9,32.2$ for $N=2,3,5,8$ --- growing, but modest), and the dense\n"
  "$S^{-1/2}$ built from that increasingly non-orthogonal Gram inflates\n"
  "$\\lambda$.  The inflation is driven by the \\emph{density} of $S^{-1/2}$ and the\n"
  "spread it induces in the coefficient distribution, not by numerical\n"
  "instability:\\ a condition number of $32$ is not ill-conditioned."),
 ("the ill-conditioned matrix of Sec.~\\ref{sec:obstruction}.",
  "the non-orthogonal matrix of Sec.~\\ref{sec:obstruction}."),
 ("The\nill-conditioning of Sec.~\\ref{sec:obstruction} is not repaired, it is\n\\emph{never present}.",
  "The\nnon-orthogonality of Sec.~\\ref{sec:obstruction} is not repaired, it is\n\\emph{never present}."),
 ("error)---where the $L^2$ overlap's intra-center block is ill-conditioned (cond\n$5.8$, itself growing with basis)",
  "error)---where the $L^2$ overlap's intra-center block is merely non-orthogonal\n"
  "(cond $5.8$, growing with basis)"),
], "abstract + 4 body loci")

patch(SYN, [
 ("inflates the block-encoding $1$-norm (the ill-conditioned overlap must be\nL\\\"owdin-inverted)",
  "inflates the block-encoding $1$-norm (the non-orthogonal overlap must be\nL\\\"owdin-inverted, and $S^{-1/2}$ is dense)"),
], "synthesis")

print("\nCLASS 3 -- the s-sector value is not the s-sector limit:")
patch(PAP, [
 ("the $s$, $+p$ and $spdf$ sectors---past the $s$-sector limit $-2.873$~Ha to",
  "the $s$, $+p$ and $spdf$ sectors---past the locked $s$-sector value $-2.873$~Ha to"),
], "abstract")

print("\nCLASS 4 -- the Acknowledgments underclaim:")
patch(PAP, [
 ("The generalized-Sturmian formulation used here is that of J.~E.~Avery and\n"
  "J.~S.~Avery; this paper adds only its reading as a quantum-encoding target and\n"
  "the resource measurements.",
  "The generalized-Sturmian formulation used here is that of J.~E.~Avery and\n"
  "J.~S.~Avery, as are the Shibuya--Wulfman integrals and the split-shell account\n"
  "of why the helium ground state in particular resists the method.  To these this\n"
  "paper adds its reading as a quantum-encoding target, the resource measurements,\n"
  "and four results we claim as our own:\\ the scale-lock identity\n"
  "Eq.~\\eqref{eq:scale_lock} together with the one-body Coulomb metric\n"
  "Eq.~\\eqref{eq:W_diagonal} it rests on;\\ the interlacing bound\n"
  "Eq.~\\eqref{eq:no_selection};\\ the identification of the Shibuya--Wulfman matrix\n"
  "\\emph{as} the $V_0$-weighted overlap in Eq.~\\eqref{eq:general_v0};\\ and the\n"
  "variational bound for the fixed-scale metric-free problem."),
], "four theorems reclaimed")
print("\ndone.")
