# Paper 60 — two named attributions with no bibitem: source verification

**Date:** 2026-09-12
**Scope:** read-only literature verification. Nothing under `papers/`, `geovac/`, `tests/`,
`CLAUDE.md`, `CHANGELOG.md` was modified.
**Trigger:** C20 defects from the citation dimension of a `/qa` run — two works named in
Paper 60 prose with no resolvable bibliography entry.

| Target | Verdict |
|:---|:---|
| 1. "Bernstein $\Theta(\kappa)$ floor" | **VERIFIED** — with a correction: the published statement of the floor is **already in Paper 60's bibliography** (`gslw2019`, Theorem 73); the word "Bernstein" names the corpus's *own* derivation tool and needs its own classical citation. `cks2017` is **not** a source for it. |
| 2. Goscinski | **VERIFIED** — with a page-range correction (51–**56**, not 51–85) and a flagged discrepancy about which volume carries the 1968 report. |

---

## Target 1 — the "Bernstein $\Theta(\kappa)$ floor"

### 1.1 What Paper 60 actually claims

Two loci (current line numbers in
`papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex`):

- **L1221–1225** (after the band-Toeplitz preconditioner result): *"This escapes the Bernstein
  $\Theta(\kappa)$ floor quoted above rather than contradicting it: that floor constrains
  polynomial approximation of $x^{-1/2}$ on $[\kappa^{-1},1]$, and preconditioning changes the
  operator rather than the polynomial…"*
- **L1426–1429** (spectrum-aware inversion probe): *"…the Bernstein $\Theta(\kappa)$ floor for
  **both** problems bounds it there (an argument, not a measurement: the floors hold for the
  interpolation and interval families alike…)"*

"quoted above" points at **L1125**, the resource-model sentence
`d_inv ~ \kappa\ln(\kappa/\epsilon)`~\cite{gslw2019,cks2017} — i.e. the paper currently cites
GSLW+CKS for the *upper* bound and then refers to a *lower* bound under the name "Bernstein"
with no citation at all. That is the C20 defect.

### 1.2 Provenance inside the repo (found first, and it settles what "Bernstein" means)

`debug/probeA_findings.md` §"Table 5 — rigorous Bernstein lower bounds" states the argument
verbatim:

> For `|p| <= 1` on `[-1,1]`, Bernstein gives `|p'(x)| <= d/sqrt(1-x^2)`, hence between any two
> accuracy points `|p(x_b) - p(x_a)| <= d*|arcsin(x_b) - arcsin(x_a)|`. … `aware/p60 floor ->
> h*kappa/6`, `generic/p60 floor -> h*kappa/2`: **both `Theta(kappa)`**.

So **"Bernstein" is the classical polynomial-derivative inequality**, used as the corpus's own
lower-bound tool — not a cited quantum result. That is why no bibitem exists, and it is why the
right repair is *two* citations, not one (see §1.7).

### 1.3 The published statement of the same floor — GSLW, already cited as `gslw2019`

Gilyén–Su–Low–Wiebe (arXiv:1806.01838; STOC 2019) contains **both** halves, and the second
half is exactly Paper 60's floor. Verified by extracting the arXiv PDF text locally (67 pp.):

**Corollary 67 (Polynomial approximations of negative power functions)** — verbatim:

> Let $\delta,\varepsilon \in (0,\tfrac12]$, $c>0$ and let $f(x) := \frac{\delta^c}{2}x^{-c}$,
> then there exist even/odd polynomials $P,P' \in \mathbb{R}[x]$, such that
> $\|P-f\|_{[\delta,1]}\le\varepsilon$, $\|P\|_{[-1,1]}\le 1$ and similarly
> $\|P'-f\|_{[\delta,1]}\le\varepsilon$, $\|P'\|_{[-1,1]}\le 1$, moreover the degree of the
> polynomials are $O\!\left(\frac{\max[1,c]}{\delta}\log\!\left(\frac1\varepsilon\right)\right)$.

With $c=\tfrac12$ and $\delta=\kappa^{-1}$ this **is** Paper 60's $x^{-1/2}$-on-$[\kappa^{-1},1]$
problem, and the degree is $O(\kappa\log(1/\varepsilon))$.

The matching lower bound, in the same proof paragraph (p. 55 of the arXiv v-text), verbatim:

> Since the derivative of the function $\frac{\delta^c}{2}x^{-c}$ at $x=\delta$ is
> $-\frac{c}{2\delta}$, we get by **Theorem 73** that the $\delta$ and $c$ dependence of the
> complexity of this procedure is **optimal**.

**Theorem 73 (Lower bound for eigenvalue transformation)** — verbatim:

> Let $I\subseteq[-1,1]$, $a\ge1$ and suppose $U$ is a $(1,a,0)$-block-encoding of an unknown
> Hermitian matrix $H$ with the only promise that the spectrum of $H$ lies in $I$. Let
> $f:I\to\mathbb{R}$, and suppose that we have a quantum circuit $V$ that implements a
> $(1,b,\varepsilon)$-block-encoding of $f(H)$ using $T$ applications of $U$, for all $U$
> fulfilling the promise. Then for all $x\ne y\in I\cap[-\tfrac12,\tfrac12]$ we have that
> $T=\Omega\!\left(\frac{|f(x)-f(y)|-2\varepsilon}{|x-y|}\right)$.

and the motivating sentence immediately before it:

> Intuitively speaking if a function has derivative $d$ on the domain of interest then we need to
> use the block-encoding $\Omega(d)$-times in order to implement the eigenvalue transformation
> corresponding to $f$.

$\delta=\kappa^{-1}\Rightarrow |f'(\delta)| = c\kappa/2 \Rightarrow T=\Omega(\kappa)$. **That is
the $\Theta(\kappa)$ floor, stated for $x^{-c}$ — hence for $x^{-1/2}$ — by the paper Paper 60
already cites.** Note the precise-form denominator of GSLW Eq. (68),
$\sqrt{2}\sqrt{1-xy-\sqrt{(1-x^2)(1-y^2)}}$, is the same arcsine-type metric the corpus's own
Bernstein integration produces; the two arguments are the same geometry in different models
(GSLW bounds *queries* $T$; Bernstein bounds *degree* $d$; in QSVT they coincide).

GSLW itself cites the classical degree-lower-bound literature as
`[SV14] Sushant Sachdeva and Nisheeth K. Vishnoi. Faster algorithms via approximation theory.
Found. Trends Theor. Comput. Sci., 9(2):125–210, 2014.` (verified verbatim in GSLW's
bibliography; Crossref: doi:10.1561/0400000065).

### 1.4 `cks2017` is NOT a source for the floor — do not cite it there

Childs–Kothari–Somma (arXiv:1511.02306) was extracted and searched: it contains **no**
polynomial-degree lower bound. Its only optimality remark is about *query* complexity and is
attributed elsewhere, verbatim:

> Ambainis later improved the $\kappa$-dependence of the HHL algorithm from quadratic to nearly
> linear [Amb12], which is essentially optimal since the dependence on $\kappa$ cannot be made
> sublinear [HHL09].

That is the Harrow–Hassidim–Lloyd query lower bound, a different statement from a polynomial
approximation floor.

### 1.5 The sharp (constant-resolved) result exists — but only for $x^{-1}$

For the **two-sided** domain $S(a)=[-1,-a]\cup[a,1]$ with $a=1/\kappa$, the minimax error of the
best odd polynomial approximation to $1/x$ is known **in closed form**:

- **I. A. Privalov, "Approximation of $1/x$ by polynomials on $[-1,-a]\cup[a,1]$", Mat. Zametki
  81(3), 472–473 (2007); Math. Notes 81(3), 415–416 (2007).**
  Verified at Crossref (doi:10.1134/S0001434607030157; Pleiades; vol. 81, issue 3–4, pp. 415–416)
  and at mathnet.ru (paper `mzm3688`, doi:10.4213/mzm3688; Russian pp. 472–473, Brief
  Communications; keywords: *polynomial approximation, extremal polynomial, polynomial of least
  deviation, uniform metric*).
- Restated in QSVT language, with the constant measured at **1.00**, by
  **Sünderhauf, Németh, Walayat, Patterson & Berntson, "Matrix inversion polynomials for the
  quantum singular value transformation", arXiv:2507.15537**, whose Theorem 1 gives
  $\varepsilon_{2n-1}(a) = (1-a)^n / (a(1+a)^{n-1})$, hence
  $d = 2n-1 \approx 1.00\,\kappa\log(\kappa/\varepsilon)+1$ — an **equality**, so both ceiling and
  floor. Their Theorem 1 is explicitly "based on the theorem stated in [10]" = Privalov.

**This is for $x^{-1}$, not $x^{-1/2}$.** No sharp-constant $x^{-1/2}$ analogue was found. The
answer to the task's explicit question:

| statement | function | domain | status |
|:---|:---|:---|:---|
| GSLW Cor. 67 (upper) + Thm 73 (lower) | $x^{-c}$, any $c>0$ — **includes $c=1/2$** | $[\delta,1]$, bounded on $[-1,1]$, definite parity | published, order-only ($\Theta(\kappa)$) |
| Privalov 2007 / Sünderhauf et al. 2025 | $x^{-1}$ **only** | $[-1,-a]\cup[a,1]$, odd | published, **exact constant** |
| probeA Table 5 (corpus's own) | $x^{-1/2}$, and the spectrum-aware *interpolation* family | $[\kappa^{-1},1]$, $|p|\le1$ on $[-1,1]$ | internal; Bernstein inequality; the L1426 locus |

### 1.6 Caveat worth keeping in the paper (it is already scoped correctly, don't lose it)

The floor is $\Theta(\kappa)$ **because of the QSVT normalisation constraints**, not because of
$x^{-1/2}$ itself. Unconstrained one-sided approximation of $x^{-1/2}$ on $[\kappa^{-1},1]$ is a
$\sqrt{\kappa}$ problem (Bernstein-ellipse parameter $\rho\approx1+2/\sqrt\kappa$); it is
$\|P\|_{[-1,1]}\le1$ **plus definite parity** — equivalently, the two-sided domain — that removes
the square root. Paper 60's own text already registers the other side of this (the PSD-shift
variant at $\approx2.8\sqrt\kappa$, amplification-gated), and `debug/probeA_findings.md` measures
both arms ($\kappa^{0.996}$ vs $\kappa^{0.503}$). Any citation added should not be allowed to
read as "approximating $x^{-1/2}$ needs degree $\kappa$" without the convention.

### 1.7 Recommended remediation

1. **L1125** (upper bound, already cited): keep `\cite{gslw2019,cks2017}`, and add the pointer
   that makes the floor citable: *"…and the $\kappa$-dependence is optimal
   (\cite[Theorem~73]{gslw2019})."*
2. **L1221** ("the Bernstein $\Theta(\kappa)$ floor quoted above"): point at GSLW Thm 73, and
   carry `\cite{bernstein1912}` for the inequality the corpus's own version of the argument uses.
3. **L1426** (the spectrum-aware *interpolation* family): this one is genuinely the corpus's own
   extension — GSLW Thm 73 does not cover the interpolation family — so cite
   `\cite{bernstein1912}` and keep the paper's existing honest hedge ("an argument, not a
   measurement").
4. Do **not** attach the floor to `cks2017` (§1.4).
5. Optional, if the paper wants the sharp constant on record: add `privalov2007` +
   `sunderhauf2025` with the explicit note that they are the $x^{-1}$ problem.

### 1.8 Ready-to-paste bibitems (Paper 60 style: `\textit{Journal}\ \textbf{vol}, page (year).`)

```latex
\bibitem{bernstein1912}
S.~N.~Bernstein, ``Sur l'ordre de la meilleure approximation des fonctions
continues par des polyn\^omes de degr\'e donn\'e,''
\textit{M\'em.\ Cl.\ Sci.\ Acad.\ Roy.\ Belg.}\ (Collection in-4$^{\circ}$)
\textbf{4}, 1--104 (1912).  The form used here --- $|p'(y)|\le
n\|p\|_{[-1,1]}/\sqrt{1-y^{2}}$ for $p$ of degree $n$ and $y\in(-1,1)$ --- is
stated as the Bernstein Inequality in P.~Borwein and T.~Erd\'elyi,
\textit{Polynomials and Polynomial Inequalities}, Graduate Texts in
Mathematics Vol.~161 (Springer, New York, 1995).

\bibitem{privalov2007}
I.~A.~Privalov, ``Approximation of $1/x$ by polynomials on
$[-1,-a]\cup[a,1]$,'' \textit{Math.\ Notes}\ \textbf{81}, 415 (2007)
[\textit{Mat.\ Zametki}\ \textbf{81}, 472 (2007)];\
doi:10.1134/S0001434607030157.

\bibitem{sunderhauf2025}
C.~S\"underhauf, Z.~N\'emeth, A.~Walayat, A.~Patterson, and B.~K.~Berntson,
``Matrix inversion polynomials for the quantum singular value
transformation,'' arXiv:2507.15537 (2025).
```

`gslw2019` needs no new entry — only the theorem pointer.

### 1.9 Verification log (Target 1)

| Claim | How verified |
|:---|:---|
| GSLW Cor. 67 text, Thm 73 text, "optimal" remark | arXiv PDF `1806.01838` downloaded, text extracted with `pypdf`, quoted verbatim from the extraction |
| GSLW has no "Bernstein" and no degree-lower-bound section | `grep -i` over the same extraction: 0 hits for "bernstein"; all 15 "lower bound" hits inspected |
| CKS2017 has no degree lower bound | arXiv PDF `1511.02306v2` extracted; only optimality remark quoted in §1.4 |
| Privalov record | Crossref work record (DOI, title, author, container, vol. 81, iss. 3–4, pp. 415–416, 2007, Pleiades) **and** mathnet.ru paper page `mzm3688` (Russian pp. 472–473, doi:10.4213/mzm3688) |
| Privalov is the source of the exact minimax | Sünderhauf et al. reference [10] + "our proof is based on the theorem stated in [10]", read in the extracted PDF |
| Sünderhauf et al. record + Theorem 1 + $d\approx1.00\kappa\log(\kappa/\varepsilon)+1$ | arXiv PDF `2507.15537` extracted and quoted |
| Bernstein 1912 memoir | Persée record `marb_0365-0952_1912_num_4_1_3751`, citation line verbatim: *"Bernstein Serge. Sur l'ordre de la meilleure approximation des fonctions continues par des polynômes de degré donné. In: Mémoires de la Classe des sciences. Académie royale de Belgique. Collection in-4°. Tome 4, 1912. 1912-1922. pp. 1-104."* |
| Bernstein inequality statement + attribution | T. Erdélyi, "Markov- and Bernstein-type inequalities for polynomials with restricted coefficients" (author PDF), which states it verbatim and cites Bernstein [1], DeVore–Lorentz, Borwein–Erdélyi |
| Borwein–Erdélyi book | Crossref book record doi:10.1007/978-1-4612-0793-1 (Springer New York, 1995, GTM, ISBN 9780387945095); GTM number 161 from the trade listings carrying the same ISBN |
| Sachdeva–Vishnoi (GSLW's [SV14]) | GSLW bibliography verbatim + Crossref doi:10.1561/0400000065 (Found. Trends Theor. Comput. Sci. **9**(2), 125–210, 2014) |

**Routes that failed:** Springer article/book pages (303 redirect to an auth endpoint);
zbMATH (403); ScienceDirect (403); ADS abstract page (405 to WebFetch); `ftp.math.utah.edu`
(connection refused); Google Books API (429, repeatedly).

---

## Target 2 — Goscinski

### 2.1 What Paper 60 attributes

L147: *"The generalized-Sturmian method of Avery and
Avery~\cite{avery1989,avery2006,averyphd,averymsc}, **built on Goscinski's weighted wave
equation** and the momentum-space Fock projection…"* — and "Goscinskian configurations" then
carries §§ around L300, L375, L415, L538, L607, L1497, L1508, L1603. No bibitem exists.

### 2.2 Verdict: **VERIFIED**, with one correction and one flagged discrepancy

**Record A — the one to cite (fully verified):**

> O. Goscinski, "Conjugate eigenvalue problems and generalized Sturmians,"
> *Advances in Quantum Chemistry* **41**, 51–56 (2002); doi:10.1016/S0065-3276(02)41046-5
> (Elsevier; ISSN 0065-3276).

- **Correction to the flagged record:** the pages are **51–56**, *not* 51–85. Crossref and
  Semantic Scholar both return `51-56`; ADS bibcode `2002AdQC...41...51G` gives first page 51.
  Volume **41** and year **2002** are confirmed.
- The publisher-supplied abstract (read through a Crossref abstract mirror) states, in substance:
  the work presents the unpublished **Research Report No. 217** of the **Quantum Chemistry Group,
  Uppsala University**, titled **"Conjugate Eigenvalue Problems and the Theory of Upper and Lower
  Bounds"**; it notes the attention the theory has received in generalized and molecular
  Sturmians, references contemporary work by John Avery and colleagues, and says the report
  appears as an appendix, in a volume honouring Per-Olov Löwdin. **The abstract does not carry a
  year for the report.**

**Record B — the 1968 report itself (the thing being reprinted):**

> O. Goscinski, *Conjugate Eigenvalue Problems and the Theory of Upper and Lower Bounds*,
> Preliminary Research Report No. 217, Quantum Chemistry Group, Uppsala University (1968),
> unpublished.

- The **report number (217), the institution, and the title** are confirmed by Record A's own
  abstract.
- The **year 1968** is confirmed from published secondary literature, not from the report (which
  is unpublished and unreachable). Verbatim, reference [5] of **M. J. Ambrosio, J. A. Del Punta,
  K. V. Rodriguez, G. Gasaneo and L. U. Ancarani, "Mathematical properties of generalized
  Sturmian functions", J. Phys. A: Math. Theor. 45, 015201 (2012)**:
  > `[5] Goscinski O 1968 Preliminary Research Report no 217 Quantum Chemistry Group, Uppsala
  > University (originally unpublished but included as an appendix in 2003 Adv. Quantum Chem. 43 207)`
  and its body text:
  > *"In 1968, Goscinsky [5] presented, in an Uppsala University internal report, a rigorous
  > mathematical generalization of this basis set. This report was unknown to the atomic physics
  > community until recently when Goscinsky and Avery presented it as an appendix in [6].
  > Originally, Goscinsky regarded Sturmian functions as solutions to the Schrödinger equation
  > with a constant and externally defined energy, and considered the magnitude of the potential
  > (the 'charges') as the eigenvalues."*

  That last sentence is, to the letter, the "weighted wave equation" Paper 60 attributes to him.

**Flagged discrepancy (report it, do not silently resolve).** Two published sources disagree on
*which* Adv. Quantum Chem. item carries the 1968 report as an appendix:

- Record A's own abstract (vol. **41**, 51–56, 2002) says the report appears as an appendix;
- Ambrosio et al. (2012) ref. [5] says it was included as an appendix in **Adv. Quantum Chem.
  43, 207 (2003)** — i.e. in **J. Avery, J. Avery and O. Goscinski, "Natural Orbitals from
  Generalized Sturmian Calculations", Adv. Quantum Chem. 43, 207–216 (2003),
  doi:10.1016/S0065-3276(03)43006-2** (Crossref-verified).

Neither article's full text was reachable (ScienceDirect 403), so this is not resolved here. It
does **not** affect Record A's own bibliographic data, which is verified independently three ways.

**Third-party citation form, for calibration.** Granados-Castro & Ancarani, Adv. Quantum Chem.
73, 3–57 (2016), reference list, verbatim:
> `Goscinski, O., Conjugate Eigenvalue Problems and Generalized Sturmians (2002) Adv. Quantum
> Chem., 41, p. 51. , Preliminary research unpublished. Included as an appendix`

So the community form is "Adv. Quantum Chem. 41, 51 (2002)" with a note that it is the
previously-unpublished preliminary research.

### 2.3 Recommended bibitem

Cite Record A (resolvable DOI, title matching the attributed concept) and name the 1968 report
inside the same entry — the form Avery's literature uses. If the paper also wants the explicit
weighted-potential/isoenergetic formulation on record, the 2003 Avery–Avery–Goscinski entry below
is the tightest match to Paper 60's own construction; its abstract reads: *"The configurations in
the basis set are solutions to an approximate Schrödinger equation with a weighted potential
$\beta_\nu V_0(x)$, the weighting factors $\beta_\nu$ being chosen in such a way as to make the
set of solutions isoenergetic."*

```latex
\bibitem{goscinski2002}
O.~Goscinski, ``Conjugate eigenvalue problems and generalized Sturmians,''
\textit{Adv.\ Quantum Chem.}\ \textbf{41}, 51 (2002);\
doi:10.1016/S0065-3276(02)41046-5.  A presentation of the previously
unpublished \textit{Conjugate Eigenvalue Problems and the Theory of Upper and
Lower Bounds}, Preliminary Research Report No.~217, Quantum Chemistry Group,
Uppsala University (1968).

\bibitem{avery2003natorb}
J.~Avery, J.~Avery, and O.~Goscinski, ``Natural orbitals from generalized
Sturmian calculations,'' \textit{Adv.\ Quantum Chem.}\ \textbf{43}, 207
(2003);\ doi:10.1016/S0065-3276(03)43006-2.
```

(The second is optional — add it only if a locus actually rests on the weighted-potential /
isoenergetic statement rather than on the 1968 priority.)

### 2.4 Verification log (Target 2)

| Claim | How verified |
|:---|:---|
| DOI, title, author, container, pages 51–56, 2002, Elsevier, ISSN 0065-3276 | Crossref work record `10.1016/s0065-3276(02)41046-5` (queried directly and via bibliographic search) |
| Volume = **41**, pages 51–56 | Semantic Scholar graph API record for the same DOI (journal: name/volume/pages) |
| Independent third confirmation of vol. 41 / p. 51 | ADS bibcode `2002AdQC...41...51G`; Granados-Castro & Ancarani (Adv. Quantum Chem. 73) reference list |
| Report No. 217, Uppsala Quantum Chemistry Group, report title, "appears as an appendix" | Publisher abstract for the DOI, read via a Crossref abstract mirror (colab.ws) |
| Report year 1968 + what the report introduced | Ambrosio et al., J. Phys. A 45, 015201 (2012): PDF downloaded from the CONICET repository, text extracted with `pypdf`, ref. [5] and body sentence quoted verbatim |
| Avery–Avery–Goscinski 2003 record + abstract | Crossref work record `10.1016/s0065-3276(03)43006-2`; abstract via the same Crossref mirror |
| "51–85" is wrong | Two aggregators return 51–56; no source anywhere returned 85 |

**Routes that failed:** ScienceDirect article pages (403); Springer (auth redirect); vdoc.pub
copy of the Adv. Quantum Chem. tribute volume (403); ADS abstract page direct fetch (405);
Google Books API (429). The *full text* of neither Adv. Quantum Chem. article was reached, which
is why the "which volume carries the appendix" question is flagged rather than answered.

*Aside, corrected while checking:* Adv. Quantum Chem. **47** (2004) — not 41 — is the volume
titled *A Tribute Volume in Honour of Professor Osvaldo Goscinski*. Volume 41 is a different
2002 volume (the abstract associates it with a Per-Olov Löwdin honour volume).

---

## What this memo does not claim

- It does not assert that any locus of Paper 60 is *wrong*. Both attributions are substantively
  correct; the defects are missing/underspecified citations and one page range.
- It does not resolve which Adv. Quantum Chem. volume physically reprints the 1968 report.
- The $x^{-1/2}$ Bernstein floor as Paper 60 states it (the interpolation family at L1426) has no
  published owner and is the corpus's own argument; the interval family at L1221 does
  (GSLW Thm 73). Those two should not be given the same citation.
