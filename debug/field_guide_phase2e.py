"""Phase 2e corrections to the field guide (stale claims + rate display)."""
import io

P = r'papers/synthesis/geovac_field_guide.tex'
s = io.open(P, encoding='utf-8').read()
fails = []


def rep(old, new, tag):
    global s
    if s.count(old) == 1:
        s = s.replace(old, new)
    else:
        fails.append((tag, s.count(old)))


rep(r'''This claim was checked, not assumed.  Working hypothesis WH1 (in the
project's internal register, CLAUDE.md~\S~1.7) was elevated to
PROVEN status in May 2026 via a five-lemma Gromov--Hausdorff
convergence theorem on truncated Camporesi--Higuchi spectral
triples~\cite{loutey_paper38} (Paper~38, released to Zenodo
2026-05-07).  The proof rests on the Latr\'emoli\`ere propinquity
framework~\cite{latremoliere2017} for quantifying when two spectral
triples are close in the Connes--Hausdorff sense, with the asymptotic
rate constant
\begin{equation}
\label{eq:propinquity_rate}
\Lambda_{n_{\max}} \le C_3 \cdot \gamma_{n_{\max}},
\qquad \gamma_{n_{\max}} \xrightarrow{n_{\max} \to \infty}
\frac{4}{\pi},
\end{equation}
universal across compact connected Lie groups~\cite{loutey_paper40}.''',
    r'''This claim was checked, not assumed --- and then re-checked.  A
five-lemma Gromov--Hausdorff convergence theorem on truncated
Camporesi--Higuchi spectral triples~\cite{loutey_paper38} (Paper~38)
gives the quantitative statement, in van~Suijlekom's state-space
Gromov--Hausdorff framework (adjacent to, but distinct from, the
Latr\'emoli\`ere propinquity~\cite{latremoliere2017}), with rate
\begin{equation}
\label{eq:propinquity_rate}
\Lambda_{n_{\max}} \le C_3 \cdot \gamma_{n_{\max}},
\qquad \gamma_{n_{\max}} \sim
\frac{4}{\pi}\,\frac{\log n_{\max}}{n_{\max}}
\xrightarrow{n_{\max} \to \infty} 0,
\end{equation}
with the constant $4/\pi$ universal across compact connected Lie
groups~\cite{loutey_paper40}.  A June 2026 adversarial verification
pass left the theorem standing \emph{conditionally}, with two named
gaps (a dual-direction reach estimate, and a kernel condition tied to
the choice of Dirac substrate;\ see Paper~38~\S~1 ``Named gaps and
history'' and the claims register,
\texttt{docs/claims\_register.md}).''', 'wh1-passage')

rep(r'''Following Paper~38, the project produced eleven math.OA standalone
papers (Papers~38--49) extending the framework along the Lorentzian
axis~\cite{loutey_paper42, loutey_paper43, loutey_paper44,
loutey_paper45, loutey_paper46, loutey_paper47, loutey_paper48,
loutey_paper49}.  Sub-results include the first published Lorentzian
propinquity convergence theorem on truncated Krein spectral triples
(Paper~45)~\cite{loutey_paper45}, an explicit strong-form extension
(Paper~46)~\cite{loutey_paper46}, and the operator-system Lorentzian
pre-length space (OSLPLS) category that closes a published-open
question of Connes--van Suijlekom on Lorentzian extensions
(Paper~49)~\cite{loutey_paper49}.  These results were checked
independently against the four published-open questions Connes--van
Suijlekom 2021~\cite{connes_vs2021} explicitly deferred to ``elsewhere.''''',
    r'''Following Paper~38, the project produced eleven math.OA standalone
papers (Papers~38--49) extending the framework along the Lorentzian
axis~\cite{loutey_paper42, loutey_paper43, loutey_paper44,
loutey_paper45, loutey_paper46, loutey_paper47, loutey_paper48,
loutey_paper49}.  The most instructive outcome of that axis is
\emph{negative}:\ the natural device for a ``Lorentzian propinquity''
--- compressing the Krein Dirac to the Krein-positive subspace ---
annihilates the spatial Dirac, so the compressed Lipschitz seminorm
vanishes identically (Paper~45's degeneracy
theorem~\cite{loutey_paper45}, with a bit-exact frozen falsifier).
An earlier convergence claim through that device, together with the
strong-form and bridge extensions built on it (Papers~46, 48, 49),
was retracted in June 2026;\ each paper carries a Status note, and
\texttt{docs/claims\_register.md} records the retraction.  What
survives:\ the operator-system substrate (Paper~44), the four-witness
modular-Hamiltonian results (Papers~42/43), the norm-resolvent
de-compactification arrow (Paper~47), and the cocycle-deficit algebra
of Paper~49.  The repaired Lorentzian target --- a Toeplitz temporal
multiplier algebra and a non-compressing Krein-compatible seminorm
--- is stated as an open problem in Paper~45.''', 'lorentzian-passage')

rep(r'''  \item \textbf{Group~1 --- Operator algebras / NCG} ($14$ papers).
    Math.OA arc.  Headline:\ WH1 PROVEN (Paper~38) with five-lemma
    Gromov--Hausdorff convergence;\ Lorentzian propinquity
    first-in-literature (Paper~45);\ Mondino--S\"amann bridge
    (Paper~48);\ OSLPLS strong-form bridge (Paper~49).  Closes
    four published-open questions of Connes--van Suijlekom
    2021~\cite{connes_vs2021}.''',
    r'''  \item \textbf{Group~1 --- Operator algebras / NCG} ($14$ papers).
    Math.OA arc.  Headline:\ conditional SU(2) state-space
    Gromov--Hausdorff convergence with explicit $4/\pi$ rate
    (Paper~38, two named gaps);\ the Lorentzian degeneracy theorem
    (Paper~45) closing the natural Lorentzian route and stating the
    repaired target;\ four-witness modular-Hamiltonian
    identification (Papers~42/43).  Papers~46--49 carry Status notes
    for descoped metric-level claims.''', 'group1-summary')

rep(r'''  \item \textbf{The mathematician (NCG, math.OA, spectral triples)}.
    Read Papers~38, 32, 45, 46, 47, 48, and 49.  Paper~38 is WH1
    PROVEN (SU(2) propinquity convergence on the Camporesi--Higuchi
    spectral triple);\ Paper~32 is the explicit spectral triple
    construction with Connes axiom audit and the \S~VIII C-arc
    closure of eight theorem-grade non-selection results;\
    Papers~45 and 46 are the first published Lorentzian propinquity
    convergence theorems (K$^+$-weak-form and strong-form,
    respectively) in the math.OA literature;\ Paper~47 closes the
    de-compactification limit at norm-resolvent level;\ Papers~48
    and 49 are the Mondino--S\"amann bridge and the OSLPLS
    strong-form bridge.''',
    r'''  \item \textbf{The mathematician (NCG, math.OA, spectral triples)}.
    Start with the three-page front-door note
    (\texttt{docs/outreach/note\_n1\_su2\_truncations.pdf}), then
    Papers~38, 45, and 32.  Paper~38 is the conditional SU(2)
    state-space Gromov--Hausdorff convergence theorem (two named
    gaps);\ Paper~45 is the Lorentzian degeneracy theorem (an
    earlier convergence claim is retracted;\ the History remark
    there explains);\ Paper~32 is the explicit spectral triple
    construction with Connes axiom audit and the \S~VIII C-arc
    closure of eight theorem-grade non-selection results.
    Papers~42/43 close the four-witness modular-Hamiltonian
    identification;\ Paper~47's norm-resolvent arrow survives the
    descope;\ Papers~46, 48, and 49 carry Status notes.''', 'readers-map')

rep(r'''\bibitem{loutey_paper45}
J.~Loutey, ``Paper 45: K$^+$-restricted weak-form Lorentzian
propinquity convergence,'' GeoVac corpus (2026).''',
    r'''\bibitem{loutey_paper45}
J.~Loutey, ``Paper 45: Lorentzian propinquity on truncated
SU(2)$\times$U(1)$_T$ Krein spectral triples:\ a degeneracy
theorem,'' GeoVac corpus (2026).''', 'bibitem45')

io.open(P, 'w', encoding='utf-8').write(s)
print('FAILS:', fails if fails else 'none')
