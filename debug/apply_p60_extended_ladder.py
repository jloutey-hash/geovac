"""Fold the extended ladder (K -> 452), the model-selection caveat, and the
measured state-preparation overlaps into Paper 60.

Three substantive changes:
 1. Headline comparison moves from K=202 to K=452 -- better converged, and the
    ladder now runs to n_max=16.
 2. A caveat the paper did not carry: the fitted floor is WINDOW-stable (1.3% /
    0.4% drift over 21 windows, drifting UP) but MODEL-family-sensitive -- a
    rejected c + b/ln K form puts the 2^1S floor BELOW chemical accuracy at
    0.82x. So "2^1S is above chemical accuracy" is a model-selection
    conclusion, and the paper must quote the measured endpoint, not the fit.
 3. "The unmeasured cost of an interior root" is now measured.
"""
import io

PAP = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
SYN = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"
REG = "debug/qa/numeric_registry.py"
RET = "debug/qa/check_retracted_terms.py"

# ----------------------------------------------------------------- registry
reg = io.open(REG, encoding="utf-8").read()

OLD_EXC = '''    "p60_exc_gap_k202": dict(
        value=1.786, convention="constant: mHa above the exact He 2^1S energy "
                                "(-2.145974046 Ha) reached by the METRIC-FREE "
                                "isoenergetic posing, full s+p+d+f, K=202. "
                                "MEASURED ladder value, not an extrapolated floor",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_excited_ladder.py. Companion "
                   "ground-state value at the SAME K and the same ||M||_1=192.9 is "
                   "7.158 mHa. Ladder: 1.972 (K=74), 1.899 (100), 1.849 (130), "
                   "1.813 (164), 1.786 (202). A free-floor fit over that window "
                   "returns 1.65 mHa but is window-sensitive, so the measured "
                   "endpoint is what is cited in the paper.",
        aliases={1.972: "K=74", 1.813: "K=164"}),
'''
NEW_EXC = '''    "p60_exc_gap_k452": dict(
        value=1.7163, convention="constant: mHa above the exact He 2^1S energy "
                                 "(-2.145974046 Ha) reached by the METRIC-FREE "
                                 "isoenergetic posing, full s+p+d+f, at the largest "
                                 "computed basis K=452 (n_max=16). A MEASURED "
                                 "ladder endpoint, deliberately NOT an extrapolated "
                                 "floor -- see the caveat below",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_excited_ladder.py, extended to "
                   "n_max=16. Companion ground-state value at the SAME K and the "
                   "same ||M||_1 is p60_gnd_gap_k452 = 6.8196 mHa. Ladder: 1.972 "
                   "(K=74), 1.813 (164), 1.786 (202), 1.7485 (290), 1.7163 (452). "
                   "CAVEAT ON THE FLOOR: a free-floor fit is WINDOW-stable (0.4% "
                   "drift over 21 windows, drifting UP, i.e. approaching from "
                   "below) and Shanks brackets it from above at [1.647, 1.676], "
                   "but it is MODEL-family sensitive -- a two-parameter c + b/lnK "
                   "form, rejected at 340x worse RMS, puts the floor at 0.82x "
                   "chemical accuracy, i.e. BELOW it. The paper therefore cites "
                   "this measured endpoint and states that 'above chemical "
                   "accuracy' is a model-selection conclusion.",
        aliases={1.786: "K=202", 1.647: "Shanks lower bracket on the floor"}),
'''
assert OLD_EXC in reg
reg = reg.replace(OLD_EXC, NEW_EXC, 1)

OLD_GND = '''    "p60_gnd_gap_k202": dict(
        value=7.158, convention="constant: mHa above the exact He ground state "
                                "(-2.903724377 Ha) reached by the METRIC-FREE "
                                "isoenergetic posing, full s+p+d+f, K=202",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_excited_ladder.py. Paired with "
                   "p60_exc_gap_k202 at the SAME K and the same ||M||_1 = 192.9 -- "
                   "the pair is the paper's state-dependence headline, so both are "
                   "registered and both ratios are DERIVED from them.",
        aliases={8.036: "K=74", 7.289: "K=164"}),
'''
NEW_GND = '''    "p60_gnd_gap_k452": dict(
        value=6.8196, convention="constant: mHa above the exact He ground state "
                                 "(-2.903724377 Ha) reached by the METRIC-FREE "
                                 "isoenergetic posing, full s+p+d+f, at the largest "
                                 "computed basis K=452 (n_max=16)",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_excited_ladder.py. Paired with "
                   "p60_exc_gap_k452 at the SAME K and the same ||M||_1 -- the pair "
                   "is the paper's state-dependence headline, so both are "
                   "registered and both ratios are DERIVED from them rather than "
                   "typed. Shanks brackets this floor from above at [6.47, 6.62].",
        aliases={8.036: "K=74", 7.158: "K=202", 6.9143: "K=340"}),
    "p60_stateprep_overlap_exc": dict(
        value=0.798, convention="constant: L2-metric overlap between the normalized "
                                "dominant single configuration and the true 2^1S "
                                "root, K=164 -- the state-preparation cost driver "
                                "for an interior root",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_stateprep_overlap.py. The four "
                   "lowest roots give 0.992, 0.798, 0.864, 0.889 -- 2^1S is the "
                   "HARDEST of the four, not the deepest, so the driver is mixing "
                   "at the bottom of the Rydberg series and not spectral depth. "
                   "Costs 1.25x the ground state in rotations, 1.54x in "
                   "repetitions; two configurations reach 0.99.",
        aliases={0.992: "ground state", 0.6371: "overlap squared, 2^1S"}),
'''
assert OLD_GND in reg
reg = reg.replace(OLD_GND, NEW_GND, 1)
reg = reg.replace('"p60_exc_ratio_k202":   ("p60_exc_gap_k202 / chem_accuracy_mha", 0.005,\n'
                  '                             "He 2^1S error at K=202, in units of chemical accuracy"),',
                  '"p60_exc_ratio_k452":   ("p60_exc_gap_k452 / chem_accuracy_mha", 0.005,\n'
                  '                             "He 2^1S error at K=452, in units of chemical accuracy"),', 1)
reg = reg.replace('"p60_gnd_ratio_k202":   ("p60_gnd_gap_k202 / chem_accuracy_mha", 0.005,\n'
                  '                             "He ground-state error at K=202, same units, same K, same ||M||_1"),',
                  '"p60_gnd_ratio_k452":   ("p60_gnd_gap_k452 / chem_accuracy_mha", 0.005,\n'
                  '                             "He ground-state error at K=452, same units, same K, same ||M||_1"),', 1)
io.open(REG, "w", encoding="utf-8").write(reg)
print("registry: k202 -> k452 swap; state-prep entry added")

# --------------------------------------------------------------------- paper
tex = io.open(PAP, encoding="utf-8").read()

A_OLD = r"""the same to encode, and at $K=202$ the ground state sits
$\gvq{p60_gnd_ratio_k202}{4.49}\times$ above chemical accuracy while
$2\,^{1}S$ sits $\gvq{p60_exc_ratio_k202}{1.12}\times$ above it, the posing"""
A_NEW = r"""the same to encode, and at the largest computed basis $K=452$ the
ground state sits $\gvq{p60_gnd_ratio_k452}{4.28}\times$ above chemical accuracy
while $2\,^{1}S$ sits $\gvq{p60_exc_ratio_k452}{1.08}\times$ above it, the posing"""
assert A_OLD in tex, "abstract K=202 locus not found"
tex = tex.replace(A_OLD, A_NEW, 1)

S_OLD = r"""At $K=202$ ($\|M\|_1=192.9$) the ground-state error is $7.16$~mHa ($4.49\times$
chemical accuracy) against $\gvq{p60_exc_gap_k202}{1.786}$~mHa ($1.12\times$)
for $2\,^{1}S$."""
S_NEW = r"""At the largest computed basis, $K=452$, the ground-state error is
$\gvq{p60_gnd_gap_k452}{6.82}$~mHa ($4.28\times$ chemical accuracy) against
$\gvq{p60_exc_gap_k452}{1.716}$~mHa ($1.08\times$) for $2\,^{1}S$."""
assert S_OLD in tex, "Sec.4 K=202 locus not found"
tex = tex.replace(S_OLD, S_NEW, 1)

O_OLD = r"""the query count $\|M\|_1/\epsilon_p\approx2\times10^{5}$ is state-independent
through $4\,^{1}S$.  The unmeasured cost of an interior root is the overlap of a
cheap trial state, not spectral resolution."""
O_NEW = r"""the query count $\|M\|_1/\epsilon_p\approx2\times10^{5}$ is state-independent
through $4\,^{1}S$.  \textbf{[MEASURED]} The remaining cost of an interior root
is state preparation, and it is modest.  In the $L^{2}$ metric the overlap of the
dominant single configuration with the true root runs $0.992$, $
\gvq{p60_stateprep_overlap_exc}{0.798}$, $0.864$ and $0.889$ across the first
four roots, so $2\,^{1}S$ --- not the deeper ones --- is the hardest, at
$1.25\times$ the ground state in rotations and $1.54\times$ in repetitions;  two
configurations already reach $0.99$.  Interior roots become \emph{easier} again,
so the driver is mixing at the bottom of the Rydberg series rather than spectral
depth.  One practical caveat:\ the dominant configuration does not track the
spectroscopic label --- the root at He $3\,^{1}S$ is dominated by
$(l,n_a,n_b)=(0,1,4)$ --- so a preparation heuristic keyed to the physical
principal quantum number selects the wrong configuration.

\textbf{[OBSERVATION]} \emph{How well the floors themselves are known.}  They
are bracketed rather than pinned, and the dominant uncertainty is not the one we
first looked for.  Across $21$ fit windows the free-floor value drifts by only
$1.3\%$ (ground) and $0.4\%$ ($2\,^{1}S$), and drifts \emph{upward} --- the fit
approaches from below --- while a model-free Shanks extrapolation descends from
above, bracketing $[6.47,6.62]$ and $[1.647,1.676]$~mHa.  Window drift is
therefore not the exposure;\ \emph{model family} is.  A two-parameter
$c+b/\ln K$ form, rejected at $340\times$ worse RMS, would place the $2\,^{1}S$
floor at $0.82\times$ chemical accuracy --- \emph{below} it.  We therefore quote
measured ladder endpoints rather than fitted asymptotes throughout, and record
that ``$2\,^{1}S$ saturates above chemical accuracy'' is a model-selection
conclusion, not a measurement."""
assert O_OLD in tex, "Sec.4 Rydberg-gap locus not found"
tex = tex.replace(O_OLD, O_NEW, 1)

C_OLD = r"""a \emph{ground-state} pathology --- $6.4$~mHa there, but ${\sim}1.8$~mHa for
$2\,^{1}S$ at identical encoding cost, and falling by $4.3\times$,"""
C_NEW = r"""a \emph{ground-state} pathology --- $6.8$~mHa there at $K=452$, but
$1.72$~mHa for $2\,^{1}S$ at identical encoding cost, and falling by $4.3\times$,"""
assert C_OLD in tex, "conclusion locus not found"
tex = tex.replace(C_OLD, C_NEW, 1)
io.open(PAP, "w", encoding="utf-8").write(tex)
print("paper: abstract + Sec.4 (2 blocks) + conclusion updated to K=452")

# ----------------------------------------------------------------- synthesis
syn = io.open(SYN, encoding="utf-8").read()
Y_OLD = r"""encode --- and at $K=202$ the ground state sits $4.49\times$ above chemical
accuracy while $2\,^{1}S$ sits $1.12\times$, the posing cost falling by"""
Y_NEW = r"""encode --- and at $K=452$ the ground state sits $4.28\times$ above chemical
accuracy while $2\,^{1}S$ sits $1.08\times$, the posing cost falling by"""
assert Y_OLD in syn, "synthesis K locus not found"
syn = syn.replace(Y_OLD, Y_NEW, 1)
io.open(SYN, "w", encoding="utf-8").write(syn)
print("synthesis: updated to K=452")

# --------------------------------------------- C16 note: better l-tail number
ret = io.open(RET, encoding="utf-8").read()
OLD_T = '"l>=4 partial-wave tail is 0.37-0.53 mHa, an order below the "\n                "6.4 mHa floor, so angular truncation cannot be the mechanism.  "'
NEW_T = '"l>=4 partial-wave tail is 0.187 mHa (ground) and 0.008 mHa "\n                "(2^1S) -- 3% and 0.5% of the respective floors -- measured from "\n                "the l-increments themselves, so angular truncation cannot be "\n                "the mechanism.  "'
if OLD_T in ret:
    ret = ret.replace(OLD_T, NEW_T, 1)
    io.open(RET, "w", encoding="utf-8").write(ret)
    print("C16: l-tail estimate replaced by the measured increments")
else:
    print("C16: l-tail string not matched -- left as is (check manually)")
