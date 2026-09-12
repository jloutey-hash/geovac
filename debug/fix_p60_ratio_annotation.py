"""C21 caught a mis-keyed annotation: the abstract cited a RATIO (1.12x chemical
accuracy) against a registry entry holding an ENERGY (1.786 mHa).  Register the
comparison properly instead of loosening the annotation:

  chem_accuracy_mha   - the conversion constant, registered once
  p60_gnd_gap_k202    - ground-state error at K=202 (the other half of the pair)
  p60_exc_ratio_k202  - DERIVED, so the ratio and the energy cannot drift apart
  p60_gnd_ratio_k202  - DERIVED, likewise

Both ratios become recomputed quantities, so C21 will re-derive them rather than
trust a typed number.
"""
import io

REG = "debug/qa/numeric_registry.py"
PAP = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"

src = io.open(REG, encoding="utf-8").read()

NEW = '''    "chem_accuracy_mha": dict(
        value=1.5936014616, convention="constant: chemical accuracy = 1 kcal/mol "
                                       "expressed in mHa (4.184 kJ/mol divided by "
                                       "2625.4996 kJ/mol per Hartree)",
        q=None,
        provenance="DEFINITION, CODATA-consistent unit conversion. Registered "
                   "2026-09-08 because Paper 60 states three accuracies as "
                   "MULTIPLES of it, and those multiples must be derived from the "
                   "energies rather than typed independently.",
        aliases={1.594: "3 dp, as printed in Sec.4"}),
    "p60_gnd_gap_k202": dict(
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
anchor = '    "p60_cond_S_converged": dict('
assert anchor in src and "chem_accuracy_mha" not in src
src = src.replace(anchor, NEW + anchor, 1)

DERIV_ANCHOR = 'DERIVED = {\n'
assert DERIV_ANCHOR in src
src = src.replace(DERIV_ANCHOR, DERIV_ANCHOR +
    '    "p60_exc_ratio_k202":   ("p60_exc_gap_k202 / chem_accuracy_mha", 0.005,\n'
    '                             "He 2^1S error at K=202, in units of chemical accuracy"),\n'
    '    "p60_gnd_ratio_k202":   ("p60_gnd_gap_k202 / chem_accuracy_mha", 0.005,\n'
    '                             "He ground-state error at K=202, same units, same K, same ||M||_1"),\n\n', 1)
io.open(REG, "w", encoding="utf-8").write(src)
print("registry: chem_accuracy_mha + p60_gnd_gap_k202 added; 2 derived ratios added")

# re-key the abstract to the ratio entries
tex = io.open(PAP, encoding="utf-8").read()
OLD = (r"and at $K=202$ the ground state sits $4.49\times$ above" "\n"
       r"chemical accuracy while $2\,^{1}S$ sits $\gvq{p60_exc_gap_k202}{1.12}\times$" "\n"
       r"above it, the posing cost falling threefold to fourfold per rung up the $^{1}S$")
NEW_T = (r"and at $K=202$ the ground state sits" "\n"
         r"$\gvq{p60_gnd_ratio_k202}{4.49}\times$ above chemical accuracy while"  "\n"
         r"$2\,^{1}S$ sits $\gvq{p60_exc_ratio_k202}{1.12}\times$ above it, the posing"  "\n"
         r"cost falling threefold to fourfold per rung up the $^{1}S$")
assert OLD in tex, "abstract ratio locus not found"
tex = tex.replace(OLD, NEW_T, 1)
io.open(PAP, "w", encoding="utf-8").write(tex)
print("paper: abstract ratios re-keyed to the derived entries")
