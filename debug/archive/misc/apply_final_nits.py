"""Remaining NITs from the DELTA code review, plus the inline test citation.

- geovac/prolate_general_m.py: a Ytab dict built and immediately discarded, a
  dead `legendre_deriv_poly_safe` returning a constant, and a
  `P.polynegative if False else ...` dead conditional.
- Paper 12 cites only tests/test_neumann_vee.py inline; the new section's
  backing test is uncited (C13 reports matrix coverage, but a reader of the
  paper cannot find it).

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import re
import sys

MOD = "geovac/prolate_general_m.py"
P12 = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"

with io.open(MOD, encoding="utf-8") as fh:
    m = fh.read()

changed = []

# 1. the discarded Ytab block
dead_block = """    # eta moments: Y[(m,s,Q)]
    Ytab = {}
    for m in m_set:
        pm = legendre_deriv_poly_safe(m)
        for s in s_set:
            yp = poly_1meta2(s)
            for Qq in range(q_max + 1):
                Ytab[(m, s, Qq)] = None  # filled per l below
    # eta needs l too (P_l^m), so tabulate on (l,m,s,Q)
    Ytab = {}
"""
if dead_block in m:
    m = m.replace(dead_block,
                  "    # eta moments, tabulated on (l, m, s, Q): the kernel's\n"
                  "    # P_l^m depends on l, so the table cannot be keyed on m alone.\n"
                  "    Ytab = {}\n", 1)
    changed.append("removed the discarded Ytab block")

# 2. dead helper
dead_fn = '''def legendre_deriv_poly_safe(m):
    return np.array([1.0])


'''
if dead_fn in m:
    m = m.replace(dead_fn, "", 1)
    changed.append("removed dead legendre_deriv_poly_safe")

# 3. dead conditional
if "t = P.polynegative if False else shift(" in m:
    m = m.replace("t = P.polynegative if False else shift(",
                  "t = shift(", 1)
    changed.append("removed the dead `if False` conditional")

with io.open(MOD, "w", encoding="utf-8") as fh:
    fh.write(m)

for c in changed:
    print("  +", c)

# ---- inline test citation in Paper 12
with io.open(P12, encoding="utf-8") as fh:
    t = fh.read()

OLD = r"""\textbf{[SCOPE]} Two limits are worth stating plainly."""
NEW = (r"""Backing:\ \texttt{tests/test\_paper12\_azimuthal\_channels.py},"""
       "\n"
       r"""whose guards are fire-tested against the specific wrong answers they"""
       "\n"
       r"""exclude (a killed azimuthal coupling, the superseded kernel prefactor,"""
       "\n"
       r"""a disabled canonical orthogonalization, and a truncation cap that would"""
       "\n"
       r"""hide the selection rule)."""
       "\n\n"
       r"""\textbf{[SCOPE]} Two limits are worth stating plainly.""")

if OLD in t:
    t = t.replace(OLD, NEW, 1)
    with io.open(P12, "w", encoding="utf-8") as fh:
        fh.write(t)
    print("  + Paper 12: inline backing-test citation added")
else:
    print("  ! Paper 12 scope anchor not found")
    sys.exit(1)
