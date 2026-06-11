"""Scan papers for unescaped underscores inside \texttt{...} commands.

Used as a one-off diagnostic in cleanup track A.
"""
import re
import sys

PATHS = [
    'papers/group1_operator_algebras/paper_32_spectral_triple.tex',
    'papers/group6_precision_observations/paper_34_projection_taxonomy.tex',
    'papers/group2_quantum_chemistry/paper_17_composed_geometries.tex',
    'papers/group3_foundations/paper_18_exchange_constants.tex',
    'papers/group2_quantum_chemistry/paper_19_coupled_composition.tex',
]

BACKSLASH = chr(92)

def scan(path):
    hits = 0
    with open(path, 'r', encoding='utf-8') as f:
        for lineno, line in enumerate(f, 1):
            for m in re.finditer(r'\\texttt\{([^{}]*)\}', line):
                arg = m.group(1)
                for j, ch in enumerate(arg):
                    if ch == '_' and (j == 0 or arg[j-1] != BACKSLASH):
                        print(f'{path}:{lineno}: unescaped _ in texttt: {arg[:80]}')
                        hits += 1
                        break
    return hits

total = sum(scan(p) for p in PATHS)
print(f'Total hits: {total}')
sys.exit(0 if total == 0 else 1)
