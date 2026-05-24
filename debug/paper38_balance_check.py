import re
text = open('papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex', encoding='utf-8').read()
begins = re.findall(r'\\begin\{([\w*]+)\}', text)
ends = re.findall(r'\\end\{([\w*]+)\}', text)
from collections import Counter
b = Counter(begins); e = Counter(ends)
mismatched = {k: (b[k], e[k]) for k in set(begins)|set(ends) if b[k] != e[k]}
print('mismatched:', mismatched)
print('total begins/ends:', sum(b.values()), '/', sum(e.values()))
print('lines:', len(text.splitlines()))
# Label collision check
labels = re.findall(r'\\label\{([^}]+)\}', text)
dupes = [l for l, c in Counter(labels).items() if c > 1]
print('duplicate labels:', dupes)
# refs to nonexistent labels
refs = set(re.findall(r'\\(?:eq)?ref\{([^}]+)\}', text))
defined = set(labels)
missing = refs - defined
print('missing refs:', sorted(missing))
