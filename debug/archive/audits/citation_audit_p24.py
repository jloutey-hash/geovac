import re
with open('papers/group3_foundations/paper_24_bargmann_segal.tex', 'r') as f:
    text = f.read()

bib_keys = set(re.findall(r'\\bibitem\{(\w+)\}', text))
cited = set()
for m in re.finditer(r'\\cite\{([^}]+)\}', text):
    for k in m.group(1).split(','):
        cited.add(k.strip())

unused = bib_keys - cited
missing = cited - bib_keys
print('Bib keys:', sorted(bib_keys))
print('Cited keys:', sorted(cited))
print('Defined but unused:', sorted(unused))
print('Cited but missing:', sorted(missing))
