import re
with open('papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex', encoding='utf-8') as f:
    text = f.read()
opens = text.count('{')
closes = text.count('}')
begins = len(re.findall(r'\\begin\{', text))
ends = len(re.findall(r'\\end\{', text))
print(f'braces: {opens} open / {closes} close (diff {opens-closes})')
print(f'environments: {begins} begin / {ends} end (diff {begins-ends})')
print(f'bibitem count: {len(re.findall(r"^\\\\bibitem", text, re.MULTILINE))}')
for key in ['farsi_latremoliere2024', 'farsi_latremoliere2025', 'hekkelman_mcdonald2024b', 'toyota2023']:
    print(f'{key}: {text.count(key)}')
