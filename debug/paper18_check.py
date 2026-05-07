import re
with open('papers/core/paper_18_exchange_constants.tex', encoding='utf-8') as f:
    text = f.read()
opens = text.count('{')
closes = text.count('}')
begins = len(re.findall(r'\\begin\{', text))
ends = len(re.findall(r'\\end\{', text))
print(f'braces: {opens} open / {closes} close (diff {opens-closes})')
print(f'environments: {begins} begin / {ends} end (diff {begins-ends})')
print(f'eta_trivialization: {text.count("eta_trivialization")}')
print(f'ac_factorization: {text.count("ac_factorization")}')
print(f'krajewski1998: {text.count("krajewski1998")}')
print(f'paschke_sitarz2000: {text.count("paschke_sitarz2000")}')
print(f'loutey_paper31: {text.count("loutey_paper31")}')
