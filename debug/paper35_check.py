import re
with open('papers/group6_precision_observations/paper_35_time_as_projection.tex', encoding='utf-8') as f:
    text = f.read()
opens = text.count('{')
closes = text.count('}')
begins = len(re.findall(r'\\begin\{', text))
ends = len(re.findall(r'\\end\{', text))
print(f'braces: {opens} open / {closes} close (diff {opens-closes})')
print(f'environments: {begins} begin / {ends} end (diff {begins-ends})')
print(f'Refined Prediction 1: {text.count("Refined Prediction")}')
print(f'LS-8a mentions: {text.count("LS-8a")}')
