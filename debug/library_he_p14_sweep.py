"""Sweep the He exponent literals in Paper 14 lines >= 2700 (He atomic-path
prose only; the composed-LiH '3.77' slopes live below line 1200)."""
p = r"papers/group4_quantum_computing/paper_14_qubit_encoding.tex"
lines = open(p, encoding="utf-8").read().split("\n")
REPS = [
    (r"Q^{1.77}", r"Q^{1.78}"),
    (r"O(Q^{4.01})", r"O(Q^{4.05})"),
    (r"its exponent ($4.01$) exceeding the term exponent ($3.77$)",
     r"its exponent ($4.05$) exceeding the term exponent ($3.79$)"),
    (r"$\alpha_{\mathrm{terms}} = 3.77$", r"$\alpha_{\mathrm{terms}} = 3.79$"),
    (r"$\alpha_{\lambda} = 1.77$", r"$\alpha_{\lambda} = 1.78$"),
    (r"scaling exponent $\alpha = 1.47$ versus~$1.77$",
     r"scaling exponent $\alpha = 1.47$ versus~$1.78$"),
]
count = 0
for i in range(2700, len(lines)):
    for a, b in REPS:
        if a in lines[i]:
            lines[i] = lines[i].replace(a, b)
            count += 1
            print(i + 1, b)
open(p, "w", encoding="utf-8").write("\n".join(lines))
print("replacements", count)
