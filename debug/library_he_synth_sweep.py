"""Exact-string replacements in the group4 synthesis (He atomic-path numbers,
physical h1 re-measure 2026-09-25).  Each anchor must match exactly once."""
import sys

p = r"papers/synthesis/group4_quantum_computing_synthesis.tex"
s = open(p, encoding="utf-8").read()
REPS = [
    (r"sub-quadratic Pauli $1$-norm $\lambda = O(Q^{1.77})$ ($R^2=0.9997$;",
     r"sub-quadratic Pauli $1$-norm $\lambda = O(Q^{1.78})$ ($R^2=0.9997$;"),
    (r"equal-qubit 1-norm advantage ($3.8\times$ at $Q=10$, $7.1\times$ at",
     r"equal-qubit 1-norm advantage ($4.0\times$ at $Q=10$, $7.3\times$ at"),
    (r"exponent (which moved $3.36 \to 4.01$, above the term exponent, withdrawing",
     r"exponent (which moved $3.36 \to 4.05$, above the term exponent, withdrawing"),
    (r"the qubit number $Q$ (fit $3.78$; corrected 2026-08-29 from the retired",
     r"the qubit number $Q$ (fit $3.79$; corrected 2026-08-29 from the retired"),
    (r"algorithms---scales as $O(Q^{4.01})$, which is \emph{steeper} than the",
     r"algorithms---scales as $O(Q^{4.05})$, which is \emph{steeper} than the"),
    (r"Pauli term exponent ($3.77$) and lies inside the Gaussian",
     r"Pauli term exponent ($3.79$) and lies inside the Gaussian"),
    (r"counts and qubitization query complexity, scales as $O(Q^{1.77})$ with",
     r"counts and qubitization query complexity, scales as $O(Q^{1.78})$ with"),
    (r"within-molecule $2.5 \to 3.17$, $1$-norm $1.69 \to 1.774$, QWC",
     r"within-molecule $2.5 \to 3.17$, $1$-norm $1.69 \to \gvq{exp_lambda_4pt}{1.784}$, QWC"),
    (r"$3.36 \to \gvq{exp_qwc_3pt}{4.013}$ --- and the QWC exponent moved far enough to overturn its",
     r"$3.36 \to \gvq{exp_qwc_3pt}{4.049}$ --- and the QWC exponent moved far enough to overturn its"),
    (r"own conclusion, since $\gvq{exp_qwc_3pt}{4.013}$ exceeds the Pauli exponent and lies inside",
     r"own conclusion, since $\gvq{exp_qwc_3pt}{4.049}$ exceeds the Pauli exponent and lies inside"),
]
for a, b in REPS:
    n = s.count(a)
    if n != 1:
        sys.exit(f"ABORT {n}x: {a[:60]}")
    s = s.replace(a, b)
open(p, "w", encoding="utf-8").write(s)
print("ok", len(REPS))
