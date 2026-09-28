"""Exponent fits + Trotter columns for the physical-h1 He series.
Reads debug/data/library_he_fix.jsonl (mode='exact')."""
import json
import math
import os

import numpy as np

rows = {}
with open(os.path.join(os.path.dirname(__file__), "data", "library_he_fix.jsonl")) as fh:
    for line in fh:
        r = json.loads(line)
        if r["mode"] == "exact":
            rows[r["n_max"]] = {**rows.get(r["n_max"], {}), **r}


def fit(ns, key):
    Q = np.array([rows[n]["Q"] for n in ns], float)
    y = np.array([rows[n][key] for n in ns], float)
    a, b = np.polyfit(np.log(Q), np.log(y), 1)
    pred = a * np.log(Q) + b
    ss = np.sum((np.log(y) - pred) ** 2)
    st = np.sum((np.log(y) - np.log(y).mean()) ** 2)
    return a, 1 - ss / st, np.abs(np.log(y) - pred).max()


for ns, lab in (((2, 3, 4, 5), "4pt"), ((2, 3, 4), "3pt")):
    for key in ("n_pauli_nonid", "lambda_nonid"):
        a, r2, mr = fit(ns, key)
        print(lab, key, "alpha=%.4f R2=%.5f maxres=%.4f" % (a, r2, mr))
a, r2, mr = fit((2, 3, 4), "qwc_groups")
print("3pt qwc alpha=%.4f R2=%.5f maxres=%.4f" % (a, r2, mr))
for n in sorted(rows):
    r = rows[n]
    lam = r["lambda_nonid"]
    print(n, r["Q"], r["n_pauli_nonid"], round(lam, 4), "lam/Q", round(lam / r["Q"], 4),
          "r1e-3", math.ceil(lam / math.sqrt(2e-3)), "r1e-6", math.ceil(lam / math.sqrt(2e-6)),
          "qwc", r.get("qwc_groups"),
          "g/t", round(r["qwc_groups"] / r["n_pauli_nonid"], 4) if "qwc_groups" in r else None)
