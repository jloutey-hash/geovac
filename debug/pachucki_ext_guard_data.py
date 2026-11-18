"""Build tests/data/paper58_pachucki_ext_2026_11_18.json from the derivation's outputs in debug/data
(the guard in tests/ must not read debug/).  Re-run after any of the inputs changes.
  sym    <- pachucki_ext_symbolic_1s.json, pachucki_ext_sym_<case>.json   (closed_form_coefficients_t)
  limit  <- pachucki_ext_limit.json, pachucki_ext_limit_<case>.json       (limit)
  denoms <- pachucki_ext_denoms_<case>.json                               (factors, keys)
Usage: python debug/pachucki_ext_guard_data.py
"""
import json
import os

D = os.path.join("debug", "data")
CASES = ("1s", "2psig", "3dsig", "ppsig", "pppi0")
out = {"provenance": "Exact outputs of the derivation debug/pachucki_extension_memo.md (2026-11-18), Q(t) runs with "
                     "a = 1, b = t. sym: coefficients of K on the nine functions, {name: {power of R: rational "
                     "function of t}}. limit: the eps^0 term of K at t = 1 + eps in e2 = e^{-2R}, Lg = gamma + ln R, "
                     "X2 = E1(2R), X4 = e^{2R} E1(4R). denoms: irreducible factors (largest multiplicity) of the "
                     "coefficient denominators of every intermediate f_alpha, F_alpha, d_alpha p^2, and every "
                     "function key. Recomputed by the slow tier of tests/test_paper58_pachucki_extension.py.",
       "sym": {}, "limit": {}, "denoms": {}}
for c in CASES:
    f = "pachucki_ext_symbolic_1s.json" if c == "1s" else f"pachucki_ext_sym_{c}.json"
    out["sym"][c] = json.load(open(os.path.join(D, f)))["closed_form_coefficients_t"]
    f = "pachucki_ext_limit.json" if c == "1s" else f"pachucki_ext_limit_{c}.json"
    out["limit"][c] = json.load(open(os.path.join(D, f)))["limit"]
    f = os.path.join(D, f"pachucki_ext_denoms_{c}.json")
    if os.path.exists(f):
        J = json.load(open(f))
        assert J["verdict"] == "PASS", (c, J["verdict"])
        keys = [{k: v for k, v in row.items() if k in ("kind", "lam", "mu", "log_power")} for row in J["keys"]]
        out["denoms"][c] = {"multi_indices": J["multi_indices"],
                            "denominator_factors_max_multiplicity": J["denominator_factors_max_multiplicity"],
                            "keys": keys}
    else:
        print("MISSING denoms for", c)
os.makedirs(os.path.join("tests", "data"), exist_ok=True)
with open(os.path.join("tests", "data", "paper58_pachucki_ext_2026_11_18.json"), "w", newline=chr(10)) as fh:
    json.dump(out, fh, indent=1, sort_keys=True)
    fh.write(chr(10))
print("wrote", {k: sorted(v) for k, v in out.items() if isinstance(v, dict)})
