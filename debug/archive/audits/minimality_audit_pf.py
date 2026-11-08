"""Minimality audit, step 4: eq:pf on the actual function N(D) (no Mellin route involved).
L(D) = int_1^inf e^{-Dx} dx / sqrt((x^2-1)(rho x^2 + 1 - rho)); L^(n) = (-1)^n int x^n (...).
Residual of D rho L'''' + 2 rho L''' + D(1-2rho) L'' + (1-2rho) L' - D(1-rho) L, relative to the
sum of |terms|, at three D and two rho, 30 digits.  Substitution x = cosh(v); v truncated at 12 (integrand below e^-24000)."""
import json
from pathlib import Path
import mpmath as mp
mp.mp.dps = 30
rec = {}
for rho in [mp.mpf(1) / 5, mp.mpf(1) / 3]:
    for D in ["0.3", "1.1", "2.7"]:
        D = mp.mpf(D)
        Ln = [(-1) ** n * mp.quad(lambda v: mp.exp(-D * mp.cosh(v)) * mp.cosh(v) ** n /
                                  mp.sqrt(rho * mp.cosh(v) ** 2 + 1 - rho), [0, 1, 3, 6, 12]) for n in range(5)]
        terms = [D * rho * Ln[4], 2 * rho * Ln[3], D * (1 - 2 * rho) * Ln[2], (1 - 2 * rho) * Ln[1], -D * (1 - rho) * Ln[0]]
        r = abs(sum(terms)) / sum(abs(x) for x in terms)
        rec[f"rho={mp.nstr(rho,4)},D={mp.nstr(D,3)}"] = mp.nstr(r, 5)
        print(mp.nstr(rho, 4), mp.nstr(D, 3), "relative residual", mp.nstr(r, 5))
Path(__file__).parent.joinpath("data", "minimality_audit_pf.json").write_text(json.dumps(rec, indent=1))
