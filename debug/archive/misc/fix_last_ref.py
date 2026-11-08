"""Last invented reference: sec:mu_level3 -> sec:mu."""
import io
import sys

P = "papers/group3_foundations/paper_18_exchange_constants.tex"

with io.open(P, encoding="utf-8") as fh:
    s = fh.read()

if "sec:mu_level3" not in s:
    print("already fixed")
    sys.exit(0)

s = s.replace(r"Sec.~\ref{sec:mu_level3}", r"Sec.~\ref{sec:mu}")

with io.open(P, "w", encoding="utf-8") as fh:
    fh.write(s)

print("fixed; remaining occurrences:", s.count("sec:mu_level3"))
