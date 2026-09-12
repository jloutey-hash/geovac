"""A second retired value alive in production code.

`geovac/sturmian_secular.py`'s docstring claim (iii) still asserts cond(S)
"climbs from ~4 into the thousands (the paper's '4 -> 3673')".  The /qa run of
2026-09-07 graded that a LARGE: 3673 is entirely a radial-box artifact -- it
switches on exactly where n_max^2 first exceeds R_MAX, and every box agrees to
4 digits below that point.  Converged cond(S) is 16.0 at K=100 and grows mildly
(~0.12*K), an ordinary Gram matrix.  It was withdrawn from the paper and the
test replaced, but the docstring was never swept -- the owner-corrected /
citer-stale pattern again, this time into code, where no paper gate looks.

Also records what the scale-lock work adds: the ill-conditioning was never the
reason to avoid the metric.
"""
import io

SEC = "geovac/sturmian_secular.py"
sec = io.open(SEC, encoding="utf-8").read()

OLD = """  (iii) restoring the L2 overlap metric S (generalized eigenproblem ``M B = p S B``)
        is ill-conditioned -- ``cond(S)`` climbs from ~4 into the thousands (the
        paper's "4 -> 3673");"""

NEW = """  (iii) restoring the L2 overlap metric S turns the problem into a generalized
        eigenproblem ``M B = p S B``.  NOTE ``cond(S)`` is MILD, not explosive:
        16.0 at K=100 on a converged radial box, growing ~0.12*K -- an ordinary
        Gram matrix.  (``4 -> 3673`` is RETIRED:  a pure radial-box artifact that
        switches on exactly where n_max^2 first exceeds R_MAX, every box agreeing
        to 4 digits below that point.)  The reason to keep the metric-free form
        is therefore NOT S's conditioning;  it is that the metric-free form
        exists only at the locked scale lambda = p_kappa, which costs the
        accuracy floor -- Paper 60 ``eq:scale_lock``;"""

assert OLD in sec, "docstring claim (iii) not found verbatim"
sec = sec.replace(OLD, NEW, 1)
io.open(SEC, "w", encoding="utf-8").write(sec)
print("geovac/sturmian_secular.py: retired cond(S) '4 -> 3673' corrected in docstring")
