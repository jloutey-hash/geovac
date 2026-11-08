"""Revisit the Paper-60 done-record, then stamp it in C16's `cited_by`.

The record's "What the paper now says" block ratified the floor claim that the
scale-lock diagnostic retires -- the same "a .done.md ratifying a retired claim"
class the 2026-09-08 /qa paper_61 run named.  Revisit first, stamp second;
stamping without revisiting is exactly what the cited_by rule forbids.
"""
import io

DONE = "docs/qa/paper_60.done.md"
RET = "debug/qa/check_retracted_terms.py"

done = io.open(DONE, encoding="utf-8").read()

OLD = """> And the cheap rule is the one that stops converging:\\ it saturates at
> **6.44 mHa, 4.0× chemical accuracy**, at any basis size. Cost growth and
> attainable accuracy are one fact."""

NEW = """> And the cheap rule carries an accuracy floor of **6.44 mHa**.
>
> **Superseded 2026-09-08 — the floor's mechanism and universality.** This
> record previously ratified "4.0× chemical accuracy, at any basis size. Cost
> growth and attainable accuracy are one fact." The floor VALUE stands; its
> attribution and its universality do not. (i) The mechanism is the **scale
> lock**, not `l_max`: metric-free holds iff `E = -λ²/2`, hence
> `λ = p_κ`, and that is not the variational optimum — freeing `λ`
> over the *identical* span reaches 1.28 mHa at K=130 against 7.46 locked, and
> 0.15 mHa of the independently known s-limit at K=136. The He `l≥4`
> partial-wave tail is 0.37–0.53 mHa, an order below the floor, so angular
> truncation cannot be it. (ii) The floor is **ground-state specific**: at
> K=202, with `‖M‖₁` identical because it does not depend on which root
> is extracted, the ground state sits 4.49× above chemical accuracy and
> 2¹S sits 1.12×; the posing cost falls 3–4× per rung up the ¹S
> ladder. (iii) The price of freeing the scale is the whole encoding advantage,
> `‖·‖₁` from `K^0.72` to `K^2.75`. Backing:
> `tests/test_paper60_scale_lock.py`; drivers `debug/p60_{variational_probe,
> scale_scan,freescale_resource,posing_cost_by_state,excited_ladder}.py`."""

assert OLD in done, "done-record locus not found verbatim"
done = done.replace(OLD, NEW, 1)
io.open(DONE, "w", encoding="utf-8").write(done)
print("done-record: floor block revisited and superseded-note added")

ret = io.open(RET, encoding="utf-8").read()
OLD_STAMP = '            "docs/qa/paper_60.done.md": None,'
NEW_STAMP = '            "docs/qa/paper_60.done.md": "reviewed 2026-09-08",'
assert OLD_STAMP in ret, "cited_by stamp locus not found"
ret = ret.replace(OLD_STAMP, NEW_STAMP, 1)
io.open(RET, "w", encoding="utf-8").write(ret)
print("C16: docs/qa/paper_60.done.md stamped reviewed 2026-09-08")
