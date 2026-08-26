"""Category 1 -- the collinear three-centre observable T2 (Paper 59).

QUOTED, NOT RECOMPUTED.  The 66-digit value is the product of a parallel
high-precision campaign whose full accounting lives in
``debug/beta2_track_a_findings.md`` (Sec. 4, Sec. 9, Sec. 10).  Reproducing it
here would mean re-running the whole parallel arbitrary-precision campaign;
this module therefore
transcribes the value together with the certification summary, and points at the
drivers that produced it.  Nothing here is re-derived, so nothing here can
silently drift away from the campaign's own numbers.

Object (Paper 59, eq:kw and the collinear observable of sec:modular):

    T2 = (8/pi) int_0^1 int_0^1 J(s,t) ds dt ,
    J(s,t) = int_0^inf j0(k(s+t)) P(s,k) P(t,k) dk ,
    P(x,k) = c e^{-D} (D^-3 + 3 D^-4 + 3 D^-5) ,  c = x(1-x) , D = sqrt(c k^2 + 1) .
"""
from __future__ import annotations

from typing import Any, Dict, List

from ._common import entry

# 66 cross-validated significant digits.  Transcribed from
# debug/beta2_track_a_findings.md Sec. 10 (runs hi1 / hi2).
T2_VALUE_66 = ("0.395355765901713964325229296804847564260563977867"
               "082108935234265469")

# The superseded corpus anchor, kept so the correction is visible in the table.
T2_ANCHOR_SUPERSEDED = "0.3953557659017139641"


def build(mode: str = "full") -> List[Dict[str, Any]]:
    del mode  # nothing here is recomputed, so the mode does not matter
    rows: List[Dict[str, Any]] = []

    rows.append(entry(
        "T2.collinear.66",
        "T2",
        "T2, the Paper 59 collinear three-centre observable",
        T2_VALUE_66,
        66,
        "Exact (k,w) factorisation (Paper 59 eq:kw): the j0 kernel is written as "
        "int_0^1 cos(z w) dw, which turns the coupling through b = s+t into a "
        "phase that factorises, so the outer (s,t) double integral collapses to "
        "the square of a one-dimensional integral: "
        "T2 = (8/pi) int_0^inf dk int_0^1 dw cos(kw) R(k,w)^2 with "
        "R(k,w) = int_0^1 cos(kw(s-1/2)) P(s,k) ds.  The k-integral is truncated "
        "at K and the remainder int_K^inf is supplied in closed form by a Watson "
        "(large-k) expansion whose terms reduce to k^-p times {1, cos k, sin k, "
        "cos 2k, sin 2k}, each integrable to incomplete Gamma functions.  Runs "
        "executed in cost-equalised parallel chunks.  Drivers: "
        "debug/beta2_t2_kw_core.py, debug/beta2_t2_tail.py, "
        "debug/beta2_t2_kw_chunk.py, debug/beta2_assemble.py.",
        "Two independent parallel configurations agree to 1.68e-67, i.e. 66 "
        "significant digits: hi1 (working precision 82 digits, K=320, panel width 4, "
        "76 nodes/panel, tail order 90) and hi2 (86 digits, K=280, panel width 5, "
        "96 nodes/panel, tail order 110), with different node-count laws, different "
        "graded s-maps and different chunk boundaries.  Both reproduce, digit for "
        "digit, the 50 printed digits of three further parameter-disjoint runs "
        "(K = 160 / 180 / 190, panel widths 4 / 3 / 6, 64 / 68 / 96 nodes per panel), "
        "which are bit-identical to each other.  Six runs in total span "
        "K in {160,180,190,200,250,280,320}.  The certification is DECOMPOSED, not "
        "two-complete-pipelines: (i) the (k,w) representation itself is checked "
        "against a direct two-dimensional quadrature of J using j0 itself -- no "
        "cosine representation, no w-integral, no symmetry folding -- agreeing to "
        "69-96 digits at k = 2, 50, 150, 300; (ii) the analytic tail is checked "
        "four ways, the sharpest being K-independence across the six runs, which "
        "bounds its relative error below 4.4e-33; (iii) the outer quadrature is "
        "checked by the six-run parameter disjointness above.  A fully independent "
        "end-to-end pipeline (the (s,t)-outer route with an analytic fibre tail) "
        "currently confirms only ~19 digits; that ceiling is a property of the old "
        "frame's oscillatory corners, not of this value.  Full accounting: "
        "debug/beta2_track_a_findings.md Sec. 3, 4, 9, 10.",
        value_kind="decimal",
        transcendence_class="open -- not a low-height element of any tested "
                            "period ring (see pslq_status)",
        pslq_status=(
            "Guarded, decoy-calibrated PSLQ at 64 digits and two working "
            "precisions: DECISIVE-NEGATIVE in every ring of dimension <= 20 at "
            "its calibrated height budget, including the corrected weight-<=3 "
            "ring {pi, K(1/2)^{+-1}, G} (dimension 20, height <= 10) and eight "
            "dedicated Catalan probes (height <= 1e12).  No candidate relation "
            "in any ring at any precision.  The negative is height-bounded: it "
            "excludes a clean closed form in those rings, not a large-height one."
        ),
        supersedes=T2_ANCHOR_SUPERSEDED,
        supersedes_note=(
            "The frozen corpus anchor 0.3953557659017139641 is correct to 18 "
            "significant digits and wrong in the 19th: the correct continuation "
            "is ...39643252... .  Paper 59 sec:modular still carries the old "
            "anchor; correcting it is a named follow-on and is NOT done here "
            "(this task makes no paper edits)."
        ),
        drivers=[
            "debug/beta2_t2_kw_core.py",
            "debug/beta2_t2_tail.py",
            "debug/beta2_t2_kw_chunk.py",
            "debug/beta2_assemble.py",
            "debug/beta2_t2_direct2d.py",
            "debug/beta2_t2_st_route.py",
        ],
        source_memo="debug/beta2_track_a_findings.md",
    ))

    rows.append(entry(
        "T2.collinear.120.pending",
        "T2",
        "T2 at ~120 digits (PENDING -- placeholder row)",
        "PENDING",
        0,
        "Same (k,w) factorisation, pushed to working precision 145 digits, "
        "K=580, tail order 220, in 12 cost-equalised parallel chunks "
        "(configuration 'u1'; see debug/beta2_u1_README.md).  A second "
        "configuration 'u2' (150 digits, K=560, panel width 5, 112 nodes/panel) "
        "is required before any digit beyond 66 may be claimed.",
        "NOT YET CERTIFIED -- and therefore claiming zero digits.  Run u1 was in "
        "flight when this table was generated; its cross-validation partner u2 "
        "has not been run at all.  The purpose of the extension is a "
        "height-<=1e4 PSLQ verdict on the pre-registered dimension-20 ring, "
        "which needs roughly 120 digits: at 64 digits that ring only supports a "
        "height-10 verdict.  Until u1 and u2 both land and agree, the certified "
        "value of T2 remains the 66-digit entry above.",
        value_kind="decimal",
        source_memo="debug/beta2_u1_README.md",
        sort_index=1,
    ))

    return rows
