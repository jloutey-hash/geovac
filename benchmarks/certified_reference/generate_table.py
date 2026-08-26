"""Generate the GeoVac certified reference-value table.

    python -m benchmarks.certified_reference.generate_table          # full run
    python -m benchmarks.certified_reference.generate_table --fast   # skip the
                                                                     # slow
                                                                     # independent
                                                                     # quadratures

Writes two artefacts:

    benchmarks/certified_reference/certified_reference_values.json   (machine)
    docs/certified_reference_values.md                               (human)

The run is deterministic: every number comes from exact symbolic evaluation,
exact rational arithmetic, arbitrary-precision arithmetic at a fixed working
precision, or deterministic adaptive quadrature at fixed tolerances.  There is
no sampling, no fitting, and no random seed anywhere in the pipeline.

CERTIFICATION DISCIPLINE.  Each entry carries four mandatory fields --
``value``, ``digits_claimed``, ``method``, ``evidence`` -- and the rule is that
``digits_claimed`` may never exceed what ``evidence`` actually establishes.
Where a value is exact (a rational, an algebraic number, an integer tuple)
``digits_claimed`` is the string ``"exact"`` rather than a number.  Where a value
is quoted from a corpus campaign rather than recomputed here, the evidence field
says so and cites the campaign's own accounting.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import date
from pathlib import Path
from typing import Any, Dict, List

REPO_ROOT = Path(__file__).resolve().parents[2]
JSON_PATH = REPO_ROOT / "benchmarks" / "certified_reference" / \
    "certified_reference_values.json"
MD_PATH = REPO_ROOT / "docs" / "certified_reference_values.md"

CATEGORY_TITLES = {
    "T2": "1. The collinear three-centre observable T2",
    "two_center_eri": "2. Two-centre hydrogenic electron-repulsion integrals",
    "slater_rational": "3. One-centre Slater repulsion integrals (exact rationals)",
    "qfd_diatomic": "4. Certified minimal-basis diatomic total energies",
    "helium_ci": "5. Graph-native helium CI ground states (fixed basis)",
    "resurgent": "6. Resurgent and connection data",
    "anchor": "7. Anchor constants",
}

CATEGORY_BLURBS = {
    "T2": (
        "One number, and the single hardest to produce in this table.  T2 is the "
        "collinear limit of the two-electron three-centre integral that blocks "
        "polyatomic closed-form evaluation in the GeoVac framework (Paper 59).  "
        "It has no known closed form; what is offered here is a certified "
        "numerical value, which is what a period-recognition search (PSLQ or "
        "similar) needs as input.  The value is QUOTED from the campaign that "
        "produced it, not recomputed by this generator: recomputing it would "
        "mean re-running that whole parallel arbitrary-precision campaign."
    ),
    "two_center_eri": (
        "Closed-form values for the four two-centre electron-repulsion integral "
        "classes over hydrogenic (Slater-type) orbitals.  These are the entries "
        "most directly useful to someone validating an integral code: each is an "
        "exact symbolic expression evaluated to 50 digits, and each is separately "
        "checked against a numerical quadrature that shares no code with it.  "
        "The transcendence class column records which special functions the "
        "closed form actually needs -- this is a structural property of the "
        "class, not of the particular numbers.\n\n"
        "**Domain restriction, measured while building this table.**  For "
        "`l_a, l_b > 0` the hybrid closed form goes through the shell "
        "reformulation, and that route requires `Z_B < Z_A` strictly.  At "
        "`Z_B = Z_A` it returns NaN -- a removable coincidence of decay rates, "
        "since approaching `Z_B -> Z_A` from below converges to the quadrature "
        "value -- and at `Z_B > Z_A` it raises `AssertionError: Ei reached step "
        "3`, because the exponential-integral argument turns negative and that "
        "branch is not implemented.  The `l > 0` hybrid rows therefore use "
        "`(Z_A, Z_B)` in `{(3,1), (4,2)}` rather than `{(1,1), (3,1)}`.  This is "
        "a coverage limit of the evaluator, not a wrong value: no closed form "
        "and quadrature anywhere in this table disagree."
    ),
    "slater_rational": (
        "The cleanest entries here.  One-centre Slater radial repulsion integrals "
        "at unit orbital exponent are RATIONAL NUMBERS, computed in exact "
        "arithmetic.  They carry no digit count because they are not "
        "approximations.  Scale to arbitrary nuclear charge by R^k(Z) = Z R^k(1)."
    ),
    "qfd_diatomic": (
        "Total energies of completely specified minimal-basis diatomic models, "
        "assembled with NO quadrature anywhere on the production path: every "
        "one-electron matrix element and every two-electron class is an exact "
        "symbolic expression, so the whole energy can be evaluated to any "
        "requested precision.  These rows are end-to-end reference data -- an "
        "error anywhere in an integral code moves the energy -- and they are "
        "the only entries in the table that exercise the full pipeline "
        "(integrals, Loewdin orthogonalisation, full CI) rather than a single "
        "integral.\n\n"
        "**What they are not.**  A two- or three-function s-only basis with "
        "unoptimised hydrogenic exponents is far from the exact energy, and no "
        "accuracy claim is made for any of these numbers.  What is certified is "
        "the value of the stated model.  Each entry carries an "
        "`honest_scope` field in the machine-readable table saying so for that "
        "system, with the corresponding textbook or exact figure for context "
        "where one exists.\n\n"
        "**Why the digit counts differ so much.**  The exchange class is built "
        "from a Neumann expansion in an index `tau`, and that expansion "
        "TERMINATES exactly when the two centres of a charge density carry the "
        "same orbital exponent (`q = (alpha - beta) R / 2 = 0`).  Homonuclear "
        "systems at equal exponent (`H2+`, `He2^2+`, and the companion `H2`) "
        "therefore have a FINITE sum and no truncation error at all: their "
        "digits are limited only by working precision.  Heteronuclear systems "
        "(`HeH+`, `BeH+`, and the companion `LiH`) have an infinite sum "
        "truncated at a per-quartet `tau_max`, and their claims are stated NET "
        "OF an explicit geometric tail bound propagated to the energy through "
        "the two-particle density matrix and the Loewdin transform, with the "
        "amplification factor recomputed for each system."
    ),
    "helium_ci": (
        "Ground-state energies of the GeoVac graph-native two-electron CI "
        "matrix for helium at fixed basis truncations `n_max`, in the singlet "
        "`M_L = 0` sector.\n\n"
        "**These certify the assembly, not the accuracy.**  Each value is the "
        "exact lowest eigenvalue of one specific finite matrix.  The physical "
        "non-relativistic infinite-mass helium ground state is "
        "`-2.903724377034119598` Ha (Pekeris/Drake), and the graph-native CI "
        "approaches it only slowly from above -- 0.19 per cent at `n_max = 7` "
        "(Paper 13).  So thirty-five certified digits here are thirty-five "
        "digits of a truncated-basis eigenvalue whose first two digits already "
        "differ from helium.  They are reference values for someone "
        "reimplementing the construction, in the same spirit as a published "
        "FCI energy in a stated finite basis -- not accuracy claims.\n\n"
        "What makes them certifiable is that the construction contains no "
        "quadrature.  The one-body diagonal is `-Z^2/(2n^2)`; the one-body "
        "off-diagonal is `kappa * (-A_ij) = +1/16` on the edges of the binary "
        "S^3 lattice; the radial Slater integrals are exact rationals; the "
        "orbital-exponent scaling is the integer `Z`; and each Gaunt angular "
        "factor is `(rational) * sqrt(rational)`.  Every matrix entry is "
        "therefore an exact algebraic number, and so is the eigenvalue.  The "
        "generator assembles the matrix in exact arithmetic, rounds it to "
        "`mpf` once at a chosen working precision, and reports both the "
        "two-precision agreement and the rigorous symmetric residual bound "
        "`|lambda_min - v^T H v| <= ||H v - lambda v||`.  At `n_max = 1` the "
        "sector holds one configuration and the answer is the rational "
        "`-11/4` outright."
    ),
    "resurgent": (
        "Connection data for the corpus's divergent-series objects.  A divergent "
        "asymptotic series still determines its function once one knows where its "
        "Borel transform is singular and with what amplitude; those amplitudes "
        "are the Stokes constants.  The corpus result is that these come out "
        "algebraic up to a power of pi fixed by the singularity type, with the "
        "transcendental content displaced to a single boundary period.  These "
        "entries are the numbers in that statement."
    ),
    "anchor": (
        "The fixed constants the rest of the table is stated against, at a "
        "uniform 50 digits, each with the relation that defines it.  Included so "
        "that a reader reproducing an entry does not have to guess a convention "
        "-- in particular whether an elliptic-integral argument is the parameter "
        "m or the modulus k."
    ),
}


def build_all(mode: str = "full") -> List[Dict[str, Any]]:
    """Build every entry.  ``mode`` is 'full' or 'fast'."""
    from . import (entries_anchors, entries_helium, entries_qfd_diatomic,
                   entries_resurgent, entries_slater, entries_t2,
                   entries_two_center)
    from ._common import sort_entries

    rows: List[Dict[str, Any]] = []
    rows += entries_t2.build(mode)
    rows += entries_two_center.build(mode)
    rows += entries_slater.build(mode)
    rows += entries_qfd_diatomic.build(mode)
    rows += entries_helium.build(mode)
    rows += entries_resurgent.build(mode)
    rows += entries_anchors.build(mode)
    return sort_entries(rows)


# ---------------------------------------------------------------------------
# Fast subset used by tests/test_certified_reference_values.py
# ---------------------------------------------------------------------------
FAST_SUBSET_IDS = [
    # category 2 -- one entry per closed-form class
    "eri.aabb.Z11.R1.4",
    "eri.hybrid_s.Z31.R2.0",
    "eri.hybrid_p.Z31.R2.0",
    "eri.exchange.Z11.R1.4",
    # category 3 -- exact rationals, both the float path and the exact path
    "slater.R0.1010_1010",
    "slater.R2.2121_2121",
    "slater.R0.5050_5050",
    # category 5 -- graph-native helium CI: the exact-rational row and the
    # cheapest certified row
    "he.gnci.Z2.nmax1",
    "he.gnci.Z2.nmax2",
    # category 7 -- anchors
    "anchor.K_half",
    "anchor.collapse_pi2_24",
    "anchor.gerade_constant",
]


def fast_subset(dps: int = 30) -> Dict[str, str]:
    """Recompute a cheap subset of the table at reduced precision.

    Returns ``{entry_id: value_string}``.  Exact-rational entries come back as
    the fraction string; everything else as a decimal printed to ``dps``
    significant digits.  The point is that this is cheap to run while
    still exercising the same code paths as the full generator.
    """
    from fractions import Fraction

    import mpmath as mp
    import sympy as sp

    from geovac.hypergeometric_slater import compute_rk_algebraic
    from geovac.two_center_eri import (aabb_closed_form, hybrid_closed_form,
                                       ordered_xi_closed)

    from ._common import frac_str

    out: Dict[str, str] = {}

    def dec(expr: sp.Expr) -> str:
        with mp.workdps(dps + 10):
            return mp.nstr(mp.mpf(str(sp.N(sp.re(expr), dps + 10))), dps,
                           strip_zeros=False)

    R14 = sp.Rational(7, 5)
    R20 = sp.Integer(2)

    out["eri.aabb.Z11.R1.4"] = dec(aabb_closed_form(
        Fraction(1), (1, 0, 0), (1, 0, 0), Fraction(1), (1, 0, 0), (1, 0, 0), R14))
    out["eri.hybrid_s.Z31.R2.0"] = dec(hybrid_closed_form(
        Fraction(3), (1, 0, 0), (1, 0, 0), (1, 0, 0), Fraction(1), (1, 0, 0), R20))
    out["eri.hybrid_p.Z31.R2.0"] = dec(hybrid_closed_form(
        Fraction(3), (2, 1, 0), (2, 1, 0), (1, 0, 0), Fraction(1), (1, 0, 0), R20))
    out["eri.exchange.Z11.R1.4"] = dec(ordered_xi_closed(R14, R14))

    out["slater.R0.1010_1010"] = frac_str(
        compute_rk_algebraic(1, 0, 1, 0, 1, 0, 1, 0, 0))
    out["slater.R2.2121_2121"] = frac_str(
        compute_rk_algebraic(2, 1, 2, 1, 2, 1, 2, 1, 2))
    out["slater.R0.5050_5050"] = frac_str(
        compute_rk_algebraic(5, 0, 5, 0, 5, 0, 5, 0, 0))

    from .entries_helium import (_exact_rational_ground_state,
                                 build_matrix_exact, ground_state_mp,
                                 realise_mp)

    out["he.gnci.Z2.nmax1"] = frac_str(_exact_rational_ground_state(2))
    he_entries, he_cfg = build_matrix_exact(2, 2)
    with mp.workdps(dps + 10):
        he_lam, _r, _s = ground_state_mp(realise_mp(he_entries, len(he_cfg)))
        out["he.gnci.Z2.nmax2"] = mp.nstr(he_lam, dps, strip_zeros=False)

    with mp.workdps(dps + 10):
        out["anchor.K_half"] = mp.nstr(mp.ellipk(mp.mpf(1) / 2), dps,
                                       strip_zeros=False)
        out["anchor.collapse_pi2_24"] = mp.nstr(mp.pi ** 2 / 24, dps,
                                                strip_zeros=False)
        x_star = mp.findroot(lambda x: mp.tan(x) - x, mp.mpf("4.4934"))
        out["anchor.gerade_constant"] = mp.nstr(
            2 / (1 + mp.sin(x_star) / x_star), dps, strip_zeros=False)

    assert set(out) == set(FAST_SUBSET_IDS), "fast subset drifted from its id list"
    return out


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
PREAMBLE = """# GeoVac certified reference values

High-precision values for integrals and connection data that, as far as the
project is aware, no other implementation currently produces -- offered as
reference data for people validating electron-repulsion integral codes and for
the precision / Bessel-moment community.

**This file is generated.**  Do not edit it by hand; edit the generators under
`benchmarks/certified_reference/` and re-run

```
python -m benchmarks.certified_reference.generate_table
```

A machine-readable copy of the same table lives at
`benchmarks/certified_reference/certified_reference_values.json`.

## The certification discipline

The one rule this table exists to enforce is: **never claim more digits than the
cross-validation supports.**  Every entry therefore carries four fields.

| field | meaning |
|:--|:--|
| `value` | the number, printed to at most `digits_claimed` significant digits |
| `digits_claimed` | how many of those digits are certified; the string `exact` when the value is a rational, an algebraic number or an integer tuple, and therefore has no digit count at all |
| `method` | how the number was produced -- which closed form, which representation, which route |
| `evidence` | what certifies `digits_claimed`, stated with the measured agreements, not asserted |

Three things this discipline deliberately keeps apart:

1. **Internal precision is not verification.**  Evaluating one exact expression
   at two working precisions shows that the printed digits are the digits of
   *that expression*.  It says nothing about whether the expression is the
   integral.  Where an entry rests on such a check, the evidence field says so
   and reports the independent-route agreement separately.
2. **An independent route caps what may be claimed as a verified integral.**  For
   the two-centre integrals the independent route is a float64 quadrature good
   to roughly 1e-10..1e-15; that number appears in every such entry as
   `independent_route_rel_agreement`, and is not silently rounded up into the
   50-digit claim.
3. **Decomposed certification is labelled as such.**  The T2 entry does not have
   two complete independent pipelines agreeing to 66 digits -- it has one
   pipeline in six parameter-disjoint configurations, plus separate
   high-precision closure of each failure mode the two-pipeline criterion is a
   proxy for.  The entry states this in full, along with the digit count
   (about 19) that a fully independent end-to-end route currently reaches.

Exact entries are the strongest rows in the table: a rational number has
infinitely many correct digits and no convergence question.  They are marked
`exact` rather than given a large digit count.

## Conventions

* Atomic units throughout: lengths in bohr, energies in hartree.
* Hydrogenic (Slater-type) orbitals `chi_{nlm}` at nuclear charge `Z`, with the
  standard normalisation `R_{nl}(r) = N_{nl} (2Zr/n)^l e^{-Zr/n} L_{n-l-1}^{2l+1}(2Zr/n)`.
* Electron-repulsion integrals in chemists' notation `(ab|cd)`, with `a, b` on
  the first-named centre.
* `K(m)` is the complete elliptic integral of the first kind in the **parameter**
  convention, `K(m) = int_0^{pi/2} dtheta / sqrt(1 - m sin^2 theta)` -- so
  `K(1/2)`, not `K(k = 1/sqrt(2))`, is the lemniscatic value quoted below.
* `E_1(z) = int_z^inf e^{-t}/t dt` is the exponential integral; `gamma` is
  Euler's constant.
* Two centres are placed at the origin and at `R zhat`.

## Known limits, and what changed

Three things a reader should know before using the table.

1. **One published value is corrected here.**  The project's own frozen anchor
   for T2, `0.3953557659017139641`, is right to 18 significant digits and wrong
   in the 19th; the correct continuation is `...39643252...`.  The 19-digit
   anchor was over-claimed relative to what its cross-validation supported.  The
   table publishes the corrected value with its full accounting.
2. **One evaluator has a domain restriction.**  The hybrid two-centre class with
   angular momentum on the one-centre pair requires `Z_B < Z_A`; see the note in
   category 2.  Values inside that domain are unaffected.
3. **One row is deliberately empty.**  A ~120-digit extension of T2 was still
   running when this table was generated.  It appears as a placeholder claiming
   zero digits rather than as a partial number, because a run without its
   cross-validation partner certifies nothing.

No entry in this table has a closed form and an independent quadrature
disagreeing beyond the tolerance the entry advertises.

## A few terms, in plain form

Some of the language below comes from the theory of periods rather than from
quantum chemistry; the entries are usable without it, but here is what it means.

* **Period.**  A number obtained by integrating a rational function over a region
  cut out by polynomial inequalities.  `pi`, `log 2` and the values of elliptic
  integrals are periods; they form a countable ring, and asking which period a
  computed number is amounts to asking what kind of geometry produced it.
* **Weight.**  A grading on periods that counts, roughly, how many nested
  integrations are needed.  `log` is weight one, the dilogarithm and `zeta(2)`
  weight two, and so on.  The statement that a class of integrals "closes at
  weight one" means it needs nothing beyond logarithms and exponential
  integrals -- no dilogarithm.
* **Height, and what a negative result means.**  A search for a closed form (by
  integer-relation algorithms such as PSLQ) asks whether the number is a
  rational combination of a fixed list of periods.  Such a search can only rule
  out combinations whose integer coefficients are smaller than some bound, the
  *height*, and that bound is set by how many digits you have and how long the
  list is.  So "no closed form at height 10" is a real statement with a real
  limit, not a proof of impossibility -- which is exactly why more digits are
  worth producing.
* **CM point / CM fibre.**  A member of a family of elliptic curves with extra
  symmetry (complex multiplication).  At such a point the period simplifies to a
  ratio of Gamma-function values -- which is why `K(1/2) = Gamma(1/4)^2 /
  (4 sqrt(pi))` appears in the anchors.

## What "transcendence class" means

Each closed form needs a specific, finite set of special functions -- its seeds.
The class is a structural fact about the integral class, and it is read off the
built expression rather than asserted: `{exp}` means the closed form is
elementary, `{exp, E_1, ln}` means it needs the exponential integral and a
logarithm, and so on.  For the integral classes below the seed set grows with
difficulty: `(AA|BB)` is elementary at any angular momentum, the hybrid class
picks up `E_1` and a logarithm once the one-centre pair has `l > 0`, and the
exchange class additionally carries Euler's `gamma`.  All of them close at
weight one -- no dilogarithm appears anywhere.
"""


def _fmt_digits(d: Any) -> str:
    return "exact" if d == "exact" else str(d)


def render_markdown(rows: List[Dict[str, Any]], meta: Dict[str, Any]) -> str:
    out: List[str] = [PREAMBLE]

    counts: Dict[str, int] = {}
    for r in rows:
        counts[r["category"]] = counts.get(r["category"], 0) + 1
    out.append("\n## Contents\n")
    out.append("| category | entries |")
    out.append("|:--|--:|")
    for cat, title in CATEGORY_TITLES.items():
        out.append(f"| {title} | {counts.get(cat, 0)} |")
    out.append(f"| **total** | **{len(rows)}** |")
    out.append("")
    out.append(f"Generated {meta['generated']} from commit-time corpus state; "
               f"generator mode `{meta['mode']}`, "
               f"Python {meta['python']}, mpmath {meta['mpmath']}, "
               f"sympy {meta['sympy']}.")
    out.append("")

    for cat, title in CATEGORY_TITLES.items():
        sub = [r for r in rows if r["category"] == cat]
        if not sub:
            continue
        out.append(f"\n---\n\n## {title}\n")
        out.append(CATEGORY_BLURBS[cat])
        out.append("")
        out.append("| entry | value | digits | class |")
        out.append("|:--|:--|:--|:--|")
        for r in sub:
            val = r["value"]
            short = val if len(val) <= 56 else val[:53] + "..."
            cls = r.get("transcendence_class", "")
            out.append(f"| `{r['id']}` | `{short}` | "
                       f"{_fmt_digits(r['digits_claimed'])} | {cls} |")
        out.append("")
        out.append("### Entry detail")
        out.append("")
        for r in sub:
            out.append(f"#### `{r['id']}`")
            out.append("")
            out.append(f"**{r['label']}**")
            out.append("")
            out.append("```")
            out.append(r["value"])
            out.append("```")
            out.append("")
            out.append(f"* **digits claimed:** {_fmt_digits(r['digits_claimed'])}")
            if r.get("defining_relation"):
                out.append(f"* **defining relation:** `{r['defining_relation']}`")
            if r.get("decimal_value"):
                out.append(f"* **decimal form:** `{r['decimal_value']}`")
            if r.get("transcendence_class"):
                out.append(f"* **transcendence class:** {r['transcendence_class']}")
            if r.get("geovac_entry_point"):
                out.append(f"* **entry point:** `{r['geovac_entry_point']}`")
            if r.get("backing_test"):
                out.append(f"* **backing test:** `{r['backing_test']}`")
            if r.get("source_memo"):
                out.append(f"* **source record:** `{r['source_memo']}`")
            if r.get("supersedes"):
                out.append(f"* **supersedes:** `{r['supersedes']}` -- "
                           f"{r.get('supersedes_note', '')}")
            if r.get("pslq_status"):
                out.append(f"* **period-recognition status:** {r['pslq_status']}")
            out.append("")
            out.append(f"**Method.** {r['method']}")
            out.append("")
            out.append(f"**Evidence.** {r['evidence']}")
            out.append("")
    return "\n".join(out) + "\n"


def main(argv: List[str] | None = None) -> int:
    import mpmath
    import sympy

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fast", action="store_true",
                        help="skip the slow independent quadratures")
    parser.add_argument("--check", action="store_true",
                        help="regenerate and compare against the stored JSON "
                             "without writing")
    args = parser.parse_args(argv)

    mode = "fast" if args.fast else "full"
    rows = build_all(mode)

    meta = {
        "generated": date.today().isoformat(),
        "mode": mode,
        "python": platform.python_version(),
        "mpmath": mpmath.__version__,
        "sympy": sympy.__version__,
        "n_entries": len(rows),
    }
    payload = {"meta": meta, "entries": rows}

    if args.check:
        stored = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        old = {e["id"]: e["value"] for e in stored["entries"]}
        new = {e["id"]: e["value"] for e in rows}
        bad = [k for k in new if k in old and old[k] != new[k]]
        missing = sorted(set(old) - set(new))
        print(f"{len(new)} entries; {len(bad)} value mismatches; "
              f"{len(missing)} missing")
        for k in bad:
            print(f"  MISMATCH {k}\n    stored {old[k]}\n    fresh  {new[k]}")
        return 1 if (bad or missing) else 0

    JSON_PATH.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n",
                         encoding="utf-8")
    MD_PATH.write_text(render_markdown(rows, meta), encoding="utf-8")
    print(f"wrote {JSON_PATH.relative_to(REPO_ROOT)} "
          f"and {MD_PATH.relative_to(REPO_ROOT)}  ({len(rows)} entries, "
          f"mode={mode})")
    for cat, title in CATEGORY_TITLES.items():
        n = sum(1 for r in rows if r["category"] == cat)
        print(f"  {title}: {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
