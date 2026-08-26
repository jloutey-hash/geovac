"""Backing test for the certified reference-value table.

The table (``docs/certified_reference_values.md`` and its machine-readable twin
``benchmarks/certified_reference/certified_reference_values.json``) publishes
high-precision values with an explicit certification discipline.  This test
guards two things:

  1. **The numbers still reproduce.**  A cheap subset spanning categories 2
     (two-centre electron-repulsion integrals), 3 (exact-rational Slater
     integrals), 5 (graph-native helium CI) and 7 (anchor constants) is
     regenerated from the live code at reduced precision and matched against
     the stored table.  If a refactor
     changes what ``geovac`` computes, this fails.
  2. **The discipline is still honoured.**  Every stored entry must carry the
     four mandatory fields, ``digits_claimed`` must be an integer or the string
     ``"exact"``, printed decimals must not exceed their claimed digit count,
     and every two-centre entry must still record an independent-route
     agreement inside the tolerance the table advertises.

Deliberately NOT re-run here: the full generator's independent quadratures
(the expensive part of a full run), and the T2 campaign, which the generator
quotes rather than recomputes in the first place.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import mpmath as mp
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = (REPO_ROOT / "benchmarks" / "certified_reference"
             / "certified_reference_values.json")
MD_PATH = REPO_ROOT / "docs" / "certified_reference_values.md"

FAST_DPS = 30

#: the tolerance the table advertises for the independent float64 quadrature
#: routes.  Loosening this is a decision, not a maintenance chore.
INDEPENDENT_ROUTE_TOL = 1e-9


@pytest.fixture(scope="module")
def table():
    assert JSON_PATH.exists(), f"missing generated table {JSON_PATH}"
    return json.loads(JSON_PATH.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def by_id(table):
    return {e["id"]: e for e in table["entries"]}


# ---------------------------------------------------------------------------
# 1. schema / discipline
# ---------------------------------------------------------------------------
def test_artifacts_exist(table):
    assert MD_PATH.exists(), f"missing human-readable table {MD_PATH}"
    assert table["meta"]["n_entries"] == len(table["entries"])
    assert len(table["entries"]) >= 25, "the table has shrunk below its gate"


def test_every_entry_has_the_four_mandatory_fields(table):
    """value / digits_claimed / method / evidence -- no entry may skip one."""
    for e in table["entries"]:
        for field in ("id", "category", "label", "value", "digits_claimed",
                      "method", "evidence"):
            assert field in e, f"{e.get('id', '?')} missing field {field}"
            assert e[field] not in (None, ""), \
                f"{e.get('id', '?')} has empty field {field}"
        assert len(e["method"]) > 40, f"{e['id']}: method field is a stub"
        assert len(e["evidence"]) > 40, f"{e['id']}: evidence field is a stub"


def test_ids_are_unique(table):
    ids = [e["id"] for e in table["entries"]]
    assert len(ids) == len(set(ids))


def test_digits_claimed_is_an_int_or_exact(table):
    for e in table["entries"]:
        d = e["digits_claimed"]
        assert d == "exact" or isinstance(d, int), \
            f"{e['id']}: digits_claimed = {d!r}"
        if isinstance(d, int):
            assert 0 <= d <= 200, f"{e['id']}: implausible digit count {d}"


def test_printed_digits_never_exceed_the_claim(table):
    """The core discipline: a value may not display more digits than it claims."""
    for e in table["entries"]:
        d = e["digits_claimed"]
        if not isinstance(d, int) or d == 0:
            continue
        if e["value_kind"] not in ("decimal",):
            continue
        digits = re.sub(r"[^0-9]", "", e["value"]).lstrip("0")
        assert len(digits) <= d, (
            f"{e['id']} prints {len(digits)} significant digits but claims {d}")


def test_two_center_entries_record_an_independent_route(by_id):
    """Every closed-form ERI entry must carry a measured quadrature agreement
    inside the advertised tolerance.  This is the table's STOP gate, frozen."""
    checked = 0
    for e in by_id.values():
        if e["category"] != "two_center_eri":
            continue
        rel = e.get("independent_route_rel_agreement")
        assert rel is not None, \
            f"{e['id']}: no independent-route agreement recorded"
        assert rel < INDEPENDENT_ROUTE_TOL, (
            f"{e['id']}: closed form and quadrature disagree at {rel:.2e}, "
            f"above the advertised tolerance {INDEPENDENT_ROUTE_TOL:.0e}")
        checked += 1
    assert checked >= 15, "the two-centre block has shrunk unexpectedly"


def test_slater_entries_are_exact_rationals(by_id):
    from fractions import Fraction
    checked = 0
    for e in by_id.values():
        if e["category"] != "slater_rational":
            continue
        assert e["digits_claimed"] == "exact"
        f = Fraction(e["value"])          # parses "p/q"; raises if it is not one
        assert f > 0
        assert abs(float(f) - float(e["decimal_value"])) < 1e-15 * abs(float(f))
        assert e["quadrature_abs_agreement"] < 1e-12, (
            f"{e['id']}: exact rational and defining-integral quadrature "
            f"disagree at {e['quadrature_abs_agreement']:.2e}")
        checked += 1
    assert checked >= 6


def test_t2_entry_states_its_decomposed_certification(by_id):
    """The 66-digit T2 value must keep its honest caveat attached."""
    e = by_id["T2.collinear.66"]
    assert e["digits_claimed"] == 66
    digits = e["value"].split(".")[1]
    assert len(digits) == 66, f"T2 prints {len(digits)} digits, claims 66"
    ev = e["evidence"].lower()
    assert "decomposed" in ev, "T2 evidence dropped the decomposed-certification caveat"
    assert "19 digits" in ev, \
        "T2 evidence dropped the fully-independent-pipeline digit count"
    # the superseded corpus anchor and the new value must agree to 18 digits and
    # differ in the 19th -- the correction this table publishes
    with mp.workdps(80):
        old, new = mp.mpf(e["supersedes"]), mp.mpf(e["value"])
        agree = int(mp.floor(-mp.log10(abs(old - new) / abs(new))))
    assert agree == 18, f"anchor correction moved: agreement is now {agree} digits"


def test_pending_entry_claims_zero_digits(by_id):
    """A run still in flight must not be able to lend the table any digits."""
    e = by_id["T2.collinear.120.pending"]
    assert e["digits_claimed"] == 0
    assert e["value"] == "PENDING"


# ---------------------------------------------------------------------------
# 2. regeneration of a fast subset
# ---------------------------------------------------------------------------
def test_fast_subset_reproduces_the_stored_table(by_id):
    """Recompute a cheap subset from live code and match the stored values.

    Spans categories 2, 3, 5 and 7: all four two-centre closed-form classes,
    three exact-rational Slater integrals (including one that goes through the
    exact-Fraction dispatch path), two graph-native helium CI energies (the
    exact rational at n_max = 1 and the certified value at n_max = 2), and
    three anchor constants.
    """
    from benchmarks.certified_reference.generate_table import (FAST_SUBSET_IDS,
                                                               fast_subset)

    fresh = fast_subset(FAST_DPS)
    assert len(fresh) >= 6
    cats = {by_id[k]["category"] for k in fresh}
    assert {"two_center_eri", "slater_rational", "helium_ci", "anchor"} <= cats

    tol = mp.mpf(10) ** (-(FAST_DPS - 2))
    for entry_id in FAST_SUBSET_IDS:
        stored = by_id[entry_id]
        got = fresh[entry_id]
        if stored["value_kind"] == "rational":
            assert got == stored["value"], \
                f"{entry_id}: exact rational changed, {got} vs {stored['value']}"
            continue
        with mp.workdps(60):
            a, b = mp.mpf(got), mp.mpf(stored["value"])
            rel = abs(a - b) / abs(b)
        assert rel < tol, (
            f"{entry_id}: regenerated {got} disagrees with stored "
            f"{stored['value']} at relative {mp.nstr(rel, 3)} "
            f"(reduced precision {FAST_DPS})")


# ---------------------------------------------------------------------------
# 3. graph-native helium CI entries
# ---------------------------------------------------------------------------
#: the production float64 pipeline must reproduce a certified helium value to
#: at least this, or the two routes have genuinely diverged.  Loosening this is
#: a decision, not a maintenance chore.
HELIUM_FLOAT_PIPELINE_TOL = 1e-12

#: physical non-relativistic infinite-mass helium ground state (Pekeris/Drake)
E_HELIUM_PHYSICAL = -2.903724377034119598


def test_helium_entries_are_certified_and_honestly_scoped(by_id):
    """Each graph-native helium row must carry its evidence AND its caveat.

    The caveat is the load-bearing part: these are exact eigenvalues of a
    truncated matrix, not accurate helium energies, and the entry has to say so
    where a reader will see it.
    """
    from fractions import Fraction

    rows = [e for e in by_id.values() if e["category"] == "helium_ci"]
    assert len(rows) >= 4, "the helium block has shrunk unexpectedly"

    certified = 0
    for e in rows:
        assert "certifies the ASSEMBLY" in e["evidence"],             f"{e['id']}: dropped the assembly-not-accuracy caveat"
        assert "2.903724377034119598" in e["evidence"],             f"{e['id']}: dropped the physical-helium comparison"
        assert e["basis_truncation"] >= 1
        assert e["ci_dimension"] >= 1

        if e["digits_claimed"] == "exact":
            # n_max = 1 is one configuration: E = -Z^2 + 5Z/8 = -11/4 at Z = 2
            assert Fraction(e["value"]) == Fraction(-11, 4),                 f"{e['id']}: exact rational moved off -11/4"
            continue

        certified += 1
        assert e["digits_claimed"] <= 35,             f"{e['id']}: claims more than the table's helium cap"
        assert e["two_precision_agree_digits"] >= e["digits_claimed"],             (f"{e['id']}: claims {e['digits_claimed']} digits on only "
             f"{e['two_precision_agree_digits']} digits of two-precision "
             f"agreement")

        with mp.workdps(150):
            claim_scale = mp.mpf(10) ** (-e["digits_claimed"])
            resid = mp.mpf(e["residual_bound"])
            rounding = mp.mpf(e["matrix_rounding_bound"])
        assert resid < claim_scale,             f"{e['id']}: residual bound {e['residual_bound']} exceeds the claim"
        assert rounding < claim_scale, (
            f"{e['id']}: matrix-rounding perturbation "
            f"{e['matrix_rounding_bound']} exceeds the claim")

        rel = e.get("float_pipeline_rel_agreement")
        if rel is not None:
            assert rel < HELIUM_FLOAT_PIPELINE_TOL, (
                f"{e['id']}: the production float pipeline and the exact route "
                f"disagree at relative {rel:.2e}, above the advertised "
                f"{HELIUM_FLOAT_PIPELINE_TOL:.0e} -- localise before relaxing")

    assert certified >= 3, "fewer than three certified helium rows"


def test_helium_series_falls_monotonically_and_stays_above_helium(by_id):
    """Structural sanity on the tabulated block, not a variational theorem.

    Enlarging the basis must lower the CI ground state, and across the
    tabulated range every value sits above the physical helium energy.  A row
    that broke either would indicate a broken assembly rather than a better
    energy -- the graph-native one-body operator is not the true kinetic
    operator, so being below helium would be a red flag, not a triumph.
    """
    rows = sorted((e for e in by_id.values() if e["category"] == "helium_ci"),
                  key=lambda e: e["basis_truncation"])
    vals = [float(mp.mpf(Fraction_or_decimal(e))) for e in rows]
    for a, b in zip(vals, vals[1:]):
        assert b < a, f"CI energy rose with basis size: {a} -> {b}"
    for e, v in zip(rows, vals):
        assert v > E_HELIUM_PHYSICAL, (
            f"{e['id']} = {v} lies BELOW the physical helium energy "
            f"{E_HELIUM_PHYSICAL} -- assembly error, not accuracy")


def Fraction_or_decimal(e):
    """Entry value as something mp.mpf can read (rationals arrive as 'p/q')."""
    from fractions import Fraction
    if e["value_kind"] == "rational":
        f = Fraction(e["value"])
        return mp.mpf(f.numerator) / mp.mpf(f.denominator)
    return mp.mpf(e["value"])
