"""Category 6 -- certified minimal-basis diatomic total energies (Paper 58).

The quadrature-free diatomic (QFD) build assembles a complete diatomic
electronic-structure calculation with NO quadrature anywhere on the production
path: every one-electron matrix element (overlap, kinetic, both
nuclear-attraction kernels) and every two-electron class (one-centre,
``(AA|BB)``, hybrid, exchange) is an exact symbolic expression, so the whole
energy can be evaluated at any requested precision.  ``geovac/qfd_core.py``
holds the closed forms; ``geovac/qfd_assemble.py`` does the Loewdin
orthogonalisation and the full CI in arbitrary-precision arithmetic.

WHAT THESE ROWS ARE, AND ARE NOT.  They are certified values of a completely
specified minimal-basis model: s-type hydrogenic orbitals with decay rate
``a = Z_orbital / n``, a stated internuclear separation, and full CI in that
basis.  Anyone who implements the same model must reproduce these digits.  They
are NOT accurate chemistry -- a two- or three-function s-only basis with
unoptimised exponents is far from the exact energy, and no accuracy claim is
made.  Their use is as reference data for validating an integral code end to
end: an error anywhere in the one- or two-electron assembly moves the energy.

THE STRUCTURAL AXIS THAT SETS THE DIGIT COUNTS.  The exchange class is built
from a Neumann (Legendre) expansion in an index ``tau``.  That expansion
TERMINATES exactly when the two centres of a charge density carry the same
orbital exponent -- the condition ``q = (alpha - beta) R / 2 = 0``.  So:

  * ``H2+`` and ``He2^2+`` are homonuclear at equal exponent: the tau sum is
    FINITE, every higher term is a symbolic zero, and the model energy has no
    truncation error at all.  Its digit count is limited only by the working
    precision of the arithmetic.
  * ``HeH+`` and ``BeH+`` are heteronuclear: the tau sum is infinite and is
    truncated at a per-quartet ``tau_max``.  These rows therefore claim digits
    NET OF AN EXPLICIT TAIL BOUND, propagated to the energy through the
    two-particle density matrix and the Loewdin transform with an amplification
    factor recomputed for that system (never inherited from another one).

Rows for the homonuclear systems are RECOMPUTED by this generator, because they
are cheap.  The two heteronuclear rows are QUOTED from the campaign that
produced them (``debug/qfd_table_ext.py``, findings in
``debug/qfd_table_findings.md``): their exchange assembly runs to ``tau = 22``
and takes far too long to sit inside a table generator.  Nothing is re-derived
here for those two, so nothing here can silently drift from the campaign.

The two systems the QFD build was demonstrated on are included below as
quoted rows: H2 (homonuclear, terminating; claim capped at 60 digits -- the
weaker run's requested precision, per this artifact's convention) and LiH
(heteronuclear, 30 digits net of the tau-tail bound).  Campaign record:
``debug/data/qfd_h2_certified.json`` / ``qfd_lih_certified.json``.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Any, Dict, List, Tuple

import mpmath as mp
import sympy as sp

from ._common import entry

CATEGORY = "qfd_diatomic"

# --------------------------------------------------------------------------
# Quoted heteronuclear rows (see the module docstring for why).
# Transcribed from debug/data/qfd_table_ext_certified.json, produced by
# debug/qfd_table_ext.py; accounting in debug/qfd_table_findings.md.
# --------------------------------------------------------------------------
_METHOD_HETERO = (
    "Quadrature-free assembly (geovac.qfd_assemble.total_energy): the "
    "one-electron matrix elements come from the Mulliken auxiliary integrals "
    "A_m(p), B_n(q) in prolate spheroidal coordinates; the two-electron tensor "
    "is built class by class -- one-centre Slater R^0, the (AA|BB) Coulomb "
    "class and the hybrid class in closed form, and the exchange class from the "
    "closed-form ordered-xi / eta factorisation term by term in the Neumann "
    "index tau.  The AO tensors are then Loewdin-orthogonalised exactly and a "
    "full CI is diagonalised, all in arbitrary-precision arithmetic.  Because "
    "the two centres carry different orbital exponents the tau sum does not "
    "terminate; each term is still a closed form, and the truncation is what "
    "the digit claim is stated net of.")

QUOTED: Dict[str, Dict[str, Any]] = {
    "qfd.h2.R1.4": {
        "label": "H2 total energy, minimal 1s/1s (zeta = 1) hydrogenic basis, "
                 "R = 1.4 bohr -- the truncation-free quadrature-free assembly",
        "value": "-1.10655660609135850801940122475993772281832890689690719549979",
        "digits": 60,
        "electronic": "-1.82084232037707279373368693904565200853261462118262148121408",
        "vnn": "0.714285714285714285714285714285714285714285714285714285714286",
        "n_electrons": 2,
        "basis": "s-only hydrogenic, a = Z_orbital / n: 1s_A a=1, 1s_B a=1",
        "transcendence_class": "{exp, E_1, ln, gamma}",
        "tau_note": "terminates symbolically (tau = 1 and all tau > 2 are exact "
                    "zeros) -- no truncation of any kind",
        "backing_test": "tests/test_paper58_qfd.py::test_h2_certified_energy_and_tau_termination (live re-assembly to ~40 digits; the 60-digit claim is campaign-level)",
        "sort_index": 0,
        "method": _METHOD_HETERO,
        "evidence":
            "Whole pipeline evaluated at working precisions 60 and 90: they "
            "agree to 6e-85 (guard digits carry the internal computations "
            "further), but the claim is CAPPED AT 60, the weaker run's "
            "requested precision, per this artifact's convention.  No "
            "truncation exists anywhere: the homonuclear exchange tau-series "
            "terminates symbolically (pinned bit-identical at tau_max 2 vs 5).  "
            "INDEPENDENT ROUTES: the exchange integral matches Sugiura's 1927 "
            "closed form to 3.1e-61 (independent {E_1, ln, gamma} content; "
            "in-suite at 1e-38, tests/test_paper58_qfd.py::"
            "test_exchange_matches_sugiura_1927); all integrals match direct "
            "quadrature at 1e-21..1e-23; the two closed-form routes to h "
            "(radial Laplacian vs hydrogenic eigen-trick) differ by exactly "
            "zero symbolically.  Backing test: tests/test_paper58_qfd.py "
            "(certified-energy leg at ~40 digits in-suite; campaign record "
            "debug/data/qfd_h2_certified.json).",
        "honest_scope": "Minimal unoptimised basis (zeta = 1): the certified "
                        "value sits 0.068 Ha above the exact H2 energy "
                        "(-1.174476 Ha).  A certified value of a fully "
                        "specified model, not an accurate H2 energy.",
    },
    "qfd.lih.R3.015": {
        "label": "LiH total energy, minimal s-only Li 1s,2s (Z_orbital = 3) / "
                 "H 1s hydrogenic basis, R = 3.015 bohr",
        "value": "-7.87261979241561217317558085835",
        "digits": 30,
        "electronic": "-8.86764466803750272043926245039",
        "vnn": "0.995024875621890547263681592040",
        "n_electrons": 4,
        "basis": "s-only hydrogenic: Li 1s a=3, Li 2s a=3/2, H 1s a=1",
        "transcendence_class": "{exp, E_1, ln, gamma}",
        "tau_note": "infinite (q = 3.015 and 0.754); truncated at tau_max = "
                    "20/16/14 per quartet, tail bounded and netted out",
        "backing_test": "tests/test_paper58_qfd.py::test_heteronuclear_exchange_tau_series + ::test_exchange_hp_matches_symbolic (one exchange quartet and the numeric accumulator; the assembled 30-digit total is campaign-level)",
        "sort_index": 1,
        "method": _METHOD_HETERO,
        "evidence":
            "Working precisions 60 vs 90 agree far beyond the claim; the "
            "BINDING constraint is the tau tail: per-quartet monotone ratio "
            "bounds sum to 1.95e-32 in the integrals, amplified through the "
            "2-RDM weight (6) times ||S^-1/2||^4 = 2.830 (lambda_min(S) = "
            "0.594448), rounded up to 64, giving an energy tail bound of "
            "1.24e-30 Ha -- hence 30 digits claimed (linear-algebra and "
            "accumulator limits sit at 55 and 44 digits).  INDEPENDENT ROUTES: "
            "one-electron integrals vs quadrature at ~1e-23; one-centre ERIs "
            "vs geovac.hypergeometric_slater at <=4.4e-16 (its float ceiling); "
            "two-electron classes vs quadrature 1e-21..1e-23 at matched "
            "truncation; exchange_hp accumulator vs the symbolic route pinned "
            "in-suite (tests/test_paper58_qfd.py::"
            "test_exchange_hp_matches_symbolic).  Campaign record "
            "debug/data/qfd_lih_certified.json + debug/qfd_track1_findings.md.",
        "honest_scope": "Minimal s-only three-function basis: far above the "
                        "exact LiH energy (-8.0705 Ha).  A certified value of "
                        "a fully specified model, not an accurate LiH energy.",
    },
    "qfd.hehplus.R1.46": {
        "label": "HeH+ total energy, minimal He 1s (Z_orbital = 2) / H 1s "
                 "(Z_orbital = 1) hydrogenic basis, R = 1.46 bohr",
        "value": "-2.895950902302325174310463610145323586925",
        "digits": 40,
        "electronic": "-4.265813916000955311296764980008337285555",
        "vnn": "1.369863013698630136986301369863013698630",
        "n_electrons": 2,
        "basis": "s-only hydrogenic, a = Z_orbital / n: He 1s a=2, H 1s a=1",
        "transcendence_class": "{exp, E_1, ln, gamma}",
        "tau_note": "infinite (q = 0.73); truncated at tau_max = 16, tail "
                    "bounded and netted out of the claim",
        "sort_index": 3,
        "method": _METHOD_HETERO,
        "evidence":
            "Whole pipeline re-run at working precisions 40 and 60: they agree "
            "to 2.53e-65 relative, so the printed digits are the digits of the "
            "assembled expression; the claim is capped at 40, the weaker run.  "
            "TRUNCATION, netted out rather than ignored: the single exchange "
            "quartet is summed to tau_max = 16, its per-tau magnitudes are "
            "measured to decrease monotonically through the tail (last ratio "
            "6.01e-4), and the geometric bound |a_taumax| r/(1-r) gives a tail "
            "below 1.19e-42 in the integral.  Propagating that to the energy "
            "with an amplification recomputed for THIS system -- two-particle "
            "density-matrix weight N(N-1)/2 = 1 times ||S^-1/2||^4 = "
            "lambda_min(S)^-2 = 3.764, lambda_min(S) = 0.5154, rounded up to 16 "
            "for margin -- bounds the energy error at 1.90e-41 Ha, i.e. 41 "
            "digits, so the tail is not the binding constraint.  The exchange "
            "accumulator's own arithmetic was checked by re-running it at "
            "working precisions 30 and 50 at fixed tau: 1.23e-44 (43 digits).  "
            "INDEPENDENT ROUTES, all sharing no code with the closed forms: raw "
            "prolate-spheroidal quadrature of the one-electron elements agrees "
            "to 1.7e-22 (overlap), 3.9e-22 (kinetic) and 5.8e-22 / 6.5e-22 (the "
            "two nuclear-attraction kernels); Newton-potential quadrature of the "
            "one-centre, (AA|BB) and hybrid two-electron classes agrees to "
            "1.7e-21, 4.7e-22 and 6.3e-22; and a fully numeric Neumann "
            "evaluation of the exchange class -- eta half and ordered-xi half "
            "both by quadrature -- agrees to 3.4e-21 when compared at matched "
            "truncation (tau <= 3), which is the accuracy of that quadrature "
            "route.  The two independent closed-form routes to h (explicit "
            "radial Laplacian vs the hydrogenic eigen-trick) differ by exactly "
            "zero, symbolically.  Backing test: tests/test_paper58_qfd.py.",
        "honest_scope":
            "Minimal unoptimised basis: a certified value of a fully specified "
            "model, not an accurate HeH+ energy.  The exact non-relativistic "
            "HeH+ ground state is about -2.978 Ha near R = 1.46 bohr; two "
            "s functions with hydrogenic exponents recover -2.896 Ha.  The "
            "variational ordering is the expected one (the model sits above "
            "the exact energy), and that is all the comparison establishes.",
        "extra": {"exchange_tau_max": 16,
                  "tail_bound_integral": "1.19e-42",
                  "energy_tail_bound_Ha": "1.90e-41",
                  "amplification_factor": 16,
                  "lambda_min_S": "0.515422776536"},
    },
    "qfd.behplus.R2.5": {
        "label": "BeH+ total energy, minimal Be 1s,2s (Z_orbital = 4) / H 1s "
                 "(Z_orbital = 1) hydrogenic basis, R = 2.5 bohr, 4 electrons",
        "value": "-14.69882770595311489573042272632",
        "digits": 31,
        "electronic": "-16.29882770595311489573042272632",
        "vnn": "1.600000000000000000000000000000",
        "n_electrons": 4,
        "basis": "s-only hydrogenic, a = Z_orbital / n: Be 1s a=4, Be 2s a=2, "
                 "H 1s a=1",
        "transcendence_class": "{exp, E_1, ln, gamma}",
        "tau_note": "infinite (q = 3.75 and 1.25); truncated per quartet at "
                    "tau_max = 22 / 18 / 16, tails bounded and netted out",
        "sort_index": 5,
        "method": _METHOD_HETERO + "  Four electrons over three spatial "
                  "orbitals: the CI is over all C(6,4) = 15 determinants.  The "
                  "Be 2s function is a genuine hydrogenic 2s, radial node "
                  "included, so the basis is not a set of simple Slater 1s "
                  "functions.",
        "evidence":
            "Whole pipeline re-run at working precisions 40 and 60: they agree "
            "to 1.76e-66 relative, so the linear algebra alone would support 40 "
            "digits.  It is not the binding constraint.  TRUNCATION, which is: "
            "the three exchange quartets are summed to tau_max = 22, 18 and 16 "
            "(the mismatch parameter q = (alpha-beta)R/2 is 3.75 for the "
            "Be1s/H1s density and 1.25 for Be2s/H1s, and larger q means slower "
            "Neumann convergence).  Per-tau magnitudes decrease monotonically "
            "through every tail; the geometric bounds |a_taumax| r/(1-r) with "
            "the last observed ratios 7.95e-3, 4.36e-3 and 2.06e-3 give "
            "2.64e-35, 6.46e-33 and 3.79e-33, summing to 1.03e-32 at the "
            "integral level.  Propagating with an amplification recomputed for "
            "THIS system -- two-particle density-matrix weight N(N-1)/2 = 6 "
            "times ||S^-1/2||^4 = lambda_min(S)^-2 = 2.975, lambda_min(S) = "
            "0.5797, rounded up to 128 for margin -- bounds the energy error at "
            "1.31e-30 Ha, i.e. 31 digits, which is the claim.  The exchange "
            "accumulator's own arithmetic was checked at working precisions 30 "
            "and 50 at fixed tau: 1.49e-43 (40 digits).  INDEPENDENT ROUTES, "
            "sharing no code with the closed forms: raw prolate-spheroidal "
            "quadrature of the one-electron elements agrees to 5.6e-22 "
            "(overlap), 6.1e-23 (kinetic) and 1.5e-22 (nuclear attraction); "
            "Newton-potential quadrature of all 18 one-centre, (AA|BB) and "
            "hybrid two-electron integrals agrees to between 0 and 7.0e-21; and "
            "a fully numeric Neumann evaluation of each of the three exchange "
            "quartets, compared at matched truncation (tau <= 3), agrees to "
            "2.7e-20, 4.6e-20 and 2.1e-20 -- the accuracy of that quadrature "
            "route.  The two independent closed-form routes to h (explicit "
            "radial Laplacian vs the hydrogenic eigen-trick) differ by exactly "
            "zero, symbolically.  Backing test: tests/test_paper58_qfd.py.",
        "honest_scope":
            "Minimal unoptimised basis: a certified value of a fully specified "
            "model, not an accurate BeH+ energy, and three s functions cannot "
            "represent the Be valence region's p character at all.  A "
            "self-contained check of how far off it is: the true BeH+ ground "
            "state must lie below its own Be+ + H dissociation asymptote, "
            "-14.325 - 0.5 = -14.825 Ha, whereas this model gives -14.699 Ha at "
            "R = 2.5 bohr -- above the asymptote, so the model does not bind.  "
            "The certified digits are digits of the model.",
        "extra": {"exchange_tau_max": [22, 18, 16],
                  "tail_bound_integral": "1.03e-32",
                  "energy_tail_bound_Ha": "1.31e-30",
                  "amplification_factor": 128,
                  "lambda_min_S": "0.579723581346",
                  "fci_dimension": 15},
    },
}


def _recompute(orbs, ZA, ZB, R, n_elec, tau_max,
               dps_lo: int = 40, dps_hi: int = 60) -> Tuple[str, str, str, int]:
    """Rebuild one homonuclear system from the live closed forms.

    Returns (E_total, E_electronic, V_NN, agreeing_digits) with the three
    numbers printed to the number of digits the two working precisions agree
    on, capped at ``dps_lo`` (a claim may never exceed the weaker run).
    """
    from geovac import qfd_assemble as AS

    S, h = AS.build_S_h(orbs, ZA, ZB, R)
    gmap, _canon = AS.build_g(orbs, R, tau_max=tau_max)
    out = {}
    for dps in (dps_lo, dps_hi):
        with mp.workdps(dps + 40):
            out[dps] = tuple(mp.mpf(x) for x in
                             AS.total_energy(S, h, gmap, orbs, n_elec,
                                             ZA, ZB, R, dps))
    with mp.workdps(dps_hi + 60):
        lo, hi = out[dps_lo][0], out[dps_hi][0]
        agree = (10 ** 6 if lo == hi
                 else int(mp.floor(-mp.log10(abs(lo - hi) / abs(hi)))))
        agree = min(agree, dps_lo)
        printed = tuple(mp.nstr(v, agree, strip_zeros=False)
                        for v in out[dps_hi])
    return printed[0], printed[1], printed[2], agree


def _tau_terminates(ZA_orb, na, ZB_orb, nb, R, tau_verify: int) -> Tuple[int, bool]:
    """Confirm the exchange Neumann sum is a FINITE sum for this pair.

    Returns (last nonzero tau, all higher terms are symbolic zeros).
    """
    from geovac import qfd_core as Q

    terms = Q.exchange_closed_form(ZA_orb, (na, 0, 0), (nb, 0, 0),
                                   ZB_orb, (na, 0, 0), (nb, 0, 0), R,
                                   tau_max=tau_verify, return_terms=True)
    simp = [sp.simplify(t) for t in terms]
    nz = [i for i, t in enumerate(simp) if t != 0]
    last = max(nz) if nz else -1
    return last, all(simp[i] == 0 for i in range(last + 1, len(simp)))


# --------------------------------------------------------------------------
def build(mode: str = "full") -> List[Dict[str, Any]]:
    del mode          # every row here is either cheap to rebuild or quoted
    rows: List[Dict[str, Any]] = []

    # ---- H2+, one electron, homonuclear, terminating ---------------------
    h2p_orbs = [("A", Fraction(1), 1), ("B", Fraction(1), 1)]
    for rlabel, R, sidx in (("2.0", sp.Integer(2), 0),
                            ("1.4", sp.Rational(7, 5), 1)):
        tot, ele, vnn, agree = _recompute(h2p_orbs, 1, 1, R, 1, tau_max=2)
        last_nz, terminates = _tau_terminates(Fraction(1), 1, Fraction(1), 1,
                                              R, 8)
        rows.append(entry(
            "qfd.h2plus.R" + rlabel,
            CATEGORY,
            "H2+ total energy, minimal 1s/1s hydrogenic basis (zeta = 1), "
            "R = " + rlabel + " bohr",
            tot, agree,
            "Quadrature-free assembly (geovac.qfd_assemble): the overlap, "
            "kinetic and both nuclear-attraction matrix elements come from the "
            "Mulliken auxiliary integrals A_m(p), B_n(q) in prolate spheroidal "
            "coordinates, in exact symbolic form; the ground state is then the "
            "one-electron full CI over the Loewdin-orthogonalised pair, "
            "computed in mpmath.  With one electron this is the 2x2 "
            "generalised secular problem (h, S), and the value is reproduced "
            "independently by the lowest generalised eigenvalue and by the "
            "closed LCAO sigma_g expression (h_AA + h_AB) / (1 + S_AB).  The "
            "transcendence class is elementary for exactly that reason: the "
            "two-electron tensor is assembled and available, and its exchange "
            "member does carry {exp, E_1, ln, gamma}, but with one electron no "
            "two-electron integral can enter the energy.",
            "Whole pipeline re-run at working precisions 40 and 60: first " +
            str(agree) + " significant digits identical (claim capped at the "
            "weaker run).  There is NO truncation error to net out: the "
            "exchange Neumann tau sum terminates at tau = " + str(last_nz) +
            " and every term through tau = 8 above it is a symbolic zero (" +
            str(terminates) + ") -- and with a single electron no two-electron "
            "integral can contribute at all.  Independent-route checks at "
            "R = 2.0 bohr: raw prolate-spheroidal quadrature of the overlap, "
            "kinetic and nuclear-attraction elements agrees to 5.8e-21, "
            "1.3e-23 and 1.7e-22 absolute; the classical literature closed "
            "forms for the 1s two-centre integrals (S_AB, T_AB, V^B_AA, "
            "V^A_AB) agree to better than 1e-40.  The two independent "
            "closed-form routes to h (explicit radial Laplacian vs the "
            "hydrogenic eigen-trick) differ by exactly zero, symbolically.  "
            "Backing test: tests/test_paper58_qfd.py.",
            transcendence_class="elementary {exp}",
            defining_relation="E_total = min spec(S^-1/2 h S^-1/2) + Z_A Z_B / R",
            electronic_energy=ele,
            nuclear_repulsion=vnn,
            n_electrons=1,
            basis="s-only hydrogenic, a = Z_orbital / n: 1s_A a=1, 1s_B a=1",
            exchange_tau_series="terminates (q = 0); no truncation error",
            geovac_entry_point="geovac.qfd_assemble.total_energy",
            backing_test="tests/test_paper58_qfd.py",
            source_memo="debug/qfd_table_findings.md",
            sort_index=sidx,
            honest_scope="Minimal unoptimised basis: this is a certified value "
                         "of a fully specified model, not an accurate H2+ "
                         "energy.  Context: the zeta = 1 LCAO H2+ curve "
                         "minimises at R = 2.49 bohr with E = -0.5648 Ha, while "
                         "the exact Born-Oppenheimer curve minimises near "
                         "R = 2.00 bohr at about -0.6026 Ha.",
        ))

    # ---- He2^2+, two electrons, homonuclear, terminating -----------------
    he2_orbs = [("A", Fraction(2), 1), ("B", Fraction(2), 1)]
    R13 = sp.Rational(13, 10)
    tot, ele, vnn, agree = _recompute(he2_orbs, 2, 2, R13, 2, tau_max=2)
    last_nz, terminates = _tau_terminates(Fraction(2), 1, Fraction(2), 1, R13, 8)
    rows.append(entry(
        "qfd.he2_2plus.R1.3",
        CATEGORY,
        "He2^2+ total energy, minimal 1s/1s hydrogenic basis (Z_orbital = 2), "
        "R = 1.3 bohr",
        tot, agree,
        "Quadrature-free assembly (geovac.qfd_assemble): closed-form one- and "
        "two-electron integrals -- one-centre Slater R^0, the (AA|BB) Coulomb "
        "class, the hybrid class and the exchange class -- followed by exact "
        "Loewdin orthogonalisation and a two-electron full CI, all in mpmath.  "
        "The exchange class is assembled fully symbolically here: at equal "
        "orbital exponents its Neumann tau sum is finite, so no numerical "
        "accumulation is needed.",
        "Whole pipeline re-run at working precisions 40 and 60: first " +
        str(agree) + " significant digits identical (claim capped at the weaker "
        "run).  There is NO truncation error to net out: the exchange tau sum "
        "terminates at tau = " + str(last_nz) + " and every term through "
        "tau = 8 above it is a symbolic zero (" + str(terminates) + ") -- the "
        "same q = 0 criterion that makes H2 exact.  Independent-route checks: "
        "raw prolate-spheroidal quadrature of the one-electron elements and a "
        "fully numeric Neumann evaluation of the exchange integral (which "
        "shares no code with the closed-form ordered-xi route) both agree at "
        "the 1e-20 level of those float/quadrature routes; the two independent "
        "closed-form routes to h differ by exactly zero, symbolically.  "
        "Backing test: tests/test_paper58_qfd.py.",
        transcendence_class="{exp, E_1, ln, gamma}",
        defining_relation="E_total = E_FCI(S^-1/2 h S^-1/2, S^-1/2 g S^-1/2) "
                          "+ Z_A Z_B / R",
        electronic_energy=ele,
        nuclear_repulsion=vnn,
        n_electrons=2,
        basis="s-only hydrogenic, a = Z_orbital / n: 1s_A a=2, 1s_B a=2",
        exchange_tau_series="terminates (q = 0); no truncation error",
        geovac_entry_point="geovac.qfd_assemble.total_energy",
        backing_test="tests/test_paper58_qfd.py",
        source_memo="debug/qfd_table_findings.md",
        sort_index=2,
        honest_scope="Minimal unoptimised basis; a certified value of a fully "
                     "specified model, not an accurate He2^2+ energy.  In this "
                     "basis the He+ + He+ dissociation limit is exactly -4 Ha "
                     "(two 1s functions at Z_orbital = 2, each -Z^2/2), and the "
                     "value at R = 1.3 bohr lies ABOVE it -- the two-function "
                     "model does not bind the real barrier-protected He2^2+ "
                     "minimum.  That is a statement about the basis, not about "
                     "the certified number.",
    ))

    # ---- quoted heteronuclear rows ---------------------------------------
    for eid, q in QUOTED.items():
        rows.append(entry(
            eid, CATEGORY, q["label"], q["value"], q["digits"],
            q["method"], q["evidence"],
            transcendence_class=q["transcendence_class"],
            defining_relation="E_total = E_FCI(S^-1/2 h S^-1/2, "
                              "S^-1/2 g S^-1/2) + Z_A Z_B / R",
            electronic_energy=q["electronic"],
            nuclear_repulsion=q["vnn"],
            n_electrons=q["n_electrons"],
            basis=q["basis"],
            exchange_tau_series=q["tau_note"],
            geovac_entry_point="geovac.qfd_assemble.total_energy",
            backing_test=q.get(
                "backing_test",
                "tests/test_paper58_qfd.py (assembly machinery only; this system "
                "is not itself exercised in-suite -- the value is quoted from the "
                "campaign, see evidence)"),
            source_memo="debug/qfd_table_findings.md",
            sort_index=q["sort_index"],
            honest_scope=q["honest_scope"],
            **q.get("extra", {}),
        ))


    # ---- the closed-form H2 PES: certified equilibrium constants ---------
    # E(R) is ONE symbolic expression (atoms {exp, E_1, log, EulerGamma});
    # constants are roots of its exact derivatives by closed-form Newton
    # (geovac.qfd_assemble.h2_pes_certify; dps 45 vs 60 Newton agree 9.5e-46;
    # dissociation to exactly -1 Ha verified at 3e-49).
    pes = {
        "qfd.h2.pes.Req": ("H2 closed-form PES: equilibrium separation R_eq (bohr)",
                           "1.66799996697274872492704641030726451334980414"),
        "qfd.h2.pes.De": ("H2 closed-form PES: D_e = -1 - E(R_eq) (Ha)",
                          "0.118650362098295311853497764333884671691909683"),
        "qfd.h2.pes.k": ("H2 closed-form PES: force constant E''(R_eq) (Ha/bohr^2)",
                         "0.254703934307969981152727935724355672678251789"),
    }
    for i, (eid, (label, val)) in enumerate(pes.items()):
        rows.append(entry(
            eid, CATEGORY, label, val, 45,
            method="root of the exact derivative of the single closed-form "
                   "expression E(R) (2x2 singlet CI over closed-form MO "
                   "integrals; minimal 1s/1s basis, zeta = 1), by Newton "
                   "iteration on symbolic dE/dR at two precisions",
            evidence="Newton at dps 45 vs 60 agree to 9.5e-46; E(R) equals the "
                     "84-digit certified FCI value at R = 1.4 (1e-42) and "
                     "dissociates to exactly -1 Ha (3e-49), so D_e is itself "
                     "closed-form; pinned in tests/test_paper58_qfd.py::"
                     "test_closed_form_pes_equilibrium_constants",
            transcendence_class="root of an {exp, E_1, log, EulerGamma} expression",
            geovac_entry_point="geovac.qfd_assemble.h2_pes_certify",
            backing_test="tests/test_paper58_qfd.py",
            source_memo="debug/data/h2_closed_pes.json",
            sort_index=90 + i,
            honest_scope="Minimal-basis (zeta = 1) model constants: exact H2 has "
                         "R_eq = 1.401 bohr; the certified numbers are exact "
                         "properties of the closed-form model curve, not of the "
                         "physical molecule.",
        ))
    return rows
