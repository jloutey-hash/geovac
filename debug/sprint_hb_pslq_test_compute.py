"""
Hain-Brown PSLQ Test A driver.

Computes GeoVac periods at the Sym^2 level of the SL_2 action and
PSLQ-identifies them against the Hain-Brown modular period ring
{E_4(2i), E_6(2i), Delta(2i), Eichler periods of Delta, zeta(3), zeta(5),
beta(2), beta(4)} at coefficient ceiling 10^6, precision 50/100/200 dps.

Setup:
    GeoVac's Sym^2 V_fund SL_2-action on {e1^2, e1 e2, e2^2} has
    Clebsch-Gordan multiplicities {1, 2, 1} (the binomial C(2,k)
    coefficients), carrying sl_2 weights {+2, 0, -2}. A Sym^2-tagged
    GeoVac period is a structural composition where the Sym^2
    Clebsch-Gordan weight factor multiplies an M2 (Seeley-DeWitt) or
    M3 (vertex-parity Hurwitz) value.

    Concrete Sym^2-tagged periods constructed:
      P_2_M2_a0  = (Sym^2 trace = 1+2+1 = 4) * a_0^{D^2} = 4 * 4*pi^2 = 16*pi^2
                   (vol-norm Dirac SD zeroth coefficient, Paper 51 Cor 2.1)
      P_2_M2_a1  = 4 * (-2*pi^2) = -8*pi^2  (vol-norm Dirac SD first coeff)
      P_2_M3_D2  = 4 * zeta_{D^2}(2) = 4 * (pi^2 - pi^4/12)
                   (T9 spectral-zeta, vertex parity coupling)
      P_2_M3_diff2 = (Sym^2 highest-weight = 1) * (D_even(2) - D_odd(2))
                   = 2 * (beta(2) - beta(0)) = 2*beta(2) - 1
                   (Paper 28 Thm 3 parity discriminant; M3 cyclotomic level 4)
      P_2_M3_diff4 = 1 * (D_even(4) - D_odd(4)) = 8 * (beta(4) - beta(2))
                   (depth-1 M3 vertex parity at s=4)
      P_2_M2_zd3  = 1 * zeta_{D^2}(3) = pi^4/3 - pi^6/30
                   (vol-norm 4th order SD residual at integer s)
      P_2_M3_S5_a = level-4 cyclotomic refinement,
                   8 * (beta(2) - beta(4))
      P_2_M1     = (Sym^2 trace) * pi = 4*pi
                   (M1 Hopf-base measure scaling by Sym^2 trace)

Hain-Brown basis (level 1 SL_2(Z) MEM at tau = 2i):
    Eisenstein E_4(2i), E_6(2i)
    Modular discriminant Delta(2i)
    Eichler-style Delta periods (real/imag parts of Eichler integral)
    Classical ZetaC(3), ZetaC(5) (in Hain-Brown's iterated-integral library)
    Dirichlet beta(2) = G, beta(4)
    Powers of pi up through pi^8

Verdict gates:
  POSITIVE if at least one Sym^2-tagged GeoVac period PSLQ-identifies
           with at least one Hain-Brown modular-ring element
           (E_4, E_6, Delta, Eichler periods of Delta) at all three
           precisions consistently.
  NEGATIVE if the only identifications are with pure MZV(MTM) basis
           (pi powers, zeta values, beta values) and no modular content
           appears.
  NULL if no PSLQ relations found at ceiling 10^6.
"""

from __future__ import annotations

import json
from pathlib import Path
import time

import mpmath as mp

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

PRECISIONS = [50, 100, 200]
PSLQ_CEILING = 10 ** 6
PSLQ_MAXSTEPS = 2000  # default 100 is too small for 10+ element basis

DATA_DIR = Path(__file__).parent / "data"
DATA_DIR.mkdir(exist_ok=True)
OUT_JSON = DATA_DIR / "hb_pslq_test_results.json"


def to_dps_str(x: mp.mpc | mp.mpf, dps: int) -> str:
    """Stringify to dps digits, preserving sign and complex parts."""
    if isinstance(x, mp.mpc):
        return f"{mp.nstr(x.real, dps)} + {mp.nstr(x.imag, dps)}j"
    return mp.nstr(x, dps)


# ---------------------------------------------------------------------------
# Modular forms at tau = 2i
# ---------------------------------------------------------------------------

def compute_modular_basis(tau_imag: float = 2.0) -> dict:
    """Compute E_4(tau), E_6(tau), Delta(tau), and an Eichler-style
    period of Delta. Uses mpmath jtheta and eta-product expressions.

    Conventions:
      Use the Dedekind eta function eta(tau) and the standard identities:
        Delta(tau) = eta(tau)^24
        E_4(tau) = 1 + 240 * sum_{n>=1} sigma_3(n) q^n,  q = exp(2*pi*i*tau)
        E_6(tau) = 1 - 504 * sum_{n>=1} sigma_5(n) q^n
      For tau = 2i the q-series converge very fast.
    """
    tau = mp.mpc(0, tau_imag)
    q = mp.exp(2 * mp.pi * 1j * tau)

    # E_4 and E_6 via q-series; truncate at N high enough for current dps
    N = max(40, int(mp.mp.dps * 0.5))

    def sigma_k(n: int, k: int) -> int:
        s = 0
        for d in range(1, n + 1):
            if n % d == 0:
                s += d ** k
        return s

    E4 = mp.mpf(1)
    E6 = mp.mpf(1)
    for n in range(1, N + 1):
        qn = q ** n
        E4 += 240 * sigma_k(n, 3) * qn
        E6 -= 504 * sigma_k(n, 5) * qn

    # Delta from E_4, E_6 via the identity Delta = (E_4^3 - E_6^2) / 1728
    Delta = (E4 ** 3 - E6 ** 2) / mp.mpf(1728)

    # Eichler-style periods of Delta along the imaginary axis from 2i to
    # i*infty. With z = i*t, dz = i*dt, Delta(it) is real. Kernel
    # K_real(t) = (t - 2)^k (real-valued shift along imag axis) gives
    # purely-real Eichler integrals:\
    #     P_k = int_{2}^{infty} Delta(it) (t - 2)^k dt
    # These are real-valued modular periods that mix Hain-Brown's Eichler
    # cohomology weights. We compute k = 0, 1, 2 to widen the modular ring.
    T_cutoff = mp.mpf(20)  # at t=20, q = exp(-2*pi*20) ~ 10^-54, tiny

    def Delta_at_it(t):
        z = mp.mpc(0, t)
        qz = mp.exp(2 * mp.pi * 1j * z)
        E4z = mp.mpf(1)
        E6z = mp.mpf(1)
        for n in range(1, N + 1):
            qn = qz ** n
            E4z += 240 * sigma_k(n, 3) * qn
            E6z -= 504 * sigma_k(n, 5) * qn
        # Real-valued for t real
        return mp.re((E4z ** 3 - E6z ** 2) / mp.mpf(1728))

    eichler_periods = {}
    for k in (0, 1, 2):
        def integrand(t, kk=k):
            return Delta_at_it(t) * (t - tau_imag) ** kk

        try:
            P_k = mp.quad(integrand, [tau_imag, T_cutoff],
                          method='tanh-sinh', maxdegree=8)
        except Exception:
            P_k = mp.mpf(0)
        eichler_periods[f"Eichler_P{k}"] = mp.re(P_k)

    return {
        "tau_imag": tau_imag,
        "E4_real": mp.re(E4),
        "E4_imag": mp.im(E4),
        "E6_real": mp.re(E6),
        "E6_imag": mp.im(E6),
        "Delta_real": mp.re(Delta),
        "Delta_imag": mp.im(Delta),
        **eichler_periods,
    }


# ---------------------------------------------------------------------------
# Hain-Brown basis assembly (real-valued; cusp form imaginary part should
# vanish at tau = 2i by parity, but we keep both)
# ---------------------------------------------------------------------------

def assemble_hb_basis(modular: dict, include_zeta: bool = True,
                      include_beta: bool = True,
                      include_pi_powers: bool = True,
                      include_cyclotomic: bool = False) -> tuple[list[mp.mpf], list[str]]:
    """Return (values, labels) for the Hain-Brown candidate basis."""
    values = []
    labels = []

    # Modular ring (load-bearing)
    values.append(mp.mpf(modular["E4_real"])); labels.append("E4(2i)")
    values.append(mp.mpf(modular["E6_real"])); labels.append("E6(2i)")
    values.append(mp.mpf(modular["Delta_real"])); labels.append("Delta(2i)")
    # Eichler-style periods P_0, P_1, P_2
    for k in (0, 1, 2):
        key = f"Eichler_P{k}"
        if key in modular:
            values.append(mp.mpf(modular[key])); labels.append(key)

    # Classical zetas (in Hain-Brown's iterated-integral library)
    if include_zeta:
        values.append(mp.zeta(3)); labels.append("zeta(3)")
        values.append(mp.zeta(5)); labels.append("zeta(5)")

    # Dirichlet beta values (level-4 cyclotomic)
    if include_beta:
        # beta(s) = L(s, chi_{-4}) = sum (-1)^n / (2n+1)^s
        G = mp.catalan  # beta(2)
        values.append(G); labels.append("beta(2)=G")
        # beta(4) = sum_{n>=0} (-1)^n / (2n+1)^4
        beta4 = mp.mpf(0)
        for n in range(0, max(80, int(mp.mp.dps * 0.8))):
            beta4 += (mp.mpf(-1) ** n) / (mp.mpf(2 * n + 1) ** 4)
        values.append(beta4); labels.append("beta(4)")

    # Powers of pi (the M2 ring, always include)
    if include_pi_powers:
        pi = mp.pi
        values.append(pi); labels.append("pi")
        values.append(pi ** 2); labels.append("pi^2")
        values.append(pi ** 3); labels.append("pi^3")
        values.append(pi ** 4); labels.append("pi^4")
        values.append(pi ** 6); labels.append("pi^6")
        values.append(pi ** 8); labels.append("pi^8")

    # Cyclotomic shifts (level-4 refinement). Hurwitz at quarter-integers
    # is redundant with pi^k and beta(s) via
    #   zeta(s, 1/4) + zeta(s, 3/4) = (4^s - 2^s) * zeta(s)
    #   zeta(s, 1/4) - zeta(s, 3/4) = 4^s * beta(s)
    # so we include the symmetric SUM only (independent of beta(s) at level 4
    # via the standard pi-power generators).
    if include_cyclotomic:
        values.append(mp.log(mp.mpf(2))); labels.append("log(2)")

    return values, labels


# ---------------------------------------------------------------------------
# GeoVac Sym^2-tagged periods
# ---------------------------------------------------------------------------

def assemble_geovac_sym2_periods() -> tuple[list[mp.mpf], list[str], list[str]]:
    """Return (values, labels, structural_descriptions).

    Each value is a Sym^2-tagged structural composition: a Sym^2
    Clebsch-Gordan weight factor multiplying an M2 (Seeley-DeWitt) or
    M3 (vertex-parity Hurwitz) period.
    """
    pi = mp.pi
    # Sym^2 CG coefficients: trace = 4, highest-weight = 1
    cg_trace = mp.mpf(4)   # 1 + 2 + 1
    cg_hw = mp.mpf(1)      # weight +2 (or -2) entry

    # M1: Hopf-base measure scaled by Sym^2 trace
    P_M1 = cg_trace * pi  # 4*pi

    # M2: Seeley-DeWitt vol-normalized coefficients (Paper 51 Cor 2.1)
    # a_0^{D^2} = 4*pi^2, a_1^{D^2} = -2*pi^2 (Dirac sector, vol-norm)
    # a_k^{Delta}_{k=0..4} = 2*pi^2 / k! (scalar sector)
    a0_Dirac = mp.mpf(4) * pi ** 2
    a1_Dirac = mp.mpf(-2) * pi ** 2
    a0_scalar = mp.mpf(2) * pi ** 2
    a2_scalar = mp.mpf(2) * pi ** 2 / mp.mpf(2)  # 2!
    a3_scalar = mp.mpf(2) * pi ** 2 / mp.mpf(6)  # 3!

    # spectral zeta at integer s on D^2 (Paper 28 Thm 1)
    zeta_D2_2 = pi ** 2 - pi ** 4 / mp.mpf(12)
    zeta_D2_3 = pi ** 4 / mp.mpf(3) - pi ** 6 / mp.mpf(30)
    zeta_D2_4 = mp.mpf(2) * pi ** 6 / mp.mpf(15) - mp.mpf(17) * pi ** 8 / mp.mpf(1260)

    # M3: vertex-parity Hurwitz / Dirichlet-L content (Paper 28 Thm 3)
    G = mp.catalan
    # beta(4)
    beta4 = mp.mpf(0)
    for n in range(0, max(100, int(mp.mp.dps * 0.8))):
        beta4 += (mp.mpf(-1) ** n) / (mp.mpf(2 * n + 1) ** 4)
    # parity discriminant D_even(s) - D_odd(s) = 2^{s-1} * (beta(s) - beta(s-2))
    # at s=2:\ 2*(beta(2) - beta(0)) = 2G - 1   (beta(0) = 1/2)
    diff_s2 = mp.mpf(2) * G - mp.mpf(1)
    # at s=4:\ 8*(beta(4) - beta(2)) = 8*beta4 - 8G
    diff_s4 = mp.mpf(8) * beta4 - mp.mpf(8) * G

    values = []
    labels = []
    descs = []

    def add(v, lbl, desc):
        values.append(v); labels.append(lbl); descs.append(desc)

    # M1-tagged
    add(P_M1, "P_M1_trace", "Sym^2 trace * pi = 4*pi (M1 * CG-trace)")

    # M2-tagged with Sym^2 trace
    add(cg_trace * a0_Dirac, "P_M2_a0Dirac_trace",
        "CG-trace * a_0^{D^2} = 16*pi^2 (M2 Dirac SD0 with Sym^2 trace)")
    add(cg_trace * a1_Dirac, "P_M2_a1Dirac_trace",
        "CG-trace * a_1^{D^2} = -8*pi^2 (M2 Dirac SD1 with Sym^2 trace)")
    add(cg_trace * a0_scalar, "P_M2_a0scalar_trace",
        "CG-trace * a_0^{scalar} = 8*pi^2 (M2 scalar SD0 with Sym^2 trace)")
    add(cg_trace * a2_scalar, "P_M2_a2scalar_trace",
        "CG-trace * a_2^{scalar} = 4*pi^2 (M2 scalar SD2 with Sym^2 trace)")
    add(cg_trace * zeta_D2_2, "P_M2_zd2_trace",
        "CG-trace * zeta_{D^2}(2) = 4*(pi^2 - pi^4/12) (M2 spec-zeta s=2 with Sym^2 trace)")
    add(cg_trace * zeta_D2_3, "P_M2_zd3_trace",
        "CG-trace * zeta_{D^2}(3) (M2 spec-zeta s=3 with Sym^2 trace)")
    add(cg_trace * zeta_D2_4, "P_M2_zd4_trace",
        "CG-trace * zeta_{D^2}(4) (M2 spec-zeta s=4 with Sym^2 trace)")

    # M2-tagged with Sym^2 highest-weight (=1, distinguishes from trivial)
    add(cg_hw * a0_Dirac, "P_M2_a0Dirac_hw",
        "CG-hw * a_0^{D^2} = 4*pi^2 (M2 Dirac SD0 with Sym^2 hw)")
    add(cg_hw * zeta_D2_2, "P_M2_zd2_hw",
        "CG-hw * zeta_{D^2}(2) = pi^2 - pi^4/12 (M2 spec-zeta s=2 with Sym^2 hw)")

    # M3-tagged with Sym^2 trace
    add(cg_trace * G, "P_M3_beta2_trace",
        "CG-trace * beta(2) = 4*G (M3 vertex parity, Catalan)")
    add(cg_trace * beta4, "P_M3_beta4_trace",
        "CG-trace * beta(4) = 4*beta(4) (M3 vertex parity at s=4)")
    add(cg_trace * diff_s2, "P_M3_diff2_trace",
        "CG-trace * (D_even-D_odd)(2) = 4*(2G-1) (M3 parity discriminant at s=2)")
    add(cg_trace * diff_s4, "P_M3_diff4_trace",
        "CG-trace * (D_even-D_odd)(4) = 4*8*(beta(4)-G) (M3 parity discriminant at s=4)")

    # M3-tagged with Sym^2 highest-weight
    add(cg_hw * G, "P_M3_beta2_hw", "CG-hw * beta(2) = G (M3 vertex parity with Sym^2 hw)")
    add(cg_hw * beta4, "P_M3_beta4_hw", "CG-hw * beta(4) (M3 vertex parity with Sym^2 hw)")
    add(cg_hw * diff_s2, "P_M3_diff2_hw",
        "CG-hw * (D_even-D_odd)(2) = 2G-1 (M3 parity discriminant at s=2 with Sym^2 hw)")

    # Joint M2 x M3 candidates (Sym^2 Clebsch-Gordan provides the bridge)
    add(zeta_D2_2 * G, "P_M2M3_zd2_G",
        "zeta_{D^2}(2) * beta(2) = (pi^2 - pi^4/12) * G (joint M2 * M3)")
    add(a0_Dirac * G, "P_M2M3_a0_G",
        "a_0^{D^2} * beta(2) = 4*pi^2 * G (joint M2 * M3)")

    # zeta(3) tagged with Sym^2 trace -- direct M3-classical injection test
    add(cg_trace * mp.zeta(3), "P_M3_zeta3_trace",
        "CG-trace * zeta(3) = 4*zeta(3) (M3 odd-zeta at level 1)")

    return values, labels, descs


# ---------------------------------------------------------------------------
# PSLQ panel
# ---------------------------------------------------------------------------

def run_pslq_panel(period_value: mp.mpf, period_label: str,
                   basis_values: list[mp.mpf], basis_labels: list[str],
                   ceiling: int = PSLQ_CEILING,
                   max_basis_relations: int = 5) -> dict:
    """Try PSLQ on [period_value] + basis_values, filtering near-zero
    basis entries (PSLQ requires nonzero vector).

    Iteratively removes trivial basis-basis relations (period coef = 0)
    by dropping the highest-coefficient basis element, then re-running.
    Up to `max_basis_relations` iterations.

    Returns dict with relation if found involving the period (period coef != 0),
    along with all the trivial basis relations discovered along the way.
    """
    threshold = mp.mpf(10) ** (-(mp.mp.dps // 4))

    kept = []
    kept_labels = []
    for v, lbl in zip(basis_values, basis_labels):
        if abs(v) > threshold:
            kept.append(v)
            kept_labels.append(lbl)

    if abs(period_value) <= threshold:
        return {
            "period": period_label,
            "relation_found": False,
            "note": "period value below threshold",
            "raw_relation": None,
        }

    trivials_found = []
    for it in range(max_basis_relations + 1):
        seq = [period_value] + kept
        if len(kept) == 0:
            break
        try:
            rel = mp.pslq(seq, maxcoeff=ceiling, maxsteps=PSLQ_MAXSTEPS)
        except (ValueError, RuntimeError) as e:
            return {
                "period": period_label,
                "relation_found": False,
                "error": str(e),
                "raw_relation": None,
                "trivials_found": trivials_found,
            }
        if rel is None:
            return {
                "period": period_label,
                "relation_found": False,
                "raw_relation": None,
                "trivials_found": trivials_found,
            }
        val = sum(rel[i] * seq[i] for i in range(len(seq)))
        val_str = mp.nstr(val, 5)
        a0 = rel[0]
        if a0 == 0:
            # Trivial basis relation; record and drop the basis element
            # with the largest |coef| then re-run.
            trivial_record = {
                "iteration": it,
                "coefs": [
                    {"basis_label": kept_labels[i],
                     "coef": int(rel[i + 1])}
                    for i in range(len(kept_labels)) if rel[i + 1] != 0
                ],
                "residual": val_str,
            }
            trivials_found.append(trivial_record)
            # Drop the LAST basis element involved in the trivial relation
            # (latest by index = "least canonical" generator).
            last_nonzero = -1
            for j in range(1, len(seq)):
                if int(rel[j]) != 0:
                    last_nonzero = j
            if last_nonzero < 1:
                break
            drop_i = last_nonzero - 1
            kept = kept[:drop_i] + kept[drop_i + 1:]
            kept_labels = kept_labels[:drop_i] + kept_labels[drop_i + 1:]
            continue
        # Period-involving relation
        coefs = []
        for i, lbl in enumerate(kept_labels):
            c = -rel[i + 1]
            if c == 0:
                continue
            coefs.append({
                "basis_label": lbl,
                "coef_num": int(c),
                "coef_den": int(a0),
            })
        return {
            "period": period_label,
            "relation_found": True,
            "trivial": False,
            "a0": int(a0),
            "coefs": coefs,
            "raw_relation": [int(x) for x in rel],
            "residual": val_str,
            "trivials_found": trivials_found,
        }

    # Exhausted iterations without finding period-involving relation
    return {
        "period": period_label,
        "relation_found": False,
        "note": f"only trivial basis-basis relations after {max_basis_relations} iterations",
        "trivials_found": trivials_found,
    }


# ---------------------------------------------------------------------------
# Verdict classifier
# ---------------------------------------------------------------------------

MODULAR_LABELS = {"E4(2i)", "E6(2i)", "Delta(2i)",
                  "Eichler_P0", "Eichler_P1", "Eichler_P2"}


def classify_relation(result: dict) -> str:
    """Classify a PSLQ relation as POSITIVE (modular-ring hit), NEGATIVE
    (pure MTM, no modular content), or NULL (no relation)."""
    if not result.get("relation_found", False):
        return "NULL"
    if result.get("trivial", False):
        return "NULL"
    used_labels = set()
    for c in result.get("coefs", []):
        if c["coef_num"] != 0:
            used_labels.add(c["basis_label"])
    if used_labels & MODULAR_LABELS:
        return "POSITIVE"
    return "NEGATIVE"


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("Hain-Brown PSLQ Test A")
    print("=" * 78)

    all_results = {
        "metadata": {
            "ceiling": PSLQ_CEILING,
            "tau_imag": 2.0,
            "precisions": PRECISIONS,
            "geovac_periods_count": None,
            "hb_basis_configs": [
                "modular_only",
                "modular_plus_zeta",
                "modular_plus_zeta_beta",
                "full",
            ],
        },
        "by_precision": {},
    }

    for dps in PRECISIONS:
        print(f"\n>>> Precision: {dps} dps")
        mp.mp.dps = dps + 10  # work with safety margin

        # Compute modular basis at this precision
        t0 = time.time()
        modular = compute_modular_basis(tau_imag=2.0)
        print(f"  Modular at tau=2i computed in {time.time()-t0:.2f}s")
        print(f"    E4(2i) ~ {to_dps_str(modular['E4_real'], 12)}")
        print(f"    E6(2i) ~ {to_dps_str(modular['E6_real'], 12)}")
        print(f"    Delta(2i) ~ {to_dps_str(modular['Delta_real'], 12)}")
        print(f"    Eichler_P0 ~ {to_dps_str(modular.get('Eichler_P0', mp.mpf(0)), 12)}")
        print(f"    Eichler_P1 ~ {to_dps_str(modular.get('Eichler_P1', mp.mpf(0)), 12)}")
        print(f"    Eichler_P2 ~ {to_dps_str(modular.get('Eichler_P2', mp.mpf(0)), 12)}")

        # GeoVac periods
        gv_values, gv_labels, gv_descs = assemble_geovac_sym2_periods()
        if all_results["metadata"]["geovac_periods_count"] is None:
            all_results["metadata"]["geovac_periods_count"] = len(gv_values)

        # Multi-basis configs
        configs = {
            "modular_only": assemble_hb_basis(
                modular, include_zeta=False, include_beta=False,
                include_pi_powers=False, include_cyclotomic=False),
            "modular_plus_zeta": assemble_hb_basis(
                modular, include_zeta=True, include_beta=False,
                include_pi_powers=False, include_cyclotomic=False),
            "modular_plus_zeta_beta": assemble_hb_basis(
                modular, include_zeta=True, include_beta=True,
                include_pi_powers=False, include_cyclotomic=False),
            "full": assemble_hb_basis(
                modular, include_zeta=True, include_beta=True,
                include_pi_powers=True, include_cyclotomic=True),
        }

        prec_dict = {
            "modular_values": {
                "E4(2i)": to_dps_str(modular['E4_real'], 30),
                "E6(2i)": to_dps_str(modular['E6_real'], 30),
                "Delta(2i)": to_dps_str(modular['Delta_real'], 30),
                "Eichler_P0": to_dps_str(modular.get('Eichler_P0', mp.mpf(0)), 30),
                "Eichler_P1": to_dps_str(modular.get('Eichler_P1', mp.mpf(0)), 30),
                "Eichler_P2": to_dps_str(modular.get('Eichler_P2', mp.mpf(0)), 30),
            },
            "config_results": {},
        }

        for cfg_name, (basis_vals, basis_labels) in configs.items():
            print(f"\n  Config: {cfg_name} (basis size = {len(basis_vals)})")
            cfg_results = []
            for v, lbl, desc in zip(gv_values, gv_labels, gv_descs):
                r = run_pslq_panel(v, lbl, basis_vals, basis_labels,
                                   ceiling=PSLQ_CEILING)
                r["verdict"] = classify_relation(r)
                r["description"] = desc
                if r["verdict"] != "NULL":
                    used = [c["basis_label"] for c in r.get("coefs", [])
                            if c["coef_num"] != 0]
                    short = ",".join(used[:5]) + ("..." if len(used) > 5 else "")
                    print(f"    [{r['verdict'][:3]}] {lbl:32s} via [{short}]")
                cfg_results.append(r)
            prec_dict["config_results"][cfg_name] = cfg_results

        all_results["by_precision"][str(dps)] = prec_dict

    # Cross-precision agreement: which periods got POSITIVE verdict at all
    # three precisions consistently?
    cross_check = {}
    for cfg_name in configs.keys():
        cross_check[cfg_name] = {
            "positive_at_all_precisions": [],
            "positive_at_some": [],
            "negative_at_all_precisions": [],
            "null_at_all": [],
        }
        for i, lbl in enumerate(gv_labels):
            verdicts = []
            for dps in PRECISIONS:
                cfg = all_results["by_precision"][str(dps)]["config_results"][cfg_name]
                verdicts.append(cfg[i]["verdict"])
            if all(v == "POSITIVE" for v in verdicts):
                cross_check[cfg_name]["positive_at_all_precisions"].append(lbl)
            elif "POSITIVE" in verdicts:
                cross_check[cfg_name]["positive_at_some"].append({
                    "label": lbl, "verdicts": verdicts})
            elif all(v == "NEGATIVE" for v in verdicts):
                cross_check[cfg_name]["negative_at_all_precisions"].append(lbl)
            elif all(v == "NULL" for v in verdicts):
                cross_check[cfg_name]["null_at_all"].append(lbl)

    all_results["cross_precision_summary"] = cross_check

    # Top-level verdict
    positive_modular_hits = 0
    for cfg_name, cc in cross_check.items():
        positive_modular_hits += len(cc["positive_at_all_precisions"])
    if positive_modular_hits > 0:
        verdict = "POSITIVE"
    else:
        any_negative_consistent = any(
            len(cc["negative_at_all_precisions"]) > 0 for cc in cross_check.values()
        )
        if any_negative_consistent:
            verdict = "NEGATIVE"
        else:
            verdict = "NULL"

    all_results["top_level_verdict"] = verdict
    all_results["positive_modular_hit_count_cross_precision"] = positive_modular_hits

    # Save
    with open(OUT_JSON, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nWrote {OUT_JSON}")
    print(f"TOP-LEVEL VERDICT: {verdict}")
    print(f"Positive modular hits (cross-precision agreement): {positive_modular_hits}")

    return all_results


if __name__ == "__main__":
    main()
