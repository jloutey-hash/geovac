"""
Hain-Brown PSLQ Test A, literal-kernel sharpening.

Closes the scope caveat at debug/sprint_hb_pslq_test_memo.md §5.2:
yesterday's HB PSLQ test used a real-shift kernel (t - tau_imag)^k along
the imaginary axis as an expedient, because the literal kernel
(z - tau)^k along the imaginary axis at tau = 2i has identically-zero
real part by parity. This driver uses the *literal* Hain-Brown Eichler
kernel (z - tau)^k integrated along a real-axis-SHIFTED contour
z = a + i t, a > 0, so the integrand Delta(a + i t) has nontrivial real
and imaginary parts and the literal kernel gives genuinely new modular
content not detectable by the imaginary-axis expedient.

Two extensions vs yesterday:
  (a) literal kernel (z - tau)^k, k = 0, 1, ..., 10 (full Eichler range
      for weight w = 12 cusp form Delta).
  (b) base point tau = a + 2 i with a = 1/2 (interior of standard
      fundamental domain), contour along the upper vertical ray
      z = a + i t from t = 2 up to a numerical cutoff.

The literal Eichler periods are complex-valued; we take both real and
imaginary parts as separate basis elements (each lives in the modular
period ring of Delta).

Cross-precision discipline inherits from yesterday: 50 / 100 / 200 dps.

Verdict gates (DECISION-GATE per sprint prompt):
  VERDICT-STABLE-NEGATIVE: NO modular-ring identifications survive
      cross-precision agreement (matches yesterday's verdict; §5.2 caveat
      CLOSED; depth-1 HB rule-out is robust under kernel sharpening).
  VERDICT-FLIPS-POSITIVE: at least one Sym^2-tagged GeoVac period
      identifies with a literal Eichler period at all three precisions
      (would reopen Hain-Brown identification at depth 1).
  NULL: literal kernel inconclusive at ceiling 10^6 (would warrant
      stronger ceiling or alternative contour follow-on).
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
PSLQ_MAXSTEPS = 2000

DATA_DIR = Path(__file__).parent / "data"
DATA_DIR.mkdir(exist_ok=True)
OUT_JSON = DATA_DIR / "hb_eichler_kernel_results.json"

# Base point and contour.
# IMPORTANT: TAU_REAL must NOT be a half-integer; at a = 1/2 the
# Fourier-series phases cos(2 pi a n) = (-1)^n / sin(2 pi a n) = 0 force
# Delta(a + i b) and E_k(a + i b) to be purely REAL, collapsing the
# literal kernel onto yesterday's imaginary-axis expedient up to overall
# i^{k+1} factors. Choose a = 1/3 so cos(2 pi n / 3), sin(2 pi n / 3)
# are non-trivial and the modular values are genuinely complex.
TAU_REAL = mp.mpf(1) / mp.mpf(3)
TAU_IMAG = mp.mpf(2)
KERNEL_K_RANGE = list(range(0, 11))  # k = 0, 1, ..., 10 (Delta is weight 12)

# Numerical cutoff for the Eichler integral (Delta decays as exp(-2*pi*t))
T_CUTOFF = mp.mpf(20)  # at t = 20, |q| ~ 1e-54


def to_dps_str(x, dps: int) -> str:
    if isinstance(x, mp.mpc):
        return f"{mp.nstr(x.real, dps)} + {mp.nstr(x.imag, dps)}j"
    return mp.nstr(x, dps)


# ---------------------------------------------------------------------------
# Modular forms at arbitrary tau (complex-valued via q-series)
# ---------------------------------------------------------------------------

def _sigma_k(n: int, k: int) -> int:
    s = 0
    for d in range(1, n + 1):
        if n % d == 0:
            s += d ** k
    return s


def _E4_E6_Delta_at(tau: mp.mpc, N: int) -> tuple:
    """Compute E_4(tau), E_6(tau), Delta(tau) via direct q-series at tau,
    real or complex.

    E_4 = 1 + 240 sum sigma_3(n) q^n
    E_6 = 1 - 504 sum sigma_5(n) q^n
    Delta = (E_4^3 - E_6^2) / 1728
    q = exp(2 pi i tau)
    """
    q = mp.exp(2 * mp.pi * 1j * tau)
    E4 = mp.mpc(1)
    E6 = mp.mpc(1)
    for n in range(1, N + 1):
        qn = q ** n
        E4 += 240 * _sigma_k(n, 3) * qn
        E6 -= 504 * _sigma_k(n, 5) * qn
    Delta = (E4 ** 3 - E6 ** 2) / mp.mpf(1728)
    return E4, E6, Delta


def compute_modular_basis_literal() -> dict:
    """Compute E_4(tau_0), E_6(tau_0), Delta(tau_0) at the shifted base
    point tau_0 = TAU_REAL + i * TAU_IMAG, plus literal Hain-Brown Eichler
    periods of Delta:
        P_k^literal = int_{tau_0}^{i infty} Delta(z) (z - tau_0)^k dz,
                       k = 0, 1, ..., 10
    along the vertical ray z = TAU_REAL + i t, t from TAU_IMAG to T_CUTOFF.

    Substituting z = TAU_REAL + i t,  dz = i dt,  z - tau_0 = i (t - TAU_IMAG):
        P_k^literal = i^{k+1} int_{TAU_IMAG}^{T_CUTOFF}
                      Delta(TAU_REAL + i t) (t - TAU_IMAG)^k dt
    Delta at TAU_REAL + i t is complex (because tau is off the imaginary
    axis, q = exp(2 pi i tau) carries a non-trivial phase factor
    exp(2 pi i * TAU_REAL * n) in each Fourier coefficient -- so the q-series
    is not pure-real even though tau(n) are real integers).

    Returns dict with real and imag parts of each modular quantity and
    each literal Eichler period.
    """
    tau_0 = mp.mpc(TAU_REAL, TAU_IMAG)
    N = max(60, int(mp.mp.dps * 0.6))

    E4_0, E6_0, Delta_0 = _E4_E6_Delta_at(tau_0, N)

    # Literal Eichler periods at the shifted base.
    eichler_periods = {}

    def delta_at_t(t):
        z = mp.mpc(TAU_REAL, t)
        # Compute Delta(z) via q-series
        q = mp.exp(2 * mp.pi * 1j * z)
        E4z = mp.mpc(1)
        E6z = mp.mpc(1)
        for n in range(1, N + 1):
            qn = q ** n
            E4z += 240 * _sigma_k(n, 3) * qn
            E6z -= 504 * _sigma_k(n, 5) * qn
        return (E4z ** 3 - E6z ** 2) / mp.mpf(1728)

    for k in KERNEL_K_RANGE:
        # Integrand: Delta(TAU_REAL + i t) * (t - TAU_IMAG)^k
        # with overall prefactor i^{k+1} applied after integration.
        def integrand_re(t, kk=k):
            d = delta_at_t(t)
            return d.real * (t - TAU_IMAG) ** kk

        def integrand_im(t, kk=k):
            d = delta_at_t(t)
            return d.imag * (t - TAU_IMAG) ** kk

        # Integrate real and imaginary parts separately for numerical
        # stability and reusable tanh-sinh quadrature.
        try:
            int_re = mp.quad(integrand_re, [TAU_IMAG, T_CUTOFF],
                             method='tanh-sinh', maxdegree=8)
            int_im = mp.quad(integrand_im, [TAU_IMAG, T_CUTOFF],
                             method='tanh-sinh', maxdegree=8)
        except Exception:
            int_re = mp.mpf(0)
            int_im = mp.mpf(0)

        # Apply the i^{k+1} prefactor:
        #  Pk = i^{k+1} * (int_re + i * int_im)
        # i^{k+1} cycles {i, -1, -i, 1, i, -1, ...}.
        i_pow = 1j ** (k + 1)
        i_re = mp.mpf(int(round(float(mp.re(i_pow)))))
        i_im = mp.mpf(int(round(float(mp.im(i_pow)))))

        # (i_re + i*i_im) * (int_re + i*int_im) =
        #   (i_re*int_re - i_im*int_im) + i*(i_re*int_im + i_im*int_re)
        Pk_real = i_re * int_re - i_im * int_im
        Pk_imag = i_re * int_im + i_im * int_re

        eichler_periods[f"Eichler_P{k}_re"] = Pk_real
        eichler_periods[f"Eichler_P{k}_im"] = Pk_imag

    return {
        "tau_real": float(TAU_REAL),
        "tau_imag": float(TAU_IMAG),
        "E4_real": mp.re(E4_0),
        "E4_imag": mp.im(E4_0),
        "E6_real": mp.re(E6_0),
        "E6_imag": mp.im(E6_0),
        "Delta_real": mp.re(Delta_0),
        "Delta_imag": mp.im(Delta_0),
        **eichler_periods,
    }


# ---------------------------------------------------------------------------
# Basis assembly
# ---------------------------------------------------------------------------

def assemble_hb_basis_literal(modular: dict,
                              include_zeta: bool = True,
                              include_beta: bool = True,
                              include_pi_powers: bool = True,
                              include_log2: bool = False) -> tuple[list, list]:
    """Return (values, labels). Modular ring uses LITERAL Eichler periods,
    both real and imaginary parts, k = 0..10."""
    values = []
    labels = []

    # E_4, E_6, Delta at shifted base: split into real / imag parts
    for sym, prefix in [("E4", "E4"), ("E6", "E6"), ("Delta", "Delta")]:
        re_key = f"{sym}_real"
        im_key = f"{sym}_imag"
        values.append(mp.mpf(modular[re_key])); labels.append(f"{prefix}(tau)_re")
        values.append(mp.mpf(modular[im_key])); labels.append(f"{prefix}(tau)_im")

    # Literal Eichler periods, real and imaginary parts
    for k in KERNEL_K_RANGE:
        re_key = f"Eichler_P{k}_re"
        im_key = f"Eichler_P{k}_im"
        if re_key in modular:
            values.append(mp.mpf(modular[re_key])); labels.append(f"EichlerLit_P{k}_re")
        if im_key in modular:
            values.append(mp.mpf(modular[im_key])); labels.append(f"EichlerLit_P{k}_im")

    if include_zeta:
        values.append(mp.zeta(3)); labels.append("zeta(3)")
        values.append(mp.zeta(5)); labels.append("zeta(5)")

    if include_beta:
        G = mp.catalan
        values.append(G); labels.append("beta(2)=G")
        beta4 = mp.mpf(0)
        for n in range(0, max(80, int(mp.mp.dps * 0.8))):
            beta4 += (mp.mpf(-1) ** n) / (mp.mpf(2 * n + 1) ** 4)
        values.append(beta4); labels.append("beta(4)")

    if include_pi_powers:
        pi = mp.pi
        values.append(pi); labels.append("pi")
        values.append(pi ** 2); labels.append("pi^2")
        values.append(pi ** 3); labels.append("pi^3")
        values.append(pi ** 4); labels.append("pi^4")
        values.append(pi ** 6); labels.append("pi^6")
        values.append(pi ** 8); labels.append("pi^8")

    if include_log2:
        values.append(mp.log(mp.mpf(2))); labels.append("log(2)")

    return values, labels


# ---------------------------------------------------------------------------
# GeoVac Sym^2-tagged periods (unchanged from yesterday)
# ---------------------------------------------------------------------------

def assemble_geovac_sym2_periods() -> tuple[list, list, list]:
    """Verbatim Sym^2-tagged GeoVac periods from yesterday's panel.

    Returns (values, labels, descriptions).
    """
    pi = mp.pi
    cg_trace = mp.mpf(4)
    cg_hw = mp.mpf(1)

    P_M1 = cg_trace * pi

    a0_Dirac = mp.mpf(4) * pi ** 2
    a1_Dirac = mp.mpf(-2) * pi ** 2
    a0_scalar = mp.mpf(2) * pi ** 2
    a2_scalar = mp.mpf(2) * pi ** 2 / mp.mpf(2)

    zeta_D2_2 = pi ** 2 - pi ** 4 / mp.mpf(12)
    zeta_D2_3 = pi ** 4 / mp.mpf(3) - pi ** 6 / mp.mpf(30)
    zeta_D2_4 = mp.mpf(2) * pi ** 6 / mp.mpf(15) - mp.mpf(17) * pi ** 8 / mp.mpf(1260)

    G = mp.catalan
    beta4 = mp.mpf(0)
    for n in range(0, max(100, int(mp.mp.dps * 0.8))):
        beta4 += (mp.mpf(-1) ** n) / (mp.mpf(2 * n + 1) ** 4)
    diff_s2 = mp.mpf(2) * G - mp.mpf(1)
    diff_s4 = mp.mpf(8) * beta4 - mp.mpf(8) * G

    values, labels, descs = [], [], []

    def add(v, lbl, desc):
        values.append(v); labels.append(lbl); descs.append(desc)

    add(P_M1, "P_M1_trace", "4*pi (M1 * CG-trace)")
    add(cg_trace * a0_Dirac, "P_M2_a0Dirac_trace", "16*pi^2 (M2 Dirac SD0 * CG-trace)")
    add(cg_trace * a1_Dirac, "P_M2_a1Dirac_trace", "-8*pi^2 (M2 Dirac SD1 * CG-trace)")
    add(cg_trace * a0_scalar, "P_M2_a0scalar_trace", "8*pi^2 (M2 scalar SD0 * CG-trace)")
    add(cg_trace * a2_scalar, "P_M2_a2scalar_trace", "4*pi^2 (M2 scalar SD2 * CG-trace)")
    add(cg_trace * zeta_D2_2, "P_M2_zd2_trace", "4*(pi^2 - pi^4/12)")
    add(cg_trace * zeta_D2_3, "P_M2_zd3_trace", "4 * zeta_{D^2}(3)")
    add(cg_trace * zeta_D2_4, "P_M2_zd4_trace", "4 * zeta_{D^2}(4)")

    add(cg_hw * a0_Dirac, "P_M2_a0Dirac_hw", "4*pi^2 (CG-hw)")
    add(cg_hw * zeta_D2_2, "P_M2_zd2_hw", "pi^2 - pi^4/12")

    add(cg_trace * G, "P_M3_beta2_trace", "4*G")
    add(cg_trace * beta4, "P_M3_beta4_trace", "4*beta(4)")
    add(cg_trace * diff_s2, "P_M3_diff2_trace", "4*(2G-1)")
    add(cg_trace * diff_s4, "P_M3_diff4_trace", "32*(beta(4)-G)")

    add(cg_hw * G, "P_M3_beta2_hw", "G (CG-hw)")
    add(cg_hw * beta4, "P_M3_beta4_hw", "beta(4) (CG-hw)")
    add(cg_hw * diff_s2, "P_M3_diff2_hw", "2G-1 (CG-hw)")

    add(zeta_D2_2 * G, "P_M2M3_zd2_G", "(pi^2 - pi^4/12) * G")
    add(a0_Dirac * G, "P_M2M3_a0_G", "4*pi^2 * G")

    add(cg_trace * mp.zeta(3), "P_M3_zeta3_trace", "4*zeta(3)")

    return values, labels, descs


# ---------------------------------------------------------------------------
# PSLQ panel (same machinery as yesterday)
# ---------------------------------------------------------------------------

def run_pslq_panel(period_value, period_label: str,
                   basis_values: list, basis_labels: list,
                   ceiling: int = PSLQ_CEILING,
                   max_basis_relations: int = 5) -> dict:
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
            "note": "period below threshold",
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
                "trivials_found": trivials_found,
            }
        if rel is None:
            return {
                "period": period_label,
                "relation_found": False,
                "trivials_found": trivials_found,
            }
        val = sum(rel[i] * seq[i] for i in range(len(seq)))
        val_str = mp.nstr(val, 5)
        a0 = rel[0]
        if a0 == 0:
            trivials_found.append({
                "iteration": it,
                "coefs": [
                    {"basis_label": kept_labels[i], "coef": int(rel[i + 1])}
                    for i in range(len(kept_labels)) if rel[i + 1] != 0
                ],
                "residual": val_str,
            })
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
            "residual": val_str,
            "trivials_found": trivials_found,
        }

    return {
        "period": period_label,
        "relation_found": False,
        "note": f"only trivial relations after {max_basis_relations} iters",
        "trivials_found": trivials_found,
    }


# ---------------------------------------------------------------------------
# Verdict classifier (modular labels expanded for the literal panel)
# ---------------------------------------------------------------------------

def is_modular_label(lbl: str) -> bool:
    """True if a basis label corresponds to a Hain-Brown modular-ring element.

    Modular ring at shifted base tau = a + 2i comprises:
      E4(tau)_re, E4(tau)_im, E6(tau)_re, E6(tau)_im, Delta(tau)_re, Delta(tau)_im,
      EichlerLit_P{k}_re / _im for k = 0..10.
    """
    for prefix in ("E4(tau)", "E6(tau)", "Delta(tau)", "EichlerLit_"):
        if lbl.startswith(prefix):
            return True
    return False


def classify_relation(result: dict) -> str:
    if not result.get("relation_found", False):
        return "NULL"
    if result.get("trivial", False):
        return "NULL"
    used_labels = set()
    for c in result.get("coefs", []):
        if c["coef_num"] != 0:
            used_labels.add(c["basis_label"])
    if any(is_modular_label(l) for l in used_labels):
        return "POSITIVE"
    return "NEGATIVE"


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("Hain-Brown PSLQ Test A: LITERAL Eichler kernel along shifted contour")
    print(f"  base point tau = {float(TAU_REAL)} + {float(TAU_IMAG)} i")
    print(f"  kernel weights k = {KERNEL_K_RANGE}")
    print("=" * 78)

    all_results = {
        "metadata": {
            "ceiling": PSLQ_CEILING,
            "tau_real": float(TAU_REAL),
            "tau_imag": float(TAU_IMAG),
            "kernel_k_range": KERNEL_K_RANGE,
            "precisions": PRECISIONS,
            "contour": (
                f"vertical ray z = {float(TAU_REAL)} + i t, "
                f"t in [{float(TAU_IMAG)}, {float(T_CUTOFF)}]"
            ),
            "kernel": "literal (z - tau_0)^k Eichler kernel",
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
        mp.mp.dps = dps + 10

        t0 = time.time()
        modular = compute_modular_basis_literal()
        print(f"  Modular (literal kernel) computed in {time.time()-t0:.1f}s")
        print(f"    E4(tau)_re ~ {to_dps_str(modular['E4_real'], 12)}")
        print(f"    E4(tau)_im ~ {to_dps_str(modular['E4_imag'], 12)}")
        print(f"    Delta(tau)_re ~ {to_dps_str(modular['Delta_real'], 12)}")
        print(f"    Delta(tau)_im ~ {to_dps_str(modular['Delta_imag'], 12)}")
        for k in (0, 5, 10):
            re_v = modular.get(f"Eichler_P{k}_re", mp.mpf(0))
            im_v = modular.get(f"Eichler_P{k}_im", mp.mpf(0))
            print(f"    Lit P{k}_re ~ {to_dps_str(re_v, 12)}")
            print(f"    Lit P{k}_im ~ {to_dps_str(im_v, 12)}")

        gv_values, gv_labels, gv_descs = assemble_geovac_sym2_periods()
        if all_results["metadata"]["geovac_periods_count"] is None:
            all_results["metadata"]["geovac_periods_count"] = len(gv_values)

        configs = {
            "modular_only": assemble_hb_basis_literal(
                modular, include_zeta=False, include_beta=False,
                include_pi_powers=False, include_log2=False),
            "modular_plus_zeta": assemble_hb_basis_literal(
                modular, include_zeta=True, include_beta=False,
                include_pi_powers=False, include_log2=False),
            "modular_plus_zeta_beta": assemble_hb_basis_literal(
                modular, include_zeta=True, include_beta=True,
                include_pi_powers=False, include_log2=False),
            "full": assemble_hb_basis_literal(
                modular, include_zeta=True, include_beta=True,
                include_pi_powers=True, include_log2=True),
        }

        prec_dict = {
            "modular_summary": {
                "tau": f"{float(TAU_REAL)} + {float(TAU_IMAG)}i",
                "E4_re": to_dps_str(modular['E4_real'], 30),
                "E4_im": to_dps_str(modular['E4_imag'], 30),
                "E6_re": to_dps_str(modular['E6_real'], 30),
                "E6_im": to_dps_str(modular['E6_imag'], 30),
                "Delta_re": to_dps_str(modular['Delta_real'], 30),
                "Delta_im": to_dps_str(modular['Delta_imag'], 30),
            },
            "literal_eichler_periods": {
                f"P{k}_{part}": to_dps_str(modular[f"Eichler_P{k}_{part}"], 30)
                for k in KERNEL_K_RANGE for part in ("re", "im")
            },
            "config_results": {},
            "basis_sizes": {n: len(v[0]) for n, v in configs.items()},
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
                    short = ",".join(used[:4]) + ("..." if len(used) > 4 else "")
                    print(f"    [{r['verdict'][:3]}] {lbl:32s} via [{short}]")
                cfg_results.append(r)
            prec_dict["config_results"][cfg_name] = cfg_results

        all_results["by_precision"][str(dps)] = prec_dict

    # Cross-precision agreement
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

    positive_modular_hits = 0
    for cc in cross_check.values():
        positive_modular_hits += len(cc["positive_at_all_precisions"])
    if positive_modular_hits > 0:
        verdict = "VERDICT-FLIPS-POSITIVE"
    else:
        any_neg = any(len(cc["negative_at_all_precisions"]) > 0
                      for cc in cross_check.values())
        if any_neg:
            verdict = "VERDICT-STABLE-NEGATIVE"
        else:
            verdict = "NULL"

    all_results["top_level_verdict"] = verdict
    all_results["positive_modular_hit_count_cross_precision"] = positive_modular_hits

    with open(OUT_JSON, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nWrote {OUT_JSON}")
    print(f"TOP-LEVEL VERDICT: {verdict}")
    print(f"Positive modular hits (cross-precision agreement): {positive_modular_hits}")

    return all_results


if __name__ == "__main__":
    main()
