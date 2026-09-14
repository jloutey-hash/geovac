"""group2 Batch-3 remediation: the FCI-atoms LARGE -- Table I, the
convergence-detail table, the abstract and the summary sentence all carry
PRE-ERI-FIX (2026-08-29 wrong-sign-q) energies.  Re-measured 2026-09-13 (this
session); N_SD unchanged (basis size), energies lower (restored correlation,
all still above the variational bound), NNZ higher (~3x, the restored
m-changing multipoles -- same census that moved 130->214 cross-block ERIs),
timings left as-is (hardware-dependent, not a physics claim).

He uses hybrid-h1 + full analytical Slater via direct CI (N_SD=5995 at n5);
Li and Be use exact-h1.  Be n4 shift (0.027 Ha) independently cross-checked:
production grid-quad R_k vs exact hypergeometric R_k agree to 1.5e-4, and the
exact route reproduces F^0(1s,1s)=5Z/8=2.500000 for Be exactly -- the shift is
physical, not an evaluator artifact.

Values (Z, n_max): E_old -> E_new / err_old -> err_new / NNZ_old -> NNZ_new
  He 2: -2.887550->-2.889303 0.56->0.50  161->203
  He 3: -2.890601->-2.892932 0.45->0.37  3394->6066
  He 4: -2.892642->-2.895227 0.38->0.29  35422->87550
  He 5: -2.893582->-2.896298 0.35->0.26  244587->777019
  Li 2: -7.101226->-7.104490 5.04->5.00  352->604
  Li 3: -7.392086->-7.394597 1.15->1.12  46200->111192
  Li 4: -7.395921->-7.398718 1.10->1.06  1455648->4385956
  Li 5: -7.397751->-7.400705 1.07->1.03  21519736->78025400
  Be 3: -14.531256->-14.558052 0.93->0.75  ---->1224507
  Be 4: -14.535460->-14.562659 0.90->0.71  37486171->118132111

Write-first; LaTeX raw strings; idempotent (marker = new triple / new phrase).
"""
from __future__ import annotations
import sys

P = "papers/group2_quantum_chemistry/paper_fci_atoms.tex"

EDITS = [
    # ---- convergence-detail table: NNZ & E & err triples (unique per row) ----
    ("conv-he2", r"203 & -2.889303 & 0.50",
     r"161 & -2.887550 & 0.56", r"203 & -2.889303 & 0.50"),
    ("conv-he3", r"6\,066 & -2.892932 & 0.37",
     r"3\,394 & -2.890601 & 0.45", r"6\,066 & -2.892932 & 0.37"),
    ("conv-he4", r"87\,550 & -2.895227 & 0.29",
     r"35\,422 & -2.892642 & 0.38", r"87\,550 & -2.895227 & 0.29"),
    ("conv-he5", r"777\,019 & -2.896298 & 0.26",
     r"244\,587 & -2.893582 & 0.35", r"777\,019 & -2.896298 & 0.26"),
    ("conv-li2", r"604 & -7.104490 & 5.00",
     r"352 & -7.101226 & 5.04", r"604 & -7.104490 & 5.00"),
    ("conv-li3", r"111\,192 & -7.394597 & 1.12",
     r"46\,200 & -7.392086 & 1.15", r"111\,192 & -7.394597 & 1.12"),
    ("conv-li4", r"4\,385\,956 & -7.398718 & 1.06",
     r"1\,455\,648 & -7.395921 & 1.10", r"4\,385\,956 & -7.398718 & 1.06"),
    ("conv-li5", r"78\,025\,400 & -7.400705 & 1.03",
     r"21\,519\,736 & -7.397751 & 1.07", r"78\,025\,400 & -7.400705 & 1.03"),
    ("conv-be3", r"1\,224\,507 & -14.558052 & 0.75",
     r"--- & -14.531256 & 0.93", r"1\,224\,507 & -14.558052 & 0.75"),
    ("conv-be4", r"118\,132\,111 & -14.562659 & 0.71",
     r"37\,486\,171 & -14.535460 & 0.90", r"118\,132\,111 & -14.562659 & 0.71"),

    # ---- Table I (tab:comparison): E & err (unique via bare-$ form) ----
    ("tab-he5", r"$-2.8963$ & 0.26 &",
     r"$-2.8936$ & 0.35 &", r"$-2.8963$ & 0.26 &"),
    ("tab-li4", r"$-7.3987$ & 1.06 &",
     r"$-7.3959$ & 1.10 &", r"$-7.3987$ & 1.06 &"),
    ("tab-li5", r"$-7.4007$ & 1.03 &",
     r"$-7.3978$ & 1.07 &", r"$-7.4007$ & 1.03 &"),
    ("tab-be3", r"$-14.558$ & 0.75 &",
     r"$-14.531$ & 0.93 &", r"$-14.558$ & 0.75 &"),
    ("tab-be4", r"$-14.563$ & 0.71 &",
     r"$-14.536$ & 0.90 &", r"$-14.563$ & 0.71 &"),

    # ---- Abstract ----
    ("abs-he", r"$-2.8963\,\text{Ha}$ (0.26\% error versus exact)",
     r"$-2.8936\,\text{Ha}$ (0.35\% error versus exact)",
     r"$-2.8963\,\text{Ha}$ (0.26\% error versus exact)"),
    ("abs-li", r"$-7.4007\,\text{Ha}$ (1.03\% error) at $n_{\max}=5$",
     r"$-7.3978\,\text{Ha}$ (1.07\% error) at $n_{\max}=5$",
     r"$-7.4007\,\text{Ha}$ (1.03\% error) at $n_{\max}=5$"),
    ("abs-be-e", r"$E = -14.563\,\text{Ha}$",
     r"$E = -14.536\,\text{Ha}$", r"$E = -14.563\,\text{Ha}$"),
    ("abs-be-err", r"(0.71\% error) at $n_{\max}=4$ with $487{,}635$",
     r"(0.90\% error) at $n_{\max}=4$ with $487{,}635$",
     r"(0.71\% error) at $n_{\max}=4$ with $487{,}635$"),

    # ---- Summary sentence (body) ----
    ("body-heli", r"achieves 0.26\% error for helium, 1.03\% for lithium",
     r"achieves 0.35\% error for helium, 1.07\% for lithium",
     r"achieves 0.26\% error for helium, 1.03\% for lithium"),
    ("body-be", r"and 0.71\% for beryllium (at $n_{\max}=4$)",
     r"and 0.90\% for beryllium (at $n_{\max}=4$)",
     r"and 0.71\% for beryllium (at $n_{\max}=4$)"),
]


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    applied, skipped, missed = [], [], []
    for name, marker, old, new in EDITS:
        if marker in t:
            skipped.append(name); continue
        if t.count(old) != 1:
            missed.append((name, t.count(old))); continue
        t = t.replace(old, new); applied.append(name)
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    for n in applied: print(f"  ok    {n}")
    for n in skipped: print(f"  skip  {n} (already applied)")
    for n, c in missed: print(f"  MISS  {n}: anchor count={c}")
    print(f"applied {len(applied)}, skipped {len(skipped)}, MISSED {len(missed)}")
    return 3 if missed else 0


if __name__ == "__main__":
    sys.exit(main())
