# MPO bond-rank support drivers (Paper 14, Thm 3.2.A.unified)

Permanent home of the Sprint S2-v2 drivers + panel data backing Paper 14's
MPO bond-rank section: the exact closed form chi^h1 = 2*rank(h_cross) + 1_LL
+ 1_RR (statement A, 29/29 LiH cuts), the (B)/(C) subadditivity inequalities,
the universal interior chi_H profile [4,16,16,9,9,9,6,3,3,2], the composed
chi=2 sub-block boundaries, and the balanced two-value quantization
chi in {9,16} = 2 + 7*N_cross across LiH/BeH2/H2O (12/12 boundaries).

Migrated from `debug/` (the prune-by-design clean-room) on 2026-07-01, 6th
group4 cert follow-on, PI-directed — same durability pattern as
`tests/wh7_support/` and `tests/rank2_rate_support/`: these files are
imported/pinned by `tests/test_paper14_mpo_rank.py`, so they must not be
pruned.

| File | Role |
|---|---|
| `sprint_s2_v2_unified_panel.py` | the Thm 3.2.A.unified LiH panel (importable: `operator_schmidt_rank`, `h1_to_qubit`, `h1_cross_cut_rank`, `build_Vee_total_qubit`, `build_Vee_per_L_qubit`) |
| `sprint_s2_v2_h1_contribution.py` | statement-A h1 closed-form study |
| `sprint_s2_v2_Vee_per_L.py` | per-multipole-L Vee rank study |
| `sprint_s2_v2_balanced_library_panel.py` | balanced-coupled 12-boundary panel (LiH/BeH2/H2O) |
| `data/*.json` | shipped panel outputs (drivers rewrite in place: `data/` relative to `__file__`) |

Consumer test: `tests/test_paper14_mpo_rank.py`. Paper cites: Paper 14 §MPO
bond rank. Sprint chronicle: CHANGELOG v3.54.0; memos remain in `debug/`
(`sprint_s2_v2_memo.md`, `sprint_s2_v2_balanced_library_memo.md`) per the
clean-room policy and are no longer cited by the paper.
