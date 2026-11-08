# 120-digit T2 run "u1" (launched 2026-08-21, v4.104.0 follow-on 3)

Target: h<=1e4 PSLQ verdict on the pre-registered dim-20 ring needs ~120 digits.
Config u1: dps=145, K=580, panel width 4, 96 nodes/panel, Ns=330+1.05k, Nw=80+0.80k, s-map p=6, 12 chunks.

- Launch (idempotent-ish: chunks that already have beta2_chunk_u1_<idx>.json are DONE;
  re-launching overwrites, so to RESUME after an interruption, edit the launcher's loop
  or just relaunch missing chunk indices manually with beta2_t2_kw_chunk.py using the
  printed chunk edges):
    python debug/beta2_launch_hi.py 145 580 4 96 330 1.05 80 0.80 6 u1 12
- Monitor:  ls debug/data/beta2_chunk_u1_*.json ; tail debug/data/beta2_chunk_u1_0.out
- Assemble when all 12 json files exist:
    python debug/beta2_assemble.py u1 145 580 220
- Cross-validation config u2 (run AFTER u1; different params, required before claiming digits):
    python debug/beta2_launch_hi.py 150 560 5 112 360 1.15 90 0.85 4 u2 12
    python debug/beta2_assemble.py u2 150 560 240
- Then re-run the guarded PSLQ (debug/beta2_t2_pslq.py) at the new precision.
