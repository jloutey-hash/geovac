#!/bin/sh
# GA3 heavy chain (sequential; each stage uses Pool(15)): exp-shape mod-p control at gamma = 2 validated against the
# CF5 float series; lin mod-p series N = 520 over 12 primes; lin float series (M 1600, 600 digits, 800 nodes, rho 2);
# the Y_lin w-grid at 32 digits.  Logs: debug/data/ga3_*.log
cd "$(dirname "$0")"
D=data
python ga3_5d_modp_series.py exp 400 P0 60 $D/cf_g2_series_g2_M1600_d600_rho2.json > $D/ga3_5d_phi_modp_exp_N400_P0.log 2>&1
for i in P0 0 1 2 3 4 5 6 7 8 9 10; do
  python ga3_5d_modp_series.py lin 520 $i > $D/ga3_5d_phi_modp_lin_N520_p$i.log 2>&1
done
python ga3_4_series.py 2 1600 600 800 2 lin > $D/ga3_series_main.log 2>&1
python ga3_3_grid.py 2 32 0.0625 > $D/ga3_3_grid.log 2>&1
echo HEAVY_DONE >> $D/ga3_3_grid.log
