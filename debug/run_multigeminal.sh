#!/bin/sh
# Launch parallel multi-geminal LM seeds + Pade runs; each writes its own log + exit code.
cd /c/Users/jlout/Desktop/Project_Geometric/debug || exit 1
L=data/mg_logs
mkdir -p $L
for s in 1 2 3 4 5 6 7 8; do
  ( python lih_multigeminal_vmc.py --mode lm --seed $s --nwalk 1000 --nsweep 9000 --nburn 4000 > $L/lm_s$s.log 2>&1; echo "exit=$?" >> $L/lm_s$s.log ) &
done
for s in 11 12 13 14; do
  ( python lih_multigeminal_vmc.py --mode pade --gamma 1.5 --seed $s --nwalk 1000 --nsweep 9000 --nburn 4000 > $L/pade15_s$s.log 2>&1; echo "exit=$?" >> $L/pade15_s$s.log ) &
done
( python lih_multigeminal_vmc.py --mode pade --gamma 1.0 --seed 21 --nwalk 1000 --nsweep 9000 --nburn 4000 > $L/pade10_s21.log 2>&1; echo "exit=$?" >> $L/pade10_s21.log ) &
( python lih_multigeminal_vmc.py --mode pade --gamma 2.0 --seed 22 --nwalk 1000 --nsweep 9000 --nburn 4000 > $L/pade20_s22.log 2>&1; echo "exit=$?" >> $L/pade20_s22.log ) &
wait
echo ALLDONE
