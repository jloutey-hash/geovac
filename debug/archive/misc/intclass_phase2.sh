#!/bin/sh
# Record IC Phase 2 (docs/predictions/2026-11-01_integrated_class.md, committed 40093e7c before this ran).
# Controls first; then two waves of one-thread runs; then the scorer.  Log: debug/data/intclass_phase2.log
cd "$(dirname "$0")/.." || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
L=debug/data/intclass_phase2.log
: > "$L"
P1=1073741789
P2=1073741833
for p in $P1 $P2; do
  python debug/intclass_controls.py $p > debug/data/intclass_controls_run_$p.log 2>&1
  echo "controls $p exit $?" >> "$L"
done
python debug/intclass_score.py run W1 $P1 > debug/data/intclass_run_W1_$P1.out 2>&1 &
a=$!
python debug/intclass_score.py run W1 $P2 > debug/data/intclass_run_W1_$P2.out 2>&1 &
b=$!
for c in N1 N2; do for p in $P1 $P2; do
  python debug/intclass_score.py na $c $p > debug/data/intclass_na_${c}_$p.out 2>&1
  echo "na $c $p exit $?" >> "$L"
done; done
wait $a; echo "W1 $P1 exit $?" >> "$L"
wait $b; echo "W1 $P2 exit $?" >> "$L"
python debug/intclass_score.py run W2 $P1 > debug/data/intclass_run_W2_$P1.out 2>&1 &
a=$!
python debug/intclass_score.py run W2 $P2 > debug/data/intclass_run_W2_$P2.out 2>&1 &
b=$!
wait $a; echo "W2 $P1 exit $?" >> "$L"
wait $b; echo "W2 $P2 exit $?" >> "$L"
python debug/intclass_score.py score > debug/data/intclass_score.out 2>&1
echo "score exit $?" >> "$L"
