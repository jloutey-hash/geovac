#!/bin/sh
# Record EX Phase 2 (launched after registration commit 8676692c)
set -x
for p in 67108859 33554393 67108837 67108819; do
  for pr in "1 1" "1 3/2" "1 2"; do python debug/exch_run.py run $pr $p --suffix EXS; echo "rc=$?"; done
done
for pr in "1 1" "1 3/2" "1 2"; do python debug/exch_run.py heldout $pr --suffix EXS; echo "rc=$?"; done
for pr in "1 1" "1 3/2" "1 2"; do python debug/exch_run.py numeric $pr --suffix EXS; echo "rc=$?"; done
echo PHASE2_DONE
