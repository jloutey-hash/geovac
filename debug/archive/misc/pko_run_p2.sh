#!/bin/bash
# PK stage 1: the kite at P2 = (1, 3/2, 7/3, 8/3), primes by index, one after another, Pool(6) each.
cd "$(dirname "$0")"
read -a PS < data/ib_primes.txt
for i in "$@"; do
  p=${PS[$i]}
  s=$(date +%s)
  IB_PRIME=$p python ib_kite_de.py pko_p2_p$i 1 3/2 7/3 8/3 5 1 80 > data/pko_p2_p${i}_run.log 2>&1
  echo "pko_p2_p$i exit=$? wall=$(( $(date +%s) - s ))s" >> data/pko_p2_status.log
done
