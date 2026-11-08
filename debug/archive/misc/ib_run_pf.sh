#!/bin/bash
# IB4: the PF kite (1, 2; 3/2, 3) at the primes listed in data/ib_primes.txt, indices from $1 to $2 (0-based)
cd "$(dirname "$0")"
read -a PS < data/ib_primes.txt
for ((i=$1; i<=$2; i++)); do
  p=${PS[$i]}
  IB_PRIME=$p python ib_kite_de.py pf_p$i 1 2 3/2 3 5 1 80 > data/ib_pf_p${i}_run.log 2>&1
  echo "pf_p$i exit=$?" >> data/ib_pf_status.log
done
