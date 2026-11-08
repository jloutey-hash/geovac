#!/bin/bash
# IB2 value by value: gamma = 2, 3, 1/2, 1 on the slice
cd "$(dirname "$0")"
for g in 2 3 1/2 1; do
  tag="g$(echo $g | tr '/' '_')"
  python ib_kite_de.py $tag 1 1 $g $g 5 1 80 > data/ib_${tag}_run.log 2>&1
  echo "$tag exit=$?" >> data/ib_gammas_status.log
done
echo ALLDONE >> data/ib_gammas_status.log
