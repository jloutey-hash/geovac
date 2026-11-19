#!/bin/bash
# usage: hybrid_ratio_batch.sh tag case region [region ...]  -- one subprocess per region, 20 min cap each
cd /c/Users/jlout/Desktop/Project_Geometric
tag=$1; case=$2; shift 2
for reg in "$@"; do
  timeout 1200 python debug/hybrid_ratio_analyse.py "$case" "$reg" > "debug/data/hybrid_ratio_${tag}_${reg}.log" 2>&1
  echo "$tag $reg exit $?"
done
