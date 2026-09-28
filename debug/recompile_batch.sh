#!/usr/bin/env bash
# Batch three-pass compile of every paper touched since the last commit (v5.16.0 prep).
# Writes one summary line per file to debug/data/recompile_batch.log.
cd /c/Users/jlout/Desktop/Project_Geometric
out=debug/data/recompile_batch.log
: > "$out"
export MIKTEX_AUTOINSTALL=1
while read -r f; do
  [ -z "$f" ] && continue
  d=$(dirname "$f"); b=$(basename "$f" .tex)
  ( cd "$d" && latexmk -pdf -interaction=nonstopmode -halt-on-error "$b.tex" > "../../debug/data/recompile_${b}.log" 2>&1 )
  rc=$?
  undef=$(grep -ciE "undefined (reference|citation)" "debug/data/recompile_${b}.log")
  echo "rc=$rc undef_lines=$undef $f" >> "$out"
done < debug/data/recompile_list.txt
echo "BATCH_DONE" >> "$out"
