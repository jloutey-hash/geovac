#!/usr/bin/env bash
# Three-pass pdflatex + gate check in ONE tool call.
# Usage: bash debug/compile_3pass.sh papers/group1_operator_algebras/paper_NN.tex
# Prints ERRORS/UNDEF/MULTIDEF counts, the output line, and GATE: PASS|FAIL.
set -u
f="$1"
d=$(dirname "$f")
b=$(basename "$f" .tex)
cd "$d" || { echo "GATE: FAIL (bad path)"; exit 2; }
for i in 1 2 3; do
    pdflatex -interaction=nonstopmode "$b.tex" > /dev/null 2>&1
done
errs=$(grep -c "^!" "$b.log" 2>/dev/null || true)
undef=$(grep -c "There were undefined" "$b.log" 2>/dev/null || true)
multi=$(grep -c "multiply-defined" "$b.log" 2>/dev/null || true)
pages=$(grep -o "Output written.*" "$b.log" | head -1)
echo "ERRORS=$errs UNDEF_BLOCKS=$undef MULTIDEF=$multi"
echo "$pages"
if [ "${errs:-1}" -eq 0 ] && [ "${undef:-1}" -eq 0 ] && [ "${multi:-1}" -eq 0 ]; then
    echo "GATE: PASS"
else
    echo "GATE: FAIL"
fi
