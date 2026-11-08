# Paper 34 support

`tx_a_signature_arithmetic.py` encodes projections 1-15 of Paper 34's dictionary
as graded maps and computes the TX-A theorems of `sec:axiomatization`
(dimensional consistency on six worked chains, order-independence, and the
15x15 composition table). Its outputs are `data/tx_a_projection_table.json`
and `data/tx_a_composition_table.json`.

Moved here from `debug/` on 2026-11-08 so the paper's citation points at a
permanent home and a backing test can import it (C22 forbids paper-backing
tests that depend on `debug/`). The same day its metadata was brought into line
with the paper: the Hopf bundle adds no variable, projection 13 is the
Drake-Swainson normalization of Bethe-log sums, and neither the hydrogen Lamb
shift (textbook assembly) nor the Uehling row (entered flat-space Pi) is a
worked chain.
