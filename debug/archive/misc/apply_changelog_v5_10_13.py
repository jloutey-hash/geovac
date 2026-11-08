"""Prepend the v5.10.13 CHANGELOG entry, add the Sec.2 one-liner, bump the cursor."""
import io

ENTRY = """## [v5.10.13] - 2026-09-08

**The general-`V₀` secular equation, the `-2.90250` attribution, and a named wall on the repair route.** Follow-up to v5.10.12, prompted by a second PI-relayed consultation of Avery's canon whose load-bearing leg was independently re-derived and numerically verified here rather than taken on report.

### `eq:general_v0` — and what the molecular metric actually is

Projecting `<Φ_μ|H−E|Φ_ν> = 0` and substituting the Sturmian equation on the ket gives

    V C = V₀ B C ,    B = diag(β_ν)

with **every L² overlap cancelling identically** — for any local `V₀`, orthonormal configurations or not. Derived two ways (relayed, and re-derived here) and verified by its atomic specialisation: `V₀B = −p_κ·I` to 6e-10, and the reconstruction `Z·diag(R_ν) − G` **bit-identical** to the assembled `M`.

The consequence is structural. Molecularly `V₀` is not diagonal, so the matrix on the right-hand side is **the weighting potential itself — not the identity, and not the L² overlap**. That *identifies* the Shibuya–Wulfman matrix rather than merely naming it: for `V₀ = −Σ_A Z_A/|r−R_A|` its off-diagonal entries are the cross-center nuclear-attraction integrals, which is literally what `geovac/shibuya_wulfman.py` computes, per that module's own docstring. **The molecular metric is not extra structure the method acquires under generalization; it *is* `V₀`.** This closes a gap Paper 60 explicitly named as the piece it could not reconstruct.

### The `-2.90250` withdrawal gains an attribution — and it confirms v5.10.12

The relay states that in that calculation the scaling parameter was **scanned as a free variational parameter and minimized**, unlinked from the output eigenvalue — i.e. our *scale-optimized* posing, not the locked one. Three internal lines already required exactly that: `eq:no_selection` makes 1.2 mHa unreachable in the locked posing at any K or selection; our free-scale ladder independently sits at 1.64 mHa at comparable K; and the relay reproduces our mechanism (1s² forcing both electrons to `Q_ν = p_κ/√2`) unprompted. Recorded in the paper as **secondary-source**, not a quotation; a primary-source check remains the one thing that could overturn the section.

**The irony is worth stating plainly:** Avery's own published practice is the *accurate* posing, and the metric-free one the quantum-computing case needs is the one that floors. They are different methods, and Paper 60 had been citing a number from the first to characterize the second.

**One phrase in the relay is not usable and is not cited:** it describes the calculation as evaluating "the metric-free matrix at each step" *while* unlinking `p_κ` from `E`. By our own `eq:scale_lock` those cannot both hold — unlinking the scale is precisely what resurrects the L² metric. The operational claim is confirmed three ways; that phrase is loose paraphrase.

### The `[OPEN]` V₀ question is now a wall with a mechanism

No weighting potential of different radial shape is used for atoms anywhere in the book or theses; all practical work is bare Coulomb. The stated cost is concrete: **losing `1/r` breaks the mapping to hyperspherical harmonics under the Fock projection**, and with it the closed-form multi-center and inter-electron integrals — the property the whole encoding rests on. That converges with what we derived from the other side (parameter-freeness dies when `V₀` carries an intrinsic length scale). The open question is therefore sharper than "does a better `V₀` exist": it is whether **any** `V₀` of different radial shape *preserves the Fock mapping*. If none does, the trade is a wall rather than unexplored ground.

Also relayed honestly and therefore **not** used: no convergence tables and no quantification of the excited-state advantage are accessible to that source, so our measured 4× state-dependence remains ours alone, uncorroborated either way.

### Backing

`eq:general_v0` backed by its own guard, written as a separate adversarial pass per §9 and fire-tested both ways, including the Z-independence leg that catches a coincidence at Z=2. Gates: C10 / C21 / C16 / C17 / C14 PASS. `memory/avery_method_and_prior_art_gaps.md` updated — the item recorded there since v5.10.12 as "the one that could falsify us" is resolved, and the molecular equation and the `V₀` obstruction are recorded as reference facts.

"""

ch = io.open("CHANGELOG.md", encoding="utf-8").read()
anchor = "## [v5.10.12] - 2026-09-08"
assert anchor in ch and "v5.10.13" not in ch
ch = ch.replace(anchor, ENTRY + anchor, 1)
io.open("CHANGELOG.md", "w", encoding="utf-8").write(ch)
print("CHANGELOG: v5.10.13 prepended")

cm = io.open("CLAUDE.md", encoding="utf-8").read()
BULLET = ("- **General-V0 secular equation + the -2.90250 attribution (2026-09-08, v5.10.13):** "
          "V C = V_0 B C, L2 cancels for ANY local V_0; the molecular metric IS V_0 (= Shibuya-Wulfman). "
          "Avery's cited number used the SCALE-OPTIMIZED posing. See CHANGELOG v5.10.13.\n")
a2 = "- **P60 floor is the scale lock, not l_max (2026-09-08, v5.10.12):**"
assert a2 in cm and "General-V0 secular equation" not in cm
cm = cm.replace(a2, BULLET + a2, 1)
cm = cm.replace("**Version:** v5.10.12 (September 8, 2026)",
                "**Version:** v5.10.13 (September 8, 2026)", 1)
io.open("CLAUDE.md", "w", encoding="utf-8").write(cm)
print("CLAUDE.md: Sec.2 one-liner + cursor -> v5.10.13")
