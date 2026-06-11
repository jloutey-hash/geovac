# L3e-P3 Phase A.2 — Operator-algebraic ε-net definition

**Date:** 2026-05-23 (parallel with A.1 literature review running in background).
**Sprint position:** Phase A.2 of Sprint L3e-P3 (scoping at `debug/sprint_l3e_p3_detailed_scoping_memo.md`). The diagnostic phase. Defines the operator-algebraic analog of the Mondino-Sämann synthetic Lorentzian causal-diamond ε-net structure, in advance of A.3 (correspondence theorem) and A.4 (apply to GeoVac wedge).
**Honest scope:** This is a definitional / structural memo, not a theorem proof. The synthetic-side definitions are drawn from the May-17 L3 literature audit (`debug/l3_literature_audit_memo.md`); the operator-algebraic-side definitions are constructed here for the first time. Phase A.3 will verify the correspondence rigorously; Phase A.4 will instantiate at the GeoVac wedge.

---

## §1. Synthetic-side ingredients (recap)

From `debug/l3_literature_audit_memo.md` §1 + Mondino-Sämann lineage (refresh pending from A.1 sub-agent):

A **Lorentzian length space** is a quintuple $(X, d, \le, \ll, \tau)$ where:
- $X$ is a topological space
- $d$ is a metric (background Riemannian-style distance, for the topology)
- $\le$ is a partial order on $X$ (causal precedence)
- $\ll \subseteq \le$ is a strict partial order (chronological precedence)
- $\tau : X \times X \to [0, \infty]$ is the **time-separation function** satisfying:
  - $\tau(x, y) = 0$ unless $x \le y$
  - $\tau(x, z) \ge \tau(x, y) + \tau(y, z)$ when $x \le y \le z$ (reverse triangle inequality)
  - Lower semicontinuity

The **causal diamond** $J(x, y) := \{z : x \le z \le y\}$ for $x \le y$ is the elementary region of the causal structure.

A **causal-diamond ε-net** on $(X, \le, \tau)$ at scale ε is a finite (or discrete) subset $X_\varepsilon \subset X$ such that:
- Every causal diamond $J(x, y)$ with $\tau(x, y) \le \varepsilon$ contains at most one point of $X_\varepsilon$
- The collection $\{J(p, q) : p, q \in X_\varepsilon, \tau(p, q) \le 2\varepsilon\}$ covers $X$
- (Cardinality bound: roughly $|X_\varepsilon| \sim \mathrm{Vol}(X)/\varepsilon^d$ in $d$ dimensions)

**GH-convergence** in Mondino-Sämann sense: $(X_n, \le_n, \tau_n) \to (X_\infty, \le_\infty, \tau_\infty)$ if there exist scale-$\varepsilon_n \to 0$ ε-nets $X_n^{(\varepsilon_n)}$ in each space that align via almost-isometries of the time-separation function (modulo $\varepsilon_n$ error).

---

## §2. Operator-algebraic ε-net definition

We construct the operator-algebraic analog by replacing each synthetic-side ingredient with a Krein-spectral-triple-side counterpart:

### Definition 2.1 (Operator-algebraic Lorentzian length space)

An **operator-algebraic Lorentzian length space** is a five-tuple
$$\bigl(\Tcal^L, \mathcal{S}^+(\Tcal^L), \preceq, \prec, \tau_{\mathrm{mod}}\bigr)$$
where:
- $\Tcal^L = (\Acal, \Krein, D^L, J^L)$ is a Lorentzian Krein spectral triple (Paper 43 §3 / Paper 44 substrate)
- $\mathcal{S}^+(\Tcal^L)$ is the set of Krein-positive states (probabilistic states on the algebra restricted to the $K^+$ cone, Paper 44 §3)
- $\preceq$ is the **modular precedence relation**: $\omega \preceq \omega'$ iff there exists $t \ge 0$ such that $\omega' = \omega \circ \sigma_t^{\omega^*}$ for the modular flow $\sigma_t^{\omega^*}$ of a chosen reference state $\omega^*$
- $\prec$ is the strict version of $\preceq$ ($t > 0$)
- $\tau_{\mathrm{mod}}(\omega, \omega') := \inf\{t \ge 0 : \omega' \in \overline{\{\omega \circ \sigma_s : s \in [0, t]\}}^{\mathrm{w}^*}\}$ — the minimal modular-flow time to reach $\omega'$ from $\omega$, with $\tau_{\mathrm{mod}}(\omega, \omega') = \infty$ if no such $t$ exists.

**Key idea:** the modular flow $\sigma_t$ of the wedge KMS state $\omega^* = \omega_W^L$ (Paper 43 §4.2) plays the role of the synthetic time-separation function. The wedge geometry provides the causal structure intrinsically via the Bisognano-Wichmann theorem (modular flow IS the Lorentz boost orbit on the wedge).

### Definition 2.2 (Modular causal diamond)

The **modular causal diamond** $J_{\mathrm{mod}}(\omega, \omega')$ for $\omega \preceq \omega'$ with $\tau_{\mathrm{mod}}(\omega, \omega') = t$ is
$$J_{\mathrm{mod}}(\omega, \omega') := \{\omega'' \in \mathcal{S}^+(\Tcal^L) : \omega \preceq \omega'' \preceq \omega'\} = \{\omega \circ \sigma_s^{\omega^*} : 0 \le s \le t\}.$$
This is the modular-flow orbit segment from $\omega$ to $\omega'$.

### Definition 2.3 (Operator-algebraic ε-net)

Given a Lorentzian Krein spectral triple $\Tcal^L$ (continuum) and a scale $\varepsilon > 0$, an **operator-algebraic ε-net at scale ε** on $\Tcal^L$ is a tuple
$$\bigl(\mathcal{E}^L_\varepsilon, P_\varepsilon, \omega^*_\varepsilon, D^L_\varepsilon\bigr)$$
where:
- $\mathcal{E}^L_\varepsilon \subset \mathcal{S}^+(\Tcal^L)$ is a finite (or discrete) subset of Krein-positive states
- $P_\varepsilon : \Acal \to \Acal_\varepsilon \subset \Acal$ is a finite-rank projection onto a finite-dimensional operator subsystem (the truncation map)
- $\omega^*_\varepsilon \in \mathcal{E}^L_\varepsilon$ is a distinguished base state (the canonical observation point, typically the BW wedge vacuum)
- $D^L_\varepsilon := P_\varepsilon D^L P_\varepsilon$ is the truncated Lorentzian Dirac

satisfying:
- (Discreteness) $|\mathcal{E}^L_\varepsilon| < \infty$
- (Modular-diamond covering) For every $\omega \in \mathcal{S}^+(\Tcal^L)$, there exist $\omega', \omega'' \in \mathcal{E}^L_\varepsilon$ with $\omega \in J_{\mathrm{mod}}(\omega', \omega'')$ and $\tau_{\mathrm{mod}}(\omega', \omega'') \le 2\varepsilon$
- (Truncation compatibility) The Lipschitz seminorm $L^L_\varepsilon(a) := \opnorm{[D^L_\varepsilon, a]}$ on $\Acal_\varepsilon$ satisfies
$|L^L_\varepsilon(P_\varepsilon a) - L^L(a)| \le \varepsilon \norm{a}_{\mathrm{op}}$ for $a \in \Acal$ (asymptotic compatibility with the continuum Lipschitz seminorm)

### Definition 2.4 (Convergence of operator-algebraic ε-nets)

A sequence of operator-algebraic ε-nets $\{\mathcal{E}^L_{\varepsilon_n}\}$ on (possibly different) truncated Lorentzian Krein spectral triples $\Tcal^L_n$ **converges to a continuum ε-net** on $\Tcal^L_\infty$ as $\varepsilon_n \to 0$ if:
- The associated truncation projections $P_{\varepsilon_n}$ have $\norm{P_{\varepsilon_n} - I}_{\mathrm{op}} \to 0$ on a dense subspace of $\Acal$
- The base states $\omega^*_{\varepsilon_n} \to \omega^*_\infty$ in the weak-* topology of $\mathcal{S}^+(\Tcal^L_\infty)$
- The truncated Lipschitz seminorms $L^L_{\varepsilon_n}$ converge to $L^L_\infty$ in some appropriate Lipschitz-uniform topology (a la Latrémolière)

---

## §3. How the GeoVac construction fits

### Identification at finite cutoff $(n_{\max}, N_t)$

For the GeoVac truncated Lorentzian Krein spectral triple $\Tcal^L_{n_{\max}, N_t}$ (Paper 43 §3 + Paper 44 substrate):

| Operator-algebraic ε-net ingredient | GeoVac realization |
|:------------------------------------|:-------------------|
| $\mathcal{E}^L_\varepsilon$ (finite state set) | Finite-dimensional $\mathcal{S}^+(\Tcal^L_{n_{\max}, N_t})$ on the truncated wedge $W_L$ (Paper 43 Def 4.1, $P_{W_L} = P_W^{\mathrm{spatial}} \otimes P_{t \ge 0}$); $\dim \mathcal{S}^+ = $ trace-class density matrices on $\dim W_L = (\dim \Hilb_{\mathrm{GV}}/2) \times N_{t,+}$ |
| $P_\varepsilon$ (truncation projection) | $P_{n_{\max}, N_t} = P_{n_{\max}} \otimes P_{N_t}$ — Paper 44 spatial + Paper 43 temporal cutoffs |
| $\omega^*_\varepsilon$ (base state) | $\omega_W^L = \omega_{\rho_W^L}$ where $\rho_W^L = e^{-K_L^{\alpha,W}}/Z$ is the BW wedge KMS state (Paper 43 §4.2) |
| $D^L_\varepsilon$ (truncated Dirac) | $\DL = i(\gamma^0 \otimes \partial_t + \DGV \otimes I_{N_t})$ restricted to $W_L$ (Paper 43 §3 + Paper 44 Lorentzian Dirac construction) |

### Scale $\varepsilon$ identification

The natural scale at finite cutoff $(n_{\max}, N_t)$ is
$$\varepsilon(n_{\max}, N_t) := \max\bigl(\gamma^{\SU(2)}_{n_{\max}}, \gamma^{U(1)}_{N_t, T}\bigr) = O\bigl(\log n_{\max}/n_{\max} + T/N_t\bigr)$$
— the Paper 45 / Paper 46 joint γ-rate. This is exactly the propinquity-decay rate at the cutoff, and matches the synthetic-side ε-net scale via the master Mellin engine M1 mechanism (the $1/\pi^2$ Hopf-base measure signature).

### Modular causal diamond at GeoVac wedge

For $\omega, \omega' \in \mathcal{S}^+(\Tcal^L_{n_{\max}, N_t})$ both supported on the wedge $W_L$, the modular flow $\sigma_t^{\omega_W^L}(\cdot) = e^{i t K_L^{\alpha, W}} \cdot e^{-i t K_L^{\alpha, W}}$ is the **polar rotation generator** action on observables. The time-separation $\tau_{\mathrm{mod}}(\omega, \omega')$ is the minimal $t$ such that $\omega'$ lies on the modular-flow orbit of $\omega$. By the integer-spectrum property of $K_L^{\alpha, W}$ (Paper 42 §5 / Paper 43 §4.4), the modular flow has period $2\pi$ exactly, so $\tau_{\mathrm{mod}}$ is bounded by $2\pi$ on each modular-flow orbit.

**This is the substantive new structural identification of A.2:** the modular flow on the wedge IS the operator-algebraic analog of the Mondino-Sämann time-separation function, with the BW canonical period $\beta = 2\pi$ playing the role of "diameter" of the synthetic Lorentzian length space.

---

## §4. Verification of synthetic-side axioms (informal)

To call the operator-algebraic structure of §2 a "Lorentzian length space" in the Mondino-Sämann sense, we need to verify the five axioms:

### Axiom 1: Topological structure

The state space $\mathcal{S}^+(\Tcal^L)$ has a natural weak-* topology inherited from $\Acal$. Topological structure is automatic.

### Axiom 2: Partial order (causal precedence)

$\omega \preceq \omega'$ iff $\omega' \in \overline{\{\omega \circ \sigma_t^{\omega^*} : t \ge 0\}}^{\mathrm{w}^*}$. Reflexivity ($\omega \preceq \omega$ via $t = 0$) and transitivity ($\omega \preceq \omega' \preceq \omega''$ implies $\omega \preceq \omega''$ via concatenation of modular-flow paths) follow from the semigroup structure of the modular flow.

Antisymmetry needs care: $\omega \preceq \omega'$ and $\omega' \preceq \omega$ means $\omega' = \omega \circ \sigma_t$ for some $t$ AND $\omega = \omega' \circ \sigma_s$ for some $s$. By the integer-spectrum property at GeoVac wedge, this forces $t + s \in 2\pi\Z$, so $\omega = \omega'$ (modulo trivial period). **Antisymmetry holds up to the $2\pi$ periodicity** — this is the operator-algebraic analog of cyclic time on the wedge KMS state.

### Axiom 3: Strict order

$\omega \prec \omega'$ iff $\omega \preceq \omega'$ and $\omega \neq \omega'$ (modulo periodicity). Follows from §3.

### Axiom 4: Time-separation function $\tau_{\mathrm{mod}}$

- $\tau_{\mathrm{mod}}(\omega, \omega') = 0$ unless $\omega \preceq \omega'$: by definition (Definition 2.1).
- $\tau_{\mathrm{mod}}(\omega, \omega'') \ge \tau_{\mathrm{mod}}(\omega, \omega') + \tau_{\mathrm{mod}}(\omega', \omega'')$ when $\omega \preceq \omega' \preceq \omega''$: follows from triangle inequality on modular-flow times. Actually wait — this is the REVERSE triangle inequality $\tau(x, z) \ge \tau(x, y) + \tau(y, z)$ which is characteristic of Lorentzian geometry. Modular flow times: $\tau(\omega, \omega'') = $ minimal $t$ such that $\omega''$ is reached from $\omega$. If $\omega \to \omega'$ in time $t_1$ and $\omega' \to \omega''$ in time $t_2$, then $\omega \to \omega''$ in time $t_1 + t_2$, so $\tau(\omega, \omega'') \le t_1 + t_2 = \tau(\omega, \omega') + \tau(\omega', \omega'')$. **This is the FORWARD triangle inequality, not reverse.** ⚠️
- Lower semicontinuity: holds because modular flow is weak-*-continuous in $t$.

**ISSUE IDENTIFIED (substantive A.2 finding):** the modular flow gives FORWARD triangle inequality, not the REVERSE inequality that characterizes Lorentzian time-separation. This is a **structural mismatch** between the operator-algebraic modular structure (which is "thermal" / forward-causal) and the synthetic Lorentzian structure (which is "geometric" / reverse-triangle).

This is the kind of structural finding that Phase A is designed to surface. **It does not necessarily kill the program** — it means the operator-algebraic ε-net is in a DIFFERENT category from the synthetic Lorentzian length space, with a different triangle direction.

### Axiom 5: ε-net cardinality

The truncated operator-system $\Op^L_{n_{\max}, N_t}$ has finite dimension $\dim \Op^L = \dim_{\mathrm{Weyl}}^2 \cdot N_t$ (Paper 44 §3). So $|\mathcal{E}^L_\varepsilon|$ is finite at finite cutoff, with the right scaling-with-ε behavior.

---

## §5. Substantive findings of Phase A.2 (preliminary)

Three substantive findings emerge at this scoping level (subject to Phase A.3 verification + A.1 literature deep-dive):

### F1 — Modular flow as operator-algebraic time-separation

The synthetic-side time-separation function $\tau$ has a natural operator-algebraic counterpart in the modular flow time $\tau_{\mathrm{mod}}$ of a BW-aligned wedge KMS state. This is the load-bearing structural identification of A.2 — without it, no bridge exists.

### F2 — Forward vs reverse triangle inequality (structural mismatch)

The modular flow time satisfies the FORWARD triangle inequality (additive composition), not the REVERSE triangle inequality (Lorentzian-characteristic) of the synthetic time-separation function. **This is a real structural difference.** Possible resolutions:
- (R1) Negate the time variable: use $\tau_{\mathrm{neg}}(\omega, \omega') := -\tau_{\mathrm{mod}}(\omega, \omega')$ as the "Lorentzian time-separation" — but this doesn't have the right range $[0, \infty]$
- (R2) Reinterpret synthetic-side axioms: the Mondino-Sämann definition uses the reverse-triangle inequality to encode Lorentzian causality; the operator-algebraic version might satisfy a DIFFERENT axiom set that captures the same "causal structure" content
- (R3) Identify the operator-algebraic structure with a DIFFERENT synthetic notion (perhaps Connes-Rovelli thermal time, rather than Mondino-Sämann Lorentzian time)

A.3 (correspondence theorem) and the A.1 sub-agent literature review should resolve which of R1/R2/R3 is the right path forward.

### F3 — The BW canonical period $2\pi$ as "diameter"

The integer spectrum of $K_L^{\alpha, W}$ gives period $\beta = 2\pi$ exactly on the modular flow. This is naturally the "diameter" of the operator-algebraic Lorentzian length space. The synthetic-side ε-net cardinality bound $|X_\varepsilon| \sim \mathrm{Vol}(X)/\varepsilon^d$ has its operator-algebraic counterpart in the truncated operator-system dimension scaling.

---

## §6. Implications for A.3 (correspondence theorem)

The Phase A.2 findings reshape Phase A.3 in important ways:

1. **The correspondence is NOT directly $\tau_{\mathrm{mod}} \leftrightarrow \tau_{\mathrm{synthetic}}$.** The triangle-inequality direction mismatch means the operator-algebraic structure is a *related-but-distinct* category. A.3 needs to establish the right kind of correspondence (functor-like, not isomorphism-like).

2. **Pre-compactness on synthetic side may not transport directly.** The Mondino-Sämann pre-compactness theorem uses the reverse triangle inequality + Lorentzian causality. The operator-algebraic version needs its own pre-compactness theorem, plausibly using the modular flow's $2\pi$-periodicity + finite-dimensional state-space-at-cutoff.

3. **The bridge is via the WEDGE, not the full Lorentzian length space.** The synthetic-side Lorentzian length space includes BOTH future and past causal cones; the operator-algebraic wedge $W_L$ only includes ONE causal cone (the future-directed wedge under BW choice). This is consistent with the wedge being a Rindler-style half-space, not the full Minkowski space.

4. **Gate-1 decision now depends on A.3 + A.1 outcomes simultaneously.** The Phase A scoping has two open structural questions: (a) does the forward-vs-reverse triangle inequality mismatch admit a clean resolution? (b) does the synthetic-side pre-compactness theorem have an operator-algebraic analog at all?

---

## §7. Honest scope of Phase A.2

This memo:
- IS a structural definition memo, not a theorem proof
- IDENTIFIES the key bridge ingredients (modular flow as time-separation, wedge KMS state as base point, integer spectrum as periodicity, propinquity rate as ε-scale)
- FLAGS a substantive structural mismatch (forward vs reverse triangle inequality) that A.3 needs to resolve
- DOES NOT depend on A.1 literature output — uses the May-17 audit memo for synthetic-side definitions
- DOES NOT verify the synthetic axioms rigorously — A.3 is the verification phase

When A.1 literature audit returns, the synthetic-side definitions in §1 may need refinement. The A.2 ingredients in §2 are robust to the expected refinements; the F2 finding (forward-vs-reverse triangle inequality) is a real algebraic statement and survives any literature-side reformulation.

**Confidence:** HIGH on F1 (modular flow as time-separation analog). HIGH on F3 ($2\pi$ as natural diameter). MEDIUM on the overall correspondence framework — depends on F2 resolution.

**Files:** `debug/l3e_p3_phase_a2_operator_algebraic_eps_net.md` (this memo, ~3000 words). Cross-references: `papers/group1_operator_algebras/paper_43_lorentzian_extension.tex` §4 wedge structure, §4.2 BW vacuum; `paper_44_lorentzian_operator_system.tex` substrate; `paper_45_lorentzian_propinquity.tex` joint γ rate; `paper_47_two_rate_hybrid_convergence.tex` two-rate hybrid; A.1 sub-agent's literature audit (pending).

---

## §8. Next steps (Phase A.3)

When A.1 literature audit returns, A.3 will:
1. Verify or refine the synthetic-side axioms in §1
2. Resolve the F2 forward-vs-reverse triangle inequality mismatch (R1, R2, or R3 from §5)
3. State the correspondence theorem rigorously
4. Verify at small panel $(n_{\max}, N_t) = (2, 3)$ that the operator-algebraic ε-net data structure is well-defined for GeoVac

A.4 then applies to the full GeoVac wedge at panel cells. A.5 synthesizes and makes the Gate-1 decision.
