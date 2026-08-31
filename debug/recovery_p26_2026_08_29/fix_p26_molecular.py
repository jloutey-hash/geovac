r"""Paper 26 molecular section: re-price to exact-rule values.
S_bond 0.303 -> 0.330; S_core 0.006 -> 0.008; ratio ~50 -> ~40.
R-independence itself survives EXACTLY (bit-identical blocks)."""
import io

P = "papers/group6_precision_observations/paper_26_entanglement.tex"
s = io.open(P, encoding="utf-8").read()
n = 0


def rep(old, new, label):
    global s, n
    assert s.count(old) == 1, f"{label}: {s.count(old)} matches"
    s = s.replace(old, new)
    n += 1
    print(f"  ok  {label}")


rep(r"""$R$-independent, with core/bond entropy ratios of $50\times$ at""",
    r"""$R$-independent, with core/bond entropy ratios of $40\times$ at""",
    "abstract ratio")

rep(r"""$R$-independence with its $\approx 50\times$ core/bond ratio""",
    r"""$R$-independence with its $\approx 40\times$ core/bond ratio""",
    "intro ratio")

rep(r"""  S_{\text{bond}}(R) = 0.303 \text{ nats}
  \quad \forall\, R \in [0.5, 10.0] \text{ bohr}.
  \label{eq:rindep}
\end{equation}""",
    r"""  S_{\text{bond}}(R) = 0.330 \text{ nats}
  \quad \forall\, R \in [0.5, 10.0] \text{ bohr}.
  \label{eq:rindep}
\end{equation}
(Value corrected 2026-08-29 with the exact global-$M_L$ ERI rule, from the
retired 0.303; the $R$-independence itself is unchanged and exact --- the
electronic blocks are bit-identical across $R$.)""",
    "eq:rindep value")

rep(r"""For LiH in the composed encoding, the core block
($Z_{\text{eff}} = 3$) has $S_{\text{core}} = 0.006$ nats, while
the bond block ($Z_{\text{eff}} = 1$) has
$S_{\text{bond}} = 0.303$ nats.  The ratio is
\begin{equation}
  \frac{S_{\text{bond}}}{S_{\text{core}}} \approx 50,
  \label{eq:ratio}
\end{equation}""",
    r"""For LiH in the composed encoding, the core block
($Z_{\text{eff}} = 3$) has $S_{\text{core}} = 0.008$ nats, while
the bond block ($Z_{\text{eff}} = 1$) has
$S_{\text{bond}} = 0.330$ nats.  The ratio is
\begin{equation}
  \frac{S_{\text{bond}}}{S_{\text{core}}} \approx 40,
  \label{eq:ratio}
\end{equation}
(both values corrected 2026-08-29 with the exact ERI rule; the retired
figures were $0.006$, $0.303$ and $\approx 50$)""",
    "ratio section values")

io.open(P, "w", encoding="utf-8").write(s)
print(f"\n{n} edits applied")
