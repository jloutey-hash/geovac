import { useMemo, useState } from 'react'
import { useDarkMode } from '../../lib/route'
import { fmt } from '../../lib/data'

// Qubit co-occurrence heatmap: cell (i, j) counts Pauli terms acting
// non-trivially on both qubit i and qubit j (diagonal: on qubit i at all).
function coOccurrence(terms: [string, number, number][], Q: number): number[][] {
  const M: number[][] = Array.from({ length: Q }, () => Array(Q).fill(0))
  for (const [s] of terms) {
    if (s === 'I') continue
    const qs = [...s.matchAll(/[XYZ](\d+)/g)].map((m) => Number(m[1]))
    for (const a of qs) {
      M[a][a] += 1
      for (const b of qs) if (b > a) {
        M[a][b] += 1
        M[b][a] += 1
      }
    }
  }
  return M
}

function hexLerp(a: string, b: string, t: number): string {
  const pa = [1, 3, 5].map((i) => parseInt(a.slice(i, i + 2), 16))
  const pb = [1, 3, 5].map((i) => parseInt(b.slice(i, i + 2), 16))
  const c = pa.map((v, i) => Math.round(v + (pb[i] - v) * t))
  return `#${c.map((v) => v.toString(16).padStart(2, '0')).join('')}`
}

// Sequential single-hue ramp (blue), stepped per mode: near-zero recedes
// toward the surface, maximum is the darkest (light) / lightest (dark) step.
function rampFor(dark: boolean): (t: number) => string {
  const stops = dark ? ['#222b39', '#2a6ac0', '#b7d3f6'] : ['#eaf2fd', '#5598e7', '#0d366b']
  return (t: number) => {
    const x = Math.min(1, Math.max(0, t)) * (stops.length - 1)
    const i = Math.min(stops.length - 2, Math.floor(x))
    return hexLerp(stops[i], stops[i + 1], x - i)
  }
}

export default function PauliHeatmap({
  terms,
  Q,
}: {
  terms: [string, number, number][]
  Q: number
}) {
  const dark = useDarkMode()
  const [tip, setTip] = useState<{ x: number; y: number; text: string } | null>(null)
  const M = useMemo(() => coOccurrence(terms, Q), [terms, Q])
  const max = Math.max(1, ...M.flat())
  const ramp = rampFor(dark)

  const cell = Q > 24 ? 11 : 15
  const pad = 26
  const size = pad + Q * cell + 4

  return (
    <div>
      <svg
        width="100%"
        style={{ maxWidth: size * 1.4 }}
        viewBox={`0 0 ${size} ${size}`}
        role="img"
        aria-label="qubit co-occurrence heatmap"
      >
        {M.map((row, i) =>
          row.map((v, j) => (
            <rect
              key={`${i}-${j}`}
              x={pad + j * cell}
              y={pad + i * cell}
              width={cell - 1}
              height={cell - 1}
              rx={1.5}
              fill={v === 0 ? 'var(--surface-1)' : ramp(v / max)}
              stroke="var(--grid)"
              strokeWidth={0.4}
              onMouseMove={(e) =>
                setTip({
                  x: e.clientX,
                  y: e.clientY,
                  text:
                    i === j
                      ? `qubit ${i}: acted on by ${fmt(v)} terms`
                      : `qubits ${i}·${j}: ${fmt(v)} shared terms`,
                })
              }
              onMouseLeave={() => setTip(null)}
            />
          )),
        )}
        {Array.from({ length: Q }, (_, i) => i)
          .filter((i) => i % 5 === 0)
          .map((i) => (
            <g key={i}>
              <text x={pad + i * cell + cell / 2} y={pad - 6} textAnchor="middle" fontSize={8.5} fill="var(--text-muted)">
                {i}
              </text>
              <text x={pad - 6} y={pad + i * cell + cell / 2 + 3} textAnchor="end" fontSize={8.5} fill="var(--text-muted)">
                {i}
              </text>
            </g>
          ))}
      </svg>
      {tip && (
        <div className="tooltip" style={{ left: tip.x + 12, top: tip.y + 12 }}>
          {tip.text}
        </div>
      )}
      <p className="small muted">
        Cell (i, j): number of Pauli terms acting non-trivially on both qubits.
        Block structure reflects the selection-rule sparsity of the encoding.
      </p>
    </div>
  )
}
