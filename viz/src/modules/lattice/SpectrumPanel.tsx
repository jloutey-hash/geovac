import { useMemo, useState } from 'react'
import type { SpectrumData } from '../../lib/types'
import { fmtSci } from '../../lib/data'

type View = 'energy' | 'dimensionless'

interface Level {
  value: number
  count: number
  label?: string
}

// Group near-degenerate eigenvalues into displayed levels.
function groupLevels(values: number[], tol = 1e-9): Level[] {
  const sorted = [...values].sort((a, b) => a - b)
  const out: Level[] = []
  for (const v of sorted) {
    const last = out[out.length - 1]
    if (last && Math.abs(v - last.value) <= tol * Math.max(1, Math.abs(v))) {
      last.count += 1
    } else {
      out.push({ value: v, count: 1 })
    }
  }
  return out
}

function LevelColumn({
  levels,
  x,
  width,
  yOf,
  color,
  onHover,
}: {
  levels: Level[]
  x: number
  width: number
  yOf: (v: number) => number
  color: string
  onHover: (tip: { x: number; y: number; text: string } | null) => void
}) {
  return (
    <g>
      {levels.map((lv, i) => (
        <line
          key={i}
          x1={x}
          x2={x + width}
          y1={yOf(lv.value)}
          y2={yOf(lv.value)}
          stroke={color}
          strokeWidth={2}
          style={{ cursor: 'default' }}
          onMouseMove={(e) =>
            onHover({
              x: e.clientX,
              y: e.clientY,
              text: `${lv.label ? lv.label + ' · ' : ''}${fmtSci(lv.value, 6)} (degeneracy ${lv.count})`,
            })
          }
          onMouseLeave={() => onHover(null)}
        />
      ))}
    </g>
  )
}

export default function SpectrumPanel({ spectrum, Z }: { spectrum: SpectrumData; Z: number }) {
  const [view, setView] = useState<View>('energy')
  const [tip, setTip] = useState<{ x: number; y: number; text: string } | null>(null)

  // Data is exported at Z=1; the artifact carries the exact scaling rules
  // (kinetic_scale.z_scaling = 'kappa * Z^2'; hydrogen levels ~ Z^2).
  const zz = Z * Z

  const { left, right, leftTitle, rightTitle, unit } = useMemo(() => {
    if (view === 'energy') {
      return {
        left: groupLevels(spectrum.graph_energies_ha.map((v) => v * zz)),
        right: spectrum.hydrogen_levels.map((h) => ({
          value: h.energy_ha * zz,
          count: h.degeneracy,
          label: `n=${h.n}`,
        })),
        leftTitle: 'Graph H = κZ²(D−A)',
        rightTitle: 'Hydrogen −Z²/2n²',
        unit: 'Ha',
      }
    }
    return {
      left: groupLevels(spectrum.graph_laplacian_spectrum),
      right: spectrum.s3_reference.levels.map((l) => ({
        value: l.eigenvalue,
        count: l.degeneracy,
        label: `n=${l.n}`,
      })),
      leftTitle: 'Graph Laplacian D−A',
      rightTitle: 'S³ reference n²−1',
      unit: 'dimensionless',
    }
  }, [spectrum, view, zz])

  const W = 380
  const H = 300
  const pad = { top: 26, bottom: 10, left: 46, right: 8 }
  const all = [...left.map((l) => l.value), ...right.map((l) => l.value)]
  const lo = Math.min(...all)
  const hi = Math.max(...all)
  const span = hi - lo || 1
  const yOf = (v: number) =>
    pad.top + (1 - (v - (lo - 0.04 * span)) / (span * 1.08)) * (H - pad.top - pad.bottom)

  const colW = (W - pad.left - pad.right - 30) / 2
  const x1 = pad.left
  const x2 = pad.left + colW + 30

  const ticks = 5
  const tickVals = Array.from({ length: ticks + 1 }, (_, i) => lo + (span * i) / ticks)

  return (
    <div>
      <div className="controls" style={{ marginBottom: 4 }}>
        <label>
          view
          <select value={view} onChange={(e) => setView(e.target.value as View)}>
            <option value="energy">energy (Ha)</option>
            <option value="dimensionless">dimensionless (Laplacian vs S³)</option>
          </select>
        </label>
      </div>
      <svg width="100%" viewBox={`0 0 ${W} ${H}`} role="img" aria-label="spectrum level diagram">
        {tickVals.map((v, i) => (
          <g key={i}>
            <line x1={pad.left} x2={W - pad.right} y1={yOf(v)} y2={yOf(v)} stroke="var(--grid)" strokeWidth={1} />
            <text x={pad.left - 5} y={yOf(v) + 3.5} textAnchor="end" fontSize={9.5} fill="var(--text-muted)">
              {Math.abs(v) < span * 1e-9 ? '0' : fmtSci(v, 3)}
            </text>
          </g>
        ))}
        <text x={x1} y={14} fontSize={11} fill="var(--text-secondary)">
          {leftTitle}
        </text>
        <text x={x2} y={14} fontSize={11} fill="var(--text-secondary)">
          {rightTitle}
        </text>
        <LevelColumn levels={left} x={x1} width={colW} yOf={yOf} color="var(--series-1)" onHover={setTip} />
        <LevelColumn levels={right} x={x2} width={colW} yOf={yOf} color="var(--series-3)" onHover={setTip} />
        {right.map(
          (lv, i) =>
            lv.label &&
            i < 3 && (
              <text
                key={i}
                x={x2 + colW + 3}
                y={yOf(lv.value) + 3.5}
                fontSize={9.5}
                fill="var(--text-secondary)"
              >
                {lv.label}
              </text>
            ),
        )}
        <text x={W - pad.right} y={H - 2} textAnchor="end" fontSize={9.5} fill="var(--text-muted)">
          {unit}
        </text>
      </svg>
      {tip && (
        <div className="tooltip" style={{ left: tip.x + 12, top: tip.y + 12 }}>
          {tip.text}
        </div>
      )}
      <p className="small muted" style={{ marginTop: 2 }}>
        {view === 'energy' ? spectrum.kinetic_scale.note : spectrum.s3_reference.note}
      </p>
      <details>
        <summary>table view</summary>
        <table className="data">
          <thead>
            <tr>
              <th>{leftTitle}</th>
              <th className="num">value</th>
              <th className="num">degeneracy</th>
            </tr>
          </thead>
          <tbody>
            {left.map((lv, i) => (
              <tr key={i}>
                <td className="muted">{i === 0 ? 'ground' : `level ${i}`}</td>
                <td className="num">{fmtSci(lv.value, 6)}</td>
                <td className="num">{lv.count}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </details>
    </div>
  )
}
