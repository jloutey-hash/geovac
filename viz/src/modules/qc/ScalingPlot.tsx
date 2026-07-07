import { useMemo, useState } from 'react'
import type { LibraryData, ScalingFit } from '../../lib/types'
import { fmt } from '../../lib/data'

interface Pt {
  x: number
  y: number
  label: string
  family: string
  color: string
}

interface Family {
  key: string
  color: string
  legend: string
  pts: Pt[]
  fit?: ScalingFit
  fitLabel?: string
}

function fitLabel(fit: ScalingFit | undefined): string | undefined {
  // A fit carrying a note (e.g. too few points) is not quotable: points only.
  if (!fit || fit.exponent === undefined || fit.note) return undefined
  const resid = fit.rms_log_residual !== undefined ? fit.rms_log_residual.toFixed(2) : '–'
  return `∝ Q^${fit.exponent.toFixed(2)} (${fit.n_points} pts, rms log-resid ${resid})`
}

export default function ScalingPlot({ lib }: { lib: LibraryData }) {
  const [tip, setTip] = useState<{ x: number; y: number; text: string } | null>(null)

  const families = useMemo<Family[]>(() => {
    const fams: Family[] = []
    const cs = lib.scaling.cross_system
    fams.push({
      key: 'cross',
      color: 'var(--series-1)',
      legend: 'library, cross-system (default build)',
      pts: cs.points.map((p) => ({
        x: p.Q,
        y: p.n_pauli,
        label: p.system,
        family: 'cross-system',
        color: 'var(--series-1)',
      })),
      fit: cs.fit,
      fitLabel: fitLabel(cs.fit),
    })
    const ws = lib.scaling.within_system
    if (ws) {
      fams.push({
        key: 'within',
        color: 'var(--series-2)',
        legend: `${ws.system}, basis growth (per-fragment composed build)`,
        pts: ws.points.map((p) => ({
          x: p.Q,
          y: p.N_pauli,
          label: `${ws.system} max_n=${p.max_n}`,
          family: 'basis growth',
          color: 'var(--series-2)',
        })),
        fit: ws.fit,
        fitLabel: fitLabel(ws.fit),
      })
    }
    const ab = lib.scaling.atomic_benchmark
    if (ab) {
      fams.push({
        key: 'atomic',
        color: 'var(--series-3)',
        legend: 'atomic basis growth (benchmark artifact)',
        pts: ab.points
          .filter((p) => p.type === 'geovac')
          .map((p) => ({
            x: p.Q,
            y: p.n_pauli,
            label: p.label,
            family: 'atomic benchmark',
            color: 'var(--series-3)',
          })),
        fit: ab.fit_he,
        fitLabel: fitLabel(ab.fit_he),
      })
      const gauss = ab.points.filter((p) => p.type !== 'geovac')
      if (gauss.length) {
        fams.push({
          key: 'gauss',
          color: 'var(--series-4)',
          legend: 'Gaussian-basis reference points',
          pts: gauss.map((p) => ({
            x: p.Q,
            y: p.n_pauli,
            label: p.label,
            family: 'Gaussian reference',
            color: 'var(--series-4)',
          })),
        })
      }
    }
    return fams
  }, [lib])

  const allPts = families.flatMap((f) => f.pts)
  const xs = allPts.map((p) => p.x)
  const ys = allPts.map((p) => p.y)
  const xlo = Math.min(...xs) / 1.3
  const xhi = Math.max(...xs) * 1.3
  const ylo = Math.min(...ys) / 1.6
  const yhi = Math.max(...ys) * 1.6

  const W = 640
  const H = 380
  const pad = { top: 14, right: 16, bottom: 42, left: 62 }
  const X = (q: number) =>
    pad.left + ((Math.log10(q) - Math.log10(xlo)) / (Math.log10(xhi) - Math.log10(xlo))) * (W - pad.left - pad.right)
  const Y = (n: number) =>
    pad.top + (1 - (Math.log10(n) - Math.log10(ylo)) / (Math.log10(yhi) - Math.log10(ylo))) * (H - pad.top - pad.bottom)

  const xticks = [2, 5, 10, 20, 50, 100].filter((t) => t >= xlo && t <= xhi)
  const yticks = [10, 100, 1e3, 1e4, 1e5, 1e6].filter((t) => t >= ylo && t <= yhi)

  return (
    <div>
      <div className="legend">
        {families.map((f) => (
          <span key={f.key}>
            <span className="swatch" style={{ background: f.color }} />
            {f.legend}
            {f.fitLabel ? ` — ${f.fitLabel}` : ''}
          </span>
        ))}
      </div>
      <svg width="100%" viewBox={`0 0 ${W} ${H}`} role="img" aria-label="Pauli terms vs qubits, log-log">
        {yticks.map((t) => (
          <g key={t}>
            <line x1={pad.left} x2={W - pad.right} y1={Y(t)} y2={Y(t)} stroke="var(--grid)" />
            <text x={pad.left - 6} y={Y(t) + 3.5} textAnchor="end" fontSize={10} fill="var(--text-muted)">
              {t >= 1000 ? `10^${Math.round(Math.log10(t))}` : fmt(t)}
            </text>
          </g>
        ))}
        {xticks.map((t) => (
          <g key={t}>
            <line y1={pad.top} y2={H - pad.bottom} x1={X(t)} x2={X(t)} stroke="var(--grid)" />
            <text y={H - pad.bottom + 14} x={X(t)} textAnchor="middle" fontSize={10} fill="var(--text-muted)">
              {fmt(t)}
            </text>
          </g>
        ))}
        <line x1={pad.left} x2={W - pad.right} y1={H - pad.bottom} y2={H - pad.bottom} stroke="var(--baseline)" />
        <line x1={pad.left} x2={pad.left} y1={pad.top} y2={H - pad.bottom} stroke="var(--baseline)" />
        <text x={(W + pad.left - pad.right) / 2} y={H - 6} textAnchor="middle" fontSize={11} fill="var(--text-secondary)">
          qubits Q (log)
        </text>
        <text
          transform={`translate(14 ${(H + pad.top - pad.bottom) / 2}) rotate(-90)`}
          textAnchor="middle"
          fontSize={11}
          fill="var(--text-secondary)"
        >
          Pauli terms (log)
        </text>
        {families.map(
          (f) =>
            f.fit?.exponent !== undefined &&
            f.fit.prefactor !== undefined &&
            !f.fit.note &&
            f.pts.length >= 2 && (
              <line
                key={`fit-${f.key}`}
                x1={X(Math.min(...f.pts.map((p) => p.x)))}
                y1={Y(f.fit.prefactor * Math.pow(Math.min(...f.pts.map((p) => p.x)), f.fit.exponent))}
                x2={X(Math.max(...f.pts.map((p) => p.x)))}
                y2={Y(f.fit.prefactor * Math.pow(Math.max(...f.pts.map((p) => p.x)), f.fit.exponent))}
                stroke={f.color}
                strokeWidth={1.5}
                strokeDasharray="5 4"
                opacity={0.7}
              />
            ),
        )}
        {allPts.map((p, i) => (
          <circle
            key={i}
            cx={X(p.x)}
            cy={Y(p.y)}
            r={4.5}
            fill={p.color}
            stroke="var(--surface-1)"
            strokeWidth={1.5}
            onMouseMove={(e) =>
              setTip({
                x: e.clientX,
                y: e.clientY,
                text: `${p.label} — Q=${fmt(p.x)}, ${fmt(p.y)} Pauli terms (${p.family})`,
              })
            }
            onMouseLeave={() => setTip(null)}
          />
        ))}
      </svg>
      {tip && (
        <div className="tooltip" style={{ left: tip.x + 12, top: tip.y + 12 }}>
          {tip.text}
        </div>
      )}
      <p className="small muted">
        Families are distinct datasets and are not interchangeable: the library
        series varies the molecule at a fixed basis; the basis-growth series grow
        the basis for a fixed system. Fit exponents, point counts, and residuals
        are read from the generated data files
        {lib.scaling.within_system?.fit.note ? ` — ${lib.scaling.within_system.fit.note}` : ''}.
      </p>
      <details>
        <summary>table view</summary>
        <table className="data">
          <thead>
            <tr>
              <th>point</th>
              <th>family</th>
              <th className="num">Q</th>
              <th className="num">Pauli terms</th>
            </tr>
          </thead>
          <tbody>
            {allPts.map((p, i) => (
              <tr key={i}>
                <td>{p.label}</td>
                <td className="muted">{p.family}</td>
                <td className="num">{fmt(p.x)}</td>
                <td className="num">{fmt(p.y)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </details>
    </div>
  )
}
