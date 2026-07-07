import { useState } from 'react'
import { useData, fmt } from '../../lib/data'
import { useDarkMode } from '../../lib/route'
import type { LatticeData, SpectrumData } from '../../lib/types'
import LatticeScene, { type LayoutKind } from './LatticeScene'
import SpectrumPanel from './SpectrumPanel'

const MAX_N_CHOICES = [2, 3, 4, 5, 6, 7, 8, 9, 10]

export default function LatticeExplorer() {
  const [maxN, setMaxN] = useState(5)
  const [Z, setZ] = useState(1)
  const [layout, setLayout] = useState<LayoutKind>('paraboloid')
  const [selected, setSelected] = useState<number | null>(null)
  const dark = useDarkMode()

  const { data: lattice, error: latErr } = useData<LatticeData>(
    `lattice/lattice_Z1_n${maxN}.json`,
  )
  const { data: spectrum } = useData<SpectrumData>(`spectrum/spectrum_Z1_n${maxN}.json`)

  const sel = selected !== null && lattice && selected < lattice.states.length
    ? lattice.states[selected]
    : null

  return (
    <div>
      <h1>Quantum-state lattice</h1>
      <p className="secondary" style={{ maxWidth: 760 }}>
        Each node is a hydrogen-like quantum state |n, l, m⟩. Edges are the two
        ladder families of the graph construction: angular (m ↔ m±1 within a
        subshell) and radial (n ↔ n±1 at fixed l, m). They are the graph&apos;s
        connectivity rules — not the radiative dipole selection rules of
        spectroscopy. Click a node to inspect it; drag to orbit, scroll to zoom.
      </p>
      <div className="controls">
        <label>
          max n
          <select value={maxN} onChange={(e) => { setMaxN(Number(e.target.value)); setSelected(null) }}>
            {MAX_N_CHOICES.map((n) => (
              <option key={n} value={n}>
                {n}
              </option>
            ))}
          </select>
        </label>
        <label>
          layout
          <select value={layout} onChange={(e) => setLayout(e.target.value as LayoutKind)}>
            <option value="paraboloid">energy paraboloid</option>
            <option value="shell">shells</option>
          </select>
        </label>
        <label>
          nuclear charge Z = {Z}
          <input
            type="range"
            min={1}
            max={10}
            value={Z}
            onChange={(e) => setZ(Number(e.target.value))}
          />
        </label>
      </div>
      <div className="legend">
        <span>
          <span className="swatch" style={{ background: 'var(--series-2)' }} />
          angular edge (m ↔ m±1)
        </span>
        <span>
          <span className="swatch" style={{ background: 'var(--series-1)' }} />
          radial edge (n ↔ n±1)
        </span>
        <span className="muted">node shade: binding magnitude (darker = deeper)</span>
      </div>
      <div className="grid-2" style={{ gridTemplateColumns: '3fr 2fr' }}>
        <div className="canvas-wrap">
          {latErr && <p className="muted">failed to load lattice data: {latErr}</p>}
          {lattice && (
            <LatticeScene
              data={lattice}
              layout={layout}
              dark={dark}
              selected={selected}
              onSelect={setSelected}
            />
          )}
        </div>
        <div>
          <div className="card">
            <h3>Selection</h3>
            {sel && lattice ? (
              <table className="data">
                <tbody>
                  <tr>
                    <td>state</td>
                    <td className="num">
                      |n={sel[0]}, l={sel[1]}, m={sel[2]}⟩
                    </td>
                  </tr>
                  <tr>
                    <td>
                      node weight ({lattice.node_weight_rule}, Z={Z})
                    </td>
                    <td className="num">
                      {(lattice.node_weights[selected!] * Z).toFixed(4)} Ha
                    </td>
                  </tr>
                  <tr>
                    <td>shell degeneracy (n²)</td>
                    <td className="num">
                      {lattice.shells.find((s) => s.n === sel[0])?.degeneracy}
                    </td>
                  </tr>
                </tbody>
              </table>
            ) : (
              <p className="muted small">
                click a node — {lattice ? `${fmt(lattice.num_states)} states, ${fmt(lattice.num_edges)} edges` : '…'}
              </p>
            )}
          </div>
          <div className="card">
            <h3>Spectrum</h3>
            {spectrum && <SpectrumPanel spectrum={spectrum} Z={Z} />}
          </div>
        </div>
      </div>
    </div>
  )
}
