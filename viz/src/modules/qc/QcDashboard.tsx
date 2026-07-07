import { useMemo, useState } from 'react'
import { useData, fmt } from '../../lib/data'
import type { LibraryData } from '../../lib/types'
import ScalingPlot from './ScalingPlot'

type SortKey = 'system' | 'n_qubits' | 'n_terms' | 'one_norm' | 'n_qwc_groups' | 'q_full'

interface Row {
  system: string
  n_qubits: number
  n_terms: number
  one_norm: number
  n_qwc_groups: number
  q_full: number | null
}

export default function QcDashboard() {
  const { data: lib, error } = useData<LibraryData>('qc/library.json')
  const [sort, setSort] = useState<{ key: SortKey; dir: 1 | -1 }>({ key: 'system', dir: 1 })

  const rows = useMemo<Row[]>(() => {
    if (!lib) return []
    return lib.systems.map((s) => ({
      system: s.system,
      n_qubits: s.n_qubits,
      n_terms: s.n_terms,
      one_norm: s.one_norm,
      n_qwc_groups: s.n_qwc_groups,
      q_full: s.tapering.find((t) => t.mode === 'full' && t.n_qubits !== undefined)?.n_qubits ?? null,
    }))
  }, [lib])

  const sorted = useMemo(() => {
    const r = [...rows]
    r.sort((a, b) => {
      const va = a[sort.key]
      const vb = b[sort.key]
      if (va === null) return 1
      if (vb === null) return -1
      if (typeof va === 'string' && typeof vb === 'string') return sort.dir * va.localeCompare(vb)
      return sort.dir * ((va as number) - (vb as number))
    })
    return r
  }, [rows, sort])

  const clickSort = (key: SortKey) =>
    setSort((s) => (s.key === key ? { key, dir: s.dir === 1 ? -1 : 1 } : { key, dir: key === 'system' ? 1 : -1 }))

  const arrow = (key: SortKey) => (sort.key === key ? (sort.dir === 1 ? ' ▲' : ' ▼') : '')

  if (error) return <p className="muted">failed to load library data: {error}</p>
  if (!lib) return <p className="muted">loading…</p>

  return (
    <div>
      <h1>Qubit Hamiltonian resources</h1>
      <p className="secondary" style={{ maxWidth: 760 }}>
        Resource metrics for the {fmt(lib.n_systems)}-system molecular Hamiltonian
        library (Jordan–Wigner encoding). Click a system for Pauli-term detail,
        symmetry-tapering variants, and export code. These metrics characterize
        encoding size and measurement cost — they are not chemical-accuracy claims.
      </p>
      <div className="card">
        <h3>Pauli terms vs qubits</h3>
        <ScalingPlot lib={lib} />
      </div>
      <div className="card">
        <h3>System library</h3>
        <table className="data">
          <thead>
            <tr>
              <th onClick={() => clickSort('system')}>system{arrow('system')}</th>
              <th className="num" onClick={() => clickSort('n_qubits')}>qubits{arrow('n_qubits')}</th>
              <th className="num" onClick={() => clickSort('n_terms')}>Pauli terms{arrow('n_terms')}</th>
              <th className="num" onClick={() => clickSort('one_norm')}>1-norm (Ha){arrow('one_norm')}</th>
              <th className="num" onClick={() => clickSort('n_qwc_groups')}>QWC groups{arrow('n_qwc_groups')}</th>
              <th className="num" onClick={() => clickSort('q_full')}>qubits, full taper{arrow('q_full')}</th>
            </tr>
          </thead>
          <tbody>
            {sorted.map((r) => (
              <tr
                key={r.system}
                className="clickable"
                onClick={() => (window.location.hash = `#/qc/${encodeURIComponent(r.system)}`)}
              >
                <td>
                  <a href={`#/qc/${encodeURIComponent(r.system)}`}>{r.system}</a>
                </td>
                <td className="num">{fmt(r.n_qubits)}</td>
                <td className="num">{fmt(r.n_terms)}</td>
                <td className="num">{r.one_norm.toFixed(2)}</td>
                <td className="num">{fmt(r.n_qwc_groups)}</td>
                <td className="num">{r.q_full !== null ? fmt(r.q_full) : '—'}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  )
}
