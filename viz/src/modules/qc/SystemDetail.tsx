import { useMemo, useState } from 'react'
import { useData, fmt, fmtSci } from '../../lib/data'
import type { QcSystemData } from '../../lib/types'
import PauliHeatmap from './PauliHeatmap'

const REPO_URL = 'https://github.com/jloutey-hash/geovac'
const PAGE = 100

function CopyButton({ text }: { text: string }) {
  const [done, setDone] = useState(false)
  return (
    <button
      className="copy"
      onClick={() => {
        void navigator.clipboard.writeText(text)
        setDone(true)
        setTimeout(() => setDone(false), 1200)
      }}
    >
      {done ? 'copied' : 'copy'}
    </button>
  )
}

export default function SystemDetail({ name }: { name: string }) {
  const { data, error } = useData<QcSystemData>(`qc/systems/${name}.json`)
  const [search, setSearch] = useState('')
  const [page, setPage] = useState(0)

  const filtered = useMemo(() => {
    if (!data) return []
    const q = search.trim().toUpperCase()
    if (!q) return data.pauli_terms
    return data.pauli_terms.filter(([s]) => s.toUpperCase().includes(q))
  }, [data, search])

  if (error) return <p className="muted">failed to load {name}: {error}</p>
  if (!data) return <p className="muted">loading {name}…</p>

  const pages = Math.max(1, Math.ceil(filtered.length / PAGE))
  const pageRows = filtered.slice(page * PAGE, (page + 1) * PAGE)

  const snippet = [
    `# clone + editable install of the research package`,
    `#   git clone ${REPO_URL}`,
    `#   pip install -e geovac openfermion`,
    `from geovac.ecosystem_export import hamiltonian`,
    ``,
    `H = hamiltonian('${data.system}')  # ${data.n_qubits} qubits, ${data.n_terms} Pauli terms`,
    `of_op = H.to_openfermion()`,
    `qk_op = H.to_qiskit()      # requires qiskit`,
    `pl_op = H.to_pennylane()   # requires pennylane`,
  ].join('\n')

  const taperRow = data.tapering.find((t) => t.mode === 'per_block' && t.n_qubits !== undefined)
  const taperSnippet = taperRow
    ? `H_t = hamiltonian('${data.system}', tapered='per_block')  # ${taperRow.n_qubits} qubits`
    : null

  return (
    <div>
      <p className="small">
        <a href="#/qc">← library</a>
      </p>
      <h1>{data.system}</h1>
      <div className="stat-row">
        <div className="stat-tile">
          <div className="value">{fmt(data.n_qubits)}</div>
          <div className="label">qubits</div>
        </div>
        <div className="stat-tile">
          <div className="value">{fmt(data.n_terms)}</div>
          <div className="label">Pauli terms</div>
        </div>
        <div className="stat-tile">
          <div className="value">{data.one_norm.toFixed(2)}</div>
          <div className="label">1-norm (Ha)</div>
        </div>
        <div className="stat-tile">
          <div className="value">{fmt(data.measurement.n_qwc_groups)}</div>
          <div className="label">QWC measurement groups</div>
        </div>
      </div>
      <p className="secondary small">
        {data.n_electrons !== null ? `${data.n_electrons} active electrons · ` : ''}
        {data.n_orbitals !== null ? `${data.n_orbitals} spatial orbitals · ` : ''}
        mean QWC group size {data.measurement.mean_group_size.toFixed(1)}
        {data.has_pk
          ? ' · Phillips–Kleinman pseudopotential separated for classical evaluation (electronic-only Hamiltonian shown)'
          : ''}
      </p>

      <div className="grid-2">
        <div className="card">
          <h3>Symmetry tapering</h3>
          <table className="data">
            <thead>
              <tr>
                <th>mode</th>
                <th className="num">qubits</th>
                <th className="num">Δ qubits</th>
                <th className="num">Pauli terms</th>
                <th className="num">1-norm (Ha)</th>
              </tr>
            </thead>
            <tbody>
              <tr>
                <td>none</td>
                <td className="num">{fmt(data.n_qubits)}</td>
                <td className="num">—</td>
                <td className="num">{fmt(data.n_terms)}</td>
                <td className="num">{data.one_norm.toFixed(2)}</td>
              </tr>
              {data.tapering.map((t) => (
                <tr key={t.mode}>
                  <td>{t.mode}</td>
                  {t.error ? (
                    <td colSpan={4} className="muted small">
                      not available: {t.error}
                    </td>
                  ) : (
                    <>
                      <td className="num">{fmt(t.n_qubits!)}</td>
                      <td className="num">−{fmt(data.n_qubits - t.n_qubits!)}</td>
                      <td className="num">{fmt(t.n_terms!)}</td>
                      <td className="num">{t.one_norm!.toFixed(2)}</td>
                    </>
                  )}
                </tr>
              ))}
            </tbody>
          </table>
          <p className="small muted">
            Z₂ symmetry tapering maps the Hamiltonian onto its fully symmetric
            sector; each mode uses a different stabilizer set.
          </p>
        </div>
        <div className="card">
          <h3>Qubit coupling map</h3>
          <PauliHeatmap terms={data.pauli_terms} Q={data.n_qubits} />
        </div>
      </div>

      <div className="card">
        <h3>
          Export <CopyButton text={snippet + (taperSnippet ? '\n' + taperSnippet : '')} />
        </h3>
        <pre className="snippet">{snippet}</pre>
        {taperSnippet && <pre className="snippet">{taperSnippet}</pre>}
        <p className="small muted">
          The <code>geovac-hamiltonians</code> standalone packaging of the same
          builder lives in the{' '}
          <a href={REPO_URL} target="_blank" rel="noreferrer">
            project repository
          </a>
          .
        </p>
      </div>

      <div className="card">
        <h3>Pauli terms ({fmt(data.n_terms)})</h3>
        <div className="controls">
          <label>
            filter
            <input
              type="text"
              placeholder="e.g. X0 or Z12"
              value={search}
              onChange={(e) => {
                setSearch(e.target.value)
                setPage(0)
              }}
            />
          </label>
          <span className="muted small">{fmt(filtered.length)} matching · sorted by |coefficient|</span>
        </div>
        <table className="data">
          <thead>
            <tr>
              <th className="num">#</th>
              <th>Pauli string</th>
              <th className="num">coefficient (Ha)</th>
              <th className="num">|coefficient|</th>
            </tr>
          </thead>
          <tbody>
            {pageRows.map(([s, re, im], i) => (
              <tr key={page * PAGE + i}>
                <td className="num muted">{page * PAGE + i + 1}</td>
                <td>
                  <code>{s}</code>
                </td>
                <td className="num">
                  {fmtSci(re, 5)}
                  {Math.abs(im) > 1e-12 ? ` + ${fmtSci(im, 3)}i` : ''}
                </td>
                <td className="num">{fmtSci(Math.hypot(re, im), 5)}</td>
              </tr>
            ))}
          </tbody>
        </table>
        {pages > 1 && (
          <div className="pagination">
            <button disabled={page === 0} onClick={() => setPage((p) => p - 1)}>
              ← prev
            </button>
            <span>
              page {page + 1} / {pages}
            </span>
            <button disabled={page >= pages - 1} onClick={() => setPage((p) => p + 1)}>
              next →
            </button>
          </div>
        )}
      </div>
    </div>
  )
}
