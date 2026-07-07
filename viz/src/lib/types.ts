// TypeScript mirrors of the geovac.viz_export JSON schemas (schema_version 1).

export interface Envelope {
  schema_version: number
  geovac_version: string
  generator: string
  kind: string
}

export interface Manifest extends Envelope {
  files: string[]
}

export type EdgeType = 'angular' | 'radial' | 'unknown'

export interface LatticeData extends Envelope {
  Z: number
  max_n: number
  states: [number, number, number][]
  edges: [number, number, EdgeType][]
  node_weights: number[]
  node_weight_rule: string
  edge_weighting: string
  shells: { n: number; degeneracy: number }[]
  num_states: number
  num_edges: number
}

export interface SpectrumData extends Envelope {
  Z: number
  max_n: number
  graph_laplacian_spectrum: number[]
  graph_energies_ha: number[]
  kinetic_scale: { value: number; z_scaling: string; note: string }
  s3_reference: {
    note: string
    levels: { n: number; eigenvalue: number; degeneracy: number }[]
  }
  hydrogen_levels: { n: number; energy_ha: number; degeneracy: number }[]
  series: {
    name: string
    n_lower: number
    lines: { n_upper: number; n_lower: number; delta_e_ha: number; wavelength_nm: number }[]
  }[]
}

export interface TaperingRow {
  mode: string
  n_qubits?: number
  n_terms?: number
  one_norm?: number
  metadata?: Record<string, unknown>
  error?: string
}

export interface QcSystemData extends Envelope {
  system: string
  n_qubits: number
  n_terms: number
  one_norm: number
  one_norm_full: number
  has_pk: boolean
  n_electrons: number | null
  n_orbitals: number | null
  ecore: number | null
  metadata: Record<string, unknown>
  measurement: {
    n_qwc_groups: number
    max_group_size: number
    min_group_size: number
    mean_group_size: number
  }
  tapering: TaperingRow[]
  pauli_terms: [string, number, number][]
}

export interface ScalingFit {
  exponent?: number
  prefactor?: number
  rms_log_residual?: number
  n_points: number
  note?: string
}

export interface LibraryData extends Envelope {
  n_systems: number
  systems: {
    system: string
    n_qubits: number
    n_terms: number
    one_norm: number
    n_qwc_groups: number
    tapering: TaperingRow[]
  }[]
  scaling: {
    cross_system: {
      label: string
      points: { system: string; Q: number; n_pauli: number }[]
      fit: ScalingFit
    }
    within_system?: {
      label: string
      system: string
      points: { max_n: number; M: number; Q: number; N_pauli: number }[]
      fit: ScalingFit
    }
    atomic_benchmark?: {
      label: string
      provenance: { source: string; generator: string; timestamp: string | null }
      points: { label: string; type: string; max_n: number | null; Q: number; n_pauli: number; one_norm: number | null }[]
      fit_he?: ScalingFit
    }
  }
}
