// Data loading. Single-source-of-truth rule (visualization_plan.md §2):
// every physics number rendered by the site comes from these artifacts;
// nothing is hardcoded in components.

const cache = new Map<string, unknown>()

export async function loadJson<T>(rel: string): Promise<T> {
  if (cache.has(rel)) return cache.get(rel) as T
  const url = `${import.meta.env.BASE_URL}data/${rel}`
  const res = await fetch(url)
  if (!res.ok) throw new Error(`fetch ${url}: HTTP ${res.status}`)
  const obj = (await res.json()) as T
  cache.set(rel, obj)
  return obj
}

import { useEffect, useState } from 'react'

export interface Loading<T> {
  data: T | null
  error: string | null
}

export function useData<T>(rel: string | null): Loading<T> {
  const [state, setState] = useState<Loading<T>>({ data: null, error: null })
  useEffect(() => {
    if (rel === null) return
    let alive = true
    setState({ data: null, error: null })
    loadJson<T>(rel)
      .then((data) => alive && setState({ data, error: null }))
      .catch((e) => alive && setState({ data: null, error: String(e) }))
    return () => {
      alive = false
    }
  }, [rel])
  return state
}

export function fmt(x: number, digits = 2): string {
  if (Number.isInteger(x) && Math.abs(x) < 1e15) return x.toLocaleString('en-US')
  return x.toLocaleString('en-US', {
    maximumFractionDigits: digits,
    minimumFractionDigits: 0,
  })
}

export function fmtSci(x: number, digits = 3): string {
  if (x === 0) return '0'
  const a = Math.abs(x)
  if (a >= 0.01 && a < 10000) return x.toPrecision(digits)
  return x.toExponential(digits - 1)
}
