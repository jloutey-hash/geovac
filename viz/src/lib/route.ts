// Minimal hash router: '#/', '#/lattice', '#/qc', '#/qc/<System>'.
import { useEffect, useState } from 'react'

export interface Route {
  page: 'home' | 'lattice' | 'qc' | 'qc-detail'
  param?: string
}

export function parseHash(hash: string): Route {
  const parts = hash.replace(/^#\/?/, '').split('/').filter(Boolean)
  if (parts.length === 0) return { page: 'home' }
  if (parts[0] === 'lattice') return { page: 'lattice' }
  if (parts[0] === 'qc') {
    if (parts.length > 1) return { page: 'qc-detail', param: decodeURIComponent(parts[1]) }
    return { page: 'qc' }
  }
  return { page: 'home' }
}

export function useRoute(): Route {
  const [route, setRoute] = useState<Route>(() => parseHash(window.location.hash))
  useEffect(() => {
    const onChange = () => setRoute(parseHash(window.location.hash))
    window.addEventListener('hashchange', onChange)
    return () => window.removeEventListener('hashchange', onChange)
  }, [])
  return route
}

export function useDarkMode(): boolean {
  const mq = window.matchMedia('(prefers-color-scheme: dark)')
  const [dark, setDark] = useState(mq.matches)
  useEffect(() => {
    const onChange = (e: MediaQueryListEvent) => setDark(e.matches)
    mq.addEventListener('change', onChange)
    return () => mq.removeEventListener('change', onChange)
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])
  return dark
}
