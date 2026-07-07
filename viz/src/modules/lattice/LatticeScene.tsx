import { useEffect, useLayoutEffect, useMemo, useRef } from 'react'
import { Canvas, useFrame, useThree, type ThreeEvent } from '@react-three/fiber'
import * as THREE from 'three'
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js'
import type { LatticeData } from '../../lib/types'

export type LayoutKind = 'paraboloid' | 'shell'

// Node positions are pure functions of the quantum numbers — the layout is a
// rendering choice, not physics (the physics is the edge list from the data).
function positions(data: LatticeData, layout: LayoutKind): Float32Array {
  const pos = new Float32Array(data.states.length * 3)
  data.states.forEach(([n, l, m], i) => {
    let x: number, y: number, z: number
    if (layout === 'paraboloid') {
      // height = energy-like -1/n^2; ring radius by l; angle by m
      const r = l === 0 ? 0 : l * 1.7
      const phi = (2 * Math.PI * (m + l)) / (2 * l + 1)
      x = r * Math.cos(phi)
      z = r * Math.sin(phi)
      y = -6 / (n * n) + 3
    } else {
      // concentric shells of radius ~ n
      const R = n * 1.25
      const theta = (Math.PI * (l + 0.5)) / n
      const phi = (2 * Math.PI * (m + l + 0.5)) / (2 * l + 1)
      x = R * Math.sin(theta) * Math.cos(phi)
      y = R * Math.cos(theta)
      z = R * Math.sin(theta) * Math.sin(phi)
    }
    pos[3 * i] = x
    pos[3 * i + 1] = y
    pos[3 * i + 2] = z
  })
  return pos
}

// Sequential ramp (blue, light->dark) for node magnitude |weight| ~ 1/n^2.
const RAMP = ['#86b6ef', '#3987e5', '#0d366b'].map((c) => new THREE.Color(c))
function rampColor(t: number): THREE.Color {
  const x = Math.min(1, Math.max(0, t)) * (RAMP.length - 1)
  const i = Math.min(RAMP.length - 2, Math.floor(x))
  return RAMP[i].clone().lerp(RAMP[i + 1], x - i)
}

function Controls() {
  const { camera, gl } = useThree()
  const ref = useRef<OrbitControls | null>(null)
  useEffect(() => {
    const c = new OrbitControls(camera, gl.domElement)
    c.enableDamping = true
    ref.current = c
    return () => c.dispose()
  }, [camera, gl])
  useFrame(() => ref.current?.update())
  return null
}

function Nodes({
  data,
  pos,
  selected,
  onSelect,
}: {
  data: LatticeData
  pos: Float32Array
  selected: number | null
  onSelect: (i: number | null) => void
}) {
  const mesh = useRef<THREE.InstancedMesh>(null!)
  const count = data.states.length

  useLayoutEffect(() => {
    const dummy = new THREE.Object3D()
    // shade nodes by exported binding magnitude (|node_weight|), data-bound
    const maxW = Math.max(...data.node_weights.map(Math.abs))
    for (let i = 0; i < count; i++) {
      dummy.position.set(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2])
      dummy.updateMatrix()
      mesh.current.setMatrixAt(i, dummy.matrix)
      mesh.current.setColorAt(i, rampColor(Math.abs(data.node_weights[i]) / maxW))
    }
    mesh.current.instanceMatrix.needsUpdate = true
    if (mesh.current.instanceColor) mesh.current.instanceColor.needsUpdate = true
  }, [data, pos, count])

  const onClick = (e: ThreeEvent<MouseEvent>) => {
    e.stopPropagation()
    onSelect(e.instanceId ?? null)
  }

  return (
    <group>
      <instancedMesh
        ref={mesh}
        args={[undefined, undefined, count] as unknown as [THREE.BufferGeometry, THREE.Material, number]}
        onClick={onClick}
        onPointerOver={() => (document.body.style.cursor = 'pointer')}
        onPointerOut={() => (document.body.style.cursor = '')}
      >
        <sphereGeometry args={[0.15, 20, 20]} />
        <meshStandardMaterial roughness={0.45} metalness={0.05} />
      </instancedMesh>
      {selected !== null && selected < count && (
        <mesh position={[pos[3 * selected], pos[3 * selected + 1], pos[3 * selected + 2]]}>
          <sphereGeometry args={[0.24, 20, 20]} />
          <meshBasicMaterial color="#ffffff" wireframe />
        </mesh>
      )}
    </group>
  )
}

function Edges({
  data,
  pos,
  type,
  color,
}: {
  data: LatticeData
  pos: Float32Array
  type: 'angular' | 'radial'
  color: string
}) {
  const geom = useMemo(() => {
    const es = data.edges.filter((e) => e[2] === type)
    const arr = new Float32Array(es.length * 6)
    es.forEach(([i, j], k) => {
      arr.set(pos.slice(3 * i, 3 * i + 3), 6 * k)
      arr.set(pos.slice(3 * j, 3 * j + 3), 6 * k + 3)
    })
    const g = new THREE.BufferGeometry()
    g.setAttribute('position', new THREE.BufferAttribute(arr, 3))
    return g
  }, [data, pos, type])
  useEffect(() => () => geom.dispose(), [geom])
  return (
    <lineSegments geometry={geom}>
      <lineBasicMaterial color={color} transparent opacity={0.55} />
    </lineSegments>
  )
}

export default function LatticeScene({
  data,
  layout,
  dark,
  selected,
  onSelect,
}: {
  data: LatticeData
  layout: LayoutKind
  dark: boolean
  selected: number | null
  onSelect: (i: number | null) => void
}) {
  const pos = useMemo(() => positions(data, layout), [data, layout])
  const dist = layout === 'shell' ? data.max_n * 3.2 : data.max_n * 2.2 + 6
  // categorical slots 1 (radial) and 2 (angular), mode-stepped
  const radialColor = dark ? '#3987e5' : '#2a78d6'
  const angularColor = dark ? '#199e70' : '#1baf7a'

  return (
    <Canvas
      key={`${data.max_n}-${layout}`}
      camera={{ position: [dist, dist * 0.55, dist], fov: 45 }}
      onPointerMissed={() => onSelect(null)}
    >
      <ambientLight intensity={dark ? 0.7 : 0.9} />
      <directionalLight position={[10, 12, 8]} intensity={dark ? 0.9 : 1.1} />
      <Controls />
      <Edges data={data} pos={pos} type="angular" color={angularColor} />
      <Edges data={data} pos={pos} type="radial" color={radialColor} />
      <Nodes data={data} pos={pos} selected={selected} onSelect={onSelect} />
    </Canvas>
  )
}
