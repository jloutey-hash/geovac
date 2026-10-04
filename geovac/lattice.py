"""
Geometric Paraboloid Lattice Generator

Generates the discrete (n, l, m) label lattice with edges that change m
(at fixed n, l) or n (at fixed l, m); no edge changes l.

Author: Refactored for production
Date: February 2026
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix, lil_matrix
from typing import List, Tuple, Dict


class GeometricLattice:
    """
    QUANTUM LATTICE on the (n, l, m) Gel'fand-Tsetlin labels.

    Nodes represent quantum states |n, l, m⟩ where:
    - n: principal quantum number (n ≥ 1)
    - l: angular momentum quantum number (0 ≤ l < n)
    - m: magnetic quantum number (-l ≤ m ≤ l)

    Edges (graph connectivity, not physical transitions; Paper 0):
    - Angular transitions: m ↔ m±1 (within same n, l)
    - Radial transitions: n ↔ n±1 (within same l, m)

    Node weights (diagonal attribute, not used by AtomicSolver):
    - For state (n,l,m) with nuclear charge Z: weight = -Z/n²
    - An n-dependent label weight; it is not the hydrogenic level energy
      -Z²/(2n²) and does not encode the Coulomb potential.

    SCOPE (corrected 2026-09-25, PI-directed re-pricing; see Papers 0/7
    and CLAUDE.md S1). This edge set is what AtomicSolver's legacy
    hamiltonian='graph' mode diagonalizes as H = kappa*Z^2*(D - A) (the node
    weights W are not part of it). No
    edge here ever changes l (angular edges fix (n, l), radial edges fix
    (l, m)), so D - A is the direct sum of max_n disconnected l-blocks
    and does NOT converge, under any rescaling, to the S^3
    Laplace-Beltrami operator (tests/test_graph_s3_nonconvergence.py):
    that operator's heat kernel is everywhere positive, while this
    graph's heat kernel is exactly zero between l-blocks. The continuum
    S^3 Laplace-Beltrami operator instead acts on these same (n, l, m)
    labels as the diagonal SO(4) Casimir L^2 + A^2 = n^2 - 1 -- an exact
    statement about the label set and its degeneracies, not a property
    of this graph's D - A spectrum. The matching scale kappa = -1/16
    maps this graph's extremal eigenvalue (band edge 8 = 2*d_max) to the
    hydrogen ground state -1/2 by construction; that match is an
    Observation, not a derivation, and does not extend to the
    excited-state (Rydberg) ladder -- the 'graph' Hamiltonian's excited
    states are non-physical. For the physical spectrum use
    AtomicSolver(hamiltonian='exact'), i.e. diag(-Z^2/(2 n^2)).

    Attributes:
    -----------
    max_n : int
        Maximum principal quantum number
    nucleus_position : np.ndarray
        3D position of the nucleus this lattice represents
    nuclear_charge : int
        Nuclear charge Z (proton number)
    states : List[Tuple[int, int, int]]
        List of all (n, l, m) quantum states
    adjacency : scipy.sparse.csr_matrix
        Sparse adjacency matrix encoding connectivity
    node_weights : np.ndarray
        Diagonal label weights -Z/n^2 for each node (not a potential)
    """

    def __init__(self,
                 max_n: int,
                 nucleus_position: Tuple[float, float, float] = (0.0, 0.0, 0.0),
                 nuclear_charge: int = 1,
                 topological_weights: bool = False):
        """
        Initialize the geometric lattice centered on a nucleus.

        Parameters:
        -----------
        max_n : int
            Maximum principal quantum number (n ≥ 1)
        nucleus_position : Tuple[float, float, float], optional
            3D coordinates of the nucleus (default: origin)
        nuclear_charge : int, optional
            Nuclear charge Z (default: 1 for hydrogen)
        topological_weights : bool, optional
            If True, use n-dependent edge weights (default: False)
        """
        if max_n < 1:
            raise ValueError("max_n must be at least 1")
        if nuclear_charge < 1:
            raise ValueError("nuclear_charge must be at least 1")

        self.max_n = max_n
        self.nucleus_position = np.array(nucleus_position, dtype=float)
        self.nuclear_charge = nuclear_charge
        self.topological_weights = topological_weights

        self.states: List[Tuple[int, int, int]] = []
        self._state_index: Dict[Tuple[int, int, int], int] = {}
        self.adjacency: csr_matrix = None
        self.node_weights: np.ndarray = None

        # Build lattice structure
        self._generate_states()
        self._build_adjacency_matrix()
        self._compute_node_weights()
    
    def _generate_states(self) -> None:
        """
        Generate all valid quantum states (n, l, m).
        
        CRITICAL: Loop bounds must match original physics:
        - n ranges from 1 to max_n (inclusive)
        - l ranges from 0 to n-1 (inclusive)
        - m ranges from -l to l (inclusive)
        """
        idx = 0
        for n in range(1, self.max_n + 1):
            for l in range(n):  # l = 0, 1, ..., n-1
                for m in range(-l, l + 1):  # m = -l, -l+1, ..., 0, ..., l
                    state = (n, l, m)
                    self.states.append(state)
                    self._state_index[state] = idx
                    idx += 1
    
    def _build_adjacency_matrix(self) -> None:
        """
        Build sparse adjacency matrix encoding all allowed transitions.
        
        CRITICAL: Connectivity rules from original physics:
        
        1. Angular momentum raising (L+): (n,l,m) → (n,l,m+1) if m < l
        2. Angular momentum lowering (L-): (n,l,m) → (n,l,m-1) if m > -l
        3. Radial raising (T+): (n,l,m) → (n+1,l,m) if n < max_n
        4. Radial lowering (T-): (n,l,m) → (n-1,l,m) if n > 1
        
        Adjacency matrix is symmetric (undirected graph).
        
        Edge Weighting:
        - Binary mode (topological_weights=False): weight = 1.0
        - Topological mode (topological_weights=True): weight = 1/(n₁*n₂),
          an n-dependent weighting (it does not encode the Coulomb potential)
        """
        n_states = len(self.states)
        
        # Use lil_matrix for efficient construction
        adj = lil_matrix((n_states, n_states), dtype=np.float64)
        
        for n, l, m in self.states:
            idx_from = self._state_index[(n, l, m)]
            
            # === Angular Momentum Transitions (m ↔ m±1) ===
            
            # L+ transition: m → m+1
            if m < l:
                target = (n, l, m + 1)
                if target in self._state_index:
                    idx_to = self._state_index[target]
                    weight = self._compute_edge_weight(n, n)  # Same n
                    adj[idx_from, idx_to] = weight
                    adj[idx_to, idx_from] = weight  # Symmetric
            
            # L- transition: m → m-1 (only count if not already connected)
            if m > -l:
                target = (n, l, m - 1)
                if target in self._state_index:
                    idx_to = self._state_index[target]
                    weight = self._compute_edge_weight(n, n)  # Same n
                    adj[idx_from, idx_to] = weight
                    adj[idx_to, idx_from] = weight  # Symmetric
            
            # === Radial Transitions (n ↔ n±1) ===
            
            # T+ transition: n → n+1
            if n < self.max_n:
                target = (n + 1, l, m)
                if target in self._state_index:
                    idx_to = self._state_index[target]
                    weight = self._compute_edge_weight(n, n + 1)
                    adj[idx_from, idx_to] = weight
                    adj[idx_to, idx_from] = weight  # Symmetric
            
            # T- transition: n → n-1 (only count if not already connected)
            if n > 1:
                target = (n - 1, l, m)
                if target in self._state_index:
                    idx_to = self._state_index[target]
                    weight = self._compute_edge_weight(n, n - 1)
                    adj[idx_from, idx_to] = weight
                    adj[idx_to, idx_from] = weight  # Symmetric
        
        # Convert to efficient CSR format for storage/computation
        self.adjacency = adj.tocsr()
    
    def _compute_edge_weight(self, n1: int, n2: int) -> float:
        """
        Compute edge weight based on quantum numbers.
        
        Topological Mode (pure geometric):
        ----------------------------------
        Weight = 1/(n₁*n₂) creates an n-dependent coupling, stronger at
        low n.

        This does NOT make D - A converge to a -1/n² (Rydberg-like)
        continuum spectrum: the edge set never changes l, so D - A
        splits into max_n disconnected l-blocks with spectrum algebraic
        in [0, 8] regardless of weighting mode
        (tests/test_graph_s3_nonconvergence.py). This is a topological
        edge-weighting scheme, not a derivation of the hydrogen
        spectrum from the graph.
        
        Binary Mode (hybrid):
        ---------------------
        Weight = 1.0 (all edges equal).
        
        Parameters:
        -----------
        n1, n2 : int
            Principal quantum numbers of connected states
        
        Returns:
        --------
        weight : float
            Edge weight for adjacency matrix
        """
        if self.topological_weights:
            # n-dependent weighting: weight ∝ 1/(n₁*n₂)
            return 1.0 / (n1 * n2)
        else:
            # Binary (uniform) weights
            return 1.0

    def _compute_node_weights(self) -> None:
        """
        Compute the diagonal node-weight attribute W, -Z/n² per state.

        W is an n-dependent label weight stored on the lattice.  It is not
        the hydrogenic level energy -Z²/(2n²), it is not the potential
        V = -Z/r, and AtomicSolver's legacy graph mode does not use it
        (that mode diagonalizes kappa*Z^2*(D - A)).  The D - A part does not
        converge to a continuum kinetic operator
        (tests/test_graph_s3_nonconvergence.py).
        """
        self.node_weights = np.zeros(self.num_states)

        for i, (n, l, m) in enumerate(self.states):
            # n-dependent label weight -Z/n²
            self.node_weights[i] = -self.nuclear_charge / (n**2)

    @property
    def num_states(self) -> int:
        """Return the total number of quantum states in the lattice."""
        return len(self.states)
    
    @property
    def num_edges(self) -> int:
        """Return the total number of edges (connections) in the lattice."""
        # Each edge counted once (upper triangle of symmetric matrix)
        return self.adjacency.nnz // 2
    
    def sparsity(self) -> float:
        """
        Compute the sparsity of the adjacency matrix.
        
        Returns:
        --------
        sparsity : float
            Fraction of zero entries: 1 - (nnz / n²)
        """
        n = self.num_states
        return 1.0 - (self.adjacency.nnz / (n * n))
    
    def get_neighbors(self, state: Tuple[int, int, int]) -> List[Tuple[int, int, int]]:
        """
        Get all neighboring states connected to the given state.
        
        Parameters:
        -----------
        state : Tuple[int, int, int]
            Quantum state (n, l, m)
        
        Returns:
        --------
        neighbors : List[Tuple[int, int, int]]
            List of connected quantum states
        """
        if state not in self._state_index:
            raise ValueError(f"State {state} not in lattice")
        
        idx = self._state_index[state]
        # Get non-zero column indices in row idx
        neighbor_indices = self.adjacency[idx].nonzero()[1]
        
        return [self.states[i] for i in neighbor_indices]
    
    def compute_graph_distance(self, state1: Tuple[int, int, int], 
                               state2: Tuple[int, int, int],
                               method: str = 'geodesic') -> float:
        """
        Compute graph distance between two quantum states.
        
        CRITICAL: This is the GEOMETRIC distance on the lattice,
        NOT Euclidean distance in 3D space.
        
        Methods:
        --------
        'geodesic': Shortest path length (sum of edge weights)
                    Uses Dijkstra's algorithm for weighted graphs
        'hops': Number of edges in shortest path (unweighted)
                Uses BFS for unweighted graphs
        
        Parameters:
        -----------
        state1, state2 : Tuple[int, int, int]
            Quantum states (n, l, m)
        method : str
            Distance metric: 'geodesic' or 'hops'
        
        Returns:
        --------
        distance : float
            Graph distance between states
            - 'geodesic': sum of edge weights along shortest path
            - 'hops': number of edges in shortest path
        """
        if state1 not in self._state_index or state2 not in self._state_index:
            raise ValueError(f"States {state1} or {state2} not in lattice")
        
        idx1 = self._state_index[state1]
        idx2 = self._state_index[state2]
        
        # Same node → zero distance
        if idx1 == idx2:
            return 0.0
        
        if method == 'geodesic':
            # Dijkstra's algorithm for weighted shortest path
            return self._dijkstra_distance(idx1, idx2)
        elif method == 'hops':
            # BFS for unweighted shortest path
            return self._bfs_distance(idx1, idx2)
        else:
            raise ValueError(f"Unknown method: {method}")
    
    def _dijkstra_distance(self, start_idx: int, end_idx: int) -> float:
        """
        Compute shortest path distance using Dijkstra's algorithm.
        
        For geometric interaction, edge weights represent "cost" to traverse.
        In topological mode, smaller n has larger weight (stronger coupling),
        but for distance we want INVERSE: distance = 1/weight.
        
        This way, paths through low-n states are "shorter" (geometrically closer).
        """
        import heapq
        
        n_states = self.num_states
        distances = np.full(n_states, np.inf)
        distances[start_idx] = 0.0
        visited = set()
        
        # Priority queue: (distance, node_idx)
        pq = [(0.0, start_idx)]
        
        while pq:
            current_dist, current_idx = heapq.heappop(pq)
            
            if current_idx in visited:
                continue
            
            if current_idx == end_idx:
                return current_dist
            
            visited.add(current_idx)
            
            # Get neighbors and edge weights
            row = self.adjacency.getrow(current_idx)
            for neighbor_idx, edge_weight in zip(row.indices, row.data):
                if neighbor_idx not in visited:
                    # For DISTANCE, use inverse of coupling strength
                    # Strong coupling (large weight) → small distance
                    # Weak coupling (small weight) → large distance
                    edge_distance = 1.0 / edge_weight if edge_weight > 0 else np.inf
                    
                    new_dist = current_dist + edge_distance
                    
                    if new_dist < distances[neighbor_idx]:
                        distances[neighbor_idx] = new_dist
                        heapq.heappush(pq, (new_dist, neighbor_idx))
        
        # No path found
        return np.inf
    
    def _bfs_distance(self, start_idx: int, end_idx: int) -> float:
        """
        Compute shortest path in terms of edge count (BFS).
        
        Ignores edge weights, just counts number of hops.
        """
        from collections import deque
        
        visited = set([start_idx])
        queue = deque([(start_idx, 0)])
        
        while queue:
            current_idx, dist = queue.popleft()
            
            if current_idx == end_idx:
                return float(dist)
            
            # Get neighbors
            row = self.adjacency.getrow(current_idx)
            for neighbor_idx in row.indices:
                if neighbor_idx not in visited:
                    visited.add(neighbor_idx)
                    queue.append((neighbor_idx, dist + 1))
        
        # No path found
        return np.inf
    
    def compute_all_pair_distances(self, method: str = 'geodesic') -> np.ndarray:
        """
        Compute graph distances between all pairs of states.
        
        Returns:
        --------
        distances : np.ndarray, shape (n_states, n_states)
            Distance matrix where distances[i, j] is graph distance
            from state i to state j
        """
        n = self.num_states
        distances = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i+1, n):
                dist = self._dijkstra_distance(i, j) if method == 'geodesic' else self._bfs_distance(i, j)
                distances[i, j] = dist
                distances[j, i] = dist  # Symmetric
        
        return distances
    
    def stitch_lattices(self, other_lattice: 'GeometricLattice', 
                       n_bridges: int = 16, bridge_weight: float = 1.0) -> Tuple[csr_matrix, int, int]:
        """
        Stitch two lattices together with sparse topological bridges.

        SUPERSEDED (legacy LCAO graph concatenation): kept for the archived
        v0.x demos. The FCI-M guardrail shows graph concatenation has an
        R-independent kinetic term and no equilibrium; molecular work uses
        natural geometry (Papers 11, 12, 17). Bridge counts below are model
        parameters, not validated physics.

        Connects boundary states (n_max) from two atomic lattices; bridge
        prioritization favours (l=0, m=0) states.

        **Connection Strategy:**
        States are prioritized by overlap probability:
        1. Primary: (n_max, l=0, m=0) states (σ-bond core)
        2. Secondary: (n_max, l=1, m=0) states (π contributions)
        3. Tertiary: Higher l and non-zero m states
        
        Parameters:
        -----------
        other_lattice : GeometricLattice
            Second atomic lattice to bond with
        n_bridges : int, optional
            Number of bridge edges connecting boundary states (default 16)
        bridge_weight : float, optional
            Weight for bridge edges (default: 1.0)
            Controls coupling strength across bond
        
        Returns:
        --------
        stitched_adjacency : scipy.sparse.csr_matrix
            Combined adjacency matrix (n₁ + n₂) × (n₁ + n₂)
        n_bridges_actual : int
            Actual number of bridge edges created
        n_total_states : int
            Total number of states in stitched system
        
        Example:
        --------
        >>> atom_A = GeometricLattice(max_n=5)
        >>> atom_B = GeometricLattice(max_n=5)
        >>> adj_H2, n_bridges, n_states = atom_A.stitch_lattices(atom_B, n_bridges=16)
        >>> print(f"H₂ molecule: {n_states} states, {n_bridges} bonds")
        """
        # Get boundary states from both lattices (sorted by priority)
        boundary_A = self._get_boundary_states_prioritized()
        boundary_B = other_lattice._get_boundary_states_prioritized()
        
        n_A = self.num_states
        n_B = other_lattice.num_states
        n_total = n_A + n_B
        
        # Build combined adjacency matrix
        A_stitched = lil_matrix((n_total, n_total))
        
        # Copy lattice A (top-left block)
        adj_A = self.adjacency.tocoo()
        for i, j, w in zip(adj_A.row, adj_A.col, adj_A.data):
            A_stitched[i, j] = w
        
        # Copy lattice B (bottom-right block, offset by n_A)
        adj_B = other_lattice.adjacency.tocoo()
        for i, j, w in zip(adj_B.row, adj_B.col, adj_B.data):
            A_stitched[n_A + i, n_A + j] = w
        
        # Add sparse bridge connections
        n_bridges_actual = 0
        max_possible = min(len(boundary_A), len(boundary_B), n_bridges)
        
        for i in range(max_possible):
            idx_A = boundary_A[i]
            idx_B = boundary_B[i]
            
            # Symmetric connection: A ↔ B
            A_stitched[idx_A, n_A + idx_B] = bridge_weight
            A_stitched[n_A + idx_B, idx_A] = bridge_weight
            n_bridges_actual += 1
        
        return A_stitched.tocsr(), n_bridges_actual, n_total
    
    def _get_boundary_states_prioritized(self) -> List[int]:
        """
        Get boundary state indices (n=n_max) sorted by bonding priority.
        
        Priority ranking for molecular bonding:
        1. (n_max, 0, 0) - s orbital, maximum overlap
        2. (n_max, 1, 0) - p_z orbital, σ-bond axis
        3. (n_max, 2, 0) - d_z² orbital
        4. Higher l with m=0
        5. Non-zero m states (lower priority)
        
        Returns:
        --------
        prioritized_indices : List[int]
            State indices sorted by bonding priority (highest first)
        """
        # Find all n_max states
        boundary_states = [(idx, state) for idx, state in enumerate(self.states) 
                          if state[0] == self.max_n]
        
        # Sort by priority: (l, |m|) - lower is better
        def priority_key(item):
            idx, (n, l, m) = item
            return (l, abs(m))
        
        boundary_states.sort(key=priority_key)
        
        # Return indices only
        return [idx for idx, state in boundary_states]
    
    def __repr__(self) -> str:
        return (f"GeometricLattice(max_n={self.max_n}, "
                f"states={self.num_states}, "
                f"edges={self.num_edges}, "
                f"sparsity={self.sparsity():.4f})")


if __name__ == "__main__":
    """
    Demonstration: Initialize lattice for n=5 and report statistics.
    """
    print("=" * 60)
    print("Geometric Paraboloid Lattice - Production Test")
    print("=" * 60)
    
    # Initialize lattice with max_n = 5
    lattice = GeometricLattice(max_n=5)
    
    # Report statistics
    print(f"\nLattice Properties:")
    print(f"  Max n:        {lattice.max_n}")
    print(f"  Num nodes:    {lattice.num_states}")
    print(f"  Num edges:    {lattice.num_edges}")
    print(f"  Sparsity:     {lattice.sparsity():.6f}")
    
    # Matrix statistics
    n = lattice.num_states
    print(f"\nMatrix Properties:")
    print(f"  Shape:        {n} × {n}")
    print(f"  Nonzero:      {lattice.adjacency.nnz}")
    print(f"  Density:      {1 - lattice.sparsity():.6f}")
    print(f"  Memory (MB):  {lattice.adjacency.data.nbytes / 1e6:.3f}")
    
    # Show example state connectivity
    print(f"\nExample Connectivity:")
    example_state = (2, 1, 0)  # n=2, l=1, m=0
    neighbors = lattice.get_neighbors(example_state)
    print(f"  State {example_state} connects to:")
    for neighbor in neighbors:
        print(f"    → {neighbor}")
    
    # Verify state count formula: sum_{n=1}^{max_n} n^2
    expected_states = sum(n**2 for n in range(1, lattice.max_n + 1))
    print(f"\nVerification:")
    print(f"  Expected states: {expected_states}")
    print(f"  Actual states:   {lattice.num_states}")
    print(f"  Match:           {expected_states == lattice.num_states}")
    
    print("\n" + "=" * 60)
    print("✓ Lattice construction complete")
    print("=" * 60)
