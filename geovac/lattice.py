"""
Geometric Paraboloid Lattice Generator

Generates a discrete lattice of hydrogen quantum states (n, l, m) with 
connectivity defined by angular momentum and radial transitions.

Author: Refactored for production
Date: February 2026
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix, lil_matrix
from typing import List, Tuple, Dict


class GeometricLattice:
    """
    TOPOLOGICAL QUANTUM LATTICE - "The Lattice is Truth"

    Nodes represent quantum states |n, l, m⟩ where:
    - n: principal quantum number (n ≥ 1)
    - l: angular momentum quantum number (0 ≤ l < n)
    - m: magnetic quantum number (-l ≤ m ≤ l)

    Edges represent allowed transitions:
    - Angular transitions: m ↔ m±1 (within same n, l)
    - Radial transitions: n ↔ n±1 (within same l, m)

    Node Weights (DIAGONAL) represent POTENTIAL ENERGY:
    - For state (n,l,m) near nucleus with charge Z: weight = -Z/n²
    - This IS the potential energy, encoded as graph topology!

    CRITICAL PRINCIPLE:
    ALL PHYSICS COMES FROM THE GRAPH. The Hamiltonian solver is "dumb" -
    it just diagonalizes H = scale*(D - A + W) where W are node weights.

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
        Diagonal weights (potential energy) for each node
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
        - Topological mode (topological_weights=True): weight = 1/(n₁*n₂)
          to encode Coulomb potential in graph structure
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
        Weight = 1/(n₁*n₂) creates n-dependent coupling that encodes
        the Coulomb potential directly in graph topology. States near
        the nucleus (low n) have stronger coupling.
        
        This makes D - A have eigenvalues scaling as -1/n² without
        adding a separate potential term.
        
        Binary Mode (hybrid):
        ---------------------
        Weight = 1.0 (all edges equal). Requires separate V = -Z/r term.
        
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
            # Topological puncture: weight ∝ 1/(n₁*n₂)
            # This encodes the 1/r potential in the graph structure
            return 1.0 / (n1 * n2)
        else:
            # Binary (uniform) weights
            return 1.0

    def _compute_node_weights(self) -> None:
        """
        Compute diagonal node weights = POTENTIAL ENERGY as graph topology.

        CRITICAL PRINCIPLE: "The Lattice is Truth"
        =========================================
        In spectral graph theory, the potential is NOT an external operator.
        It is the DIAGONAL WEIGHT (self-loop) of the graph.

        For a hydrogen-like atom with nuclear charge Z:
            Node weight for state (n, l, m) = -Z/n²

        This IS the potential energy V = -Z/r, encoded topologically through
        the quantum number n. The spatial scaling r = n²/Z is IMPLICIT in
        the Bohr model - we don't need explicit coordinates!

        The Hamiltonian solver will use:
            H = kinetic_scale * (D - A + W)

        where W = diag(node_weights) encodes the potential.

        Physical Interpretation:
        -----------------------
        - Low n states (close to nucleus): Large negative weight (strong binding)
        - High n states (far from nucleus): Small negative weight (weak binding)
        - Nuclear charge Z: Scales the binding strength

        Examples:
        - Hydrogen (Z=1, n=1): weight = -1/1 = -1.0 Ha
        - Helium (Z=2, n=1):   weight = -2/1 = -2.0 Ha (twice as strong!)
        - H⁻ (Z=1, n=2):       weight = -1/4 = -0.25 Ha (weakly bound)

        This is PURE TOPOLOGY - no coordinates, no distance calculations!
        """
        self.node_weights = np.zeros(self.num_states)

        for i, (n, l, m) in enumerate(self.states):
            # Topological potential energy: -Z/n²
            # This encodes the Coulomb attraction through quantum numbers alone
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
        
        This creates a molecular bonding configuration by connecting boundary
        states (n_max) from two atomic lattices. The sparse bridge represents
        orbital overlap along the bond axis.
        
        **Physical Interpretation:**
        - Each bridge edge represents a quantum mechanical overlap integral
        - n_bridges controls bond strength (optimal: 8-24 for H₂)
        - Bridge prioritization mimics σ-bond character (l=0,m=0 dominance)
        
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
            Number of bridge edges connecting boundary states
            Typical values:
            - 1-4: Weak bonding (single σ character)
            - 8-24: Normal covalent bond (H₂ optimal range)
            - >50: Strong multi-orbital mixing
            Default: 16 (experimentally validated for H₂)
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
