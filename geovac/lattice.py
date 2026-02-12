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
    Clean implementation of the hydrogen quantum state lattice.
    
    Nodes represent quantum states |n, l, m⟩ where:
    - n: principal quantum number (n ≥ 1)
    - l: angular momentum quantum number (0 ≤ l < n)
    - m: magnetic quantum number (-l ≤ m ≤ l)
    
    Edges represent allowed transitions:
    - Angular transitions: m ↔ m±1 (within same n, l)
    - Radial transitions: n ↔ n±1 (within same l, m)
    
    Attributes:
    -----------
    max_n : int
        Maximum principal quantum number
    states : List[Tuple[int, int, int]]
        List of all (n, l, m) quantum states
    adjacency : scipy.sparse.csr_matrix
        Sparse adjacency matrix encoding connectivity
    """
    
    def __init__(self, max_n: int):
        """
        Initialize the geometric lattice.
        
        Parameters:
        -----------
        max_n : int
            Maximum principal quantum number (n ≥ 1)
        """
        if max_n < 1:
            raise ValueError("max_n must be at least 1")
        
        self.max_n = max_n
        self.states: List[Tuple[int, int, int]] = []
        self._state_index: Dict[Tuple[int, int, int], int] = {}
        self.adjacency: csr_matrix = None
        
        # Build lattice structure
        self._generate_states()
        self._build_adjacency_matrix()
    
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
        Uses binary weights (1 = connected, 0 = not connected).
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
                    adj[idx_from, idx_to] = 1.0
                    adj[idx_to, idx_from] = 1.0  # Symmetric
            
            # L- transition: m → m-1 (only count if not already connected)
            if m > -l:
                target = (n, l, m - 1)
                if target in self._state_index:
                    idx_to = self._state_index[target]
                    adj[idx_from, idx_to] = 1.0
                    adj[idx_to, idx_from] = 1.0  # Symmetric
            
            # === Radial Transitions (n ↔ n±1) ===
            
            # T+ transition: n → n+1
            if n < self.max_n:
                target = (n + 1, l, m)
                if target in self._state_index:
                    idx_to = self._state_index[target]
                    adj[idx_from, idx_to] = 1.0
                    adj[idx_to, idx_from] = 1.0  # Symmetric
            
            # T- transition: n → n-1 (only count if not already connected)
            if n > 1:
                target = (n - 1, l, m)
                if target in self._state_index:
                    idx_to = self._state_index[target]
                    adj[idx_from, idx_to] = 1.0
                    adj[idx_to, idx_from] = 1.0  # Symmetric
        
        # Convert to efficient CSR format for storage/computation
        self.adjacency = adj.tocsr()
    
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
