from collections import defaultdict

cdef class DisjointSetUnion:
    """
    Disjoint Set Union (Union-Find) data structure with path compression and union by rank.
    """
    cdef public int[:] parent
    cdef public int[:] rank
    cdef public int n
    
    def __init__(self, int n):
        self.n = n
        # Use memoryviews for better performance
        import array
        self.parent = array.array('i', range(n))
        self.rank = array.array('i', [0] * n)

    cdef int _find(self, int x):
        """Internal find method with path compression."""
        # Bounds check to prevent segfaults
        if x < 0 or x >= self.n:
            raise IndexError(f"Index {x} out of bounds for size {self.n}")

        cdef int original_x = x
        cdef int root
        
        # Find root
        while self.parent[x] != x:
            x = self.parent[x]
        root = x
        
        # Path compression
        x = original_x
        while self.parent[x] != root:
            original_x = self.parent[x]
            self.parent[x] = root
            x = original_x
        
        return root

    cdef bint _union(self, int a, int b):
        """Internal union method. Returns True if merged, False if already in the same set."""
        # Bounds check to prevent segfaults
        if a < 0 or a >= self.n or b < 0 or b >= self.n:
            raise IndexError(f"Index out of bounds: a={a}, b={b} for size {self.n}")

        cdef int ra = self._find(a)
        cdef int rb = self._find(b)
        
        if ra == rb:
            return False
            
        # Union by rank
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1
        return True

    def find(self, int x):
        """Find the root of the set containing x with path compression."""
        return self._find(x)
    
    def union(self, int a, int b):
        """Union the sets containing a and b. Returns True if merged, False if
        already in the same set."""
        return self._union(a, b)

    def update(self, edges, int offset_a = 0, int offset_b = 0):
        """Union all pairs in edges."""
        cdef int u, v
        for u, v in edges:
            self._union(u + offset_a, v + offset_b)

    def get_components(self):
        """Get all components as lists of nodes."""
        components = defaultdict(list)
        cdef int i, r
        for i in range(self.n):
            r = self._find(i)
            components[r].append(i)
        return list(components.values())

    @staticmethod
    def deduplicate_edges(edges):
        """Remove duplicate edges while preserving component connectivity."""
        if not edges:
            return []

        cdef dict node_to_index = {}
        cdef list unique_nodes = []
        next_index = 0

        # Assign compact indices to each node encountered in the edges
        for edge in edges:
            for node in edge:
                if node not in node_to_index:
                    node_to_index[node] = next_index
                    unique_nodes.append(node)
                    next_index += 1

        # Construct connected components using DSU
        dsu = DisjointSetUnion(next_index)
        for u, v in edges:
            dsu._union(node_to_index[u], node_to_index[v])
        components = dsu.get_components()

        # Reconstruct deduplicated edges from components
        deduplicated_edges = []
        for component in components:
            if len(component) < 2:
                continue
            start = unique_nodes[component[0]]
            for i in range(1, len(component)):
                end = unique_nodes[component[i]]
                deduplicated_edges.append((start, end))
                start = end

        return deduplicated_edges
