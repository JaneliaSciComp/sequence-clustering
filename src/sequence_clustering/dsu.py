from collections import defaultdict


class DisjointSetUnion:
    """
    Disjoint Set Union (Union-Find) data structure with path compression and union by rank.
    """
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x: int) -> int:
        """Find the root of the set containing x with path compression."""
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int):
        """Union the sets containing a and b. Returns True if merged, False if
        already in the same set."""
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1
        return True

    def update(self, edges: list[tuple[int, int]], offset_a: int = 0, offset_b: int = 0):
        """Union all pairs in edges."""
        for u, v in edges:
            self.union(u + offset_a, v + offset_b)

    def get_components(self) -> list[list[int]]:
        """Get all components as lists of nodes."""
        components = defaultdict(list)
        for i in range(len(self.parent)):
            r = self.find(i)
            components[r].append(i)
        return list(components.values())
