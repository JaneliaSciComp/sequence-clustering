import numpy as np
import faiss

# ---------- Encoding (exact, 4-bit one-hot per base) ----------
_BASE2IDX = np.full(256, -1, dtype=np.int8)
for i, b in enumerate(b"ACGT"):
    _BASE2IDX[b] = i

def encode_onehot_4bit(seq: str) -> np.ndarray:
    """Return packed binary code (uint8 array) with 4 bits per nucleotide.
       Exactly 1 bit set per position; length = 4*L bits."""
    s = np.frombuffer(seq.strip().upper().encode("ascii"), dtype=np.uint8)
    idx = _BASE2IDX[s]
    if np.any(idx < 0):
        raise ValueError("Only A/C/G/T are allowed")
    L = len(idx)
    nbits = 4 * L
    nbytes = (nbits + 7) // 8
    code = np.zeros(nbytes, dtype=np.uint8)

    # For position p with base index k in {0,1,2,3}, set bit at bidx = 4*p + k
    bidx = 4 * np.arange(L, dtype=np.int64) + idx.astype(np.int64)

    byte_pos = bidx >> 3
    bit_pos  = bidx & 7
    code[byte_pos] |= (1 << bit_pos).astype(np.uint8)
    return code

# ---------- Exact Hamming check on nucleotides (fast, short-circuit) ----------
def hamming_leq2_nt(a: str, b: str) -> bool:
    if len(a) != len(b):
        return False
    diff = 0
    for ca, cb in zip(a, b):
        if ca != cb:
            diff += 1
            if diff > 2:
                return False
    return True

# ---------- Main dedup (order-dependent, exact ≤2 nts) ----------
class FaissUniqueRadius2:
    """
    Maintains one FAISS binary index per sequence length (dimension = 4*L bits).
    We only insert representatives. For a new sequence, we:
      1) encode to 4L-bit one-hot
      2) find nearest neighbor (k=1) in binary Hamming space
      3) if nearest distance <= 4, verify true nucleotide Hamming <= 2
         - if yes: SUPPRESS
         - if no: treat as new rep (hash collision defense)
      4) if no neighbor or distance > 4: add as new rep
    """
    def __init__(self):
        self._perL = {}   # L -> dict(index=faiss.IndexBinary*, reps=[str])

    def _get_index_for_L(self, L: int):
        entry = self._perL.get(L)
        if entry is not None:
            return entry
        d_bits = 4 * L
        # Try MIH; if unavailable, fall back to flat binary index
        try:
            index = faiss.IndexBinaryMultiHash(d_bits)  # deterministic MIH, if present
        except Exception:
            index = faiss.IndexBinaryFlat(d_bits)       # always available
        entry = {"index": index, "reps": []}
        self._perL[L] = entry
        return entry

    def process(self, seq: str) -> bool:
        """
        Process one sequence in stream order.
        Returns True if the sequence becomes a new representative, False if suppressed.
        """
        s = seq.strip().upper()
        L = len(s)
        if L == 0:
            return False
        entry = self._get_index_for_L(L)
        index = entry["index"]
        reps  = entry["reps"]

        code = encode_onehot_4bit(s)[None, :]  # shape (1, d_bits/8)

        if index.ntotal > 0:
            # k=1 nearest neighbor is sufficient: if the closest is >4, nothing is within 4.
            D, I = index.search(code, k=1)     # D: uint16 hamming distances
            d0 = int(D[0, 0])
            i0 = int(I[0, 0])
            if i0 != -1 and d0 <= 4:
                # Verify on nucleotides to guard against rare bit collisions
                if hamming_leq2_nt(s, reps[i0]):
                    return False  # SUPPRESS (within ≤2 of earlier rep)

        # New representative (first-wins)
        reps.append(s)
        index.add(code)
        return True

    def representatives(self):
        # Flatten in ascending L order for reproducibility, but original order is kept within each L
        out = []
        for L in sorted(self._perL):
            out.extend(self._perL[L]["reps"])
        return out

# ---------- Example usage ----------
if __name__ == "__main__":
    seqs = [
        "ACGTACGTAC",      # rep
        "ACGTACGTAT",      # 1 mismatch -> suppressed
        "ACGTACGTAA",      # 2 mismatches vs rep? (count and see)
        "ACGTACGCAC",      # 2 mismatches -> suppressed
        "TCGTACGTAC",      # 1 mismatch at pos0 -> suppressed
        "GGGGACGTAC",      # many mismatches -> new rep
    ]

    deduper = FaissUniqueRadius2()
    kept_flags = [deduper.process(s) for s in seqs]
    reps = deduper.representatives()

    print("Kept flags:", kept_flags)   # [True, False, ?, False, False, True] depending on exact distances
    print("Representatives:", reps)
