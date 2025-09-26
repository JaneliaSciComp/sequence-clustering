from dataclasses import dataclass


@dataclass(slots=True)
class UniqueSequence:
    """Represents a unique sequence observation."""
    sequence: str
    count: int
