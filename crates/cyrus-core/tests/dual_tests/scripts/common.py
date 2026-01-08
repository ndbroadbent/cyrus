"""Shared utilities for CYTools fixture generation."""

import json
import numpy as np

class NumpyEncoder(json.JSONEncoder):
    """JSON encoder that handles numpy types."""
    def default(self, obj):
        if isinstance(obj, (np.integer, np.int64, np.int32)):
            return int(obj)
        if isinstance(obj, (np.floating, np.float64, np.float32)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def to_list(obj):
    """Recursively convert numpy arrays and tuples to lists."""
    if isinstance(obj, np.ndarray):
        return [to_list(x) for x in obj.tolist()]
    if isinstance(obj, (list, tuple)):
        return [to_list(x) for x in obj]
    if isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    if isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    return obj


def write_fixture(data: dict, filename: str):
    """Write fixture to JSON file."""
    with open(filename, 'w') as f:
        json.dump(data, f, indent=2, cls=NumpyEncoder)
    print(f"Wrote {filename}")


# McAllister dual polytope - small test case (h11=4)
DUAL_POINTS = [
    [0, 0, 0, 0],     # 0: origin
    [-1, 2, -1, -1],  # 1
    [1, -1, 0, 0],    # 2
    [-1, -1, 1, 1],   # 3
    [-1, -1, 1, 2],   # 4
    [-1, -1, 2, 1],   # 5
    [-1, -1, 2, 3],   # 6
    [-1, -1, 3, 2],   # 7
    [-1, -1, 2, 2],   # 8
]
