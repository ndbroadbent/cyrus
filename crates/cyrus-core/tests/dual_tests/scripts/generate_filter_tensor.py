#!/usr/bin/env python3
"""Generate filter_tensor_indices fixture - tests the utility function directly."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from common import write_fixture
from cytools.utils import filter_tensor_indices

def main():
    # Test case from CYTools docstring
    tensor = {(0,1): 0, (1,1): 1, (1,2): 2, (1,3): 3, (2,3): 4}
    indices = [1, 3]

    result = filter_tensor_indices(tensor, indices)

    # Convert tuple keys to lists for JSON
    result_list = [{"indices": list(k), "value": v} for k, v in sorted(result.items())]

    fixture = {
        "description": "filter_tensor_indices utility function test",
        "test_cases": [
            {
                "input_tensor": [{"indices": list(k), "value": v} for k, v in sorted(tensor.items())],
                "indices": indices,
                "expected": result_list,
            }
        ],
    }

    # Add more test cases
    tensor2 = {(0,0,0): 1, (0,0,1): 2, (0,1,1): 3, (1,1,1): 4, (0,1,2): 5, (1,2,2): 6}
    indices2 = [0, 2]
    result2 = filter_tensor_indices(tensor2, indices2)
    fixture["test_cases"].append({
        "input_tensor": [{"indices": list(k), "value": v} for k, v in sorted(tensor2.items())],
        "indices": indices2,
        "expected": [{"indices": list(k), "value": v} for k, v in sorted(result2.items())],
    })

    write_fixture(fixture, Path(__file__).parent.parent / "fixtures" / "filter_tensor.json")

if __name__ == "__main__":
    main()
