#!/usr/bin/env python3
"""Generate full intersection numbers fixture (point-indexed)."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from common import DUAL_POINTS, write_fixture
from cytools import Polytope

def main():
    p = Polytope(DUAL_POINTS)
    tri = p.triangulate()
    cy = tri.get_cy()

    kappa = cy.intersection_numbers(in_basis=False)

    # Convert to sorted list format
    entries = []
    for (i, j, k), val in kappa.items():
        entries.append({
            "indices": sorted([i, j, k]),
            "value": int(val)
        })
    entries.sort(key=lambda x: x["indices"])

    fixture = {
        "description": "Full intersection numbers (point-indexed) from CYTools",
        "input_points": DUAL_POINTS,
        "n_entries": len(entries),
        "entries": entries,
    }

    write_fixture(fixture, Path(__file__).parent.parent / "fixtures" / "intersection_full.json")

if __name__ == "__main__":
    main()
