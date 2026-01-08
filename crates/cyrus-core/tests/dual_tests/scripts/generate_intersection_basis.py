#!/usr/bin/env python3
"""Generate intersection numbers in basis fixture."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from common import DUAL_POINTS, to_list, write_fixture
from cytools import Polytope

def main():
    p = Polytope(DUAL_POINTS)
    tri = p.triangulate()
    cy = tri.get_cy()

    basis = cy.divisor_basis()
    kappa = cy.intersection_numbers(in_basis=True)

    # Convert to sorted list format
    entries = []
    for (i, j, k), val in kappa.items():
        entries.append({
            "indices": sorted([i, j, k]),
            "value": int(val)
        })
    entries.sort(key=lambda x: x["indices"])

    fixture = {
        "description": "Intersection numbers in divisor basis from CYTools",
        "input_points": DUAL_POINTS,
        "basis": to_list(basis),
        "h11": int(cy.h11()),
        "n_entries": len(entries),
        "entries": entries,
    }

    write_fixture(fixture, Path(__file__).parent.parent / "fixtures" / "intersection_basis.json")

if __name__ == "__main__":
    main()
