#!/usr/bin/env python3
"""Generate divisor basis fixture."""

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

    fixture = {
        "description": "Divisor basis from CYTools",
        "input_points": DUAL_POINTS,
        "basis": to_list(basis),
        "h11": int(cy.h11()),
        "h21": int(cy.h21()),
    }

    write_fixture(fixture, Path(__file__).parent.parent / "fixtures" / "divisor_basis.json")

if __name__ == "__main__":
    main()
