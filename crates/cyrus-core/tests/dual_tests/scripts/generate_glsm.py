#!/usr/bin/env python3
"""Generate GLSM charge matrix fixture."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from common import DUAL_POINTS, to_list, write_fixture
from cytools import Polytope

def main():
    p = Polytope(DUAL_POINTS)
    tri = p.triangulate()
    cy = tri.get_cy()

    glsm = cy.glsm_charge_matrix()

    fixture = {
        "description": "GLSM charge matrix from CYTools",
        "input_points": DUAL_POINTS,
        "matrix": to_list(glsm),
        "shape": list(glsm.shape),
        "h11": int(cy.h11()),
    }

    write_fixture(fixture, Path(__file__).parent.parent / "fixtures" / "glsm.json")

if __name__ == "__main__":
    main()
