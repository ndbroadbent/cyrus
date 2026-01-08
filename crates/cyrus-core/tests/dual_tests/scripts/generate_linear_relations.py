#!/usr/bin/env python3
"""Generate GLSM linear relations fixture."""

import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))

from common import DUAL_POINTS, to_list, write_fixture
from cytools import Polytope

def main():
    p = Polytope(DUAL_POINTS)
    tri = p.triangulate()
    cy = tri.get_cy()

    # Linear relations with origin (used for basis transformation)
    lr_with_origin = cy.glsm_linear_relations(include_origin=True)

    # Linear relations without origin (used for intersection computation)
    lr_no_origin = cy.glsm_linear_relations(include_origin=False)

    # Verify: L @ Q^T = 0
    glsm = cy.glsm_charge_matrix()
    product = lr_with_origin @ glsm.T

    fixture = {
        "description": "GLSM linear relations from CYTools",
        "input_points": DUAL_POINTS,
        "with_origin": {
            "matrix": to_list(lr_with_origin),
            "shape": list(lr_with_origin.shape),
        },
        "without_origin": {
            "matrix": to_list(lr_no_origin),
            "shape": list(lr_no_origin.shape),
        },
        "verification": {
            "L_times_Q_transpose": to_list(product),
            "is_zero": bool(np.allclose(product, 0)),
        },
    }

    write_fixture(fixture, Path(__file__).parent.parent / "fixtures" / "linear_relations.json")

if __name__ == "__main__":
    main()
