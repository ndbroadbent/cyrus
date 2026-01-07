#!/usr/bin/env python3
"""
Detailed CYTools intersection computation - exports ALL intermediate values.

This extracts internal computations to debug exactly where our Rust implementation diverges.
"""

import json
import numpy as np
from cytools import Polytope

class NumpyEncoder(json.JSONEncoder):
    """Custom JSON encoder for numpy types."""
    def default(self, obj):
        if isinstance(obj, (np.integer, np.int64, np.int32)):
            return int(obj)
        if isinstance(obj, (np.floating, np.float64, np.float32)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)

def to_list(obj):
    """Recursively convert to list."""
    if isinstance(obj, np.ndarray):
        return [to_list(x) for x in obj.tolist()]
    if isinstance(obj, (list, tuple)):
        return [to_list(x) for x in obj]
    if isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    if isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    return obj

def main():
    # McAllister dual polytope points
    dual_points = [
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

    result = {"input_points": dual_points}

    # Create polytope
    p = Polytope(dual_points)

    # Get points that CYTools actually uses
    result["polytope_points"] = to_list(p.points())
    result["polytope_points_not_interior_facets"] = to_list(p.points_not_interior_to_facets())

    # Triangulation
    tri = p.triangulate()
    result["simplices"] = to_list(tri.simplices())

    # Get CY
    cy = tri.get_cy()

    # GLSM charge matrix (full)
    glsm = cy.glsm_charge_matrix()
    result["glsm_charge_matrix"] = to_list(glsm)

    # Linear relations
    lr = cy.glsm_linear_relations()
    result["linear_relations"] = to_list(lr)

    # Points used for intersection (NOT interior to facets)
    pts_nif = p.points_not_interior_to_facets()
    result["points_not_interior_facets_count"] = len(pts_nif)

    # Divisor basis
    basis = cy.divisor_basis()
    result["divisor_basis"] = to_list(basis)
    result["h11"] = int(cy.h11())

    # Now get the internal intersection computation
    # This extracts the 4-forms and the linear system

    # Get intersection numbers with in_basis=False to see point-indexed values
    kappa = cy.intersection_numbers(in_basis=False)
    result["intersection_numbers_point_indexed"] = {
        str(list(k)): int(v) for k, v in kappa.items()
    }

    # Get intersection numbers with in_basis=True
    kappa_basis = cy.intersection_numbers(in_basis=True)
    result["intersection_numbers_basis_indexed"] = {
        str(list(k)): int(v) for k, v in kappa_basis.items()
    }

    # Try to extract more internal data
    # The _construct_intnum_equations_4d method builds the linear system
    # Let's trace through the algorithm manually

    # Get all 4-simplices (maximal simplices in 4D)
    simplices_arr = np.array(tri.simplices())
    result["n_simplices"] = len(simplices_arr)

    # For each simplex, compute the volume contribution
    simplex_volumes = []
    for simplex in simplices_arr:
        pts = [p.points()[i] for i in simplex]
        # Volume of simplex in 4D is |det([p1-p0, p2-p0, p3-p0, p4-p0])| / 4!
        mat = np.array([pts[i] - pts[0] for i in range(1, 5)])
        vol = abs(np.linalg.det(mat)) / 24
        simplex_volumes.append(vol)
    result["simplex_volumes"] = simplex_volumes
    result["total_volume"] = sum(simplex_volumes)

    # The 4-form summation approach
    # For each pair of 2-faces, compute contribution
    # This is complex - let's just verify the key intermediate values

    # Also check the fan structure
    result["dimension"] = int(p.dim())
    result["is_reflexive"] = p.is_reflexive()

    print(json.dumps(result, indent=2, cls=NumpyEncoder))

if __name__ == "__main__":
    main()
