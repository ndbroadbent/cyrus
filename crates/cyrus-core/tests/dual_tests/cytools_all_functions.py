#!/usr/bin/env python3
"""
CYTools comprehensive extraction - ALL intermediate values for Rust comparison.

Extracts every function's output so we can verify our Rust implementations match exactly.
"""

import json
import numpy as np
from cytools import Polytope
from itertools import combinations

class NumpyEncoder(json.JSONEncoder):
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

def compute_simplex_determinant(points, simplex_indices):
    """Compute |det| for a simplex (matches CYTools internal)."""
    pts = np.array([points[i] for i in simplex_indices])
    # Augment with 1s column
    pts_aug = np.hstack([pts, np.ones((len(pts), 1))])
    return abs(int(round(np.linalg.det(pts_aug))))

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

    result = {
        "test_case": "mcallister_dual",
        "input_points": dual_points,
        "n_points": len(dual_points),
        "dim": 4,
    }

    # Create polytope
    p = Polytope(dual_points)

    # ==================== POLYTOPE FUNCTIONS ====================
    result["polytope"] = {
        "points": to_list(p.points()),
        "vertices": to_list(p.vertices()),
        "is_reflexive": p.is_reflexive(),
        "dimension": int(p.dim()),
        "points_not_interior_to_facets": to_list(p.points_not_interior_to_facets()),
    }

    # ==================== TRIANGULATION ====================
    # Default triangulation (FRST - placing heights)
    tri = p.triangulate()
    simplices = to_list(tri.simplices())

    # Compute determinants for each simplex
    simplex_dets = []
    for simp in simplices:
        det = compute_simplex_determinant(dual_points, simp)
        simplex_dets.append({"simplex": simp, "abs_det": det})

    # Get the heights used for this triangulation
    heights = to_list(tri.heights())

    result["triangulation"] = {
        "simplices": simplices,
        "n_simplices": len(simplices),
        "simplex_determinants": simplex_dets,
        "heights": heights,
    }

    # ==================== TRIANGULATION FROM CUSTOM HEIGHTS ====================
    # Test that we can reproduce the same triangulation from the heights
    tri_from_heights = p.triangulate(heights=heights)
    simplices_from_heights = to_list(tri_from_heights.simplices())

    result["triangulation_from_heights"] = {
        "input_heights": heights,
        "simplices": simplices_from_heights,
        "matches_original": simplices == simplices_from_heights,
    }

    # ==================== TORIC VARIETY / CALABI-YAU ====================
    cy = tri.get_cy()

    # GLSM charge matrix
    glsm = cy.glsm_charge_matrix()
    result["glsm_charge_matrix"] = {
        "matrix": to_list(glsm),
        "shape": list(glsm.shape),
    }

    # Linear relations (full, with origin)
    lr_full = cy.glsm_linear_relations(include_origin=True)
    result["linear_relations_with_origin"] = {
        "matrix": to_list(lr_full),
        "shape": list(lr_full.shape),
    }

    # Linear relations (without origin) - used for intersection computation
    lr_no_origin = cy.glsm_linear_relations(include_origin=False)
    result["linear_relations_no_origin"] = {
        "matrix": to_list(lr_no_origin),
        "shape": list(lr_no_origin.shape),
    }

    # Verify: L @ Q^T = 0
    product = lr_full @ glsm.T
    result["linear_relations_glsm_product"] = to_list(product)
    result["linear_relations_glsm_product_is_zero"] = np.allclose(product, 0)

    # ==================== DIVISOR BASIS ====================
    basis = cy.divisor_basis()
    result["divisor_basis"] = to_list(basis)
    result["h11"] = int(cy.h11())
    result["h21"] = int(cy.h21())

    # ==================== INTERSECTION NUMBERS ====================
    # Full intersection numbers (point-indexed)
    kappa_full = cy.intersection_numbers(in_basis=False)
    kappa_full_list = []
    for (i, j, k), val in kappa_full.items():
        indices = sorted([i, j, k])
        kappa_full_list.append({
            "indices": indices,
            "value": int(val)
        })
    kappa_full_list.sort(key=lambda x: x["indices"])
    result["intersection_numbers_full"] = kappa_full_list

    # Intersection numbers in basis
    kappa_basis = cy.intersection_numbers(in_basis=True)
    kappa_basis_list = []
    for (i, j, k), val in kappa_basis.items():
        indices = sorted([i, j, k])
        kappa_basis_list.append({
            "indices": indices,
            "value": int(val)
        })
    kappa_basis_list.sort(key=lambda x: x["indices"])
    result["intersection_numbers_basis"] = kappa_basis_list

    # ==================== DISTINCT INTERSECTION NUMBERS (4-form) ====================
    # These are the "known" values: κ^V_{ijkl} = 1/|det| for distinct non-origin indices
    # The determinant is computed from the 4x4 matrix of point coordinates
    pts = np.array(dual_points)

    distinct_4form = []
    for simp in simplices:
        # Remove origin from simplex
        non_origin = [i for i in simp if i != 0]
        if len(non_origin) >= 4:
            # All 4-combinations of non-origin indices
            for combo in combinations(non_origin, 4):
                combo_sorted = sorted(combo)
                # Compute determinant for these 4 points (4x4 matrix of coordinates)
                sub_pts = pts[list(combo_sorted)]
                det = abs(np.linalg.det(sub_pts))
                if det > 1e-10:
                    distinct_4form.append({
                        "indices": list(combo_sorted),
                        "det": det,
                        "kappa_4form": 1.0 / det
                    })

    # Remove duplicates
    seen = set()
    unique_4form = []
    for entry in distinct_4form:
        key = tuple(entry["indices"])
        if key not in seen:
            seen.add(key)
            unique_4form.append(entry)
    unique_4form.sort(key=lambda x: x["indices"])
    result["distinct_4form_values"] = unique_4form

    # ==================== 3-FACE PROBES ====================
    # All 3-tuples (with replacement) from simplices, excluding origin
    all_3tuples = set()
    for simp in simplices:
        non_origin = sorted([i for i in simp if i != 0])
        # Generate all 3-tuples with replacement (i <= j <= k)
        for i in non_origin:
            for j in non_origin:
                if j < i:
                    continue
                for k in non_origin:
                    if k < j:
                        continue
                    all_3tuples.add((i, j, k))

    result["probe_3faces"] = sorted([list(t) for t in all_3tuples])

    # ==================== SANITY CHECKS ====================
    result["sanity_checks"] = {
        "kappa_345": kappa_full.get((3, 4, 5), kappa_full.get((3, 5, 4), 0)),
        "kappa_111": kappa_full.get((1, 1, 1), 0),
        "kappa_333": kappa_full.get((3, 3, 3), 0),
        "kappa_123": kappa_full.get((1, 2, 3), 0),
        "kappa_888": kappa_full.get((8, 8, 8), 0),
    }

    # ==================== FLAT DIRECTION (Stage 4) ====================
    # Use test flux vectors
    K_flux = [-1, 1, 0, 1]   # Sample K flux
    M_flux = [1, 0, -1, 1]   # Sample M flux
    h11_dual = int(cy.h11())

    # Compute N matrix: N_ab = κ_abc M^c
    kappa_basis_dict = cy.intersection_numbers(in_basis=True)
    N_matrix = np.zeros((h11_dual, h11_dual))
    for a in range(h11_dual):
        for b in range(h11_dual):
            for c in range(h11_dual):
                key = tuple(sorted([a, b, c]))
                kappa_val = kappa_basis_dict.get(key, 0)
                N_matrix[a, b] += kappa_val * M_flux[c]

    result["flat_direction"] = {
        "K_flux": K_flux,
        "M_flux": M_flux,
        "N_matrix": to_list(N_matrix),
    }

    # Solve for p: N @ p = K
    try:
        p = np.linalg.solve(N_matrix, K_flux)
        result["flat_direction"]["p"] = to_list(p)
        result["flat_direction"]["solve_success"] = True

        # Compute κ_abc p^a p^b p^c
        kappa_ppp = 0.0
        for a in range(h11_dual):
            for b in range(h11_dual):
                for c in range(h11_dual):
                    key = tuple(sorted([a, b, c]))
                    kappa_val = kappa_basis_dict.get(key, 0)
                    kappa_ppp += kappa_val * p[a] * p[b] * p[c]
        result["flat_direction"]["kappa_ppp"] = kappa_ppp

        # e^K0 = (4/3) / κ_ppp
        if abs(kappa_ppp) > 1e-10:
            result["flat_direction"]["eK0"] = (4.0 / 3.0) / kappa_ppp
        else:
            result["flat_direction"]["eK0"] = None
    except np.linalg.LinAlgError:
        result["flat_direction"]["p"] = None
        result["flat_direction"]["solve_success"] = False
        result["flat_direction"]["kappa_ppp"] = None
        result["flat_direction"]["eK0"] = None

    print(json.dumps(result, indent=2, cls=NumpyEncoder))

if __name__ == "__main__":
    main()
