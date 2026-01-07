#!/usr/bin/env python3
"""
CYTools intersection number computation - exports all intermediate values for Rust comparison.

This script computes intersection numbers step-by-step and exports every intermediate
result as JSON so we can compare against our Rust implementation.

Usage:
    python cytools_intersection.py > intersection_fixture.json
"""

import json
import numpy as np
from cytools import Polytope

def to_json_serializable(obj):
    """Convert numpy arrays and other types to JSON-serializable format."""
    if isinstance(obj, np.ndarray):
        return to_json_serializable(obj.tolist())
    if isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    if isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    if isinstance(obj, dict):
        # Handle tuple keys by converting to string
        result = {}
        for k, v in obj.items():
            if isinstance(k, tuple):
                key = str(list(k))
            else:
                key = str(k)
            result[key] = to_json_serializable(v)
        return result
    if isinstance(obj, (list, tuple)):
        return [to_json_serializable(x) for x in obj]
    if isinstance(obj, int):
        return obj
    if isinstance(obj, float):
        return obj
    return obj

def compute_intersection_fixture():
    """Compute intersection numbers for the McAllister dual polytope."""

    # McAllister dual polytope points (from the paper)
    # This is the DUAL of polytope 4-214-647
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
        "input": {
            "points": dual_points,
            "n_points": len(dual_points),
            "dim": 4,
        }
    }

    # Create polytope
    p = Polytope(dual_points)

    # Get triangulation (use default/FRST)
    tri = p.triangulate()

    result["triangulation"] = {
        "simplices": to_json_serializable(tri.simplices()),
        "n_simplices": len(tri.simplices()),
    }

    # Get CalabiYau
    cy = tri.get_cy()

    # GLSM charge matrix
    glsm = cy.glsm_charge_matrix()
    result["glsm"] = {
        "matrix": to_json_serializable(glsm),
        "shape": list(glsm.shape),
    }

    # Linear relations (reduced form)
    # This is the key matrix for intersection computation
    lr = cy.glsm_linear_relations()
    result["linear_relations"] = {
        "matrix": to_json_serializable(lr),
        "shape": list(lr.shape),
    }

    # Divisor basis
    basis = cy.divisor_basis()
    result["divisor_basis"] = to_json_serializable(basis)

    # Full intersection numbers (point-indexed, before basis projection)
    # CYTools returns a dict with (i,j,k) tuple keys
    kappa_full = cy.intersection_numbers(in_basis=False)

    # Convert dict to list of (i,j,k,value) tuples
    kappa_entries = []
    for (i, j, k), val in kappa_full.items():
        # Ensure canonical order (i <= j <= k)
        indices = sorted([i, j, k])
        entry = (indices[0], indices[1], indices[2], int(val))
        if entry not in kappa_entries:
            kappa_entries.append(entry)

    kappa_entries.sort()
    result["kappa_full"] = kappa_entries

    # Intersection numbers in basis
    kappa_basis = cy.intersection_numbers(in_basis=True)

    kappa_basis_entries = []
    for (i, j, k), val in kappa_basis.items():
        indices = sorted([i, j, k])
        entry = (indices[0], indices[1], indices[2], int(val))
        if entry not in kappa_basis_entries:
            kappa_basis_entries.append(entry)

    kappa_basis_entries.sort()
    result["kappa_basis"] = {
        "entries": kappa_basis_entries,
        "dim": len(basis),
    }

    # Also get some specific values for quick sanity checks
    result["sanity_checks"] = {
        "kappa_345": kappa_full.get((3, 4, 5), kappa_full.get((3, 5, 4), 0)),
        "kappa_111": kappa_full.get((1, 1, 1), 0),
        "kappa_333": kappa_full.get((3, 3, 3), 0),
        "kappa_123": kappa_full.get((1, 2, 3), 0),
    }

    return result

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

def main():
    fixture = compute_intersection_fixture()
    print(json.dumps(fixture, indent=2, cls=NumpyEncoder))

if __name__ == "__main__":
    main()
