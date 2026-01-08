#!/usr/bin/env python3
"""
Generate GV invariant test fixtures from CYTools for dual testing.

This script computes:
1. Mori cone cap generators
2. Grading vector
3. GV invariants

The outputs are saved as JSON for comparison with Rust implementations.
"""

import json
import numpy as np
from pathlib import Path

# Use installed CYTools
from cytools import Polytope

DATA_DIR = Path('/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647')
OUTPUT_DIR = Path(__file__).parent / 'assertions'

def load_points():
    """Load polytope points from McAllister's data."""
    with open(DATA_DIR / 'points.dat') as f:
        lines = f.readlines()
    points = []
    for line in lines:
        if line.strip():
            coords = [int(x) for x in line.strip().split(',')]
            points.append(coords)
    return points

def load_heights():
    """Load triangulation heights from McAllister's data."""
    with open(DATA_DIR / 'heights.dat') as f:
        content = f.read()
    return [float(x) for x in content.strip().split(',')]

def main():
    print("Loading McAllister polytope data...")
    points = load_points()
    heights = load_heights()

    print(f"Creating polytope with {len(points)} points...")
    p = Polytope(points)

    print("Getting triangulation points (not interior to facets)...")
    tri_points = p.points_not_interior_to_facets()
    print(f"Triangulation uses {len(tri_points)} points")

    print("Creating triangulation from heights...")
    t = p.triangulate(heights=heights)

    print("Getting CalabiYau...")
    cy = t.get_cy()

    print(f"h11 = {cy.h11()}, h21 = {cy.h21()}")

    # Get divisor basis
    print("Getting divisor basis...")
    basis = cy.divisor_basis()
    print(f"Divisor basis shape: {basis.shape if hasattr(basis, 'shape') else len(basis)}")

    # Get curve basis (GLSM charge matrix)
    print("Getting curve basis...")
    curve_basis = cy.curve_basis(include_origin=False, as_matrix=True)
    print(f"Curve basis shape: {curve_basis.shape}")

    # Get intersection numbers
    print("Getting intersection numbers...")
    kappa = cy.intersection_numbers(in_basis=True, format="dok")
    print(f"Number of non-zero intersection numbers: {len(kappa)}")

    # Get Mori cone cap
    print("Computing Mori cone cap...")
    mori = cy.mori_cone_cap(in_basis=True)
    mori_rays = mori.rays()
    print(f"Mori cone has {len(mori_rays)} rays")

    # Get grading vector
    print("Computing grading vector...")
    grading = mori.find_grading_vector()
    print(f"Grading vector: {grading}")

    # Compute GV invariants
    print("Computing GV invariants (this may take a while)...")
    gvs = cy.compute_gvs(min_points=100)

    # Extract GV data
    gv_charges = list(gvs._charge2invariant.keys())
    gv_values = list(gvs._charge2invariant.values())
    print(f"Computed {len(gv_charges)} GV invariants")

    # Save fixtures
    print("\nSaving fixtures...")

    # Mori cone generators
    mori_fixture = {
        "description": "Mori cone cap generators for McAllister 4-214-647",
        "in_basis": True,
        "rays": mori_rays.tolist(),
        "n_rays": len(mori_rays),
        "dim": len(mori_rays[0]) if len(mori_rays) > 0 else 0
    }
    with open(OUTPUT_DIR / 'mori_cone_cap.json', 'w') as f:
        json.dump(mori_fixture, f, indent=2)
    print(f"Saved mori_cone_cap.json")

    # Grading vector
    grading_fixture = {
        "description": "Grading vector for Mori cone (interior point of dual)",
        "grading_vector": grading.tolist() if grading is not None else None,
        "dim": len(grading) if grading is not None else 0
    }
    with open(OUTPUT_DIR / 'grading_vector.json', 'w') as f:
        json.dump(grading_fixture, f, indent=2)
    print(f"Saved grading_vector.json")

    # GV invariants
    gv_fixture = {
        "description": "Gopakumar-Vafa invariants from CYTools",
        "n_invariants": len(gv_charges),
        "charges": [list(c) for c in gv_charges],
        "values": [int(v) for v in gv_values],
        "first_10": [
            {"charge": list(gv_charges[i]), "value": int(gv_values[i])}
            for i in range(min(10, len(gv_charges)))
        ]
    }
    with open(OUTPUT_DIR / 'gv_invariants_cytools.json', 'w') as f:
        json.dump(gv_fixture, f, indent=2)
    print(f"Saved gv_invariants_cytools.json")

    # Curve basis for reference
    curve_basis_fixture = {
        "description": "Curve basis (GLSM without origin) for GV computation",
        "shape": list(curve_basis.shape),
        "matrix": curve_basis.tolist()
    }
    with open(OUTPUT_DIR / 'curve_basis.json', 'w') as f:
        json.dump(curve_basis_fixture, f, indent=2)
    print(f"Saved curve_basis.json")

    # Summary of intersection numbers (first few for verification)
    kappa_items = list(kappa.items())[:20]
    kappa_fixture = {
        "description": "Sample intersection numbers (first 20) for verification",
        "in_basis": True,
        "sample": [
            {"indices": list(k), "value": int(v)}
            for k, v in kappa_items
        ]
    }
    with open(OUTPUT_DIR / 'intersection_sample.json', 'w') as f:
        json.dump(kappa_fixture, f, indent=2)
    print(f"Saved intersection_sample.json")

    print("\nDone! Fixtures saved to:", OUTPUT_DIR)

if __name__ == '__main__':
    main()
