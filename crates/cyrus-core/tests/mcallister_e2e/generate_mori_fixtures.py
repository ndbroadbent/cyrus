#!/usr/bin/env python3
"""
Generate Mori cone fixtures quickly (without full GV computation).
"""

import json
import numpy as np
from pathlib import Path

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

    print("Creating triangulation from heights...")
    t = p.triangulate(heights=heights)

    print("Getting CalabiYau...")
    cy = t.get_cy()
    print(f"h11 = {cy.h11()}, h21 = {cy.h21()}")

    # Get divisor basis
    print("Getting divisor basis...")
    basis = cy.divisor_basis()
    basis_list = basis.tolist() if hasattr(basis, 'tolist') else list(basis)
    print(f"Divisor basis: {basis_list[:10]}... (length {len(basis_list)})")

    # Get curve basis (GLSM charge matrix)
    print("Getting curve basis...")
    curve_basis = cy.curve_basis(include_origin=False, as_matrix=True)
    print(f"Curve basis shape: {curve_basis.shape}")

    # Get intersection numbers
    print("Getting intersection numbers (in_basis=True)...")
    kappa = cy.intersection_numbers(in_basis=True, format="dok")
    print(f"Number of non-zero intersection numbers: {len(kappa)}")

    # Get Mori cone cap
    print("Computing Mori cone cap (in_basis=True)...")
    mori = cy.mori_cone_cap(in_basis=True)
    mori_rays = mori.rays()
    print(f"Mori cone has {len(mori_rays)} rays")

    # Get grading vector
    print("Computing grading vector...")
    grading = mori.find_grading_vector()
    print(f"Grading vector shape: {grading.shape if hasattr(grading, 'shape') else len(grading)}")

    # Save fixtures
    print("\nSaving fixtures...")
    OUTPUT_DIR.mkdir(exist_ok=True)

    # Divisor basis
    basis_fixture = {
        "description": "Divisor basis indices for McAllister 4-214-647",
        "h11": cy.h11(),
        "h21": cy.h21(),
        "basis": basis_list
    }
    with open(OUTPUT_DIR / 'divisor_basis_cytools.json', 'w') as f:
        json.dump(basis_fixture, f, indent=2)
    print(f"Saved divisor_basis_cytools.json")

    # Curve basis
    curve_basis_fixture = {
        "description": "Curve basis (GLSM without origin) for GV computation",
        "shape": list(curve_basis.shape),
        "matrix": curve_basis.tolist()
    }
    with open(OUTPUT_DIR / 'curve_basis_cytools.json', 'w') as f:
        json.dump(curve_basis_fixture, f, indent=2)
    print(f"Saved curve_basis_cytools.json")

    # Mori cone generators (sample only - full is too large)
    mori_sample_fixture = {
        "description": "Mori cone cap sample (first 100 rays) for McAllister 4-214-647",
        "in_basis": True,
        "total_rays": len(mori_rays),
        "dim": len(mori_rays[0]) if len(mori_rays) > 0 else 0,
        "sample_rays": mori_rays[:100].tolist()
    }
    with open(OUTPUT_DIR / 'mori_cone_sample.json', 'w') as f:
        json.dump(mori_sample_fixture, f, indent=2)
    print(f"Saved mori_cone_sample.json (sample of {len(mori_rays)} rays)")

    # Grading vector
    grading_fixture = {
        "description": "Grading vector for Mori cone (interior point of dual)",
        "grading_vector": grading.tolist() if grading is not None else None,
        "dim": len(grading) if grading is not None else 0
    }
    with open(OUTPUT_DIR / 'grading_vector_cytools.json', 'w') as f:
        json.dump(grading_fixture, f, indent=2)
    print(f"Saved grading_vector_cytools.json")

    # Sample intersection numbers
    kappa_items = sorted(kappa.items())[:50]
    kappa_fixture = {
        "description": "Sample intersection numbers (first 50 sorted) for verification",
        "in_basis": True,
        "total_nonzero": len(kappa),
        "sample": [
            {"indices": list(k), "value": int(v)}
            for k, v in kappa_items
        ]
    }
    with open(OUTPUT_DIR / 'intersection_sample_cytools.json', 'w') as f:
        json.dump(kappa_fixture, f, indent=2)
    print(f"Saved intersection_sample_cytools.json")

    print("\nDone! Fixtures saved to:", OUTPUT_DIR)

if __name__ == '__main__':
    main()
