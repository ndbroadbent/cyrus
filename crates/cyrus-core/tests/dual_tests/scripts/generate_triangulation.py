#!/usr/bin/env python3
"""Generate triangulation fixture."""

import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))

from common import DUAL_POINTS, to_list, write_fixture
from cytools import Polytope

def compute_simplex_determinant(points, simplex_indices):
    """Compute |det| for a simplex (matches CYTools internal)."""
    pts = np.array([points[i] for i in simplex_indices])
    pts_aug = np.hstack([pts, np.ones((len(pts), 1))])
    return abs(int(round(np.linalg.det(pts_aug))))

def main():
    p = Polytope(DUAL_POINTS)
    tri = p.triangulate()

    simplices = to_list(tri.simplices())
    heights = to_list(tri.heights())

    # Compute determinants
    simplex_dets = []
    for simp in simplices:
        det = compute_simplex_determinant(DUAL_POINTS, simp)
        simplex_dets.append({"simplex": simp, "abs_det": det})

    fixture = {
        "description": "Triangulation from CYTools default (FRST)",
        "input_points": DUAL_POINTS,
        "simplices": simplices,
        "n_simplices": len(simplices),
        "heights": heights,
        "simplex_determinants": simplex_dets,
    }

    write_fixture(fixture, Path(__file__).parent.parent / "fixtures" / "triangulation.json")

if __name__ == "__main__":
    main()
