#!/usr/bin/env python3
"""Generate flat direction fixture."""

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

    h11 = int(cy.h11())
    kappa = cy.intersection_numbers(in_basis=True)

    # Test flux vectors
    K_flux = [-1, 1, 0, 1]
    M_flux = [1, 0, -1, 1]

    # Compute N matrix: N_ab = kappa_abc M^c
    N_matrix = np.zeros((h11, h11))
    for a in range(h11):
        for b in range(h11):
            for c in range(h11):
                key = tuple(sorted([a, b, c]))
                kappa_val = kappa.get(key, 0)
                N_matrix[a, b] += kappa_val * M_flux[c]

    result = {
        "K_flux": K_flux,
        "M_flux": M_flux,
        "N_matrix": to_list(N_matrix),
    }

    # Solve N @ p = K
    try:
        p = np.linalg.solve(N_matrix, K_flux)
        result["p"] = to_list(p)
        result["solve_success"] = True

        # Compute kappa_abc p^a p^b p^c
        kappa_ppp = 0.0
        for a in range(h11):
            for b in range(h11):
                for c in range(h11):
                    key = tuple(sorted([a, b, c]))
                    kappa_val = kappa.get(key, 0)
                    kappa_ppp += kappa_val * p[a] * p[b] * p[c]
        result["kappa_ppp"] = kappa_ppp

        # e^K0 = 1 / ((4/3) * kappa_ppp).  Negative kappa_ppp is not a
        # physically valid flat direction for the exponential Kähler factor.
        if kappa_ppp > 1e-10:
            result["eK0"] = 1.0 / ((4.0 / 3.0) * kappa_ppp)
        else:
            result["eK0"] = None
    except np.linalg.LinAlgError:
        result["p"] = None
        result["solve_success"] = False
        result["kappa_ppp"] = None
        result["eK0"] = None

    fixture = {
        "description": "Flat direction computation from CYTools intersection numbers",
        "input_points": DUAL_POINTS,
        "h11": h11,
        "result": result,
    }

    write_fixture(fixture, Path(__file__).parent.parent / "fixtures" / "flat_direction.json")

if __name__ == "__main__":
    main()
