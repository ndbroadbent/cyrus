#!/usr/bin/env python3
import json
import sys
import time
from fractions import Fraction


def load_request():
    data = sys.stdin.read()
    if not data.strip():
        raise ValueError("No input received")
    return json.loads(data)


def solve_points(req):
    try:
        from ortools.sat.python import cp_model
    except Exception as exc:
        raise RuntimeError(f"ortools import failed: {exc}") from exc

    hyperplanes = req["hyperplanes"]
    grading = req["grading_vector"]
    min_points = req.get("min_points")
    max_deg = req.get("max_deg")
    max_coord = req.get("max_coord", 1000)
    deg_window = req.get("deg_window", 0)
    max_deg_limit = req.get("max_deg_limit")
    max_time_sec = req.get("max_time_sec")
    max_solutions = req.get("max_solutions")
    num_search_workers = req.get("num_search_workers")
    strict = req.get("strict", True)
    c = req.get("c", 0)
    check_grading = req.get("check_grading", True)

    if (min_points is None) == (max_deg is None):
        raise ValueError("Either min_points or max_deg must be set (exclusively)")

    if max_coord is None:
        max_coord = cp_model.INT32_MAX - 1

    n = len(grading)
    if n == 0:
        return []

    def normalize_c_list():
        if isinstance(c, list):
            return c
        return [c] * len(hyperplanes)

    def add_hyperplane_constraints(model, xs):
        c_list = normalize_c_list()
        for row, cc in zip(hyperplanes, c_list):
            cc_rat = Fraction(cc).limit_denominator()
            denom = cc_rat.denominator
            numer = cc_rat.numerator
            model.Add(sum(row[i] * xs[i] * denom for i in range(n)) >= numer)

    start_time = time.time()

    def remaining_time():
        if max_time_sec is None:
            return None
        return max_time_sec - (time.time() - start_time)

    def build_model(deg_low=None, deg_high=None):
        model = cp_model.CpModel()
        xs = [model.NewIntVar(-max_coord, max_coord, f"x_{i}") for i in range(n)]
        add_hyperplane_constraints(model, xs)
        deg_expr = sum(grading[i] * xs[i] for i in range(n))
        if deg_low is not None:
            model.Add(deg_low <= deg_expr)
        if deg_high is not None:
            model.Add(deg_expr <= deg_high)
        return model, xs, deg_expr

    counter = {"count": 0, "hit_limit": False}

    if check_grading:
        model = cp_model.CpModel()
        xs = [model.NewIntVar(-max_coord, max_coord, f"x_{i}") for i in range(n)]
        add_hyperplane_constraints(model, xs)
        model.Add(sum(grading[i] * xs[i] for i in range(n)) <= 0)

        solver = cp_model.CpSolver()
        rem = remaining_time()
        if rem is not None:
            if rem <= 0:
                raise RuntimeError("lattice-point grading check timed out")
            solver.parameters.max_time_in_seconds = rem
        if num_search_workers is not None:
            solver.parameters.num_search_workers = num_search_workers

        class SolutionCheck(cp_model.CpSolverSolutionCallback):
            def __init__(self, counter_ref):
                super().__init__()
                self.counter = counter_ref

            def on_solution_callback(self):
                self.counter["count"] += 1
                if self.counter["count"] > 1:
                    raise RuntimeError("grading vector must be wrong")

        check_counter = {"count": 0}
        solution = SolutionCheck(check_counter)
        try:
            solver.SearchForAllSolutions(model, solution)
        except RuntimeError as exc:
            raise RuntimeError(str(exc)) from exc

    if max_deg is not None:
        model, xs, _ = build_model(None, max_deg)
        solver = cp_model.CpSolver()
        rem = remaining_time()
        if rem is not None:
            if rem <= 0:
                raise RuntimeError("lattice-point search timed out")
            solver.parameters.max_time_in_seconds = rem
        if num_search_workers is not None:
            solver.parameters.num_search_workers = num_search_workers

        class SolutionStore(cp_model.CpSolverSolutionCallback):
            def __init__(self, variables, counter_ref):
                super().__init__()
                self.vars = variables
                self.counter = counter_ref

            def on_solution_callback(self):
                pt = [self.Value(v) for v in self.vars]
                sys.stdout.write(json.dumps(pt) + "\n")
                self.counter["count"] += 1
                if max_solutions is not None and self.counter["count"] >= max_solutions:
                    self.counter["hit_limit"] = True
                    self.StopSearch()

        solution = SolutionStore(xs, counter)
        status = solver.SearchForAllSolutions(model, solution)
        if counter["hit_limit"]:
            raise RuntimeError("max_solutions reached before completion")
        if strict and status != cp_model.OPTIMAL:
            raise RuntimeError(f"search failed with status: {solver.StatusName(status)}")
    else:
        deg = 0
        while counter["count"] < min_points:
            if max_deg_limit is not None and deg > max_deg_limit:
                raise RuntimeError("max_deg_limit reached before min_points")
            model, xs, _ = build_model(deg, deg + deg_window)
            solver = cp_model.CpSolver()
            rem = remaining_time()
            if rem is not None:
                if rem <= 0:
                    raise RuntimeError("lattice-point search timed out")
                solver.parameters.max_time_in_seconds = rem
            if num_search_workers is not None:
                solver.parameters.num_search_workers = num_search_workers

            class SolutionStore(cp_model.CpSolverSolutionCallback):
                def __init__(self, variables, counter_ref):
                    super().__init__()
                    self.vars = variables
                    self.counter = counter_ref

                def on_solution_callback(self):
                    pt = [self.Value(v) for v in self.vars]
                    sys.stdout.write(json.dumps(pt) + "\n")
                    self.counter["count"] += 1
                    if max_solutions is not None and self.counter["count"] >= max_solutions:
                        self.counter["hit_limit"] = True
                        self.StopSearch()

            solution = SolutionStore(xs, counter)
            status = solver.SearchForAllSolutions(model, solution)
            if counter["hit_limit"]:
                raise RuntimeError("max_solutions reached before min_points")
            if strict and status != cp_model.OPTIMAL:
                raise RuntimeError(
                    "search failed with status: "
                    f"{solver.StatusName(status)} at deg {deg}"
                )
            deg += deg_window + 1

    return None


def main():
    req = load_request()
    solve_points(req)


if __name__ == "__main__":
    main()
