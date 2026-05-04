#!/usr/bin/env bash

resolve_cyrus_python() {
  local import_check='from ortools.sat.python import cp_model'

  if [[ -n "${CYRUS_PYTHON:-}" ]]; then
    if ! command -v "$CYRUS_PYTHON" >/dev/null 2>&1; then
      echo "CYRUS_PYTHON points to an executable that cannot be found: $CYRUS_PYTHON" >&2
      return 1
    fi
    if ! "$CYRUS_PYTHON" -c "$import_check" >/dev/null 2>&1; then
      echo "CYRUS_PYTHON cannot import OR-Tools: $CYRUS_PYTHON" >&2
      return 1
    fi
    export CYRUS_PYTHON
    return 0
  fi

  local script_dir repo_root code_root project_dir candidate
  script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  repo_root="$(cd "$script_dir/.." && pwd)"
  code_root="$(cd "$repo_root/.." && pwd)"
  project_dir="${CYRUS_MCALLISTER_PYPROJECT:-$code_root/string_theory/mcallister_2107}"

  for candidate in "$project_dir/.venv/bin/python3" "$project_dir/.venv/bin/python"; do
    if [[ -x "$candidate" ]] && "$candidate" -c "$import_check" >/dev/null 2>&1; then
      export CYRUS_PYTHON="$candidate"
      return 0
    fi
  done

  if [[ -d "$project_dir" ]] && command -v uv >/dev/null 2>&1; then
    if (cd "$project_dir" && uv run python -c "$import_check") >/dev/null 2>&1; then
      for candidate in "$project_dir/.venv/bin/python3" "$project_dir/.venv/bin/python"; do
        if [[ -x "$candidate" ]] && "$candidate" -c "$import_check" >/dev/null 2>&1; then
          export CYRUS_PYTHON="$candidate"
          return 0
        fi
      done
    fi
  fi

  echo "Unable to resolve a CYRUS_PYTHON interpreter with OR-Tools." >&2
  echo "Set CYRUS_PYTHON explicitly, or run uv sync in: $project_dir" >&2
  return 1
}
