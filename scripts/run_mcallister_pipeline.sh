#!/usr/bin/env bash
set -euo pipefail

TIMEOUT="${CYRUS_PIPELINE_TIMEOUT:-20m}"
CACHE_DIR="${CYRUS_CACHE_DIR:-crates/cyrus-core/target/cyrus-cache}"
LOCK_FILE="${CYRUS_RUN_LOCK:-/tmp/cyrus_mcallister.lock}"

if ! command -v gtimeout >/dev/null 2>&1; then
  echo "gtimeout not found. Install coreutils or set a different timeout tool." >&2
  exit 1
fi

if [[ -z "${CYRUS_PYTHON:-}" ]]; then
  echo "CYRUS_PYTHON is not set (needed for OR-Tools lattice points)." >&2
  exit 1
fi

export CYRUS_CACHE_DIR="$CACHE_DIR"
export CYRUS_LATTICE_MAX_TIME_SEC="${CYRUS_LATTICE_MAX_TIME_SEC:-60}"
export CYRUS_LATTICE_MAX_SOLUTIONS="${CYRUS_LATTICE_MAX_SOLUTIONS:-2000}"
export CYRUS_LATTICE_MAX_DEG="${CYRUS_LATTICE_MAX_DEG:-12}"
export CYRUS_LATTICE_THREADS="${CYRUS_LATTICE_THREADS:-1}"
export CYRUS_GV_THREADS="${CYRUS_GV_THREADS:-1}"
export CYRUS_GV_POOL_SIZE="${CYRUS_GV_POOL_SIZE:-1000}"

if pgrep -f 'mcallister_e2e|mcallister_racetrack|mcallister_gv|mcallister_first_principles' >/dev/null 2>&1; then
  echo "Another McAllister job is already running. Refusing to start in parallel." >&2
  exit 1
fi

LOCK_DIR="${LOCK_FILE}.d"
if ! mkdir "$LOCK_DIR" 2>/dev/null; then
  if [[ -f "$LOCK_DIR/pid" ]]; then
    old_pid="$(cat "$LOCK_DIR/pid" || true)"
    if [[ -n "$old_pid" ]] && ps -p "$old_pid" >/dev/null 2>&1; then
      echo "Run lock is held by PID $old_pid. Refusing to start." >&2
      exit 1
    fi
  fi
  echo "Stale lock detected. Removing and continuing." >&2
  rm -rf "$LOCK_DIR"
  mkdir "$LOCK_DIR"
fi

echo "$$" > "$LOCK_DIR/pid"

cleanup() {
  local code=$?
  rm -rf "$LOCK_DIR"
  if pgrep -f 'target/.*/mcallister_(e2e|racetrack|gv|first_principles)' >/dev/null 2>&1; then
    pkill -f 'target/.*/mcallister_(e2e|racetrack|gv|first_principles)' >/dev/null 2>&1 || true
  fi
  if pgrep -f 'find_lattice_points_ortools.py' >/dev/null 2>&1; then
    pkill -f 'find_lattice_points_ortools.py' >/dev/null 2>&1 || true
  fi
  exit "$code"
}
trap cleanup EXIT INT TERM

EXTRA_ARGS=()
if [[ -n "${CYRUS_MCALLISTER_DATA_DIR:-}" ]]; then
  EXTRA_ARGS+=(--data-dir "$CYRUS_MCALLISTER_DATA_DIR")
elif [[ "${CYRUS_ALLOW_FIXTURES:-}" == "1" ]]; then
  if [[ -n "${CYRUS_FIRST_PRINCIPLES:-}" ]]; then
    echo "CYRUS_ALLOW_FIXTURES cannot be used with CYRUS_FIRST_PRINCIPLES." >&2
    exit 1
  fi
  EXTRA_ARGS+=(--allow-fixtures)
else
  echo "CYRUS_MCALLISTER_DATA_DIR not set. Refusing to use fixtures." >&2
  echo "Set CYRUS_ALLOW_FIXTURES=1 to permit JSON fixture fallback." >&2
  exit 1
fi

exec gtimeout --foreground "$TIMEOUT" \
  --kill-after=10s \
  cargo run -p cyrus-core --bin mcallister_first_principles --release -- "${EXTRA_ARGS[@]}" "$@"
