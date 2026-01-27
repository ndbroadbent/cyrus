#!/usr/bin/env bash
set -euo pipefail

TIMEOUT="${CYRUS_FLAT_DIRECTION_TIMEOUT:-5m}"
LOCK_FILE="${CYRUS_RUN_LOCK:-/tmp/cyrus_mcallister.lock}"

if ! command -v gtimeout >/dev/null 2>&1; then
  echo "gtimeout not found. Install coreutils or set a different timeout tool." >&2
  exit 1
fi

if pgrep -f 'mcallister_e2e|mcallister_racetrack|mcallister_gv|mcallister_first_principles|mcallister_flat_direction' >/dev/null 2>&1; then
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
  if pgrep -f 'target/.*/mcallister_(e2e|racetrack|gv|first_principles|flat_direction)' >/dev/null 2>&1; then
    pkill -f 'target/.*/mcallister_(e2e|racetrack|gv|first_principles|flat_direction)' >/dev/null 2>&1 || true
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
  cargo run -p cyrus-core --bin mcallister_flat_direction --release -- "${EXTRA_ARGS[@]}" "$@"
