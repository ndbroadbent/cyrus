#!/bin/bash
# Check file length limits: 500 for code, 1000 for test files
# Applies to Rust source files
#
# Usage:
#   scripts/check_file_length.sh [files...]
#   If no files provided, checks all tracked files

set -e

# Exclusion patterns (paths or glob patterns to skip)
EXCLUSIONS=(
  "*.md"
  "*.json"
  "*.yml"
  "*.yaml"
  "*.toml"
  "*.txt"
  "*.lock"
  "Cargo.lock"
  "*.generated.*"
  "*/generated/*"
  # Binary file extensions
  "*.jpg"
  "*.jpeg"
  "*.png"
  "*.gif"
  "*.ico"
  "*.svg"
  "*.pdf"
  "*.zip"
  "*.tar"
  "*.gz"
)

exitcode=0

# Function to check if file should be excluded
should_exclude() {
  local file="$1"
  for pattern in "${EXCLUSIONS[@]}"; do
    # Remove leading ./ from file path for matching
    local clean_file="${file#./}"
    # shellcheck disable=SC2053
    if [[ "$clean_file" == $pattern ]]; then  # We want glob matching here
      return 0
    fi
  done
  return 1
}

# Legacy files grandfathered above the limits, capped near their size when the
# exception was added so they can only shrink. Do NOT add new entries — split
# new code into modules instead. TODO: refactor these below the normal limits.
legacy_limit() {
  case "${1#./}" in
    crates/cyrus-core/src/gv.rs)                                       echo 17450 ;;
    crates/cyrus-core/src/bin/mcallister_first_principles/main.rs)     echo 4700 ;;  # TODO: keep splitting into modules
    crates/cyrus-core/src/kklt.rs)                                     echo 4400 ;;
    crates/cyrus-core/src/triangulation/secondary.rs)                  echo 4150 ;;
    crates/cyrus-core/src/cone/mod.rs)                                 echo 1450 ;;
    crates/cyrus-core/src/intersection/complete_intersection.rs)       echo 1000 ;;
    crates/cyrus-core/src/bin/mcallister_flat_direction.rs)            echo 950 ;;
    crates/cyrus-core/src/curve_basis.rs)                              echo 950 ;;
    crates/cyrus-core/src/intersection/cytools_algorithm/ambient.rs)   echo 600 ;;
    crates/cyrus-core/src/cone/tip.rs)                                 echo 520 ;;
    crates/cyrus-core/tests/potent_ray_affine_circuits.rs)             echo 1550 ;;
    crates/cyrus-core/tests/mcallister_e2e/stage5_gv.rs)               echo 1350 ;;
    crates/cyrus-core/tests/mcallister_e2e/stage0_data_integrity.rs)   echo 1200 ;;
    *)                                                                 echo "" ;;
  esac
}

# Check a single file
check_file() {
  local file="$1"

  # Skip if file doesn't exist or is excluded
  if [[ ! -f "$file" ]] || should_exclude "$file"; then
    return 0
  fi

  # Only check Rust files
  if [[ "$file" != *.rs ]]; then
    return 0
  fi

  lines=$(wc -l < "$file" 2>/dev/null || echo "0")

  # Determine if this is a test file
  is_test=false
  if [[ "$file" == *_test.rs ]] || \
     [[ "$file" == */tests/* ]] || \
     [[ "$file" == */benches/* ]]; then
    is_test=true
  fi

  local grandfathered
  grandfathered=$(legacy_limit "$file")
  if [[ -n "$grandfathered" ]]; then
    if [ "$lines" -gt "$grandfathered" ]; then
      echo "❌ $file has $lines lines (grandfathered cap $grandfathered; shrink it, don't grow it)"
      exitcode=1
    fi
    return 0
  fi

  if [ "$is_test" = true ]; then
    # Test files can be up to 1000 lines
    if [ "$lines" -gt 1000 ]; then
      echo "❌ $file has $lines lines (max 1000 allowed for test files)"
      exitcode=1
    fi
  else
    # Regular files max 500 lines
    if [ "$lines" -gt 500 ]; then
      echo "❌ $file has $lines lines (max 500 allowed)"
      exitcode=1
    fi
  fi
}

# If files provided as arguments, check only those
if [ $# -gt 0 ]; then
  for file in "$@"; do
    check_file "$file"
  done
else
  # Otherwise check all tracked files
  while IFS= read -r file; do
    check_file "$file"
  done < <(git ls-files)
fi

if [ $exitcode -eq 0 ]; then
  echo "✅ All files within length limits"
fi

exit $exitcode
