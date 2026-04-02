#!/usr/bin/env bash
# Create a transfer archive with code/runtime assets only (no large local data).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

STAMP="$(date +%Y%m%d_%H%M%S)"
OUT_PATH="${ROOT_DIR}/pangenome_mapping_code_${STAMP}.tar.gz"
INCLUDE_TESTS="false"

usage() {
  cat <<'EOF'
Usage:
  scripts/create_cluster_archive.sh [--output /path/archive.tar.gz] [--include-tests]

Description:
  Builds a compressed tarball containing code and runtime files required to run
  the mapping pipeline on a cluster, excluding local data like cluster_test_data.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --output|-o)
      OUT_PATH="$2"
      shift 2
      ;;
    --include-tests)
      INCLUDE_TESTS="true"
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

# Paths needed for runtime execution.
INCLUDE_PATHS=(
  "cli.py"
  "pipeline.py"
  "requirements.txt"
  "run_mapping.slurm"
  "README.md"
  "src"
  "stages"
  "validation"
  "scripts"
  "benchmarks"
  "alignment_viewer"
)

if [[ "${INCLUDE_TESTS}" == "true" ]]; then
  INCLUDE_PATHS+=("tests")
fi

# Keep only paths that actually exist.
EXISTING_PATHS=()
for rel in "${INCLUDE_PATHS[@]}"; do
  if [[ -e "${ROOT_DIR}/${rel}" ]]; then
    EXISTING_PATHS+=("${rel}")
  fi
done

if [[ ${#EXISTING_PATHS[@]} -eq 0 ]]; then
  echo "No files selected for archive. Aborting." >&2
  exit 1
fi

mkdir -p "$(dirname "${OUT_PATH}")"

echo "Creating archive: ${OUT_PATH}"
(
  cd "${ROOT_DIR}"
  tar -czf "${OUT_PATH}" "${EXISTING_PATHS[@]}"
)

echo "Archive created successfully."
echo "Size:"
du -h "${OUT_PATH}"
echo "Top-level contents:"
tar -tzf "${OUT_PATH}" | awk -F/ 'NF{print $1}' | sort -u

cat <<EOF

Transfer command example:
  scp "${OUT_PATH}" <user>@<cluster-host>:<remote-dir>/
EOF
