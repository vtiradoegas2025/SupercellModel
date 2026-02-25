#!/usr/bin/env bash
# Isolated component mode smoke test for Vulkan volume rendering.
# Selects one representative field from exports and displays it
# with natural texture detail to inspect structure and camera behavior.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
INPUT_DIR="${1:-${ROOT_DIR}/data/exports}"
BIN="${ROOT_DIR}/bin/vulkan_viewer"

if [[ ! -x "${BIN}" ]]; then
  make -C "${ROOT_DIR}" vulkan
fi

STEP_DIR="$(find "${INPUT_DIR}" -maxdepth 1 -type d -name 'step_*' | sort | head -n 1)"
if [[ -z "${STEP_DIR}" ]]; then
  echo "No step_* directories found in ${INPUT_DIR}" >&2
  exit 1
fi

mapfile -t FIELDS < <(find "${STEP_DIR}" -maxdepth 1 -type f -name 'th0_*.npy' -exec basename {} \; | sed -E 's/^th0_//' | sed -E 's/\.npy$//' | sort -u)
if [[ ${#FIELDS[@]} -eq 0 ]]; then
  echo "No th0_*.npy fields found in ${STEP_DIR}" >&2
  exit 1
fi

ISOLATE=""
for candidate in w theta vorticity_z qr; do
  for field in "${FIELDS[@]}"; do
    if [[ "${field}" == "${candidate}" ]]; then
      ISOLATE="${field}"
      break
    fi
  done
  [[ -n "${ISOLATE}" ]] && break
done

if [[ -z "${ISOLATE}" ]]; then
  ISOLATE="${FIELDS[0]}"
fi

echo "[isolated] input=${INPUT_DIR} isolate=${ISOLATE}"
exec "${BIN}" \
  --window-test \
  --render-backend volume \
  --input "${INPUT_DIR}" \
  --field "${ISOLATE}" \
  --volume-mode isolated \
  --isolate-field "${ISOLATE}" \
  --texture-mode natural \
  --camera-mode orbit \
  --camera-orbit-fps 0.030 \
  --camera-distance 2.10 \
  --camera-height 0.80 \
  --camera-fov-deg 58 \
  --playback-fps 1.4 \
  --ray-steps 220 \
  --ray-threshold 0.24 \
  --ray-opacity 1.20 \
  --ray-brightness 1.12 \
  --ray-ambient 0.78 \
  --ray-anisotropy 0.56 \
  --ray-max-distance 4.9 \
  --sun-dir 0.34,0.21,0.92
