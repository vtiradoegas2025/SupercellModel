#!/usr/bin/env bash
# Composite mode smoke test for Vulkan volume rendering.
# Chooses a compact multi-field set from exported data and runs
# a blended render with realistic texture detail and orbit camera.
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

PREFERRED=(theta w qr vorticity_z)
SELECTED=()
for name in "${PREFERRED[@]}"; do
  for field in "${FIELDS[@]}"; do
    if [[ "${field}" == "${name}" ]]; then
      SELECTED+=("${field}")
      break
    fi
  done
done
for field in "${FIELDS[@]}"; do
  [[ ${#SELECTED[@]} -ge 4 ]] && break
  found=0
  for chosen in "${SELECTED[@]}"; do
    if [[ "${chosen}" == "${field}" ]]; then
      found=1
      break
    fi
  done
  [[ ${found} -eq 0 ]] && SELECTED+=("${field}")
done

if [[ ${#SELECTED[@]} -lt 2 ]]; then
  echo "Composite mode needs at least 2 fields; found ${#SELECTED[@]}" >&2
  exit 1
fi

FIELDS_CSV="$(IFS=,; echo "${SELECTED[*]}")"

echo "[composite] input=${INPUT_DIR} fields=${FIELDS_CSV}"
exec "${BIN}" \
  --window-test \
  --render-backend volume \
  --input "${INPUT_DIR}" \
  --fields "${FIELDS_CSV}" \
  --volume-mode composite \
  --texture-mode natural \
  --camera-mode orbit \
  --camera-orbit-fps 0.022 \
  --camera-distance 2.30 \
  --camera-height 0.84 \
  --camera-fov-deg 54 \
  --playback-fps 1.2 \
  --ray-steps 240 \
  --ray-threshold 0.26 \
  --ray-opacity 1.22 \
  --ray-brightness 1.05 \
  --ray-ambient 0.85 \
  --ray-anisotropy 0.55 \
  --ray-max-distance 5.4 \
  --sun-dir 0.30,0.24,0.92
