#!/usr/bin/env bash
# Cycle mode smoke test for Vulkan volume rendering.
# Detects multiple export fields, rotates through them over time, and
# runs with orbit camera settings to review temporal field transitions.
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

PREFERRED=(theta w qr qi vorticity_z)
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
  [[ ${#SELECTED[@]} -ge 5 ]] && break
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
  echo "Cycle mode needs at least 2 fields; found ${#SELECTED[@]}" >&2
  exit 1
fi

FIELDS_CSV="$(IFS=,; echo "${SELECTED[*]}")"
START_FIELD="${SELECTED[1]}"

echo "[cycle] input=${INPUT_DIR} fields=${FIELDS_CSV} start=${START_FIELD}"
exec "${BIN}" \
  --window-test \
  --render-backend volume \
  --input "${INPUT_DIR}" \
  --fields "${FIELDS_CSV}" \
  --volume-mode cycle \
  --isolate-field "${START_FIELD}" \
  --component-cycle-fps 0.75 \
  --texture-mode natural \
  --camera-mode orbit \
  --camera-orbit-fps 0.020 \
  --camera-distance 2.28 \
  --camera-height 0.86 \
  --camera-fov-deg 53 \
  --playback-fps 1.1 \
  --ray-steps 236 \
  --ray-threshold 0.25 \
  --ray-opacity 1.25 \
  --ray-brightness 1.06 \
  --ray-ambient 0.82 \
  --ray-anisotropy 0.57 \
  --ray-max-distance 5.3 \
  --sun-dir 0.32,0.23,0.92
