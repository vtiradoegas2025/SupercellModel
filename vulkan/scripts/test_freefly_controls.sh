#!/usr/bin/env bash
# Free-fly camera control smoke test for Vulkan volume rendering.
# Starts in supercell mode and enables constrained interactive camera motion,
# so you can validate yaw/pitch movement and bounded navigation behavior.
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

FIELDS_CSV="$(IFS=,; echo "${SELECTED[*]}")"

echo "[freefly] input=${INPUT_DIR} fields=${FIELDS_CSV}"
echo "Controls: WASD strafe/forward, Q/E down/up, arrows look, right-mouse look, Shift boost, R reset"

exec "${BIN}" \
  --window-test \
  --render-backend volume \
  --input "${INPUT_DIR}" \
  --fields "${FIELDS_CSV}" \
  --volume-mode supercell \
  --texture-mode natural \
  --camera-mode freefly \
  --camera-distance 2.35 \
  --camera-height 0.90 \
  --camera-fov-deg 56 \
  --playback-fps 1.0 \
  --ray-steps 256 \
  --ray-threshold 0.25 \
  --ray-opacity 1.30 \
  --ray-brightness 1.08 \
  --ray-ambient 0.80 \
  --ray-anisotropy 0.58 \
  --ray-max-distance 5.6 \
  --sun-dir 0.28,0.25,0.93
