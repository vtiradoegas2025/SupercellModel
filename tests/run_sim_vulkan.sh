#!/usr/bin/env bash
set -euo pipefail

MODE="full"  # full | sim | view
CONFIG="/tmp/full_duration_all_fields.yaml"
OUTDIR="/tmp/full_duration_all_fields_out"
COMPONENT="all"  # all or a single field (e.g. theta, vorticity_z)
STYLE="cinematic-bw"  # default | cinematic-bw

PLAYBACK_FPS="1.5"
RAY_STEPS="320"
RAY_THRESHOLD="0.22"
RAY_OPACITY="1.65"
RAY_BRIGHTNESS="1.05"
RAY_AMBIENT="0.55"
RAY_ANISOTROPY="0.62"
RAY_MAX_DISTANCE="6.2"
SUN_DIR="0.25,0.18,0.95"
WINDOW_WIDTH="1280"
WINDOW_HEIGHT="720"
WINDOW_FRAMES="0"

SKIP_BUILD=0

print_help() {
  cat <<'EOF'
Run TornadoModel simulation and Vulkan viewer with optional cinematic style.

Usage:
  ./run_sim_vulkan.sh [options]

Options:
  --mode <full|sim|view>         Run full pipeline, simulation only, or viewer only (default: full)
  --config <path>                Simulation config path (default: /tmp/full_duration_all_fields.yaml)
  --outdir <path>                Simulation output directory (default: /tmp/full_duration_all_fields_out)
  --component <all|field>        Viewer field selection; "all" auto-detects exported fields (default: all)
  --style <default|cinematic-bw> Vulkan viewer style (default: cinematic-bw)
  --window-width <n>             Vulkan window width (default: 1280)
  --window-height <n>            Vulkan window height (default: 720)
  --window-frames <n>            Frames to render; 0 means until window close (default: 0)
  --playback-fps <f>             Volume playback FPS (default: 1.5)
  --ray-steps <n>                Raymarch steps (default: 320)
  --ray-threshold <f>            Density threshold (default: 0.22)
  --ray-opacity <f>              Extinction scale (default: 1.65)
  --ray-brightness <f>           Brightness scale (default: 1.05)
  --ray-ambient <f>              Ambient scale (default: 0.55)
  --ray-anisotropy <f>           HG anisotropy (default: 0.62)
  --ray-max-distance <f>         Maximum ray distance (default: 6.2)
  --sun-dir <x,y,z>              Sun direction (default: 0.25,0.18,0.95)
  --skip-build                   Skip build steps
  -h, --help                     Show this help

Examples:
  ./run_sim_vulkan.sh
  ./run_sim_vulkan.sh --mode sim
  ./run_sim_vulkan.sh --mode view --component theta
  ./run_sim_vulkan.sh --mode view --component all --style cinematic-bw
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      MODE="${2:-}"
      shift 2
      ;;
    --config)
      CONFIG="${2:-}"
      shift 2
      ;;
    --outdir)
      OUTDIR="${2:-}"
      shift 2
      ;;
    --component)
      COMPONENT="${2:-}"
      shift 2
      ;;
    --style)
      STYLE="${2:-}"
      shift 2
      ;;
    --window-width)
      WINDOW_WIDTH="${2:-}"
      shift 2
      ;;
    --window-height)
      WINDOW_HEIGHT="${2:-}"
      shift 2
      ;;
    --window-frames)
      WINDOW_FRAMES="${2:-}"
      shift 2
      ;;
    --playback-fps)
      PLAYBACK_FPS="${2:-}"
      shift 2
      ;;
    --ray-steps)
      RAY_STEPS="${2:-}"
      shift 2
      ;;
    --ray-threshold)
      RAY_THRESHOLD="${2:-}"
      shift 2
      ;;
    --ray-opacity)
      RAY_OPACITY="${2:-}"
      shift 2
      ;;
    --ray-brightness)
      RAY_BRIGHTNESS="${2:-}"
      shift 2
      ;;
    --ray-ambient)
      RAY_AMBIENT="${2:-}"
      shift 2
      ;;
    --ray-anisotropy)
      RAY_ANISOTROPY="${2:-}"
      shift 2
      ;;
    --ray-max-distance)
      RAY_MAX_DISTANCE="${2:-}"
      shift 2
      ;;
    --sun-dir)
      SUN_DIR="${2:-}"
      shift 2
      ;;
    --skip-build)
      SKIP_BUILD=1
      shift
      ;;
    -h|--help)
      print_help
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      print_help
      exit 1
      ;;
  esac
done

if [[ "$MODE" != "full" && "$MODE" != "sim" && "$MODE" != "view" ]]; then
  echo "Invalid --mode: $MODE (expected full|sim|view)" >&2
  exit 1
fi

if [[ "$STYLE" != "default" && "$STYLE" != "cinematic-bw" ]]; then
  echo "Invalid --style: $STYLE (expected default|cinematic-bw)" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

run_simulation() {
  if [[ ! -f "$CONFIG" ]]; then
    echo "Simulation config not found: $CONFIG" >&2
    exit 1
  fi
  echo "[run] tornado_sim -> $OUTDIR"
  ./bin/tornado_sim --headless --config="$CONFIG" --outdir="$OUTDIR" --log-profile normal
}

resolve_fields_csv() {
  if [[ "$COMPONENT" != "all" ]]; then
    printf '%s' "$COMPONENT"
    return
  fi

  local step0="$OUTDIR/step_000000"
  if [[ ! -d "$step0" ]]; then
    echo "Expected step directory missing: $step0" >&2
    exit 1
  fi

  local fields
  fields="$(
    find "$step0" -maxdepth 1 -type f -name 'th0_*.npy' \
      -exec basename {} \; \
      | sed -E 's/^th0_//' \
      | sed -E 's/\.npy$//' \
      | sort -u \
      | paste -sd, -
  )"

  if [[ -z "$fields" ]]; then
    echo "No th0_*.npy fields found in $step0" >&2
    exit 1
  fi

  printf '%s' "$fields"
}

run_viewer() {
  if [[ ! -d "$OUTDIR" ]]; then
    echo "Output directory not found: $OUTDIR" >&2
    exit 1
  fi

  local fields_csv
  fields_csv="$(resolve_fields_csv)"
  echo "[run] vulkan_viewer fields=$fields_csv style=$STYLE"

  ./bin/vulkan_viewer --window-test --render-backend volume \
    --input "$OUTDIR" \
    --fields "$fields_csv" \
    --style "$STYLE" \
    --window-width "$WINDOW_WIDTH" \
    --window-height "$WINDOW_HEIGHT" \
    --window-frames "$WINDOW_FRAMES" \
    --playback-fps "$PLAYBACK_FPS" \
    --ray-steps "$RAY_STEPS" \
    --ray-threshold "$RAY_THRESHOLD" \
    --ray-opacity "$RAY_OPACITY" \
    --ray-brightness "$RAY_BRIGHTNESS" \
    --ray-ambient "$RAY_AMBIENT" \
    --ray-anisotropy "$RAY_ANISOTROPY" \
    --ray-max-distance "$RAY_MAX_DISTANCE" \
    --sun-dir "$SUN_DIR"
}

if [[ "$SKIP_BUILD" -eq 0 ]]; then
  if [[ "$MODE" == "full" || "$MODE" == "sim" ]]; then
    echo "[build] make -j4"
    make -j4
  fi
  if [[ "$MODE" == "full" || "$MODE" == "view" ]]; then
    echo "[build] make vulkan"
    make vulkan
  fi
fi

if [[ "$MODE" == "full" || "$MODE" == "sim" ]]; then
  run_simulation
fi

if [[ "$MODE" == "full" || "$MODE" == "view" ]]; then
  run_viewer
fi

