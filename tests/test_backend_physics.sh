#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR"

make bin/tornado_sim bin/field_validator

duration_s="${BACKEND_PHYS_DURATION_S:-5}"
write_every_s="${BACKEND_PHYS_WRITE_EVERY_S:-1}"
archive_root="${BACKEND_PHYS_ARCHIVE_DIR:-$ROOT_DIR/tests/artifacts/backend_physics_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$archive_root"

declare -a configs
if [[ "$#" -gt 0 ]]; then
  configs=("$@")
else
  configs=(lp hp cyclic sharpy_lp)
fi

ensure_sharpy_fixture() {
  local sounding_file="/tmp/tmv_sharpy_baseline.nc"
  if [[ -f "$sounding_file" ]]; then
    return 0
  fi

  python3 - <<'PY'
import numpy as np
import xarray as xr
from pathlib import Path

out = Path("/tmp/tmv_sharpy_baseline.nc")
z = np.array([0, 250, 500, 1000, 1500, 2000, 3000, 4000, 5500, 7000, 9000, 12000, 14000, 16000, 18000], dtype=float)
p = np.array([1000, 970, 940, 900, 860, 820, 740, 670, 590, 520, 430, 300, 220, 150, 100], dtype=float)
t = np.array([300.0, 298.8, 297.5, 294.5, 291.5, 288.5, 281.5, 274.5, 266.0, 258.0, 246.0, 223.0, 214.0, 208.0, 202.0], dtype=float)
td = np.array([295.0, 293.7, 292.4, 288.5, 284.5, 280.5, 270.0, 260.0, 247.0, 236.0, 220.0, 200.0, 190.0, 180.0, 170.0], dtype=float)
ws = np.array([4.0, 5.0, 6.0, 8.0, 11.0, 14.0, 20.0, 25.0, 30.0, 34.0, 40.0, 48.0, 55.0, 60.0, 65.0], dtype=float)
wd = np.array([160.0, 165.0, 170.0, 178.0, 186.0, 195.0, 210.0, 222.0, 235.0, 248.0, 260.0, 270.0, 278.0, 285.0, 292.0], dtype=float)
ds = xr.Dataset(
    data_vars={
        "height_m": (("level",), z),
        "pressure_hpa": (("level",), p),
        "temperature_k": (("level",), t),
        "dewpoint_k": (("level",), td),
        "wind_speed_ms": (("level",), ws),
        "wind_direction_deg": (("level",), wd),
    },
    attrs={
        "station_id": "KOUN",
        "timestamp_utc": "2026-02-20T00:00:00Z",
        "latitude_deg": 35.2,
        "longitude_deg": -97.4,
        "elevation_m": 360.0,
    },
)
ds.to_netcdf(out, engine="scipy")
print(out)
PY
}

for cfg in "${configs[@]}"; do
  if [[ "$cfg" == "sharpy_lp" ]]; then
    ensure_sharpy_fixture
    break
  fi
done

overall_status=0
for cfg in "${configs[@]}"; do
  echo "[backend-physics] config=${cfg} duration=${duration_s}s write_every=${write_every_s}s archive=${archive_root}"
  if ! python3 - "$ROOT_DIR" "$cfg" "$duration_s" "$write_every_s" "$archive_root" <<'PY'
import json
import os
import pathlib
import resource
import shutil
import subprocess
import sys
import time
import numpy as np

root = pathlib.Path(sys.argv[1])
cfg = sys.argv[2]
duration_s = sys.argv[3]
write_every_s = sys.argv[4]
archive_root = pathlib.Path(sys.argv[5])

out_dir = pathlib.Path(f"/tmp/tornado_backend_phys_{cfg}")
report_dir = pathlib.Path(f"/tmp/tornado_backend_phys_{cfg}_reports")
if out_dir.exists():
    shutil.rmtree(out_dir)
if report_dir.exists():
    shutil.rmtree(report_dir)

config_path = root / "configs" / f"{cfg}.yaml"
if not config_path.exists():
    print(f"[backend-physics] missing config: {config_path}")
    raise SystemExit(2)

sim_cmd = [
    str(root / "bin" / "tornado_sim"),
    "--headless",
    f"--config={config_path}",
    f"--duration={duration_s}",
    f"--write-every={write_every_s}",
    f"--outdir={out_dir}",
    "--guard-mode", "strict",
    "--guard-fail-on", "both",
    "--guard-scope", "exported",
    "--guard-report", str(report_dir),
    "--log-profile", "quiet",
]
val_cmd = [
    str(root / "bin" / "field_validator"),
    "--input", str(out_dir),
    "--contract", "cm1",
    "--mode", "strict",
    "--scope", "exported",
    "--json", str(out_dir / "offline_validation.json"),
]

combined_logs = []
start_sim = time.perf_counter()
sim = subprocess.run(sim_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
sim_elapsed_s = time.perf_counter() - start_sim
combined_logs.append(sim.stdout)

val_elapsed_s = 0.0
val = None
if sim.returncode == 0:
    start_val = time.perf_counter()
    val = subprocess.run(val_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    val_elapsed_s = time.perf_counter() - start_val
    combined_logs.append(val.stdout)

log_text = "".join(combined_logs)
(archive_root / f"{cfg}.log").write_text(log_text)

rss_raw = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
if rss_raw > 10_000_000:
    peak_rss_mib_est = rss_raw / (1024.0 * 1024.0)
    peak_rss_unit_guess = "bytes"
else:
    peak_rss_mib_est = rss_raw / 1024.0
    peak_rss_unit_guess = "kilobytes"

validation_path = out_dir / "offline_validation.json"
guard_summary_path = report_dir / "summary.json"
if validation_path.exists():
    shutil.copy2(validation_path, archive_root / f"{cfg}_offline_validation.json")
if guard_summary_path.exists():
    shutil.copy2(guard_summary_path, archive_root / f"{cfg}_guard_summary.json")

nonfinite = 0
bounds = 0
missing_required = 0
missing_exported = 0
known_not_implemented = 0
report_count = 0
validation_failed = True
if validation_path.exists():
    payload = json.loads(validation_path.read_text())
    reports = payload.get("reports", [])
    report_count = len(reports)
    validation_failed = bool(payload.get("failed", False))
    for report in reports:
        missing_required += len(report.get("missing_required", []))
        missing_exported += len(report.get("missing_exported", []))
        known_not_implemented = max(known_not_implemented, len(report.get("known_not_implemented", [])))
        for field in report.get("fields", []):
            stats = field.get("stats", {})
            nonfinite += int(stats.get("nan_count", 0)) + int(stats.get("inf_count", 0))
            bounds += int(stats.get("below_min_count", 0)) + int(stats.get("above_max_count", 0))

integration_success = False
integration_error = ""
integration_shape = []
if sim.returncode == 0:
    try:
        step_dirs = sorted([p for p in out_dir.iterdir() if p.is_dir() and p.name.startswith("step_")])
        if not step_dirs:
            raise RuntimeError("no exported step directories found")

        step0 = step_dirs[0]
        manifest_path = step0 / "manifest.json"
        if not manifest_path.exists():
            raise RuntimeError(f"missing manifest: {manifest_path}")

        manifest = json.loads(manifest_path.read_text())
        manifest_fields = {entry.get("field_id", "") for entry in manifest.get("fields", [])}
        required_manifest_fields = {"u", "v", "w", "rho", "theta", "radar", "temperature"}
        missing_manifest_fields = sorted(required_manifest_fields - manifest_fields)
        if missing_manifest_fields:
            raise RuntimeError(
                "manifest missing required fields: " + ",".join(missing_manifest_fields))

        theta_slice = np.load(step0 / "th0_theta.npy")
        radar_slice = np.load(step0 / "th0_radar.npy")
        u_slice = np.load(step0 / "th0_u.npy")
        v_slice = np.load(step0 / "th0_v.npy")
        w_slice = np.load(step0 / "th0_w.npy")

        if theta_slice.ndim != 2:
            raise RuntimeError(f"isolated field read failed: theta ndim={theta_slice.ndim}")
        if radar_slice.shape != theta_slice.shape:
            raise RuntimeError(
                f"renderer compatibility mismatch: radar shape={radar_slice.shape} theta shape={theta_slice.shape}")
        if u_slice.shape != v_slice.shape or u_slice.shape != w_slice.shape:
            raise RuntimeError(
                f"combined field read shape mismatch: u={u_slice.shape} v={v_slice.shape} w={w_slice.shape}")

        combined_vector = np.stack([u_slice, v_slice, w_slice], axis=0)
        if combined_vector.shape[1:] != theta_slice.shape:
            raise RuntimeError(
                f"combined field read shape mismatch against scalar grid: {combined_vector.shape} vs {theta_slice.shape}")

        integration_shape = list(theta_slice.shape)
        integration_success = True
    except Exception as exc:  # noqa: BLE001
        integration_error = str(exc)

step_time_s = sim_elapsed_s / max(1, report_count)
gate_pass = (
    sim.returncode == 0 and
    val is not None and val.returncode == 0 and
    nonfinite == 0 and
    bounds == 0 and
    missing_required == 0 and
    missing_exported == 0 and
    not validation_failed and
    integration_success
)

metrics = {
    "case": cfg,
    "sim_elapsed_s": sim_elapsed_s,
    "validator_elapsed_s": val_elapsed_s,
    "step_time_s": step_time_s,
    "peak_rss_raw": rss_raw,
    "peak_rss_unit_guess": peak_rss_unit_guess,
    "peak_rss_mib_est": peak_rss_mib_est,
    "report_count": report_count,
    "validator": {
        "nonfinite": nonfinite,
        "bounds": bounds,
        "missing_required": missing_required,
        "missing_exported": missing_exported,
        "known_not_implemented": known_not_implemented,
        "failed": validation_failed,
    },
    "integration": {
        "success": integration_success,
        "shape": integration_shape,
        "error": integration_error,
    },
    "gate": {
        "pass": gate_pass,
        "require_nonfinite_zero": nonfinite == 0,
        "require_bounds_zero": bounds == 0,
        "require_missing_required_zero": missing_required == 0,
        "require_missing_exported_zero": missing_exported == 0,
        "integration_success": integration_success,
        "validator_success": (val is not None and val.returncode == 0 and not validation_failed),
    },
}

(archive_root / f"{cfg}_metrics.json").write_text(json.dumps(metrics, indent=2))

print(f"[backend-physics] metrics {json.dumps(metrics, sort_keys=True)}")
if sim.returncode != 0:
    raise SystemExit(sim.returncode)
if val is None or val.returncode != 0:
    raise SystemExit(1 if val is None else val.returncode)
if not gate_pass:
    raise SystemExit(1)
PY
  then
    overall_status=1
  fi
done

python3 - "$archive_root" "${configs[@]}" <<'PY'
import json
import pathlib
import sys

archive_root = pathlib.Path(sys.argv[1])
configs = sys.argv[2:]
summary = {"archive_root": str(archive_root), "cases": []}
for cfg in configs:
    metrics_path = archive_root / f"{cfg}_metrics.json"
    if not metrics_path.exists():
        summary["cases"].append({"case": cfg, "missing_metrics": True})
        continue
    summary["cases"].append(json.loads(metrics_path.read_text()))

summary_path = archive_root / "summary_metrics.json"
summary_path.write_text(json.dumps(summary, indent=2))
print(f"[backend-physics] archived_summary={summary_path}")
PY

if [[ "$overall_status" -ne 0 ]]; then
  echo "Backend physics validation failed for one or more cases. Archive: ${archive_root}" >&2
  exit 1
fi

echo "Backend physics validation passed for: ${configs[*]}"
echo "Archive: ${archive_root}"
