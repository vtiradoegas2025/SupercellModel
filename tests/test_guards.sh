#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR"

make bin/tornado_sim bin/field_validator

OUT_DIR="/tmp/tornado_guard_test"
rm -rf "$OUT_DIR"

./bin/tornado_sim \
  --headless \
  --config=configs/classic.yaml \
  --duration=1 \
  --write-every=1 \
  --outdir="$OUT_DIR" \
  --guard-mode strict \
  --guard-fail-on both \
  --guard-scope exported \
  --log-profile quiet

./bin/field_validator \
  --input "$OUT_DIR" \
  --contract cm1 \
  --mode strict \
  --scope exported \
  --json "$OUT_DIR/offline_validation.json"

# Removing an exported diagnostic should fail strict exported validation.
find "$OUT_DIR" -name "*_theta_v.npy" -delete
if ./bin/field_validator \
  --input "$OUT_DIR" \
  --contract cm1 \
  --mode strict \
  --scope exported \
  --json "$OUT_DIR/offline_validation_missing_exported.json"; then
  echo "Expected strict exported validation to fail when theta_v files are missing, but it succeeded"
  exit 1
fi

# Required computed field should fail strict mode when forced out of bounds.
REQ_CFG="/tmp/classic_required_computed_fail.yaml"
cp configs/classic.yaml "$REQ_CFG"
cat >> "$REQ_CFG" <<'YAML'
validation:
  guard_mode: strict
  guard_fail_on: both
  guard_scope: required
  field_overrides:
    p_prime:
      min: 1.0
      max: 0.0
YAML

if ./bin/tornado_sim \
  --headless \
  --config="$REQ_CFG" \
  --duration=1 \
  --write-every=1 \
  --outdir="/tmp/tornado_guard_required_fail" \
  --log-profile quiet; then
  echo "Expected strict failure for required computed field (p_prime), but run succeeded"
  exit 1
fi

# Report-only computed field should not fail strict mode when out of bounds.
REPORT_CFG="/tmp/classic_report_only_computed_pass.yaml"
cp configs/classic.yaml "$REPORT_CFG"
cat >> "$REPORT_CFG" <<'YAML'
validation:
  guard_mode: strict
  guard_fail_on: both
  guard_scope: required
  field_overrides:
    dynamic_pressure:
      min: 1.0
      max: 0.0
YAML

./bin/tornado_sim \
  --headless \
  --config="$REPORT_CFG" \
  --duration=1 \
  --write-every=1 \
  --outdir="/tmp/tornado_guard_report_only_pass" \
  --log-profile quiet

echo "Guard validation smoke test passed"
