#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR"

make bin/radiation_regression_test
./bin/radiation_regression_test

echo "Radiation regression test passed"
