#!/usr/bin/env bash
set -euo pipefail
APP_DIR="${APP_DIR:-/app}"
cp "$(dirname "$0")/solve.cpp" "$APP_DIR/solve.cpp"
g++ -std=c++17 -O2 -pipe "$APP_DIR/solve.cpp" -o "$APP_DIR/solve"
