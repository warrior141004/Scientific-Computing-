#!/usr/bin/env bash
set -uo pipefail

mkdir -p /logs/verifier
printf '0\n' > /logs/verifier/reward.txt
if [ "$PWD" = "/" ]; then
    echo "Error: No working directory set." >&2
    exit 0
fi
PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 -m pytest -p no:cacheprovider -rA -v /tests/test_outputs.py
rc=$?
if [ "$rc" -eq 0 ]; then
    echo 1 > /logs/verifier/reward.txt
else
    echo 0 > /logs/verifier/reward.txt
fi
