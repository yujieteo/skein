#!/bin/sh
# Checks that each model in tests/unsupported makes skein run report the
# diagnostic `unsupported` and not a wrong answer.
#
# Needs build/skein (./build.sh).
set -eu
root=$(cd "$(dirname "$0")/.." && pwd)
status=0
for model in "$root"/tests/unsupported/*.skn; do
    if ! "$root/build/skein" run "$model" 2>&1 | grep -q 'error,unsupported,'; then
        echo "FAIL $(basename "$model"): skein run did not report unsupported"
        status=1
    fi
done
exit $status
