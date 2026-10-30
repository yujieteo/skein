#!/bin/sh
# Builds and tests the skein compiler.
#
#   ./build.sh          build build/skein
#   ./build.sh test     build, then run each src/<part>/test.goose and
#                       check the generated C (scripts/generated_c.sh)
#
# Needs the Goose compiler on PATH (or in $GOOSE), and the Goose standard
# library in $GOOSE_STDLIB when goose cannot find it beside itself.
set -eu

root=$(cd "$(dirname "$0")" && pwd)
goose=${GOOSE:-goose}
out="$root/build"
mkdir -p "$out"

build() {
    "$goose" build -q "$root/src/skein.goose" -o "$out/skein"
}

run_tests() {
    status=0
    for test in "$root"/src/*/test.goose; do
        part=$(basename "$(dirname "$test")")
        if ! "$goose" build -q "$test" -o "$out/test_$part"; then
            echo "FAIL $part: the test does not build"
            status=1
            continue
        fi
        if ! (cd "$root" && "$out/test_$part"); then
            status=1
        fi
    done
    if ! "$root/scripts/generated_c.sh"; then
        status=1
    fi
    return $status
}

case "${1:-build}" in
    build) build ;;
    test) build && run_tests ;;
    *) echo "usage: ./build.sh [build|test]" >&2; exit 2 ;;
esac
