#!/bin/sh
# Builds and tests the skein compiler.
#
#   ./build.sh          build build/skein
#   ./build.sh test     build, then build each src/<part>/test.goose with
#                       the checks of goose --debug (Power of Ten rule 10),
#                       run it, and check the generated C
#                       (scripts/generated_c.sh). That script runs each
#                       model, whose solves check rule 3, and skein run
#                       checks rule 7. Then scripts/power_of_ten.sh checks
#                       rules 1, 4 and 5.
#   ./build.sh gpu      build, then check the GPU path of the page of
#                       `skein build --target wasm` against its CPU path
#                       in Node (tests/gpu/gpu.test.mjs). Needs emcc, Node
#                       and `npm install` in tests/gpu. With no GPU adapter
#                       the checks skip, and SKEIN_GPU=1 makes that a
#                       failure.
#   ./build.sh browser  build, then solve 7 explicit models in the page in
#                       headless Chromium with WebGPU, on the GPU and on the
#                       CPU, and check their expect values
#                       (tests/gpu/browser.mjs). It fails if the adapter is
#                       not Nvidia or a GPU run uses the CPU. Needs emcc,
#                       Node 22+ and Chromium; see the script for the
#                       environment.
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
        if ! "$goose" build -q --debug "$test" -o "$out/test_$part"; then
            echo "FAIL $part: the test does not build"
            status=1
            continue
        fi
        # The solvers test keeps AMG levels in SKEIN_WARM (test_warm). Only
        # that test gets it: the runs of the other tests must not warm start.
        warm=
        if [ "$part" = solvers ]; then warm=$(mktemp -d); fi
        if ! (cd "$root" && SKEIN_WARM="$warm" "$out/test_$part"); then
            status=1
        fi
        if [ -n "$warm" ]; then rm -rf "$warm"; fi
    done
    if ! "$root/scripts/generated_c.sh"; then
        status=1
    fi
    return $status
}

case "${1:-build}" in
    build) build ;;
    test) build && run_tests ;;
    gpu) build && node --test --test-reporter=spec "$root/tests/gpu/gpu.test.mjs" ;;
    browser) build && node "$root/tests/gpu/browser.mjs" ;;
    *) echo "usage: ./build.sh [build|test|gpu|browser]" >&2; exit 2 ;;
esac
