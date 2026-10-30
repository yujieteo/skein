#!/bin/sh
# Checks or writes the C that Goose generates from each model that skein run
# can run. The C is in the repository next to its model, for audit.
#
#   scripts/generated_c.sh          fail with a diff if a committed C differs
#   scripts/generated_c.sh write    write the C again
#
# For each model, skein run writes .skein/build/<model>/main.goose, and builds
# and runs it. It must pass. Then `goose -o` writes the C of that program,
# which is the same C that the build of skein run compiles.
#
# Goose writes the absolute path of each imported file into the C. This
# script removes the repository directory from these paths, so the C is the
# same on each machine. That is the only change to the C.
#
# Needs build/skein (./build.sh), and goose on PATH or in $GOOSE.
set -eu

# The models that skein run can run. Add a model here when skein run can run it.
models="
examples/femahoy/grid_cantilever.skn
"

root=$(cd "$(dirname "$0")/.." && pwd)
goose=${GOOSE:-goose}
mode=${1:-check}
case "$mode" in
    check|write) ;;
    *) echo "usage: scripts/generated_c.sh [check|write]" >&2; exit 2 ;;
esac

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

status=0
for model in $models; do
    name=$(basename "$model" .skn)
    committed="$root/${model%.skn}.c"
    if ! (cd "$work" && "$root/build/skein" run "$root/$model" >/dev/null 2>"$work/run.err"); then
        cat "$work/run.err"
        echo "FAIL $model: skein run fails"
        status=1
        continue
    fi
    (cd "$work" && "$goose" -q -o "$work/$name.c" ".skein/build/$name/main.goose")
    sed "s|\"$root/|\"|" "$work/$name.c" > "$work/$name.normal.c"
    if [ "$mode" = write ]; then
        cp "$work/$name.normal.c" "$committed"
        echo "wrote ${model%.skn}.c"
    elif ! diff -u "$committed" "$work/$name.normal.c" > "$work/$name.diff"; then
        head -n 40 "$work/$name.diff"
        echo "FAIL $model: the committed C differs. Run scripts/generated_c.sh write."
        status=1
    fi
done
exit $status
