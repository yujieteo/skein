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
# script removes the repository directory from these paths. That is the only
# change to the C. The C is then the same on each machine with the Goose
# version in goose_version. After you change Goose, set goose_version to the
# output of `goose --version` and write the C again.
#
# Needs build/skein (./build.sh), and goose on PATH or in $GOOSE.
set -eu

# The models that skein run can run. Add a model here when skein run can run it.
models="
examples/closed_cylinder.skn
examples/elements/column_weight_tet.skn
examples/elements/conduction_slab_tet.skn
examples/elements/heated_bar_tet.skn
examples/elements/heated_bar_template_tet.skn
examples/elements/kirsch_plate_tri.skn
examples/elements/uniform_temperature_tet.skn
examples/femahoy/grid_cantilever.skn
examples/femahoy/grid_kirsch_hole.skn
examples/femahoy/scordelis_lo_roof.skn
examples/shells/thin_cylinder_shell.skn
"

# The Goose that wrote the committed C.
goose_version="goose v2026.3-120-ga108677d (runtime fba599d313fe)"

root=$(cd "$(dirname "$0")/.." && pwd -P)
goose=${GOOSE:-goose}
mode=${1:-check}
case "$mode" in
    check|write) ;;
    *) echo "usage: scripts/generated_c.sh [check|write]" >&2; exit 2 ;;
esac

local_version=$("$goose" --version)
if [ "$local_version" != "$goose_version" ]; then
    echo "FAIL: goose is \"$local_version\", the committed C needs \"$goose_version\". Build goose at that commit and set GOOSE to it."
    exit 1
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

status=0
for model in $models; do
    # skein run names the build directory after the `model <name>:` line.
    name=$(sed -n 's/^model \([A-Za-z0-9_]*\):.*/\1/p' "$root/$model" | head -n 1)
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
