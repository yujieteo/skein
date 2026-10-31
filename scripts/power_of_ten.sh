#!/bin/sh
# Checks rules of the Power of Ten (docs/tradeoffs.md §6) on generated
# programs and on each runtime file that they import:
#
#   rule 1   no `recursive fn`
#   rule 4   each function has at most 60 lines, from its `fn` line to its
#            closing brace (the splitter, src/codegen/split.goose, makes it
#            so for the generated functions)
#   rule 5   each call of tests/power_of_ten/faults.goose, with a broken
#            precondition, stops on an assert of the runtime
#
#   scripts/power_of_ten.sh main.goose...
#
# scripts/generated_c.sh runs it on the program of each model that skein
# run can run. Each of these runs also checks rule 3: each solve of the
# runtime asserts that it grows no array of the state. skein run checks
# rule 7 on each program that it writes (src/codegen/status.goose), and
# ./build.sh test builds the tests with goose --debug for rule 10.
#
# Needs goose on PATH or in $GOOSE.
set -eu

if [ $# -eq 0 ]; then
    echo "usage: scripts/power_of_ten.sh main.goose..." >&2
    exit 2
fi

root=$(cd "$(dirname "$0")/.." && pwd -P)

# The absolute path of file $1.
absolute() { echo "$(cd "$(dirname "$1")" && pwd -P)/$(basename "$1")"; }

# The arguments and the files that they import, each once. A package import
# `import name from "dir";` is the file dir/name.goose. The standard library
# is not checked.
todo=""
for f in "$@"; do todo="$todo $(absolute "$f")"; done
seen=""
while [ -n "$todo" ]; do
    set -- $todo
    f=$1
    shift
    todo="$*"
    case " $seen " in *" $f "*) continue ;; esac
    seen="$seen $f"
    dir=$(dirname "$f")
    for next in $(sed -n -e 's/^import \.\([a-z_0-9]*\);.*/\1.goose/p' \
                         -e 's/^import \([a-z_0-9]*\) from "\(.*\)";.*/\2\/\1.goose/p' "$f"); do
        case $next in /*) ;; *) next="$dir/$next" ;; esac
        todo="$todo $(absolute "$next")"
    done
done

# A function starts with `fn` at the start of a line (after its keywords)
# and ends at the first `}` at the start of a line.
status=0
awk -v root="$root/" '
    function where() { f = FILENAME; sub("^" root, "", f); return f ":" FNR }
    /^(recursive |inline |simd |checkfree )*fn / {
        if ($0 ~ /^recursive /) { print "FAIL " where() ": recursive fn (rule 1)"; bad = 1 }
        name = $0; sub(/\(.*/, "", name); sub(/.*fn /, "", name)
        start = FNR; at = where(); open = ($0 !~ /}[ \t]*$/)
        next
    }
    open && /^}/ {
        open = 0
        if (FNR - start + 1 > 60) { print "FAIL " at ": fn " name " has " FNR - start + 1 " lines, more than 60 (rule 4)"; bad = 1 }
    }
    END { exit bad }
' $seen || status=1

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
"${GOOSE:-goose}" build -q "$root/tests/power_of_ten/faults.goose" -o "$work/faults"
for call in runtime_material runtime_cells_for runtime_mesh_box runtime_solve_static \
            elements_elastic elements_fix shells_laminate_ply shells_fix solvers_multiply \
            solvers_ilu0_apply; do
    if "$work/faults" "$call" 2>"$work/err" || ! grep -q "assert failed" "$work/err"; then
        cat "$work/err"
        echo "FAIL faults.goose $call: no assert stops it (rule 5)"
        status=1
    fi
done
exit $status
