#!/bin/sh
# Checks two rules of the Power of Ten (docs/tradeoffs.md §6) on generated
# programs and on each runtime file that they import:
#
#   rule 1   no `recursive fn`
#   rule 4   each function has at most 60 lines, from its `fn` line to its
#            closing brace (the splitter, src/codegen/split.goose, makes it
#            so for the generated functions)
#
#   scripts/power_of_ten.sh main.goose...
#
# scripts/generated_c.sh runs it on the program of each model that skein
# run can run.
set -eu

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
' $seen
