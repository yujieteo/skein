#!/bin/sh
# Checks rules of the Power of Ten (docs/tradeoffs.md §6) on generated
# programs and on each runtime file that they import:
#
#   rule 1   no `recursive fn`
#   rule 4   each function has at most 60 lines, from its `fn` line to its
#            closing brace (the splitter, src/codegen/split.goose, makes it
#            so for the generated functions)
#   rule 5   each call of tests/power_of_ten/faults.goose, with a broken
#            precondition, stops on the abort of that precondition
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
# Each call of faults.goose and the message of the abort that must stop it.
while IFS='|' read -r call message; do
    if "$work/faults" "$call" 2>"$work/err" || ! grep -qF "goose runtime error: $message (" "$work/err"; then
        cat "$work/err"
        echo "FAIL faults.goose $call: no abort \"$message\" stops it (rule 5)"
        status=1
    fi
done <<'CALLS'
runtime_material|material: E > 0 and -1 < nu < 0.5
runtime_cells_for|cells_for: length > 0 and size > 0
runtime_mesh_box|build_lattice: counts > 0 and hi > lo on each axis
runtime_fix|fix: a mesh and its dofs
runtime_solve_static|solve_linear_static: a material, a mesh and its dofs
elements_elastic|set_elastic: E > 0, -1 < nu < 0.5, t > 0 and a plane assumption
elements_fix|fix_region: a physics, -1 <= comp < ncomp and the arrays of the dofs
elements_traction|add_region_load: a physics and the arrays of the dofs
elements_solve|solve_system: a physics and the arrays of the dofs
shells_laminate_ply|laminate_ply: t, e1, e2, g12, g13 and g23 > 0
shells_fix|fix: 0 <= n < node_count(), and 0 <= dof < 6
shells_force|force: 0 <= n < node_count(), and 0 <= dof < 6
shells_mesh_cylinder|mesh_cylinder: r > 0, x1 > x0, t1 > t0, nx, nt > 0, and 0 <= section < sections.len
flow_mesh|mesh_from_arrays: each vertex of a triangle < xy.len / 2
flow_fluid|fluid: density > 0, viscosity > 0 and a stabilize value
flow_fix|fix_velocity: comp 0 or 1, and clear_bcs on a mesh
flow_solve|solve_steady: clear_bcs, and rtol > 0
flow_value|value_at: a value, and a solve
grid_mesh|mesh_rect: xb > xa, yb > ya, and cells_x, cells_y >= 1
grid_equation|equation: a mesh, c_t >= 0 and d >= 0
grid_dg|dg: a mesh, 1 <= p <= 3, the flux upwind or lax_friedrichs, and sigma > 0
grid_fix|fix_selected: k >= 0, a method and select_faces
solvers_multiply|multiply: row_start has y.len + 1 entries, col and val one size, x and y one size
solvers_ilu0_apply|ilu0_apply: row_start has r.len + 1 entries, z and diag r.len, lu col.len
CALLS
exit $status
