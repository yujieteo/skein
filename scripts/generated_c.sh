#!/bin/sh
# Checks or writes the C that Goose generates from each model that skein run
# can run. The C is in the repository next to its model, for audit.
#
#   scripts/generated_c.sh          fail with a diff if a committed C differs
#   scripts/generated_c.sh write    write the C again
#
# For each model, skein run writes .skein/build/<model>/main.goose, and builds
# and runs it. It must pass. Then `goose -o` writes the C of that program,
# which is the same C that the build of skein run compiles. Then
# scripts/power_of_ten.sh checks each program and the runtime that it imports.
#
# Each solve of a model must check its memory first (src/solvers/memory.goose):
# its records must hold at least as many memory_checks as solves. With the
# abort of memory_check on a path that has no memory model, a new solution
# path cannot skip the model, if a model here runs it.
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
examples/elements/bar_contact_tri.skn
examples/elements/bar_wave_tri.skn
examples/elements/bloch_cell_tri.skn
examples/elements/buckling_column_tet.skn
examples/elements/column_weight_tet.skn
examples/elements/conduction_slab_tet.skn
examples/elements/cook_membrane_tri.skn
examples/elements/darcy_annulus_tri.skn
examples/elements/far_field_box_tet.skn
examples/elements/frequency_cantilever_tet.skn
examples/elements/harmonic_duct_tri.skn
examples/elements/heated_bar_study_tet.skn
examples/elements/heated_bar_tet.skn
examples/elements/heated_bar_template_tet.skn
examples/elements/kirsch_nurbs_tri.skn
examples/elements/kirsch_plate_tri.skn
examples/elements/pml_duct_tri.skn
examples/elements/radiation_annulus_tri.skn
examples/elements/taylor_bar_tet.skn
examples/elements/uniform_temperature_tet.skn
examples/femahoy/bar_pressure.skn
examples/femahoy/cantilever_2step.skn
examples/femahoy/contact_cattaneo.skn
examples/femahoy/contact_hertz.skn
examples/femahoy/dynamics_bar_step.skn
examples/femahoy/grid_cantilever.skn
examples/femahoy/grid_conductivity_slab.skn
examples/femahoy/grid_cooling_slab.skn
examples/femahoy/grid_coax_capacitor.skn
examples/femahoy/grid_coax_tem.skn
examples/femahoy/grid_free_block_modes_ccx.skn
examples/femahoy/grid_frequency_cantilever.skn
examples/femahoy/grid_heated_bar.skn
examples/femahoy/grid_kirsch_hole.skn
examples/femahoy/grid_lame_cylinder.skn
examples/femahoy/grid_pipe_wall.skn
examples/femahoy/grid_radiating_fin_ccx.skn
examples/femahoy/grid_radiation_slab.skn
examples/femahoy/grid_round_cantilever.skn
examples/femahoy/heat_fin.skn
examples/femahoy/hyperelastic_tube.skn
examples/femahoy/plasticity_sphere.skn
examples/femahoy/scordelis_lo_roof.skn
examples/femahoy/taylor_bar_impact.skn
examples/femahoy/topology_mbb_beam.skn
examples/femahoy/warpage_thick_copper.skn
examples/flow/cavity_coarse.skn
examples/flow/channel_sst_tri.skn
examples/flow/poiseuille_stokes.skn
examples/flow/turek_hron_fsi1_tri.skn
examples/grid/mms_dg_static.skn
examples/grid/mms_dg_transient.skn
examples/grid/mms_fv_static.skn
examples/grid/mms_fv_transient.skn
examples/lame_profile.skn
examples/meshless/mms_gfd.skn
examples/meshless/mms_rkpm.skn
examples/optim/stiffness_bar.skn
examples/shells/beam_buckling.skn
examples/shells/beam_modes.skn
examples/shells/cantilever_beam.skn
examples/shells/stiffened_barrel.skn
examples/shells/thin_cylinder_shell.skn
examples/shells/window_barrel.skn
samples/02_cavity.skn
samples/05_meshless_plate.skn
samples/06_cylinder_sweep.skn
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
programs=""
for model in $models; do
    # skein run names the build directory after the `model <name>:` line.
    name=$(sed -n 's/^model \([A-Za-z0-9_]*\):.*/\1/p' "$root/$model" | head -n 1)
    committed="$root/${model%.skn}.c"
    if ! (cd "$work" && "$root/build/skein" run "$root/$model" >"$work/run.toon" 2>"$work/run.err"); then
        cat "$work/run.err"
        echo "FAIL $model: skein run fails"
        status=1
        continue
    fi
    solves=$(sed -n 's/^solve\[\([0-9]*\)\].*/\1/p' "$work/run.toon")
    checks=$(sed -n 's/^memory_checks: //p' "$work/run.toon")
    if [ -z "$checks" ] || [ "$checks" -lt "${solves:-0}" ]; then
        echo "FAIL $model: ${solves:-0} solves, and ${checks:-no} memory checks (src/solvers/memory.goose)"
        status=1
    fi
    programs="$programs $work/.skein/build/$name/main.goose"
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
if ! "$root/scripts/power_of_ten.sh" $programs; then
    status=1
fi
exit $status
