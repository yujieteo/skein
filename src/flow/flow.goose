// The entry of the package flow: steady incompressible Navier-Stokes and
// Stokes flow in 2D on P1 or P2 triangles, with SUPG, PSPG and LSIC
// stabilization (spec §17.2). A generated program imports it with
// `import flow from "<path to src>/flow";` and calls, in this sequence:
//
//     flow::read_mesh(".skein/results/cavity", 1);
//     flow::fluid(1.0, 0.01, true, flow::SUPG_PSPG);
//     flow::clear_bcs();
//     flow::select_boundary() { x, y => abs(y - 1.0) <= 1.0e-6 };
//     flow::fix_velocity(0, 1.0);
//     flow::fix_velocity(1, 0.0);
//     flow::select_point(0.0, 0.0);
//     flow::fix_pressure(0.0);
//     let ok = flow::solve_steady(1.0e-8);
//     let u, found = flow::value_at(flow::VALUE_VX, 0.5, 0.5);
//
// skein run meshes the domain with src/mesher/ and writes the triangles in
// the results directory (gen.goose). Reason: the mesher has recursive
// functions, and a generated program and its runtime have none (Power of
// Ten rule 1, scripts/power_of_ten.sh).
//
// The package has no units: the generated program gives each value in the
// scaled units of the `scales` line (spec §21.2), and multiplies each
// result by its reference value.
namespace flow;
import .mesh;
import .cell;
import .solve;
import .post;
