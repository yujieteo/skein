// The entry of the package grid: the cell-centred finite volume (`fv`,
// spec §17.3) and the discontinuous Galerkin (`dg`, spec §17.4) methods of
// one scalar unknown u in 2D, on a grid of rectangles:
//
//     ct ddt(u) + div(a u) - d div(grad(u)) = f(x, y, t)
//
// with ct, the velocity a and the diffusivity d constant. A generated
// program imports it with `import grid from "<path to src>/grid";` and
// calls, in this sequence:
//
//     grid::mesh_rect(0.0, 0.0, 1.0, 1.0, 20, 20);
//     grid::equation(1.0, 1.0, 0.5, 0.01);
//     grid::fv(grid::LINEAR_UPWIND, grid::LEAST_SQUARES, grid::VAN_LEER);   // or grid::dg(2, grid::UPWIND, 3.0)
//     grid::clear_bcs();
//     grid::select_faces() { x, y => abs(x) <= tol };
//     grid::fix_selected(0);
//     grid::start(initial);
//     let ok = grid::solve_static(1.0e-10, source, boundary);
//     let ok = grid::solve_transient(1.0, 0.0, source, boundary);
//     let error, norm = grid::l2_error(exact);
//
// `source(x, y, t)` is f, `boundary(k, x, y, t)` is the value of the
// essential condition k on a face that fix_selected(k) marked, and
// `exact(x, y, t)` and `initial(x, y)` are fields. A face that no
// condition marks has no flux: F·n = 0 (spec §12.3).
//
// The package has no units: the program gives each value in the unit
// system of the model header.
namespace grid;
import .mesh;
import .fv;
import .dg;
import .solve;
import .post;
