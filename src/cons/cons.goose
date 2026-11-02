// The entry of the package cons: the cell-centred finite volume method for
// a system of conservation laws in 2D on triangles,
//
//     ddt(U) + div(F(U)) = S(U, x, y, t)
//
// with the explicit SSP Runge-Kutta scheme of order 3 (spec §17.3). A
// generated program imports it with `import cons from "<path to src>/cons";`
// and calls, in this sequence:
//
//     cons::read_mesh(".skein/results/wedge", 2);
//     cons::setup(4, cons::VAN_LEER);
//     cons::start(initial);
//     let ok = cons::solve_transient(2.0, 0.0006, flux, source, ghost, wall);
//     let m = cons::mean(0, p_field);
//
// `flux(u, nx, ny, f)` is F·n, `source(u, x, y, t, s)` is S,
// `ghost(mask, u, x, y, t, nx, ny, g)` is the state outside a boundary
// side, and `wall(mask, u, x, y, t, nx, ny, f)` sets the components of the
// boundary flux that a condition on a flux gives (solve.goose). The program
// writes them from the equations of the model, and no function of this
// package knows a physics.
//
// skein run meshes the domain with src/mesher/ and writes the triangles,
// the masks of the boundary sides and the flags of the regions in the
// results directory (gen.goose). Reason: the mesher has recursive functions,
// and a generated program and its runtime have none (Power of Ten rule 1,
// scripts/power_of_ten.sh).
//
// The package has no units: the program gives each value in SI.
namespace cons;
import .mesh;
import .solve;
import .post;
