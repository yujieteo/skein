// The entry of the package optim: plane linear elasticity on a grid of
// equal rectangles with bilinear quadrilaterals (Q4), the adjoint of a
// point objective, and the topology loop of femahoy gap G1 (SIMP, the
// sensitivity filter and the optimality criteria update of Andreassen et
// al. 2011). A generated program imports it with
// `import optim from "<path to src>/optim";` and calls, in this sequence:
//
//     optim::mesh_rect(0.0, 0.0, 60.0, 20.0, 60, 20);
//     optim::elastic(1.0, 0.3, 1.0, optim::PLANE_STRESS);
//     optim::clear_bcs();
//     let found = optim::fix_line(optim::X, 0.0, optim::X, 0.0);    // or fix_point(x, y, comp, v)
//     let found = optim::point_load(0.0, 20.0, 0.0, -1.0);           // or edge_load(axis, value, tx, ty)
//     optim::objective(optim::Y, 0.0, 20.0, -1.0);
//     let ok = optim::solve();                                        // or optim::topology(0.5, 3.0, 2.4, 200)
//     let j = optim::objective_value();
//     let dj = optim::input_gradient();                               // the adjoint: dJ/dE, dJ/dnu, dJ/dt, dJ/d of each load
//
// The package has no units: the program gives each value in the unit
// system of the model header.
namespace optim;
import .q4;
import .band;
import .adjoint;
import .topology;
import .post;
