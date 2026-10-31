// The entry of the package meshless: the meshless methods `rkpm` and
// `gfd` (spec §17.5) on a cloud of points in 2D. The physics is one scalar
// unknown u with
//
//     div(a u) - d div(grad(u)) = f(x, y)
//
// and a, d constant, or linear elasticity in plane stress or plane
// strain (rkpm only). A generated program imports it with
// `import meshless from "<path to src>/meshless";` and calls, in this
// sequence:
//
//     meshless::points_rect(0.0, 0.0, 1.0, 1.0, 16, 16);    // or meshless::set_cloud(points, triangles)
//     meshless::equation(1.0, 0.5, 0.01);                   // or meshless::elastic(E, nu, true, t)
//     let ok = meshless::rkpm(2.5);                         // or meshless::gfd(12)
//     meshless::clear_bcs();
//     meshless::select_faces() { x, y => abs(x) <= tol };   // or meshless::select_points(mask)
//     meshless::fix_selected(0);                            // or meshless::fix_component(c, k), meshless::load_selected(c, v)
//     let ok = meshless::solve_static(source, boundary);
//     let error, norm = meshless::l2_error(exact);
//
// `source(x, y, t)` is f, `boundary(k, x, y, t)` is the value of the
// condition k on an edge that fix_selected(k) or fix_component(c, k)
// marked, and `exact(x, y, t)` is a field. Only static steps run.
//
// The package has no units: the program gives each value in the unit
// system of the model header.
namespace meshless;
import .cloud;
import .rk;
import .model;
import .rkpm;
import .gfd;
import .solve;
import .direct;
import .post;
