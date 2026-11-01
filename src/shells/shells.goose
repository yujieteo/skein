// The entry of the package shells: shells, beams and laminates (roadmap
// lane 5). A program imports it, and the solvers of its memory check, with
// `import shells from "<path to src>/shells";` and `import solvers from
// "<path to src>/solvers";`, and calls the functions of model.goose, laminate.goose and solve.goose:
//
//     shells::clear();
//     let s = shells::isotropic_section(E, nu, h);
//     let i = shells::mesh_cylinder(R, 0.0, L, 0.0, theta, nx, nt, s, 1.0, 0.0, 0.0);
//     shells::fix(shells::grid_node(i, 0, 0), shells::UX);
//     shells::shell_load(c, 0.0, 0.0, -q, false);
//     let ok = shells::solve();
//     let uz = shells::displacement(n, shells::UZ);
//
// The runtime has no units, as src/runtime/.
namespace shells;
import .laminate;
import .model;
import .mitc4;
import .beam;
import .solve;
import .eigen;     // modes and buckling of the beams
import .mesh;
import .vtu;
import .kernels;   // shell_stiffness and beam_stiffness
