// The physics templates of the standard library (spec §4.2), in Skein
// source. The checker parses this text after the model, and a physics that
// uses a template takes its statements from it.
//
// Only `linear_elastic` has its source in the spec. `navier_stokes` is the
// source of this checker: the incompressible Newtonian equations of
// sample 02. The other templates of spec §4.2 come with their solvers.
import std;
namespace check;

let STD_PATH = "std.skn";

let STD_SOURCE = """
    model std: 3d, SI

    template physics linear_elastic(b: Vector<Force/Volume> = 0) solves u: Vector<Length> {
        assume small_strain, isotropic
        let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u))
        eq div(sigma) + b = 0
    }

    template physics navier_stokes(f: Vector<Force/Volume> = 0) solves v: Vector<Velocity>, p: Pressure {
        assume incompressible, newtonian, laminar
        let sigma: Tensor2<Pressure> = -p*I + 2*mu*sym(grad(v))
        eq rho*(ddt(v) + grad(v)*v) - div(sigma) = f
        eq div(v) = 0
    }
    """;
