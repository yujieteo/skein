// The physics templates of the standard library (spec §4.2), in Skein
// source. The checker parses this text after the model, and a physics that
// uses a template takes its statements from it.
//
// `linear_elastic`, `reissner_mindlin` and `timoshenko` have their source in
// the spec (§4.2, §4.3). `heat` is the steady conduction that src/elements
// solves. `navier_stokes` and `stokes` are the incompressible Newtonian
// equations that src/flow solves. `stokes` has no convection term.
// `k_omega_sst` is the k-ω SST model of Menter (1994) with the mean flow,
// kinematic k and w, and the production γ S² of ω. Its constant 1e-20 s^-2
// is the floor of the cross diffusion. The other templates of spec §4.2
// come with their solvers.
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

    template physics heat(q_v: Power/Volume = 0) solves T: Temperature {
        assume isotropic
        let q: Vector<HeatFlux> = fourier(k, T)
        eq div(q) = q_v
    }

    template physics navier_stokes(f: Vector<Force/Volume> = 0) solves v: Vector<Velocity>, p: Pressure {
        assume incompressible, newtonian, laminar
        let sigma: Tensor2<Pressure> = -p*I + 2*mu*sym(grad(v))
        eq rho*(ddt(v) + grad(v)*v) - div(sigma) = f
        eq div(v) = 0
    }

    template physics stokes(f: Vector<Force/Volume> = 0) solves v: Vector<Velocity>, p: Pressure {
        assume incompressible, newtonian, laminar
        let sigma: Tensor2<Pressure> = -p*I + 2*mu*sym(grad(v))
        eq rho*ddt(v) - div(sigma) = f
        eq div(v) = 0
    }

    template physics k_omega_sst(f: Vector<Force/Volume> = 0) solves v: Vector<Velocity>, p: Pressure, k: Velocity^2, w: Frequency {
        assume incompressible, newtonian, turbulent(k_omega_sst)
        let nu_m: KinViscosity = mu/rho
        let y_w: Length = wall_distance()
        let S_t: Frequency = sqrt(2*ddot(sym(grad(v)), sym(grad(v))))
        let CD_kw: Frequency^2 = max(2*0.856/w*dot(grad(k), grad(w)), 1e-20 s^-2)
        let F1: Real = tanh(min(max(sqrt(k)/(0.09*w*y_w), 500*nu_m/(y_w^2*w)), 4*0.856*k/(CD_kw*y_w^2))^4)
        let F2: Real = tanh(max(2*sqrt(k)/(0.09*w*y_w), 500*nu_m/(y_w^2*w))^2)
        let nu_t: KinViscosity = 0.31*k/max(0.31*w, S_t*F2)
        let P_k: Velocity^2*Frequency = min(nu_t*S_t^2, 10*0.09*k*w)
        let sigma: Tensor2<Pressure> = -p*I + 2*rho*(nu_m + nu_t)*sym(grad(v))
        eq rho*(ddt(v) + grad(v)*v) - div(sigma) = f
        eq div(v) = 0
        eq ddt(k) + dot(v, grad(k)) - div((nu_m + (0.85*F1 + 1 - F1)*nu_t)*grad(k)) = P_k - 0.09*k*w
        eq ddt(w) + dot(v, grad(w)) - div((nu_m + (0.5*F1 + 0.856*(1 - F1))*nu_t)*grad(w)) = (0.5532*F1 + 0.4403*(1 - F1))*S_t^2 - (0.075*F1 + 0.0828*(1 - F1))*w^2 + 2*(1 - F1)*0.856/w*dot(grad(k), grad(w))
    }

    template physics reissner_mindlin(layup: Angle = 0 deg, mirror: Bool = false, axis: Vector<Real> = (1, 0, 0), top: Vector<Real>, q: Vector<Pressure> = 0) solves w: Vector<Length>, theta: Vector<Angle> {
        assume small_strain
        let eps_m: Tensor2<Real> = tangent(strain(w), top)
        let kappa: Tensor2<Real/Length> = tangent(sym(grad(cross(theta, top))), top)
        let gamma: Vector<Real> = tangent(grad(dot(w, top)) + cross(theta, top), top)
        let Q_ply: Tensor4<Pressure> = ply_stiffness(E1, E2, G12, nu12)
        let A_lam: Tensor4<Force/Length> = laminate(Q_ply, layup, mirror, t_ply, h, axis, top, 0)
        let B_lam: Tensor4<Force> = laminate(Q_ply, layup, mirror, t_ply, h, axis, top, 1)
        let D_lam: Tensor4<Force*Length> = laminate(Q_ply, layup, mirror, t_ply, h, axis, top, 2)
        let V_s: Vector<Force/Length> = 5/6*G12*h*gamma
        let N: Tensor2<Force/Length> = ddot(A_lam, eps_m) + ddot(B_lam, kappa) + outer(top, V_s)
        let M: Tensor2<Force> = ddot(B_lam, eps_m) + ddot(D_lam, kappa)
        eq div(N) + q = 0
        eq div(M) - V_s = 0
    }

    template physics timoshenko(A: Area, I_y: Length^4, I_z: Length^4, J: Length^4, y_axis: Vector<Real>, offset: Vector<Length> = 0, q: Vector<Force/Length> = 0) solves w: Vector<Length>, theta: Vector<Angle> {
        assume small_strain
        let e_yb: Vector<Real> = y_axis - dot(y_axis, tau)*tau
        let e_zb: Vector<Real> = cross(tau, e_yb)
        let G_b: Pressure = E/(2*(1 + nu))
        let gam_b: Vector<Real> = along(w + cross(theta, offset)) + cross(tau, theta)
        let kap_b: Vector<Real/Length> = along(theta)
        let Nx: Force = E*A*dot(tau, gam_b)
        let F_b: Vector<Force> = Nx*tau + 5/6*G_b*A*(dot(e_yb, gam_b)*e_yb + dot(e_zb, gam_b)*e_zb)
        let M_b: Vector<Force*Length> = G_b*J*dot(tau, kap_b)*tau + E*I_y*dot(e_yb, kap_b)*e_yb + E*I_z*dot(e_zb, kap_b)*e_zb
        eq along(F_b) + q = 0
        eq along(M_b) + cross(tau, F_b) = 0
    }
    """;

// The text of `skein help std`: STD_SOURCE with no model line.
fn std_text() -> const u8[:] {
    var i = 0;
    while i < STD_SOURCE.len && STD_SOURCE[i] != '\n' { i++; }
    while i < STD_SOURCE.len && STD_SOURCE[i] == '\n' { i++; }
    return STD_SOURCE[i..];
}
