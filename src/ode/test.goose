import std;
import math;
import ode from "../ode";
import solvers from "../solvers";

var failures = 0;

fn fail(what: u8[:]) {
    failures++;
    print("FAIL ", what);
}

fn near(what: u8[:], got: f64, want: f64, rel: f64) {
    if !(abs(got - want) <= rel * max(abs(want), 1.0e-300)) { fail(str(what, ": ", got, ", not ", want, " within ", rel)); }
}

fn test_jets() {
    let a: f64[5] = [2.0, 0.5, 0.0, 0.0, 0.0];
    near("sin'", ode::jsin(a)[1], 0.5 * cos(2.0), 1.0e-12);
    near("cos'", ode::jcos(a)[1], -0.5 * sin(2.0), 1.0e-12);
    near("sin''/2", ode::jsin(a)[2], -0.25 * sin(2.0) / 2.0, 1.0e-12);
    near("exp'''/6", ode::jexp(a)[3], 0.125 * exp(2.0) / 6.0, 1.0e-12);
    near("sqrt'", ode::jsqrt(a)[1], 0.25 / sqrt(2.0), 1.0e-12);
    near("ln'", ode::jlog(a)[1], 0.25, 1.0e-12);
    near("1/a''/2", ode::jdiv(ode::jconst(1.0), a)[2], 0.25 / 8.0, 1.0e-12);
    near("a^3'", ode::jpow(a, 3)[1], 3.0 * 4.0 * 0.5, 1.0e-12);
    near("ddt(a^2)", ode::jddt(ode::jpow(a, 2))[0], 2.0, 1.0e-12);
    let quartic: f64[5] = [0.0, 0.0, 0.0, 0.0, 1.0];
    near("four ddt of t^4", ode::jddt(ode::jddt(ode::jddt(ode::jddt(quartic))))[0], 24.0, 1.0e-12);
}

fn decay(t: f64) { ode::resid[0] = ode::jadd(ode::jddt(ode::uj[0]), ode::uj[0])[0]; }

fn test_decay() {
    ode::begin(1, 0.001, 6, 1000);
    ode::initial(0, 1.0);
    for n in 1..1001 {
        let ok = ode::advance(n) { t => decay(t); };
        if !ok { fail("decay does not converge"); return; }
    }
    near("exp(-1)", ode::value(0), exp(-1.0), 1.0e-3);
}

fn test_second_derivative() {
    ode::begin(2, 0.001, 6, 500);
    for n in 1..501 {
        let ok = ode::advance(n) { t => residual(t); };
        if !ok { fail("second derivative does not converge"); return; }
    }
    near("x", ode::value(0), sin(1.5), 1.0e-9);
    near("a", ode::value(1), -9.0 * sin(1.5), 1.0e-6);
}

fn residual(t: f64) {
    ode::resid[0] = ode::jsub(ode::uj[0], ode::jsin(ode::jmul(ode::jconst(3.0), ode::tj)))[0];
    ode::resid[1] = ode::jsub(ode::uj[1], ode::jddt(ode::jddt(ode::uj[0])))[0];
}

fn main() {
    test_jets();
    test_decay();
    test_second_derivative();
    if failures > 0 { exit(1); }
    print("ode: all tests pass");
}
