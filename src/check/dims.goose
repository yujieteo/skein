// Dimensions and units (spec §6). A dimension is a vector of seven integer
// exponents over M, L, T, Θ, I, N and J (spec §6.1).
import std;
namespace check;

let BASE_DIMS = ["M", "L", "T", "Θ", "I", "N", "J"];

// The named dimension types of spec §6.2. `Stress` is `Pressure`.
let DIM_NAMES = [
    "Length", "Area", "Volume", "Mass", "Time", "Temperature", "Force", "Power",
    "KinViscosity", "SpecificHeat", "HeatFlux", "Current", "Velocity", "Acceleration",
    "Frequency", "Density", "Pressure", "Stress", "Energy", "DynViscosity", "Conductivity",
    "Expansion", "HeatTransfer", "Amount",
];
let DIM_VALUES: i64[7][24] = [
    [0, 1, 0, 0, 0, 0, 0], [0, 2, 0, 0, 0, 0, 0], [0, 3, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0],
    [1, 1, -2, 0, 0, 0, 0], [1, 2, -3, 0, 0, 0, 0], [0, 2, -1, 0, 0, 0, 0],
    [0, 2, -2, -1, 0, 0, 0], [1, 0, -3, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0],
    [0, 1, -1, 0, 0, 0, 0], [0, 1, -2, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0],
    [1, -3, 0, 0, 0, 0, 0], [1, -1, -2, 0, 0, 0, 0], [1, -1, -2, 0, 0, 0, 0],
    [1, 2, -2, 0, 0, 0, 0], [1, -1, -1, 0, 0, 0, 0], [1, 1, -3, -1, 0, 0, 0],
    [0, 0, 0, -1, 0, 0, 0], [1, 0, -3, -1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0],
];

let NO_DIM: i64[7] = [0, 0, 0, 0, 0, 0, 0];
// In a `nondim` model each quantity is `Real` (spec §6.5): set_nondim
// makes these dimensions empty.
var LENGTH: i64[7] = [0, 1, 0, 0, 0, 0, 0];
var TIME: i64[7] = [0, 0, 1, 0, 0, 0, 0];
var FORCE: i64[7] = [1, 1, -2, 0, 0, 0, 0];

fn set_nondim(on: bool) {
    LENGTH = if on { NO_DIM } else { [0, 1, 0, 0, 0, 0, 0] };
    TIME = if on { NO_DIM } else { [0, 0, 1, 0, 0, 0, 0] };
    FORCE = if on { NO_DIM } else { [1, 1, -2, 0, 0, 0, 0] };
}

fn is_dimless(d: i64[7]) -> bool { d == NO_DIM }

// The dimension of a named dimension type, and true, or false.
fn named_dim(name: u8[:]) -> i64[7], bool {
    for n, i in DIM_NAMES {
        if n == name { return DIM_VALUES[i], true; }
    }
    return NO_DIM, false;
}

fn dim_pow(d: i64[7], k: i64) -> i64[7] {
    var r = d;
    for i in 7 { r[i] = d[i] * k; }
    return r;
}

// Writes a dimension: its name when it has one, else a product of base
// dimensions such as `M L^-1 T^-1`, or `1`.
fn format_dim(out: u8[>..]&, d: i64[7]) {
    for n, i in DIM_NAMES {
        if DIM_VALUES[i] == d { out.append(n); return; }
    }
    if is_dimless(d) { out.append("1"); return; }
    var first = true;
    for i in 7 {
        if d[i] == 0 { continue; }
        if !first { out.push(' '); }
        first = false;
        out.append(BASE_DIMS[i]);
        if d[i] != 1 { format(out, "^", d[i]); }
    }
}

// --- units (spec §6.3) ------------------------------------------------------------

let UNIT_NAMES = ["m", "g", "s", "K", "A", "mol", "cd", "N", "Pa", "J", "W", "Hz", "C", "V",
                  "rad", "deg", "L", "t", "min", "h", "bar", "degC"];
let UNIT_DIMS: i64[7][22] = [
    [0, 1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 1], [1, 1, -2, 0, 0, 0, 0], [1, -1, -2, 0, 0, 0, 0],
    [1, 2, -2, 0, 0, 0, 0], [1, 2, -3, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0],
    [0, 0, 1, 0, 1, 0, 0], [1, 2, -3, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0], [0, 3, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0], [1, -1, -2, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0],
];
let UNIT_SCALES: f64[22] = [1.0, 0.001, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                            1.0, 1.0, 1.0, 0.017453292519943295, 0.001, 1000.0, 60.0,
                            3600.0, 100000.0, 1.0];
// The units that take an SI prefix (spec §6.3).
let PREFIXED_UNITS = ["m", "g", "s", "N", "Pa", "J", "W", "Hz", "A", "V", "L"];
let PREFIX_CHARS = "numckMG";
let PREFIX_SCALES: f64[7] = [1e-9, 1e-6, 1e-3, 1e-2, 1e3, 1e6, 1e9];

struct Unit {
    dim: i64[7],
    scale: f64,     // the factor to SI
    angle: bool,    // `rad` or `deg`
    affine: bool,   // `degC`
}

// The unit with name `name`, and true, or false.
fn find_unit(name: u8[:]) -> Unit, bool {
    for u, i in UNIT_NAMES {
        if u == name {
            return Unit { dim: UNIT_DIMS[i], scale: UNIT_SCALES[i], angle: u == "rad" || u == "deg",
                          affine: u == "degC" }, true;
        }
    }
    let none = Unit { dim: NO_DIM, scale: 1.0, angle: false, affine: false };
    guard name.len > 1 else { return none, false; }
    let p = find(PREFIX_CHARS, name[..1]);
    guard p >= 0 else { return none, false; }
    guard contains_name(PREFIXED_UNITS, name[1..]) else { return none, false; }
    for u, i in UNIT_NAMES {
        if u == name[1..] {
            return Unit { dim: UNIT_DIMS[i], scale: UNIT_SCALES[i] * PREFIX_SCALES[p],
                          angle: false, affine: false }, true;
        }
    }
    return none, false;
}

// True when `name` is one of `names`.
fn contains_name(names: const u8[:][:], name: u8[:]) -> bool {
    for n in names {
        if n == name { return true; }
    }
    return false;
}
