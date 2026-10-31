// The meaning of a type written in the source (spec §3 `type`, §5, §6.2).
import std;
namespace check;
import syntax from "../syntax";
import .dims;
import .types;
import .state;

let TYPE_NAMES = ["Real", "Int", "Angle", "Bool", "Str", "Solid", "Surface", "Curve", "Point",
                  "Mesh", "Points", "Results"];

// The type of the type node `n`. With `report`, an unknown name is an error.
recursive fn type_of(n: i64, report: bool) -> Ty {
    let k = syntax::nodes[n];
    match k.kind {
        syntax::N_TYPE_NAME => {
            let name = syntax::text(k.tok);
            if name == "Real" { return base_ty(B_REAL); }
            if name == "Int" { return base_ty(B_INT); }
            if name == "Angle" { return base_ty(B_ANGLE); }
            if name == "Bool" { return kind_ty(K_BOOL); }
            if name == "Str" { return kind_ty(K_STR); }
            if name == "Point" { return geom_ty(0); }
            if name == "Curve" { return geom_ty(1); }
            if name == "Surface" { return geom_ty(2); }
            if name == "Solid" { return geom_ty(3); }
            if name == "Mesh" { return kind_ty(K_MESH); }
            if name == "Points" { return kind_ty(K_POINTS); }
            if name == "Results" { return kind_ty(K_RESULTS); }
            let d, ok = named_dim(name);
            if ok { return num_ty(if nondim { NO_DIM } else { d }, 0); }
            if report { unknown_type(k.tok); }
            return any_ty();
        },
        syntax::N_TYPE_APPLY => {
            let name = syntax::text(k.tok);
            var t = type_of(k.first, report);
            if t.kind != K_NUM { return t; }
            if name == "Vector" || name == "Tensor2" || name == "Tensor4" {
                if t.rank != 0 && report { nested_tensor(k.tok); }
                t.rank = if name == "Vector" { 1 } else if name == "Tensor2" { 2 } else { 4 };
                return t;
            }
            t.num = if name == "Complex" { NS_COMPLEX } else if name == "Quaternion" { NS_QUATERNION } else { NS_DUAL };
            return t;
        },
        syntax::N_DIM_MUL, syntax::N_DIM_DIV => {
            let a = type_of(k.first, report);
            let b = type_of(syntax::nodes[k.first].next, report);
            guard a.kind == K_NUM && b.kind == K_NUM else { return any_ty(); }
            var t = num_ty(if k.kind == syntax::N_DIM_MUL { a.dim + b.dim } else { a.dim - b.dim }, 0);
            return t;
        },
        syntax::N_DIM_POW => {
            let a = type_of(k.first, report);
            guard a.kind == K_NUM else { return any_ty(); }
            // The base is a named type with whole exponents, so each half power is valid.
            let d, _ok = dim_pow_halves(a.dim, k.aux);
            return num_ty(d, 0);
        },
        _ => { return any_ty(); },
    }
}

fn unknown_type(tok: i64) {
    error(tok, "unknown_name", str("'", syntax::text(tok), "' is not a type"),
          "use Real, Int, Angle, Bool, Str, a geometry type or a dimension of spec §6.2");
}

fn nested_tensor(tok: i64) {
    error(tok, "type_mismatch", "a tensor type takes a scalar type", "write Vector<Length>, not Vector<Vector<Length>>");
}
