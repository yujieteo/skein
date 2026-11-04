import std;
namespace codegen;

fn shells_owns_template(word: u8[:]) -> bool { return false; }

fn shells_read_template(word: u8[:], n: i64) -> bool { return false; }

fn shells_emit_template(out: u8[>..]&) {}

fn shells_read_geometry(n: i64, value: i64) -> bool { return false; }

fn shells_surface(n: i64, value: i64, set: u8[..64]&) -> bool { return false; }

fn shells_emit_mesh(out: u8[>..]&) {}

fn shells_cell(cell: u8[:], order: i64) -> bool { return false; }
