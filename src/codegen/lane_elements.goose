import std;
namespace codegen;

fn elements_owns_template(word: u8[:]) -> bool { return word == "heat"; }

fn elements_read_template(word: u8[:], n: i64) -> bool { return false; }

fn elements_emit_template(out: u8[>..]&) {}

fn elements_cell(cell: u8[:], order: i64) -> bool { return false; }

fn elements_emit_mesh(out: u8[>..]&) -> bool { return false; }
