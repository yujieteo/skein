// Unit names (spec §6.3). The parser needs them to find the end of a unit
// expression after a literal.
import std;
namespace syntax;

let UNITS = ["m", "g", "s", "K", "A", "mol", "cd", "N", "Pa", "J", "W", "Hz", "C", "V",
             "rad", "deg", "L", "t", "min", "h", "bar", "degC"];
// The units that take an SI prefix.
let PREFIXED = ["m", "g", "s", "N", "Pa", "J", "W", "Hz", "A", "V", "L"];
let PREFIXES = "numckMG";

// True when `word` is a unit name: a unit, or an SI prefix and a unit that
// takes one.
fn is_unit_name(word: u8[:]) -> bool {
    for u in UNITS {
        if u == word { return true; }
    }
    guard word.len > 1 && contains(PREFIXES, word[0]) else { return false; }
    for u in PREFIXED {
        if u == word[1..] { return true; }
    }
    return false;
}
