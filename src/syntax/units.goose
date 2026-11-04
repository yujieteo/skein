import std;
namespace syntax;

let UNITS = ["m", "g", "s", "K", "A", "mol", "cd", "N", "Pa", "J", "W", "Hz", "C", "V",
             "rad", "deg", "L", "t", "min", "h", "bar", "degC"];
let PREFIXED = ["m", "g", "s", "N", "Pa", "J", "W", "Hz", "A", "V", "L"];
let PREFIXES = "numckMG";

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
