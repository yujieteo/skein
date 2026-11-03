#!/bin/sh
# The strictest checks that Goose offers, on all of src/: a warning of the
# type checker is an error (a redundant cast, a needless `&`, a deprecated
# extern), and each file is in the layout of goose format. pyright checks the
# Python of scripts/ in strict mode. The tests build with --debug
# (build.sh test): overflow, range and tag checks.
#
# Needs goose on PATH or in $GOOSE, and pyright (or nix).
set -eu

root=$(cd "$(dirname "$0")/.." && pwd -P)
goose=${GOOSE:-goose}
status=0

for entry in "$root"/src/skein.goose "$root"/src/*/test.goose; do
    out=$("$goose" check -q --diagnostics json "$entry" 2>&1)
    if [ "$out" != "[]" ]; then
        echo "FAIL strict: goose check reports warnings or errors for ${entry#"$root"/}"
        echo "$out" | sed -n 's/.*"text": "\([^"]*\)".*/\1/p'
        status=1
    fi
done
if ! "$goose" format --verify "$root/src" >/dev/null; then
    echo "FAIL strict: a file is not in the layout of goose format. Run goose format src."
    status=1
fi
if command -v pyright >/dev/null 2>&1; then
    pyright_run="pyright"
else
    pyright_run="nix shell nixpkgs#pyright --command pyright"
fi
if ! (cd "$root" && $pyright_run -p pyrightconfig.json >"${TMPDIR:-/tmp}/strict_pyright.$$" 2>&1); then
    echo "FAIL strict: pyright (strict mode, pyrightconfig.json) reports errors in scripts/"
    cat "${TMPDIR:-/tmp}/strict_pyright.$$"
    status=1
fi
rm -f "${TMPDIR:-/tmp}/strict_pyright.$$"
[ $status -eq 0 ] && echo "strict: ok"
exit $status
