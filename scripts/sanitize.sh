#!/bin/sh
# Builds the C of each src/<part>/test.goose with AddressSanitizer and
# UndefinedBehaviorSanitizer (-fno-sanitize-recover=all, so a report stops
# the test), and runs it, at the same time as the other tests. The runtime
# is part of the build (--standalone), so it is sanitized too. A report is a
# fault of the generator or of the runtime: fix that, not the generated C.
#
# Goose lays out a struct packed, without padding (Goose reference, structs),
# so every sanitizer report of type "misaligned address" comes from the layout
# of the generator and not from Skein. goose_packed_layout leaves out that one
# class of report. GOOSE_HOME is a new directory, so the sanitized runtime
# that the nested skein run builds is not in the cache of normal builds.
#
# Needs a C compiler (cc, or $SANITIZE_CC) with both sanitizers.
set -eu

root=$(cd "$(dirname "$0")/.." && pwd -P)
goose=${GOOSE:-goose}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

cc=${SANITIZE_CC:-cc}
goose_packed_layout=-fno-sanitize=alignment
export GOOSE_HOME="$work/goose_home"
cat >"$work/cc" <<WRAP
#!/bin/sh
exec "$cc" -fsanitize=undefined,address -fno-sanitize-recover=all -fno-omit-frame-pointer $goose_packed_layout "\$@"
WRAP
chmod +x "$work/cc"

one() {
    part=$1
    if ! GOOSE_CC="$work/cc" "$goose" build -q --debug --standalone \
        "$root/src/$part/test.goose" -o "$work/test_$part" >"$work/$part.log" 2>&1; then
        echo "FAIL sanitize $part: the test does not build"
        cat "$work/$part.log"
        return 1
    fi
    state=$(mktemp -d)
    warm=
    if [ "$part" = solvers ]; then warm=$state; fi
    if ! (cd "$root" && SKEIN_WARM="$warm" "$work/test_$part" >"$work/$part.log" 2>&1); then
        echo "FAIL sanitize $part: a test or a sanitizer report"
        tail -n 40 "$work/$part.log"
        rm -rf "$state"
        return 1
    fi
    rm -rf "$state"
}

parts=$(cd "$root/src" && ls -d */ | tr -d / | while read -r p; do [ -f "$root/src/$p/test.goose" ] && echo "$p"; done)
status=0
pids=""
for part in $parts; do
    one "$part" &
    pids="$pids $!"
done
for pid in $pids; do
    wait "$pid" || status=1
done
[ $status -eq 0 ] && echo "sanitize: ok ($(echo $parts | wc -w) test programs)"
exit $status
