#!/bin/sh
# Process replay of the outputs of skein run for each model of
# scripts/generated_c.sh, with the skein of the base
# and with the skein of this checkout, and compares the two outputs:
#
#   the generated C                  must be the same
#   the solver results (skein run)   must agree within the tolerance that the
#                                    model states (its `expect` and verify
#                                    limits). Seconds and memory peaks are
#                                    left out.
#
# The base is origin/main (main when there is no origin). The base is
# extracted with git archive, so no worktree is made. A model that is not in
# the base is new and is not compared.
#
#   scripts/replay.sh            fail on a difference; print it
#   scripts/replay.sh --accept   print the difference of the generated C and
#                                pass, once. Use it when the change must change
#                                the C. A result beyond the tolerance always
#                                fails.
#   REPLAY_ACCEPT=1 ./build.sh test   the same for the whole test run.
#
# Needs build/skein (./build.sh), and goose on PATH or in $GOOSE.
set -eu

root=$(cd "$(dirname "$0")/.." && pwd -P)
goose=${GOOSE:-goose}

prepare() {
    tree=$1 model=$2
    work=$(mktemp -d)
    for d in examples samples; do ln -s "$tree/$d" "$work/$d"; done
    name=$(sed -n 's/^model \([A-Za-z0-9_]*\):.*/\1/p' "$tree/$model" | head -n 1)
}

capture_c() {
    tree=$1 out=$2 model=$3
    prepare "$tree" "$model"
    mkdir -p "$work/shim" "$out/$(dirname "$model")"
    printf '#!/bin/sh\nexit 1\n' >"$work/shim/goose"
    chmod +x "$work/shim/goose"
    (cd "$work" && PATH="$work/shim:$PATH" "$tree/build/skein" run "$model" >/dev/null 2>&1) || true
    if (cd "$work" && "$goose" -q -o "$work/$name.c" ".skein/build/$name/main.goose") 2>"$work/err"; then
        sed -e "s|\"$tree/|\"|" -e 's|"[^"]*/stdlib/|"stdlib/|' "$work/$name.c" >"$out/${model%.skn}.c"
    else
        cp "$work/err" "$out/${model%.skn}.fail"
    fi
    rm -rf "$work"
}

capture_run() {
    tree=$1 out=$2 model=$3
    prepare "$tree" "$model"
    if (cd "$work" && "$tree/build/skein" run "$model" >"$work/run.toon" 2>"$work/run.err"); then
        cp "$work/run.toon" "$out/${model%.skn}.toon"
    else
        cp "$work/run.err" "$out/${model%.skn}.runfail"
    fi
    rm -rf "$work"
}

case "${1:-}" in
    capture_c) shift; capture_c "$@"; exit 0 ;;
    capture_run) shift; capture_run "$@"; exit 0 ;;
esac

accept=
[ -z "${REPLAY_ACCEPT:-}" ] || accept=--accept
case "${1:-}" in
    "") ;;
    --accept) accept=--accept ;;
    *) echo "usage: scripts/replay.sh [--accept]" >&2; exit 2 ;;
esac

base=origin/main
git -C "$root" rev-parse -q --verify "$base" >/dev/null || base=main
if git -C "$root" diff --quiet "$base" --; then
    echo "replay: this checkout is the base ($base): nothing to compare"
    exit 0
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
mkdir -p "$work/base/build" "$work/base_out" "$work/head_out"
git -C "$root" archive "$base" | tar -x -C "$work/base"
"$goose" build -q "$work/base/src/skein.goose" -o "$work/base/build/skein"

models=$(sed -n '/^models="$/,/^"$/p' "$root/scripts/generated_c.sh" | grep '\.skn$')
in_parallel() {
    step=$1
    shift
    for model in "$@"; do
        echo "$work/base $work/base_out $model"
        echo "$root $work/head_out $model"
    done | xargs -P "$(nproc)" -L 1 sh -c '"$0" "$1" "$2" "$3" "$4"' "$root/scripts/replay.sh" "$step" 
}

# The program of a model fixes its results, so the models run only when
# their generated C differs.
in_parallel capture_c $models
changed=$(python3 -I "$root/scripts/guardrails/replay_diff.py" --changed "$work/base_out" "$work/head_out")
if [ -n "$changed" ]; then
    # shellcheck disable=SC2086
    in_parallel capture_run $changed
fi

exec python3 -I "$root/scripts/guardrails/replay_diff.py" $accept "$work/base_out" "$work/head_out"
