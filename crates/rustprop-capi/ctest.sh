#!/usr/bin/env bash
#
# Build the C library and prove a C and a C++ program can actually use it.
#
# This is the check behind include/rustprop.h. The header is hand-written, so
# nothing guarantees it matches src/lib.rs except compiling a real translation
# unit against it and linking the real artifact — which is what happens below,
# once for the shared library and once for the static one. A signature that
# drifted shows up here as a compile or link error.
#
# Usage:  ctest.sh [FEATURES] [TARGET]
#
#   FEATURES  cargo feature list           (default: all-backends)
#   TARGET    target triple to build for   (default: the host's own)
#
# TARGET matters in CI, where the matrix builds an explicit triple and the
# artifacts land under target/<triple>/ rather than target/.
#
# Linux and macOS. Windows is covered by the equivalent MSVC steps in
# .github/workflows/release.yml — cl.exe against the import library — because
# the compiler driver and the artifact names differ enough that one script
# doing both would be less clear than two that each do one.
set -euo pipefail

FEATURES="${1:-all-backends}"
TARGET="${2:-}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CRATE="$ROOT/crates/rustprop-capi"
# An explicit --target changes where cargo puts the artifacts, so the two have
# to be decided together or the script looks in the wrong directory.
TARGET_ARG=""
TARGET_DIR=""
if [ -n "$TARGET" ]; then
    TARGET_ARG="--target $TARGET"
    TARGET_DIR="$TARGET/"
fi
# `release-capi` and not `release`: the shipped C library must unwind, so the
# catch_unwind at each entry point can turn a panic into a status instead of
# taking the host process down. See the root Cargo.toml.
PROFILE="release-capi"
OUT="$ROOT/target/${TARGET_DIR}$PROFILE"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

case "$(uname -s)" in
    Linux)  DYLIB="librustprop.so";    LD_VAR="LD_LIBRARY_PATH" ;;
    Darwin) DYLIB="librustprop.dylib"; LD_VAR="DYLD_LIBRARY_PATH" ;;
    *)      echo "ctest.sh covers Linux and macOS; on Windows see release.yml" >&2
            exit 2 ;;
esac

echo "==> building rustprop-capi (--features $FEATURES, --profile $PROFILE)"
cargo build -p rustprop-capi --features "$FEATURES" --profile "$PROFILE" $TARGET_ARG

[ -f "$OUT/$DYLIB" ]        || { echo "missing $OUT/$DYLIB" >&2; exit 1; }
[ -f "$OUT/librustprop.a" ] || { echo "missing $OUT/librustprop.a" >&2; exit 1; }

# The system libraries a Rust staticlib needs are a property of the toolchain
# and the target, not something to hardcode: ask rustc rather than guess, or
# this breaks the first time a platform differs.
echo "==> resolving native-static-libs"
NATIVE_LIBS="$(cargo rustc -p rustprop-capi --features "$FEATURES" \
    --profile "$PROFILE" $TARGET_ARG --crate-type staticlib -- \
    --print native-static-libs 2>&1 \
    | sed -n 's/^note: native-static-libs: //p' | tail -1)"
echo "    ${NATIVE_LIBS:-(none reported)}"

CC="${CC:-cc}"
CXX="${CXX:-c++}"
WARN="-Wall -Wextra -Werror"

run() {  # run <label> <binary>
    echo "--- $1"
    env "$LD_VAR=$OUT" "$2"
}

echo
echo "==> C, shared"
$CC $WARN -I "$CRATE/include" "$CRATE/examples/smoke.c" \
    -L "$OUT" -lrustprop -lm -o "$WORK/c_shared"
run "c_shared" "$WORK/c_shared"

echo
echo "==> C, static"
# RUSTPROP_STATIC suppresses the __declspec(dllimport) the header applies on
# Windows; harmless elsewhere, and defining it here keeps the two paths honest
# about which one they are.
$CC $WARN -DRUSTPROP_STATIC -I "$CRATE/include" "$CRATE/examples/smoke.c" \
    "$OUT/librustprop.a" -lm $NATIVE_LIBS -o "$WORK/c_static"
run "c_static" "$WORK/c_static"

echo
echo "==> C++, shared"
$CXX -std=c++17 $WARN -I "$CRATE/include" "$CRATE/examples/smoke.cc" \
    -L "$OUT" -lrustprop -o "$WORK/cc_shared"
run "cc_shared" "$WORK/cc_shared"

echo
echo "==> C++, static"
$CXX -std=c++17 $WARN -DRUSTPROP_STATIC -I "$CRATE/include" \
    "$CRATE/examples/smoke.cc" "$OUT/librustprop.a" $NATIVE_LIBS \
    -o "$WORK/cc_static"
run "cc_static" "$WORK/cc_static"

# A C consumer holding only the header must find every symbol it declares,
# whichever engines this build carries — that is the promise the header makes
# in "WHICH ENGINES DOES MY COPY HAVE?". Check it against the artifact rather
# than trusting the source.
echo
echo "==> every declared symbol is exported by $DYLIB"
missing=0
declared="$(sed -n 's/^RUSTPROP_API [^(]* \**\(rustprop_[a-z0-9_]*\)(.*/\1/p' \
    "$CRATE/include/rustprop.h" | sort -u)"
[ -n "$declared" ] || { echo "parsed no declarations out of the header" >&2; exit 1; }
if command -v nm >/dev/null 2>&1; then
    # `nm -D`, not `nm -g`: the release profile sets `strip = "symbols"`, which
    # empties the static symbol table and leaves only the DYNAMIC one — which
    # is the table a consumer links against anyway, so it is also the right
    # one to be asking. macOS `nm` has no -D; -gU reads the export table there.
    exported="$(nm -D --defined-only "$OUT/$DYLIB" 2>/dev/null \
        || nm -gU "$OUT/$DYLIB")"
    while read -r sym; do
        # Mach-O prefixes an underscore; match either spelling.
        if ! printf '%s\n' "$exported" | grep -qE "[[:space:]]_?${sym}$"; then
            echo "    MISSING: $sym"
            missing=1
        fi
    done <<< "$declared"
    [ "$missing" -eq 0 ] && echo "    all $(printf '%s\n' "$declared" | wc -l | tr -d ' ') declared symbols present"
else
    echo "    nm unavailable; skipped"
fi
[ "$missing" -eq 0 ] || exit 1

echo
echo "PASSED: C and C++, shared and static"
