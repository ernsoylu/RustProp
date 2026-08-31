#!/usr/bin/env bash
#
# Run every language binding in examples/bindings/ against the built library.
#
# USAGE.md shows these bindings as the way to call rustprop from Python, Go,
# Java and Fortran. Documentation that is never executed drifts, and FFI
# documentation drifts silently and expensively — a wrong `argtypes` line in a
# ctypes example does not fail loudly, it passes garbage to a solver and
# returns a plausible number. So each one is a real program that checks its
# answers against known CoolProp values, and this runs them.
#
# Toolchains that are not installed are SKIPPED, not failed: a contributor
# without a JDK should still be able to run this, and CI reports which ones
# actually ran.
#
# Usage:  bindings-test.sh [FEATURES] [TARGET]
#
#   FEATURES  cargo feature list           (default: all-backends)
#   TARGET    target triple to build for   (default: the host's own)
#
# TARGET matters in CI, where naming it explicitly keeps this, ctest.sh and
# package.sh all pointing at one build directory instead of three.
set -euo pipefail

FEATURES="${1:-all-backends}"
TARGET="${2:-}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BINDINGS="$ROOT/crates/rustprop-capi/examples/bindings"
PROFILE="release-capi"
TARGET_ARG=""
TARGET_DIR=""
if [ -n "$TARGET" ]; then
    TARGET_ARG="--target $TARGET"
    TARGET_DIR="$TARGET/"
fi
OUT="$ROOT/target/${TARGET_DIR}$PROFILE"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

case "$(uname -s)" in
    Linux)  LIB="$OUT/librustprop.so";    LD_VAR="LD_LIBRARY_PATH" ;;
    Darwin) LIB="$OUT/librustprop.dylib"; LD_VAR="DYLD_LIBRARY_PATH" ;;
    *)      echo "bindings-test.sh covers Linux and macOS" >&2; exit 2 ;;
esac

echo "==> building rustprop-capi (--features $FEATURES)"
cargo build -p rustprop-capi --features "$FEATURES" --profile "$PROFILE" $TARGET_ARG
[ -f "$LIB" ] || { echo "missing $LIB" >&2; exit 1; }

ran=0
skipped=""

have() { command -v "$1" >/dev/null 2>&1; }

run_one() {  # run_one <label> <command...>
    local label="$1"; shift
    echo
    echo "--- $label"
    if "$@"; then
        ran=$((ran + 1))
    else
        echo "!! $label FAILED" >&2
        return 1
    fi
}

# --- Python (ctypes) -------------------------------------------------------
if have python3; then
    run_one "Python (ctypes)" env "RUSTPROP_LIB=$LIB" python3 "$BINDINGS/rustprop.py"
else
    skipped="$skipped python3"
fi

# --- Go (cgo) --------------------------------------------------------------
# cgo needs a module directory, so the example is copied into a scratch one.
if have go && have cc; then
    mkdir -p "$WORK/go"
    cp "$BINDINGS/rustprop.go" "$WORK/go/main.go"
    ( cd "$WORK/go" && go mod init rustprop_bindings_test >/dev/null 2>&1 )
    run_one "Go (cgo)" env \
        "CGO_CFLAGS=-I$ROOT/crates/rustprop-capi/include" \
        "CGO_LDFLAGS=-L$OUT -lrustprop" \
        "$LD_VAR=$OUT" \
        sh -c "cd '$WORK/go' && go run main.go"
else
    skipped="$skipped go"
fi

# --- Java (Foreign Function & Memory API, 22+) -----------------------------
# `java <file>.java` compiles in memory, so no build step is needed. FFM is
# stable from 22; older JDKs would need --enable-preview and are not supported
# by this example.
if have java; then
    jver="$(java -version 2>&1 | head -1 | sed -n 's/.*"\([0-9][0-9]*\).*/\1/p')"
    if [ -n "$jver" ] && [ "$jver" -ge 22 ]; then
        run_one "Java (FFM)" env "RUSTPROP_LIB=$LIB" \
            java --enable-native-access=ALL-UNNAMED "$BINDINGS/Rustprop.java"
    else
        skipped="$skipped java(need>=22,have=${jver:-?})"
    fi
else
    skipped="$skipped java"
fi

# --- Fortran (iso_c_binding) -----------------------------------------------
if have gfortran; then
    # -J puts the generated .mod in the scratch directory; without it gfortran
    # drops one in the current working directory, which is the repo root in CI.
    gfortran -J "$WORK" "$BINDINGS/rustprop.f90" -o "$WORK/fortran_demo" \
        -L "$OUT" -lrustprop
    run_one "Fortran (iso_c_binding)" env "$LD_VAR=$OUT" "$WORK/fortran_demo"
else
    skipped="$skipped gfortran"
fi

echo
echo "==> $ran binding(s) ran and passed"
[ -n "$skipped" ] && echo "    skipped (toolchain absent):$skipped"
[ "$ran" -gt 0 ] || { echo "no binding ran at all" >&2; exit 1; }
echo "PASSED"
