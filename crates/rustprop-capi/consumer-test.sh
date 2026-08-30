#!/usr/bin/env bash
#
# Consume a packaged SDK the way a downstream project would.
#
# package.sh producing a plausible-looking directory tree proves nothing: the
# generated .pc and CMake config are the parts most likely to be wrong, and
# they are wrong in ways that only appear when something actually tries to
# build against them — a stale prefix, a missing static dependency, an
# IMPORTED target that CMake refuses. So this builds four real consumers out
# of the extracted tree, in a scratch directory, with no reference to the
# source repository at all.
#
# Usage:  consumer-test.sh <sdk-dir>
#
#   consumer-test.sh dist/rustprop-0.1.0-x86_64-unknown-linux-gnu
#
# Linux and macOS. Exits non-zero on the first failure.
set -euo pipefail

SDK="$(cd "${1:?usage: consumer-test.sh <sdk-dir>}" && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

case "$(uname -s)" in
    Linux)  LD_VAR="LD_LIBRARY_PATH" ;;
    Darwin) LD_VAR="DYLD_LIBRARY_PATH" ;;
    *)      echo "consumer-test.sh covers Linux and macOS" >&2; exit 2 ;;
esac

echo "==> consuming $SDK"
[ -f "$SDK/include/rustprop.h" ] || { echo "no header in the SDK" >&2; exit 1; }

cat > "$WORK/app.c" <<'EOF'
/* A downstream consumer: nothing but the header and the library. */
#include <stdio.h>
#include <math.h>
#include "rustprop.h"

int main(void) {
    printf("rustprop %s, CoolProp %s, backends [%s]\n", rustprop_version(),
           rustprop_upstream_version(), rustprop_backends());
    double d;
    if (rustprop_props_si("Dmolar", "T", 300, "P", 101325, "Water", &d)) {
        char msg[512];
        rustprop_last_error_message(msg, sizeof msg);
        fprintf(stderr, "FAILED: %s\n", msg);
        return 1;
    }
    if (fabs((d - 55317.35277350119) / d) > 1e-8) {
        fprintf(stderr, "FAILED: got %.17g\n", d);
        return 1;
    }
    printf("Dmolar(Water, 300 K, 101325 Pa) = %.15g mol/m^3  ok\n", d);
    return 0;
}
EOF

# --- 1. pkg-config, shared -------------------------------------------------
if command -v pkg-config >/dev/null 2>&1; then
    export PKG_CONFIG_PATH="$SDK/lib/pkgconfig"
    echo
    echo "--- pkg-config --modversion"
    pkg-config --modversion rustprop
    echo "--- pkg-config --cflags --libs"
    pkg-config --cflags --libs rustprop

    echo "--- building against pkg-config (shared)"
    # shellcheck disable=SC2046
    cc "$WORK/app.c" $(pkg-config --cflags --libs rustprop) -lm -o "$WORK/app_pc"
    env "$LD_VAR=$SDK/lib" "$WORK/app_pc"

    # The static path is the one that breaks when Libs.private is wrong, so it
    # is worth its own check rather than assuming the shared result carries.
    echo "--- building against pkg-config (static)"
    # shellcheck disable=SC2046
    cc -DRUSTPROP_STATIC "$WORK/app.c" -I "$SDK/include" \
        "$SDK/lib/librustprop.a" $(pkg-config --libs --static rustprop | sed 's/-lrustprop//') \
        -lm -o "$WORK/app_pc_static"
    "$WORK/app_pc_static"
else
    echo "!! pkg-config not installed; skipping that half" >&2
fi

# --- 2. CMake find_package -------------------------------------------------
if command -v cmake >/dev/null 2>&1; then
    mkdir -p "$WORK/cm"
    cp "$WORK/app.c" "$WORK/cm/app.c"
    cat > "$WORK/cm/CMakeLists.txt" <<'EOF'
cmake_minimum_required(VERSION 3.16)
project(rustprop_consumer C)
find_package(rustprop REQUIRED)

add_executable(app_shared app.c)
target_link_libraries(app_shared PRIVATE rustprop::rustprop m)

add_executable(app_static app.c)
target_link_libraries(app_static PRIVATE rustprop::rustprop_static m)
EOF
    echo
    echo "--- cmake configure"
    cmake -S "$WORK/cm" -B "$WORK/cm/build" \
        -DCMAKE_PREFIX_PATH="$SDK" -DCMAKE_BUILD_TYPE=Release >/dev/null
    echo "--- cmake build"
    cmake --build "$WORK/cm/build" >/dev/null
    echo "--- running the CMake consumers"
    env "$LD_VAR=$SDK/lib" "$WORK/cm/build/app_shared"
    "$WORK/cm/build/app_static"

    # A version request the package cannot satisfy must be REFUSED. Without
    # this, rustprop-config-version.cmake could be silently accepting anything
    # and nobody would notice until a real incompatibility shipped.
    echo "--- cmake version gate"
    cat > "$WORK/cm/CMakeLists.txt" <<'EOF'
cmake_minimum_required(VERSION 3.16)
project(rustprop_version_gate C)
find_package(rustprop 99.0 REQUIRED)
EOF
    if cmake -S "$WORK/cm" -B "$WORK/cm/vbuild" \
            -DCMAKE_PREFIX_PATH="$SDK" >/dev/null 2>&1; then
        echo "FAILED: find_package(rustprop 99.0) should not have succeeded" >&2
        exit 1
    fi
    echo "    find_package(rustprop 99.0) correctly refused"
else
    echo "!! cmake not installed; skipping that half" >&2
fi

# --- 3. the shipped example still builds from the SDK alone ----------------
echo
echo "--- the SDK's own example"
cc -I "$SDK/include" "$SDK/share/rustprop/examples/smoke.c" \
    -L "$SDK/lib" -lrustprop -lm -o "$WORK/smoke"
env "$LD_VAR=$SDK/lib" "$WORK/smoke" | tail -3

echo
echo "PASSED: pkg-config and CMake consumers build and run against the SDK"
