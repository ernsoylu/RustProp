#!/usr/bin/env bash
#
# Assemble the distributable C/C++ SDK for one target triple.
#
# Produces a relocatable tree — extract it anywhere and both pkg-config and
# CMake find it, because neither generated file contains an absolute path.
# That matters more than it sounds: a release artifact is unpacked by someone
# whose directory layout we cannot know, and a baked-in /home/runner/work
# prefix is the classic way a prebuilt SDK arrives broken.
#
#   rustprop-<version>-<target>/
#     include/rustprop.h
#     lib/librustprop.{so,dylib,a}          (rustprop.dll + .lib on Windows)
#     lib/pkgconfig/rustprop.pc
#     lib/cmake/rustprop/rustprop-config.cmake
#     lib/cmake/rustprop/rustprop-config-version.cmake
#     share/rustprop/examples/{smoke.c,smoke.cc}
#     BUILD-INFO.txt  LICENSE  README.md
#
# Usage:  package.sh <target-triple> <staging-dir> [features]
#
#   package.sh x86_64-unknown-linux-gnu dist all-backends
#
# The caller is responsible for having built the artifacts already; this only
# collects and describes them. `TARGET_CPU` may be set to record which
# microarchitecture baseline the binaries were compiled for.
set -euo pipefail

TARGET="${1:?usage: package.sh <target-triple> <staging-dir> [features]}"
STAGE="${2:?usage: package.sh <target-triple> <staging-dir> [features]}"
FEATURES="${3:-all-backends}"
PROFILE="release-capi"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CRATE="$ROOT/crates/rustprop-capi"
BUILT="$ROOT/target/$TARGET/$PROFILE"

VERSION="$(sed -n 's/^version *= *"\(.*\)"/\1/p' "$ROOT/Cargo.toml" | head -1)"
[ -n "$VERSION" ] || { echo "could not read version from Cargo.toml" >&2; exit 1; }

NAME="rustprop-${VERSION}-${TARGET}${TARGET_CPU:+-${TARGET_CPU}}"
DEST="$STAGE/$NAME"

rm -rf "$DEST"
mkdir -p "$DEST/include" "$DEST/lib/pkgconfig" "$DEST/lib/cmake/rustprop" \
         "$DEST/share/rustprop/examples"

# --- the artifacts ---------------------------------------------------------
# Which files exist is target-dependent, so copy what is there and insist that
# SOMETHING was: a silent empty package is the worst outcome here.
copied=0
for f in librustprop.so librustprop.dylib librustprop.a \
         rustprop.dll rustprop.dll.lib rustprop.lib; do
    if [ -f "$BUILT/$f" ]; then
        cp "$BUILT/$f" "$DEST/lib/$f"
        copied=$((copied + 1))
    fi
done
[ "$copied" -gt 0 ] || {
    echo "no library artifacts found in $BUILT" >&2
    echo "did you build with --target $TARGET --profile $PROFILE ?" >&2
    exit 1
}
# Which linkages this package actually offers. Not every target has both:
# rustc drops the `cdylib` crate type wherever the target defaults to
# `crt-static`, which is true of every musl target. Say so in BUILD-INFO.txt
# rather than letting a consumer discover a missing .so on their own. (No musl
# target is in the release matrix today — see release.yml — but this stays
# because the condition is a property of the target, not of that decision.)
LINKAGE="static and shared"
if [ ! -f "$DEST/lib/librustprop.so" ] && [ ! -f "$DEST/lib/librustprop.dylib" ] \
   && [ ! -f "$DEST/lib/rustprop.dll" ]; then
    LINKAGE="STATIC ONLY (no shared library for this target)"
fi

# A DLL belongs next to the executable, not in lib/; keep a copy in bin/ so
# the layout works either way on Windows.
if [ -f "$DEST/lib/rustprop.dll" ]; then
    mkdir -p "$DEST/bin"
    cp "$DEST/lib/rustprop.dll" "$DEST/bin/rustprop.dll"
fi

cp "$CRATE/include/rustprop.h"  "$DEST/include/"
cp "$CRATE/examples/smoke.c"    "$DEST/share/rustprop/examples/"
cp "$CRATE/examples/smoke.cc"   "$DEST/share/rustprop/examples/"
cp "$ROOT/LICENSE"              "$DEST/"

# --- pkg-config ------------------------------------------------------------
# ${pcfiledir} is what makes this relocatable: pkg-config expands it to the
# directory the .pc was found in, so the prefix follows the extracted tree.
cat > "$DEST/lib/pkgconfig/rustprop.pc" <<EOF
prefix=\${pcfiledir}/../..
exec_prefix=\${prefix}
libdir=\${prefix}/lib
includedir=\${prefix}/include

Name: rustprop
Description: Thermophysical properties (CoolProp 8.0.0 semantics), pure Rust
URL: https://github.com/ernsoylu/RustProp
Version: ${VERSION}
Cflags: -I\${includedir}
Libs: -L\${libdir} -lrustprop
Libs.private: ${RUSTPROP_NATIVE_LIBS:-}
EOF

# --- CMake package config --------------------------------------------------
# IMPORTED targets rather than plain variables, so a consumer writes
# `target_link_libraries(app PRIVATE rustprop::rustprop)` and inherits the
# include directory with it. Both linkages are offered; the shared one is the
# default alias because that is what most consumers of a prebuilt SDK want.
cat > "$DEST/lib/cmake/rustprop/rustprop-config.cmake" <<'EOF'
# rustprop — CMake package configuration.
#
#   find_package(rustprop REQUIRED)
#   target_link_libraries(myapp PRIVATE rustprop::rustprop)          # shared
#   target_link_libraries(myapp PRIVATE rustprop::rustprop_static)   # static
#
# Relocatable: the prefix is derived from this file's own location, so the
# extracted tree can live anywhere.

get_filename_component(_rustprop_prefix "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE)
set(rustprop_INCLUDE_DIRS "${_rustprop_prefix}/include")

include(CMakeFindDependencyMacro)
find_package(Threads QUIET)

# --- shared -----------------------------------------------------------------
foreach(_cand librustprop.so librustprop.dylib rustprop.dll)
  if(EXISTS "${_rustprop_prefix}/lib/${_cand}")
    set(_rustprop_shared "${_rustprop_prefix}/lib/${_cand}")
    break()
  endif()
endforeach()

if(_rustprop_shared AND NOT TARGET rustprop::rustprop)
  add_library(rustprop::rustprop SHARED IMPORTED)
  set_target_properties(rustprop::rustprop PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${rustprop_INCLUDE_DIRS}")
  if(WIN32)
    # On Windows the DLL is the runtime and the import library is what the
    # linker consumes; CMake needs both, and gets them wrong if only one is set.
    set_target_properties(rustprop::rustprop PROPERTIES
      IMPORTED_LOCATION "${_rustprop_shared}"
      IMPORTED_IMPLIB   "${_rustprop_prefix}/lib/rustprop.dll.lib")
  else()
    set_target_properties(rustprop::rustprop PROPERTIES
      IMPORTED_LOCATION "${_rustprop_shared}")
  endif()
endif()

# --- static -----------------------------------------------------------------
foreach(_cand librustprop.a rustprop.lib)
  if(EXISTS "${_rustprop_prefix}/lib/${_cand}")
    set(_rustprop_static "${_rustprop_prefix}/lib/${_cand}")
    break()
  endif()
endforeach()

if(_rustprop_static AND NOT TARGET rustprop::rustprop_static)
  add_library(rustprop::rustprop_static STATIC IMPORTED)
  set_target_properties(rustprop::rustprop_static PROPERTIES
    IMPORTED_LOCATION "${_rustprop_static}"
    INTERFACE_INCLUDE_DIRECTORIES "${rustprop_INCLUDE_DIRS}"
    # RUSTPROP_STATIC drops the __declspec(dllimport) the header would
    # otherwise apply on Windows.
    INTERFACE_COMPILE_DEFINITIONS "RUSTPROP_STATIC")
  # A Rust staticlib carries no record of the system libraries it needs, so
  # they are named here. The list came from `rustc --print native-static-libs`
  # at package time; without it, a static link fails at the very last step.
  set(_rustprop_static_deps "@RUSTPROP_STATIC_DEPS@")
  if(_rustprop_static_deps)
    set_target_properties(rustprop::rustprop_static PROPERTIES
      INTERFACE_LINK_LIBRARIES "${_rustprop_static_deps}")
  endif()
endif()

if(NOT _rustprop_shared AND NOT _rustprop_static)
  set(rustprop_FOUND FALSE)
  set(rustprop_NOT_FOUND_MESSAGE
      "no rustprop library found under ${_rustprop_prefix}/lib")
else()
  set(rustprop_FOUND TRUE)
endif()
EOF

# Substitute the static-link dependency list gathered by the caller. Semicolons
# are CMake's list separator, so the space-separated rustc output is converted.
DEPS="$(printf '%s' "${RUSTPROP_NATIVE_LIBS:-}" \
    | sed 's/^ *//; s/ *$//; s/-l//g; s/  */;/g')"
# BSD and GNU sed disagree about -i, so rewrite through a temp file.
sed "s|@RUSTPROP_STATIC_DEPS@|${DEPS}|" \
    "$DEST/lib/cmake/rustprop/rustprop-config.cmake" > "$DEST/.cfg.tmp"
mv "$DEST/.cfg.tmp" "$DEST/lib/cmake/rustprop/rustprop-config.cmake"

cat > "$DEST/lib/cmake/rustprop/rustprop-config-version.cmake" <<EOF
set(PACKAGE_VERSION "${VERSION}")

# rustprop is pre-1.0: treat the MINOR component as the compatibility
# boundary, which is what Cargo's own semver rules do for 0.x crates. Once
# 1.0 ships this becomes major-compatible.
if(PACKAGE_FIND_VERSION VERSION_GREATER PACKAGE_VERSION)
  set(PACKAGE_VERSION_COMPATIBLE FALSE)
else()
  string(REGEX MATCH "^([0-9]+)\\\\.([0-9]+)" _pv "\${PACKAGE_VERSION}")
  set(_pv_major "\${CMAKE_MATCH_1}")
  set(_pv_minor "\${CMAKE_MATCH_2}")
  string(REGEX MATCH "^([0-9]+)\\\\.([0-9]+)" _fv "\${PACKAGE_FIND_VERSION}")
  set(_fv_major "\${CMAKE_MATCH_1}")
  set(_fv_minor "\${CMAKE_MATCH_2}")
  if(_pv_major STREQUAL "0")
    if(_pv_major STREQUAL _fv_major AND _pv_minor STREQUAL _fv_minor)
      set(PACKAGE_VERSION_COMPATIBLE TRUE)
    else()
      set(PACKAGE_VERSION_COMPATIBLE FALSE)
    endif()
  elseif(_pv_major STREQUAL _fv_major)
    set(PACKAGE_VERSION_COMPATIBLE TRUE)
  else()
    set(PACKAGE_VERSION_COMPATIBLE FALSE)
  endif()
  if(PACKAGE_FIND_VERSION STREQUAL PACKAGE_VERSION)
    set(PACKAGE_VERSION_EXACT TRUE)
  endif()
endif()
EOF

# --- what this binary actually is ------------------------------------------
# A prebuilt library is opaque: engines and fluids were chosen at compile time
# and cannot be inspected from the outside without running it. This file, and
# rustprop_backends() at runtime, are the two ways to find out.
cat > "$DEST/BUILD-INFO.txt" <<EOF
rustprop ${VERSION} — C/C++ SDK
CoolProp semantics: 8.0.0

target            ${TARGET}
target-cpu        ${TARGET_CPU:-(default baseline for this target)}
cargo features    ${FEATURES}
cargo profile     ${PROFILE}  (release + panic=unwind)
linkage           ${LINKAGE}
built             $(date -u '+%Y-%m-%dT%H:%M:%SZ')
toolchain         $(rustc --version 2>/dev/null || echo unknown)

Engines and fluids are compile-time selections. Ask the library itself:

    rustprop_backends()      -> "heos,if97,..."
    rustprop_has_backend(x)  -> 1 / 0
    rustprop_fluid_count()   -> how many HEOS fluids are compiled in
    rustprop_fluid_name(i)   -> their names

Every function in include/rustprop.h is exported by every build. A call into
an engine this binary does not carry returns RUSTPROP_UNAVAILABLE (102); it
never fails to link.

Getting started: share/rustprop/examples/smoke.c is a worked example of every
call, and README.md has the three ways to build against this tree.
EOF

cp "$CRATE/README-C.md" "$DEST/README.md"

echo "$DEST"
