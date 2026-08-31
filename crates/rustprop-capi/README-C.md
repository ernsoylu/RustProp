# rustprop — C/C++ SDK

Thermophysical properties with [CoolProp 8.0.0](https://github.com/CoolProp/CoolProp)
semantics, implemented in pure Rust and validated against the upstream oracle
(41,629 records across 123 fixtures). No runtime, no data files to install, no
`REFPROP` — one library and one header.

Upstream project: <https://github.com/ernsoylu/RustProp>

## What is in this tree

```
include/rustprop.h                     the whole API, documented
lib/librustprop.{so,dylib,a}           shared and static  (rustprop.dll on Windows)
lib/pkgconfig/rustprop.pc              for pkg-config
lib/cmake/rustprop/                    for find_package(rustprop)
share/rustprop/examples/smoke.c        a worked example of every call
share/rustprop/examples/smoke.cc       the same, plus an idiomatic C++ wrapper
BUILD-INFO.txt                         which target, which engines, which toolchain
```

The tree is relocatable — extract it anywhere. Neither the `.pc` nor the CMake
config contains an absolute path.

`BUILD-INFO.txt` names the linkage this particular package offers. Some
targets are **static-only** — Rust emits no `cdylib` where the target defaults
to `crt-static`, which is true of every musl target. There is no shared
library in those packages; use `rustprop::rustprop_static` (CMake) or link
`lib/librustprop.a` directly.

## Three ways to build against it

**CMake** (recommended):

```cmake
find_package(rustprop REQUIRED)
target_link_libraries(myapp PRIVATE rustprop::rustprop)         # shared
# target_link_libraries(myapp PRIVATE rustprop::rustprop_static)  # static
```

```bash
cmake -B build -DCMAKE_PREFIX_PATH=/path/to/rustprop-0.1.0-<target>
```

**pkg-config**:

```bash
export PKG_CONFIG_PATH=/path/to/rustprop-0.1.0-<target>/lib/pkgconfig
cc myapp.c $(pkg-config --cflags --libs rustprop) -o myapp
```

**By hand**:

```bash
cc -I <sdk>/include myapp.c -L <sdk>/lib -lrustprop -o myapp
cc -I <sdk>/include -DRUSTPROP_STATIC myapp.c <sdk>/lib/librustprop.a \
   $(pkg-config --libs --static rustprop) -o myapp     # static
```

Static linking needs the system libraries the Rust runtime uses; they are in
the `.pc` file's `Libs.private` and in the CMake static target, so prefer
either of those over spelling them out.

## Hello, water

```c
#include <stdio.h>
#include "rustprop.h"

int main(void) {
    double d;
    int rc = rustprop_props_si("Dmolar", "T", 300, "P", 101325, "Water", &d);
    if (rc != RUSTPROP_OK) {
        char msg[512];
        rustprop_last_error_message(msg, sizeof msg);
        fprintf(stderr, "rustprop: %s\n", msg);
        return 1;
    }
    printf("%.15g mol/m^3\n", d);   /* 55317.3527735012 */
    return 0;
}
```

Fluid names take the upstream spellings, backend prefix included: `Water`,
`R134a`, `IF97::Water`, `SRK::Propane`, `HEOS::Methane&Ethane`.

## Things worth knowing before you start

**Engines are chosen at compile time.** rustprop exists to be small — a
build carrying only IAPWS-IF97 is 128 KB, one carrying everything is 4.2 MB.
The release artifacts named `all-backends` carry every engine and all 130
fluids; if you need a smaller one, build from source with the features you
want. Either way, ask the library what it has:

```c
printf("%s\n", rustprop_backends());          /* "heos,if97,humid-air" */
if (!rustprop_has_backend("pcsaft")) { /* ... */ }
for (size_t i = 0; i < rustprop_fluid_count(); i++)
    puts(rustprop_fluid_name(i));
```

**Every symbol exists in every build.** A call into an engine your copy does
not carry returns `RUSTPROP_UNAVAILABLE` (102). It never fails to link, so one
header and one set of call sites work against any build.

**Errors are statuses, not sentinels.** `RUSTPROP_OK` (0) means `*out` was
written; nothing else writes it. `rustprop_last_error_message()` explains what
happened, and the codes 1–11 map one-to-one onto CoolProp's exception types,
because this port reproduces upstream's refusals as well as its numbers.

**Threads are fine.** Every function is safe to call concurrently, with no
initialisation call and no lock on your side. The last-error slot is
thread-local. There is no global state to configure and nothing to shut down.

**Nothing needs freeing.** No function returns memory you own. `const char *`
returns live as long as the process; everything else goes into buffers you
supply.

**Batches beat loops.** `rustprop_props_si_many` parses the parameter names
once instead of once per state, which dominates a large sweep. A cell that
fails becomes `NaN` rather than failing the call, so clipping the phase
envelope costs you that cell and nothing else.

## Checking your copy

`share/rustprop/examples/smoke.c` exercises every entry point and verifies the
results against known CoolProp values. It is the same program CI runs against
each released artifact:

```bash
cc -I include share/rustprop/examples/smoke.c -L lib -lrustprop -lm -o smoke
LD_LIBRARY_PATH=lib ./smoke        # DYLD_LIBRARY_PATH on macOS
```

## Numbers, and which build you have

Every released artifact is checked against the CoolProp 8.0.0 oracle before it
is published, and `BUILD-INFO.txt` records the target and the microarchitecture
baseline it was compiled for. Artifacts with a `-v2`/`-v3`/`-v4` suffix are
built for newer x86-64 baselines and will not run on older processors; the
unsuffixed one runs anywhere.

## License

MIT. Derivative work of CoolProp (MIT, © 2012–2018 Ian H. Bell and other
CoolProp developers) — see `LICENSE`.
