# Using rustprop

Thermophysical properties with [CoolProp 8.0.0](https://github.com/CoolProp/CoolProp)
semantics — `PropsSI` and `HAPropsSI` — in pure Rust, validated against the
upstream oracle (41,629 records across 123 fixtures).

This guide covers every way to call it. If you only want the API reference for
C, it is [`crates/rustprop-capi/include/rustprop.h`](crates/rustprop-capi/include/rustprop.h),
which is written to be read.

---

## 1. Pick your route

| You are writing | Route | Start at |
|---|---|---|
| Rust | the crates, from crates.io | [§3](#3-rust) |
| C | the C ABI + `rustprop.h` | [§4](#4-c) |
| C++ | the same, with a small wrapper | [§5](#5-c-1) |
| Python | `ctypes` on the shared library | [§6](#6-python) |
| Go | `cgo` | [§7](#7-go) |
| Java / Kotlin / Scala | the FFM API (JDK 22+) | [§8](#8-java) |
| Fortran | `iso_c_binding` | [§9](#9-fortran) |
| JavaScript / TypeScript | the WebAssembly bundle | [§10](#10-javascript-and-typescript) |
| C#, Julia, MATLAB, R, Ruby, … | any C FFI | [§11](#11-any-other-language) |
| Nothing — I just want a number | the CLI | [§12](#12-the-command-line) |

Everything except Rust and JavaScript goes through the same C ABI, so §4 is
worth skimming even if you are here for Python: the error model, the
build-time engine selection and the threading rules are identical everywhere.

---

## 2. Getting rustprop

### The idea you need first: engines are chosen at compile time

rustprop exists to be small. A build carrying only IAPWS-IF97 is **128 KB**;
one carrying every engine and all 130 fluids is **4.2 MB**. Which engines and
which fluids a binary contains is fixed when it is compiled — there is no
"load a fluid" call, by design.

That means:

- **From source (Rust, or building the C library yourself):** you choose, with
  Cargo features. Pay for exactly what you use.
- **From a release artifact:** you get `all-backends` — everything — because a
  prebuilt binary cannot be re-specialised after the fact. Ask it what it has
  at run time (`rustprop_backends()`), or build from source for a smaller one.

### Prebuilt artifacts

Attached to every [release](https://github.com/ernsoylu/RustProp/releases):

| Artifact | For |
|---|---|
| `rustprop-<ver>-<target>[-<baseline>].tar.gz` | The C/C++ SDK: libraries, header, pkg-config, CMake, examples, CLI |
| `rustprop-cli-<ver>-<target>.tar.gz` | Just the CLI |
| `rustprop-wasm-<preset>.tar.gz` | WebAssembly, as `web` / `nodejs` / `bundler` |
| `rustprop-<ver>-rust-sources.tar.gz` | Rust sources for offline and vendored builds |

Targets: Linux x86-64 (four instruction-set baselines), Linux arm64, Linux
armv7, macOS arm64 and x86-64, Windows x86-64 and arm64.

> **No musl artifact, for now.** It was in the matrix and was removed: musl's
> libm disagrees with glibc's in the `validity` golden suite — seven parity
> failures on `PR::Propane` at T = 1e30 K, by factors of exactly 2, 4 and 8.
> Every other suite passed, so the divergence is narrow and lives in the
> extreme-value tail, but this project does not ship numbers it has not
> checked. Re-adding musl means characterising where its libm diverges first.

### Which x86-64 artifact?

Four exist. They differ **only** in which processors they run on:

| Suffix | Needs | Runs on |
|---|---|---|
| *(none)* | baseline x86-64 | anything since 2003 |
| `-x86-64-v2` | SSE4.2, POPCNT | ~2009+ (Nehalem, Bulldozer) |
| `-x86-64-v3` | AVX2, FMA, BMI2 | ~2013+ (Haswell, Excavator) |
| `-x86-64-v4` | AVX-512 | some Xeon/Core; **not** most AMD |

**The numbers are identical across all four** — verified, not assumed: the
full test suite passes at every baseline, and a 29,848-value sweep across
every engine returns byte-identical results from each. So this is purely a
reach question. When in doubt, take the unsuffixed one.

`BUILD-INFO.txt` inside each SDK records its target, its baseline, whether it
was executed or only cross-built, and whether it offers shared linkage.

---

## 3. Rust

```bash
cargo add rustprop --features if97
```

Default features are empty. Pick engines explicitly:

| Feature | Engine |
|---|---|
| `heos` | Multiparameter Helmholtz EOS — **also add fluids**, see below |
| `heos-mixtures` | HEOS mixtures (binary-pair + departure data) |
| `if97` | IAPWS-IF97 water/steam, self-contained |
| `cubics` | SRK / Peng-Robinson, 116 fluids included |
| `incompressible` | Brines and secondary fluids, 126 included |
| `pcsaft` | PC-SAFT, 180 fluids included |
| `humid-air` | `HAPropsSI` psychrometrics |
| `tabular` | TTSE / bicubic tables (low-level API; pulls `heos`) |
| `svdsbtl` | SVD-compressed lookup (low-level API) |
| `all-backends` | Everything, all 130 fluids |

`heos` selects the *engine*; `rustprop-data` selects the *fluids*, because
per-fluid data dominates binary size:

```toml
[dependencies]
rustprop = { version = "0.1", features = ["heos"] }
rustprop-data = { version = "0.1", features = ["water", "r134a"] }
```

```rust
use rustprop::{props_si, ha_props_si};

fn main() -> rustprop::Result<()> {
    // PropsSI("Dmolar", "T", 300, "P", 101325, "Water")
    let d = props_si("Dmolar", "T", 300.0, "P", 101325.0, "Water")?;
    assert!((d - 55317.35277350119).abs() / d < 1e-8);

    let h = props_si("H", "T", 300.0, "P", 101325.0, "IF97::Water")?;
    let w = ha_props_si("W", "T", 300.0, "P", 101325.0, "R", 0.5)?;
    println!("{d} {h} {w}");
    Ok(())
}
```

Errors are `rustprop::Error`, one variant per CoolProp exception type — this
port reproduces upstream's *refusals* as faithfully as its numbers.

### Offline and vendored builds

Take `rustprop-<ver>-rust-sources.tar.gz`. Its `source/` directory is the
workspace with path dependencies intact:

```toml
[dependencies]
rustprop = { path = "vendor/rustprop-0.1.0-rust-sources/source/crates/rustprop",
             features = ["if97"] }
```

```bash
CARGO_NET_OFFLINE=true cargo build
```

The `crates/` directory in the same bundle holds the published `.crate` files
for archival and registry mirroring. Those do **not** build offline on their
own — `cargo package` rewrites path dependencies into registry dependencies,
so an unpacked crate goes looking for its siblings on crates.io. Use
`source/`, or `cargo vendor` on a networked machine.

### Why there is no prebuilt Rust library

Rust has no stable ABI, and an `.rlib` also encodes the compiler version and
the resolved feature set, so one built by CI would link only against an
identical `rustc`. If you want a prebuilt binary callable from Rust, use the C
ABI below — `extern "C"` is the one interface Rust keeps stable.

---

## 4. C

### Build against the SDK

**CMake** — recommended:

```cmake
find_package(rustprop REQUIRED)
target_link_libraries(myapp PRIVATE rustprop::rustprop)          # shared
# target_link_libraries(myapp PRIVATE rustprop::rustprop_static) # static
```

```bash
cmake -B build -DCMAKE_PREFIX_PATH=/path/to/rustprop-0.1.0-<target>
```

**pkg-config**:

```bash
export PKG_CONFIG_PATH=/path/to/rustprop-0.1.0-<target>/lib/pkgconfig
cc myapp.c $(pkg-config --cflags --libs rustprop) -lm -o myapp
```

**By hand**:

```bash
cc -I <sdk>/include myapp.c -L <sdk>/lib -lrustprop -lm -o myapp
LD_LIBRARY_PATH=<sdk>/lib ./myapp          # DYLD_LIBRARY_PATH on macOS
```

Static linking needs the system libraries the Rust runtime uses. Do not guess
them — they are in the `.pc` file's `Libs.private` and in the CMake static
target:

```bash
cc -DRUSTPROP_STATIC -I <sdk>/include myapp.c <sdk>/lib/librustprop.a \
   $(pkg-config --libs --static rustprop | sed 's/-lrustprop//') -lm -o myapp
```

### Hello, water

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

### The five rules

**1. Errors are statuses.** `RUSTPROP_OK` (0) means `*out` was written.
Nothing else writes it — your variable keeps whatever it held. Codes 1–11 map
one-to-one onto CoolProp's exception types; 100+ are boundary conditions
(NULL argument, bad UTF-8, engine absent).

**2. Every symbol exists in every build.** A call into an engine your copy was
not compiled with returns `RUSTPROP_UNAVAILABLE` (102). It never fails to
link, so one header and one set of call sites work against any build. Ask what
you have:

```c
printf("%s\n", rustprop_backends());        /* "heos,if97,humid-air" */
if (!rustprop_has_backend("pcsaft")) { /* degrade gracefully */ }
for (size_t i = 0; i < rustprop_fluid_count(); i++)
    puts(rustprop_fluid_name(i));
```

**3. Threads are fine.** Every function is safe to call concurrently, with no
initialisation call and no lock on your side. The last-error slot is
thread-local, so read it on the thread that failed, before that thread makes
another call.

**4. Nothing needs freeing.** No function returns memory you own. `const char *`
returns live for the life of the process; everything else goes into buffers
you supply.

**5. Batch, don't loop.** `rustprop_props_si_many` parses the parameter names
once instead of once per state, which dominates a large sweep:

```c
double t[1000], p[1000], out[1000];
/* ... fill t and p ... */
rustprop_props_si_many("Dmolar", "T", t, "P", p, "Water", 1000, out);
/* A state that fails is NaN in its slot; the batch does not abort. */
```

### Error messages, in full

`rustprop_last_error_message` returns the length **excluding** the NUL, so a
return `>= len` means truncation. The two-call idiom gets everything:

```c
size_t need = rustprop_last_error_message(NULL, 0);
char *msg = malloc(need + 1);
rustprop_last_error_message(msg, need + 1);
```

Fluid names take the upstream spellings, backend prefix included: `Water`,
`R134a`, `IF97::Water`, `SRK::Propane`, `HEOS::Methane&Ethane`.

---

## 5. C++

The header is `extern "C"`-guarded, so it just works. A small wrapper gives
you exceptions and `std::vector`; the full version is in
[`examples/smoke.cc`](crates/rustprop-capi/examples/smoke.cc).

```cpp
#include <stdexcept>
#include <string>
#include <vector>
#include "rustprop.h"

namespace rustprop {

inline std::string last_error() {
    std::string s(rustprop_last_error_message(nullptr, 0), '\0');
    if (!s.empty()) rustprop_last_error_message(&s[0], s.size() + 1);
    return s;
}

struct Error : std::runtime_error {
    int status;
    explicit Error(int s)
        : std::runtime_error(std::string(rustprop_status_string(s)) + ": " + last_error()),
          status(s) {}
};

inline double props_si(const char* out, const char* n1, double v1,
                       const char* n2, double v2, const char* fluid) {
    double d = 0.0;
    if (int rc = rustprop_props_si(out, n1, v1, n2, v2, fluid, &d)) throw Error(rc);
    return d;
}

inline std::vector<double> props_si_many(const char* out, const char* n1,
                                         const std::vector<double>& v1, const char* n2,
                                         const std::vector<double>& v2, const char* fluid) {
    std::vector<double> r(v1.size());
    if (int rc = rustprop_props_si_many(out, n1, v1.data(), n2, v2.data(), fluid,
                                        v1.size(), r.data())) throw Error(rc);
    return r;   // failing states are NaN, not exceptions
}

}  // namespace rustprop
```

`find_package(rustprop)` works the same from a C++ project. Link
`rustprop::rustprop_static` and the `RUSTPROP_STATIC` define is applied for
you (it drops `__declspec(dllimport)` on Windows).

---

## 6. Python

No extension module, no build step — `ctypes` on the shared library. The
ready-made wrapper is
[`examples/bindings/rustprop.py`](crates/rustprop-capi/examples/bindings/rustprop.py);
copy it next to the library and import it.

```python
import rustprop

print(rustprop.version(), rustprop.upstream_version())   # 0.1.0 8.0.0
print(rustprop.backends())

d = rustprop.props_si("Dmolar", "T", 300, "P", 101325, "Water")   # 55317.3527735012
h = rustprop.props_si("H", "T", 300, "P", 101325, "IF97::Water")
w = rustprop.ha_props_si("W", "T", 300, "P", 101325, "R", 0.5)

# One output over many states; failures come back as nan.
xs = rustprop.props_si_many("Dmolar", "T", [300.0, 400.0, 500.0],
                            "P", [101325.0] * 3, "Water")

try:
    rustprop.props_si("Dmolar", "T", 300, "P", 101325, "NoSuchFluid")
except rustprop.RustpropError as e:
    print(e)     # ValueError: key [NoSuchFluid] was not found in ...
```

Point it at the library with `RUSTPROP_LIB=/path/to/librustprop.so`, or put
the library beside the module.

> **If you write your own `ctypes` layer, declare `argtypes` and `restype`.**
> Without them ctypes assumes `int` parameters and every `double` you pass
> arrives as garbage — silently, as a plausible-looking wrong number.

With NumPy, `props_si_many` takes any sequence; for large sweeps pass the
array's buffer directly to avoid a copy.

---

## 7. Go

`cgo`, against the header and library. Full file:
[`examples/bindings/rustprop.go`](crates/rustprop-capi/examples/bindings/rustprop.go).

```bash
CGO_CFLAGS="-I<sdk>/include" \
CGO_LDFLAGS="-L<sdk>/lib -lrustprop" \
LD_LIBRARY_PATH=<sdk>/lib \
go run main.go
```

```go
d, err := PropsSI("Dmolar", "T", 300, "P", 101325, "Water")
if err != nil {
    var rperr *Error
    if errors.As(err, &rperr) { log.Fatalf("status %d: %s", rperr.Status, rperr.Message) }
}

many, _ := PropsSIMany("Dmolar", "T", temps, "P", press, "Water")
```

> Read the error message on the **same goroutine**, immediately after the
> failing call. The slot is thread-local and Go may move a goroutine between
> OS threads at a call boundary. The wrapper does exactly this.

For a self-contained binary, link the static library instead:
`CGO_LDFLAGS="<sdk>/lib/librustprop.a -lm"` plus the `Libs.private` list.

---

## 8. Java

The Foreign Function & Memory API, stable since **JDK 22** — no JNI, no
wrapper library to compile. Full file:
[`examples/bindings/Rustprop.java`](crates/rustprop-capi/examples/bindings/Rustprop.java).

```bash
RUSTPROP_LIB=<sdk>/lib/librustprop.so \
  java --enable-native-access=ALL-UNNAMED Rustprop.java
```

```java
try (var rp = new Rustprop(Path.of(System.getenv("RUSTPROP_LIB")))) {
    double d = rp.propsSi("Dmolar", "T", 300, "P", 101325, "Water");
    double[] many = rp.propsSiMany("Dmolar", "T", temps, "P", press, "Water");
    if (rp.hasBackend("humid-air"))
        System.out.println(rp.haPropsSi("W", "T", 300, "P", 101325, "R", 0.5));
}
```

Works from Kotlin and Scala unchanged. Two things that bite:

- A `const char *` return arrives as a zero-length `MemorySegment`; call
  `.reinterpret(...)` before `getString(0)`.
- `Arena.ofShared()` for the library handle, so it outlives one call.

---

## 9. Fortran

`iso_c_binding`, directly against the C ABI — no shim. Full module:
[`examples/bindings/rustprop.f90`](crates/rustprop-capi/examples/bindings/rustprop.f90).

```bash
gfortran rustprop.f90 -o demo -L<sdk>/lib -lrustprop
LD_LIBRARY_PATH=<sdk>/lib ./demo
```

```fortran
use rustprop
real(c_double) :: d
integer(c_int) :: status

call props_si("Dmolar", "T", 300.0_c_double, "P", 101325.0_c_double, "Water", d, status)
if (status /= RUSTPROP_OK) then
    print *, "rustprop: ", last_error()
end if
```

The module handles NUL-termination for you (`cstr`) and reading `const char *`
returns (`from_cptr`), which are the two things that make raw
`iso_c_binding` tedious.

---

## 10. JavaScript and TypeScript

For the browser and Node, use the **WebAssembly** bundle, not the C ABI. This
is what rustprop was built for, and the bundles are small because engines and
fluids are selected per preset.

```bash
npm  # unpack rustprop-wasm-<preset>.tar.gz, or build your own:
wasm-pack build crates/rustprop-wasm --target web --features heos,water
```

Three targets ship in each preset archive: `web` (browsers, ES modules),
`nodejs` (CommonJS `require`), `bundler` (webpack/Vite/Rollup).

```js
import init, { props_si, props_si_many, ha_props_si } from "./pkg/rustprop_wasm.js";

await init();                                  // once, before any call
props_si("D", "T", 400, "P", 101325, "IF97::Water");    // 0.55492158...

// The typed fast path: strings stay on this side, numbers cross as
// Float64Array. Use it for sweeps — per-call string marshalling otherwise
// dominates your runtime long before the thermodynamics does.
const t = new Float64Array([300, 400, 500]);
const p = new Float64Array([101325, 101325, 101325]);
props_si_many("Dmolar", "T", t, "P", p, "Water");   // failures are NaN
```

Errors arrive as thrown JS exceptions carrying rustprop's message, so
`try`/`catch` behaves the way you expect.

Presets and their sizes are in [WASM-SIZES.md](WASM-SIZES.md) — IF97 alone is
128 KB, HEOS with Water 339 KB, everything 4.2 MB. Prefer a narrow preset;
that is the entire point of the project.

> Unlike the C ABI, the wasm bindings **omit** exports for engines that were
> not compiled in — a JS caller finds out by catching a `TypeError`. Check
> before calling if you built a narrow bundle.

---

## 11. Any other language

If it can call a C function, it can call rustprop. The pattern is the same
everywhere: load the shared library, declare the twelve functions from
`rustprop.h`, pass UTF-8 NUL-terminated strings and a pointer for the result,
check the `int` status.

The examples below follow that pattern but are **not machine-checked here** —
unlike §4–§10, which are real programs run by
[`bindings-test.sh`](crates/rustprop-capi/bindings-test.sh) and `ctest.sh`.

**C# / .NET** — P/Invoke:

```csharp
[DllImport("rustprop", CallingConvention = CallingConvention.Cdecl)]
static extern int rustprop_props_si(
    [MarshalAs(UnmanagedType.LPUTF8Str)] string output,
    [MarshalAs(UnmanagedType.LPUTF8Str)] string name1, double val1,
    [MarshalAs(UnmanagedType.LPUTF8Str)] string name2, double val2,
    [MarshalAs(UnmanagedType.LPUTF8Str)] string fluid, out double result);
```

**Julia** — `ccall`:

```julia
function props_si(out, n1, v1, n2, v2, fluid)
    r = Ref{Cdouble}(0)
    rc = ccall((:rustprop_props_si, "librustprop"), Cint,
               (Cstring, Cstring, Cdouble, Cstring, Cdouble, Cstring, Ref{Cdouble}),
               out, n1, v1, n2, v2, fluid, r)
    rc == 0 || error("rustprop status $rc")
    r[]
end
```

**MATLAB** — `loadlibrary`:

```matlab
loadlibrary('librustprop', 'rustprop.h');
[rc, ~, ~, ~, ~, d] = calllib('librustprop', 'rustprop_props_si', ...
                              'Dmolar', 'T', 300, 'P', 101325, 'Water', 0);
```

**Ruby** (`fiddle`), **R** (`.C`/Rcpp), **Lua** (LuaJIT `ffi`), **Zig**
(`@cImport`), **Swift** (module map), **Nim**, **D**, **Crystal** — all the
same shape.

---

## 12. The command line

Every SDK contains `bin/rustprop-cli`, and it ships standalone as
`rustprop-cli-<ver>-<target>.tar.gz`:

```bash
$ rustprop-cli props Dmolar T 300 P 101325 Water
55317.35277350119
$ rustprop-cli props H T 300 P 101325 IF97::Water
112665.04341853978
$ rustprop-cli ha H T 298.15 P 101325 R 0.5
50423.45039075701
```

One query per invocation, full round-trip precision on stdout, non-zero exit
and a message on stderr for a failure — so it composes with shell tooling.

---

## 13. Containers and cloud

Take the `x86_64-unknown-linux-gnu` or `aarch64-unknown-linux-gnu` SDK and use
a glibc-based image — `debian:stable-slim`, `ubuntu`, or a distroless glibc
base. Link statically against `librustprop.a` so the image only needs libc:

```dockerfile
FROM debian:stable-slim AS build
RUN apt-get update && apt-get install -y --no-install-recommends gcc libc6-dev \
    && rm -rf /var/lib/apt/lists/*
COPY rustprop-0.1.0-x86_64-unknown-linux-gnu/ /opt/rustprop/
COPY app.c .
RUN cc -DRUSTPROP_STATIC -I/opt/rustprop/include app.c \
       /opt/rustprop/lib/librustprop.a -lm -o /app

FROM gcr.io/distroless/base-debian12
COPY --from=build /app /app
ENTRYPOINT ["/app"]
```

There is **no musl artifact** and therefore no `scratch`-based fully static
route today — see the note in §2 for why. If you need one, build it yourself
from source and validate it against your own states; the target compiles
cleanly, it is only the extreme-value goldens that disagree.

For arm64 hosts (Graviton, Ampere, Apple silicon CI) take the
`aarch64-unknown-linux-gnu` artifact — it is built and tested on a native
arm64 runner, not cross-compiled.

There is no global state, no configuration file, no initialisation call and
nothing to shut down, so a request handler can call rustprop directly from any
thread. Memory use is bounded: the per-fluid caches are keyed by fluid, and
there is at most one entry per fluid compiled in.

---

## 14. Troubleshooting

**`undefined reference to rustprop_props_si`** — the linker found the header
but not the library. Add `-L<sdk>/lib -lrustprop`, or use pkg-config / CMake.

**`error while loading shared libraries: librustprop.so`** — it linked, but
the loader cannot find it at run time. Set `LD_LIBRARY_PATH` (Linux),
`DYLD_LIBRARY_PATH` (macOS), or put the DLL beside the executable (Windows).
For a permanent fix use an RPATH: `-Wl,-rpath,<sdk>/lib`.

**Static link fails with a wall of `undefined reference`** — you left out the
system libraries. Use `pkg-config --libs --static rustprop`, or CMake's
`rustprop::rustprop_static`.

**Status 102, `RUSTPROP_UNAVAILABLE`** — the engine is not in this binary.
`rustprop_backends()` says what is. Release artifacts are `all-backends`; a
build from source carries only the features you asked for.

**`Illegal instruction` / `SIGILL` on startup** — you took a `-v2`/`-v3`/`-v4`
artifact and this processor is older than that baseline. Use the unsuffixed
one.

**A number is `NaN` from `props_si_many`** — that state failed; the batch
does not abort on one cell. Re-run that state through the scalar call to get
the reason.

**No `.so` in the package** — check `BUILD-INFO.txt`'s `linkage` line. Some
targets are static-only (Rust does not emit a `cdylib` where the target
defaults to `crt-static`); link `librustprop.a` there.

**`ctypes` returns nonsense** — you did not set `argtypes`/`restype`. See §6.

**Results differ from CoolProp** — they should not; that is the whole point.
Known, deliberate divergences are tabulated in
[NEXT-STEPS.md](NEXT-STEPS.md#known-divergences-from-upstream), and
[CHANGELOG.md](CHANGELOG.md) lists what is not ported. If your case is in
neither, please open an issue — it is a bug.

---

## 15. Checking your copy

Every SDK carries a self-test that verifies real values against known CoolProp
results:

```bash
cc -I include share/rustprop/examples/smoke.c -L lib -lrustprop -lm -o smoke
LD_LIBRARY_PATH=lib ./smoke
```

From the repository, the same checks plus the packaging and every language
binding:

```bash
crates/rustprop-capi/ctest.sh          # C and C++, shared and static
crates/rustprop-capi/bindings-test.sh  # Python, Go, Java, Fortran
```

---

## License

MIT. Derivative work of CoolProp (MIT, © 2012–2018 Ian H. Bell and other
CoolProp developers) — see [LICENSE](LICENSE).
