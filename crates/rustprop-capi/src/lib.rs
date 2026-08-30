//! C ABI for rustprop — `PropsSI` / `HAPropsSI` for C, C++, and any host that
//! can call a C function.
//!
//! The hand-written header is `include/rustprop.h`; it is the contract, and
//! `examples/smoke.c` + `examples/smoke.cc` compile against it and link the
//! built library, which is how CI proves the two agree. Nothing here is
//! generated, so a change to a signature below is a change to the header in
//! the same edit.
//!
//! # What differs from `rustprop-wasm`
//!
//! The JS bindings `#[cfg]` their exports away when an engine is not compiled
//! in, because a JS caller finds out by catching a `TypeError`. A C caller has
//! only the header, and a missing symbol is a link failure — or, for a
//! `dlopen`ed module, a crash. So **every function below is exported by every
//! build**, and a call into an engine this binary does not carry answers
//! [`RUSTPROP_UNAVAILABLE`]. [`rustprop_has_backend`] is how a caller asks what
//! it actually got.
//!
//! # Thread safety
//!
//! Every entry point is safe to call concurrently from any number of threads.
//! This is not a promise talked into existence: `rustprop`'s per-fluid caches
//! live in `static`s, and a `static` requires `Sync`, so the compiler has
//! already proven the shared state is race-free. The per-solve derivative
//! memos are stack locals. `tests::concurrent_props_si_agrees_with_serial`
//! pins the behaviour.
//!
//! The one piece of per-thread state is the last-error slot read by
//! [`rustprop_last_error_message`], which is a `thread_local!` — an error
//! raised on one thread is never visible from another.
//!
//! # Panics do not cross the boundary
//!
//! Unwinding out of an `extern "C"` function is undefined behaviour, so every
//! entry point wraps its body in `catch_unwind` and reports
//! [`RUSTPROP_PANIC`]. This is ABI hygiene, not one of the invented guards
//! CLAUDE.md forbids: it constrains no thermodynamic path and can only fire
//! where the alternative is UB in the caller's process.
//!
//! For it to be able to fire at all the library must be built with
//! `--profile release-capi` (root `Cargo.toml`), which is `release` with
//! `panic = "unwind"`. Under the plain `release` profile — `panic = "abort"`,
//! correct for wasm, where there is no host process to protect — a panic still
//! aborts, and `catch_unwind` is dead weight that costs nothing.

// The exported symbols need `#[unsafe(no_mangle)]` (edition 2024), which the
// workspace-wide `deny(unsafe_code)` rejects. This crate and
// `tools/wasm-size-probe` are the only two exceptions in the tree, and this is
// the only one that ships: a C ABI cannot be written without it.
#![allow(unsafe_code)]

use std::cell::RefCell;
use std::ffi::{CStr, c_char, c_int};
use std::panic::{AssertUnwindSafe, catch_unwind};

// ---------------------------------------------------------------------------
// Status codes
// ---------------------------------------------------------------------------
//
// 1..=11 are one-to-one with `rustprop_core::Error`'s variants, in declaration
// order. `Error` is `#[non_exhaustive]`, so a variant added upstream lands in
// the `_` arm as RUSTPROP_ERROR rather than silently taking another code's
// number — the mapping in `status_of` is exhaustive over what exists today and
// will fail to compile as a warning-free build only if someone reorders it.
//
// 100+ are conditions that exist only at this boundary and have no `Error`
// counterpart.

/// Success.
pub const RUSTPROP_OK: c_int = 0;
/// `NotImplementedError`
pub const RUSTPROP_NOT_IMPLEMENTED: c_int = 1;
/// `SolutionError` — an iterative solver failed to converge.
pub const RUSTPROP_SOLUTION: c_int = 2;
/// `AttributeError`
pub const RUSTPROP_ATTRIBUTE: c_int = 3;
/// `OutOfRangeError` — input outside the range of validity.
pub const RUSTPROP_OUT_OF_RANGE: c_int = 4;
/// `ValueError` — invalid parameter/fluid/input names and values.
pub const RUSTPROP_VALUE: c_int = 5;
/// `WrongFluidError`
pub const RUSTPROP_WRONG_FLUID: c_int = 6;
/// `CompositionError`
pub const RUSTPROP_COMPOSITION: c_int = 7;
/// `InputError`
pub const RUSTPROP_INPUT: c_int = 8;
/// `NotAvailableError` — property not available for this fluid/model.
pub const RUSTPROP_NOT_AVAILABLE: c_int = 9;
/// `KeyError`
pub const RUSTPROP_KEY: c_int = 10;
/// `MultipleSolutionsError`
pub const RUSTPROP_MULTIPLE_SOLUTIONS: c_int = 11;
/// An `Error` variant this ABI predates. See [`rustprop_last_error_message`].
pub const RUSTPROP_ERROR: c_int = 99;

/// A required pointer argument was NULL.
pub const RUSTPROP_NULL_ARGUMENT: c_int = 100;
/// A string argument was not valid UTF-8.
pub const RUSTPROP_INVALID_UTF8: c_int = 101;
/// The engine or fluid data this call needs was not compiled into this binary.
pub const RUSTPROP_UNAVAILABLE: c_int = 102;
/// A panic was caught at the boundary. Always a bug; please report it.
pub const RUSTPROP_PANIC: c_int = 104;

// The mirror image of `unavailable`'s allow: every call site sits behind a
// backend cfg, so a build with NO engines compiles them all away.
#[allow(dead_code)]
fn status_of(e: &rustprop::Error) -> c_int {
    use rustprop::Error as E;
    match e {
        E::NotImplemented(_) => RUSTPROP_NOT_IMPLEMENTED,
        E::Solution(_) => RUSTPROP_SOLUTION,
        E::Attribute(_) => RUSTPROP_ATTRIBUTE,
        E::OutOfRange(_) => RUSTPROP_OUT_OF_RANGE,
        E::Value(_) => RUSTPROP_VALUE,
        E::WrongFluid(_) => RUSTPROP_WRONG_FLUID,
        E::Composition(_) => RUSTPROP_COMPOSITION,
        E::Input(_) => RUSTPROP_INPUT,
        E::NotAvailable(_) => RUSTPROP_NOT_AVAILABLE,
        E::Key(_) => RUSTPROP_KEY,
        E::MultipleSolutions(_) => RUSTPROP_MULTIPLE_SOLUTIONS,
        _ => RUSTPROP_ERROR,
    }
}

// ---------------------------------------------------------------------------
// Per-thread last error
// ---------------------------------------------------------------------------

thread_local! {
    /// (status, message) of the most recent call on THIS thread.
    static LAST_ERROR: RefCell<(c_int, String)> = const {
        RefCell::new((RUSTPROP_OK, String::new()))
    };
}

fn set_error(status: c_int, message: impl Into<String>) -> c_int {
    let message = message.into();
    LAST_ERROR.with(|e| *e.borrow_mut() = (status, message));
    status
}

fn clear_error() {
    LAST_ERROR.with(|e| {
        let mut e = e.borrow_mut();
        e.0 = RUSTPROP_OK;
        e.1.clear();
    });
}

/// The status of the most recent call on the calling thread.
#[unsafe(no_mangle)]
pub extern "C" fn rustprop_last_error_code() -> c_int {
    LAST_ERROR.with(|e| e.borrow().0)
}

/// Copy the most recent error message on the calling thread into `buf`.
///
/// Returns the message length in bytes, NOT counting the terminating NUL —
/// so a return value `>= len` means the message was truncated, and the call
/// can be repeated with a buffer of `returned + 1` bytes. What is written is
/// always NUL-terminated as long as `len > 0`.
///
/// Passing `buf = NULL` / `len = 0` is how you ask for the length alone.
///
/// # Safety
///
/// `buf` must be NULL, or point to at least `len` writable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rustprop_last_error_message(buf: *mut c_char, len: usize) -> usize {
    LAST_ERROR.with(|e| {
        let e = e.borrow();
        let bytes = e.1.as_bytes();
        if !buf.is_null() && len > 0 {
            // Leave room for the NUL, and never split a UTF-8 sequence: a C
            // caller that hands the buffer to a UTF-8-aware API would
            // otherwise get a byte string it cannot decode.
            let mut n = bytes.len().min(len - 1);
            while n > 0 && !e.1.is_char_boundary(n) {
                n -= 1;
            }
            // SAFETY: `n <= len - 1` writable bytes, then the NUL at `n`.
            unsafe {
                std::ptr::copy_nonoverlapping(bytes.as_ptr(), buf as *mut u8, n);
                *buf.add(n) = 0;
            }
        }
        bytes.len()
    })
}

/// A short, static, human-readable name for a status code (never NULL).
#[unsafe(no_mangle)]
pub extern "C" fn rustprop_status_string(status: c_int) -> *const c_char {
    let s: &CStr = match status {
        RUSTPROP_OK => c"ok",
        RUSTPROP_NOT_IMPLEMENTED => c"NotImplementedError",
        RUSTPROP_SOLUTION => c"SolutionError",
        RUSTPROP_ATTRIBUTE => c"AttributeError",
        RUSTPROP_OUT_OF_RANGE => c"OutOfRangeError",
        RUSTPROP_VALUE => c"ValueError",
        RUSTPROP_WRONG_FLUID => c"WrongFluidError",
        RUSTPROP_COMPOSITION => c"CompositionError",
        RUSTPROP_INPUT => c"InputError",
        RUSTPROP_NOT_AVAILABLE => c"NotAvailableError",
        RUSTPROP_KEY => c"KeyError",
        RUSTPROP_MULTIPLE_SOLUTIONS => c"MultipleSolutionsError",
        RUSTPROP_ERROR => c"Error",
        RUSTPROP_NULL_ARGUMENT => c"NullArgument",
        RUSTPROP_INVALID_UTF8 => c"InvalidUtf8",
        RUSTPROP_UNAVAILABLE => c"Unavailable",
        RUSTPROP_PANIC => c"Panic",
        _ => c"UnknownStatus",
    };
    s.as_ptr()
}

// ---------------------------------------------------------------------------
// Argument marshalling
// ---------------------------------------------------------------------------

/// Borrow a C string argument. `Err(status)` is already recorded as the
/// thread's last error, so a caller can return it directly.
///
/// # Safety
///
/// `p` must be NULL or a NUL-terminated C string valid for the call.
unsafe fn str_arg<'a>(p: *const c_char, what: &str) -> Result<&'a str, c_int> {
    if p.is_null() {
        return Err(set_error(
            RUSTPROP_NULL_ARGUMENT,
            format!("argument `{what}` was NULL"),
        ));
    }
    // SAFETY: non-NULL, and the caller guarantees NUL termination.
    unsafe { CStr::from_ptr(p) }.to_str().map_err(|_| {
        set_error(
            RUSTPROP_INVALID_UTF8,
            format!("argument `{what}` was not valid UTF-8"),
        )
    })
}

/// Borrow several C string arguments, stopping at the FIRST bad one.
///
/// Short-circuiting is the point. The earlier shape evaluated all of them and
/// chose a status afterwards, which could return one argument's status
/// alongside a different argument's message: `str_arg` records into the
/// thread's last-error slot as it goes, so the LAST failure won the message
/// while the FIRST won the code. A caller doing the normal thing — check the
/// code, then read the message — was told about two different arguments.
/// `status_and_message_describe_the_same_argument` pins it.
///
/// Expands to a `return` from the enclosing closure, so it is only usable
/// inside an entry point's `guard(|| ...)` body.
macro_rules! str_args {
    ($($p:expr => $name:literal),+ $(,)?) => {
        ($(match str_arg($p, $name) {
            Ok(v) => v,
            Err(status) => return status,
        },)+)
    };
}

/// Run an entry point's body, turning a panic into [`RUSTPROP_PANIC`] rather
/// than unwinding into C. See the module docs.
fn guard(f: impl FnOnce() -> c_int) -> c_int {
    match catch_unwind(AssertUnwindSafe(f)) {
        Ok(status) => status,
        Err(_) => set_error(
            RUSTPROP_PANIC,
            "rustprop panicked; this is a bug — please report it at \
             https://github.com/ernsoylu/RustProp/issues",
        ),
    }
}

/// The status/message for a call into an engine this binary was not built
/// with.
///
/// `dead_code`-allowed on purpose: it is reached only from the `cfg(not(...))`
/// arms, so an `all-backends` build — which is what CI lints — compiles every
/// one of its call sites away. Named so the message tells the caller what to do about it, since a
/// prebuilt binary cannot be recompiled by whoever hit this.
#[allow(dead_code)]
fn unavailable(backend: &str) -> c_int {
    set_error(
        RUSTPROP_UNAVAILABLE,
        format!(
            "this rustprop build does not carry the `{backend}` backend; \
             rebuild with --features {backend}, or use a release artifact \
             built with all-backends (see rustprop_backends())"
        ),
    )
}

// ---------------------------------------------------------------------------
// Introspection
// ---------------------------------------------------------------------------

/// This library's own version, e.g. `"0.1.0"` (never NULL, static).
#[unsafe(no_mangle)]
pub extern "C" fn rustprop_version() -> *const c_char {
    // A NUL is appended at compile time so the pointer is a valid C string
    // without allocating.
    concat!(env!("CARGO_PKG_VERSION"), "\0").as_ptr() as *const c_char
}

/// The upstream CoolProp release this port tracks, e.g. `"8.0.0"` (never NULL,
/// static).
#[unsafe(no_mangle)]
pub extern "C" fn rustprop_upstream_version() -> *const c_char {
    concat!("8.0.0", "\0").as_ptr() as *const c_char
}

/// Every backend compiled into THIS binary, comma-separated and without a
/// trailing separator, e.g. `"heos,if97,humid-air"`. The empty string if none.
/// Never NULL; the pointer stays valid for the life of the process.
///
/// Engine selection is a compile-time choice, so this is the only way a caller
/// holding a prebuilt library can find out what it may ask for.
#[unsafe(no_mangle)]
pub extern "C" fn rustprop_backends() -> *const c_char {
    use std::ffi::CString;
    use std::sync::OnceLock;
    // Built once and leaked rather than assembled from `concat!`: stable Rust
    // has no way to drop a trailing separator at compile time, and one
    // allocation over the life of the process is not worth an ugly string.
    static LIST: OnceLock<&'static CStr> = OnceLock::new();
    LIST.get_or_init(|| {
        let names: Vec<&str> = BACKEND_TABLE
            .iter()
            .filter(|(_, present)| *present)
            .map(|(name, _)| *name)
            .collect();
        let joined = CString::new(names.join(",")).expect("backend names contain no NUL");
        &*Box::leak(joined.into_boxed_c_str())
    })
    .as_ptr()
}

/// Every backend this ABI knows about, and whether THIS build carries it.
///
/// One table drives both [`rustprop_backends`] and [`rustprop_has_backend`],
/// so the two can never disagree — `introspection_describes_this_build` would
/// catch it if they did. Order matches the README's feature table.
const BACKEND_TABLE: [(&str, bool); 9] = [
    ("heos", cfg!(feature = "heos")),
    ("heos-mixtures", cfg!(feature = "heos-mixtures")),
    ("if97", cfg!(feature = "if97")),
    ("cubics", cfg!(feature = "cubics")),
    ("incompressible", cfg!(feature = "incompressible")),
    ("pcsaft", cfg!(feature = "pcsaft")),
    ("tabular", cfg!(feature = "tabular")),
    ("svdsbtl", cfg!(feature = "svdsbtl")),
    ("humid-air", cfg!(feature = "humid-air")),
];

/// Is `name` a backend this binary carries? 1 for yes, 0 for no (including an
/// unknown or NULL name).
///
/// Accepted names are the Cargo feature names: `heos`, `heos-mixtures`,
/// `if97`, `cubics`, `incompressible`, `pcsaft`, `tabular`, `svdsbtl`,
/// `humid-air`.
///
/// # Safety
///
/// `name` must be NULL or a NUL-terminated C string valid for the call.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rustprop_has_backend(name: *const c_char) -> c_int {
    if name.is_null() {
        return 0;
    }
    // SAFETY: non-NULL, and the caller guarantees NUL termination.
    let Ok(name) = unsafe { CStr::from_ptr(name) }.to_str() else {
        return 0;
    };
    // Read from the same table `rustprop_backends` prints, so the two cannot
    // drift apart. An unknown name simply is not in it.
    let present = BACKEND_TABLE
        .iter()
        .any(|(known, compiled_in)| *known == name && *compiled_in);
    c_int::from(present)
}

/// How many HEOS fluids this binary carries data for.
///
/// Fluid data is opt-in per fluid, so this is as build-dependent as
/// [`rustprop_backends`]. 0 when the `heos` backend is absent.
#[unsafe(no_mangle)]
pub extern "C" fn rustprop_fluid_count() -> usize {
    #[cfg(feature = "heos")]
    {
        rustprop_data::fluids::all().len()
    }
    #[cfg(not(feature = "heos"))]
    {
        0
    }
}

/// The name of compiled-in HEOS fluid `index`, or NULL if out of range.
///
/// The pointer is static and stays valid for the life of the process. Iterate
/// `0 .. rustprop_fluid_count()` to enumerate what this binary can answer for.
#[unsafe(no_mangle)]
pub extern "C" fn rustprop_fluid_name(index: usize) -> *const c_char {
    #[cfg(feature = "heos")]
    {
        // The generated names are plain ASCII identifiers with no interior
        // NUL, but they are Rust `&'static str` and therefore NOT
        // NUL-terminated. One CString per fluid is leaked on first call —
        // bounded by the fluid count (130 at most), built once.
        use std::collections::BTreeMap;
        use std::ffi::CString;
        use std::sync::OnceLock;
        // The generated names are `&'static str` and therefore NOT
        // NUL-terminated, so one CString per fluid is built on first call and
        // leaked. Bounded by the fluid count (130 at most) and built once, so
        // the pointers handed out stay valid for the life of the process —
        // which is what the C contract above promises.
        //
        // Sorted, so the index a caller iterates is stable across builds that
        // enable the same fluids in a different feature order.
        static TABLE: OnceLock<Vec<&'static CStr>> = OnceLock::new();
        let table = TABLE.get_or_init(|| {
            let sorted: BTreeMap<&str, ()> = rustprop_data::fluids::all()
                .iter()
                .map(|(name, _)| (*name, ()))
                .collect();
            sorted
                .into_keys()
                .map(|name| {
                    let owned = CString::new(name).expect("fluid name contains no interior NUL");
                    &*Box::leak(owned.into_boxed_c_str())
                })
                .collect()
        });
        match table.get(index) {
            Some(name) => name.as_ptr(),
            None => std::ptr::null(),
        }
    }
    #[cfg(not(feature = "heos"))]
    {
        let _ = index;
        std::ptr::null()
    }
}

// ---------------------------------------------------------------------------
// PropsSI
// ---------------------------------------------------------------------------

/// `PropsSI(output, name1, val1, name2, val2, fluid)` — writes the result to
/// `*out` and returns [`RUSTPROP_OK`], or leaves `*out` untouched and returns
/// a non-zero status whose message is [`rustprop_last_error_message`].
///
/// # Safety
///
/// Every string argument must be a NUL-terminated C string valid for the call;
/// `out` must point to a writable `double`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rustprop_props_si(
    output: *const c_char,
    name1: *const c_char,
    val1: f64,
    name2: *const c_char,
    val2: f64,
    fluid: *const c_char,
    out: *mut f64,
) -> c_int {
    guard(|| {
        clear_error();
        if out.is_null() {
            return set_error(RUSTPROP_NULL_ARGUMENT, "argument `out` was NULL");
        }
        // SAFETY: the caller guarantees these are NUL-terminated or NULL.
        let (output, name1, name2, fluid) = unsafe {
            str_args!(
                output => "output",
                name1 => "name1",
                name2 => "name2",
                fluid => "fluid",
            )
        };
        #[cfg(any(
            feature = "heos",
            feature = "if97",
            feature = "cubics",
            feature = "incompressible",
            feature = "pcsaft"
        ))]
        {
            match rustprop::props_si(output, name1, val1, name2, val2, fluid) {
                Ok(v) => {
                    // SAFETY: checked non-NULL above.
                    unsafe { *out = v };
                    RUSTPROP_OK
                }
                Err(e) => set_error(status_of(&e), e.to_string()),
            }
        }
        #[cfg(not(any(
            feature = "heos",
            feature = "if97",
            feature = "cubics",
            feature = "incompressible",
            feature = "pcsaft"
        )))]
        {
            let _ = (output, name1, val1, name2, val2, fluid);
            unavailable("heos")
        }
    })
}

/// `PropsSI` over many states at once: one output, one fluid, two input
/// vectors of `n` values each, `n` results written to `out`.
///
/// A state that fails yields NaN in its slot and does not stop the batch — a
/// sweep over a grid routinely clips the phase envelope, and one bad cell
/// should not cost the caller the other ten thousand. The return status
/// therefore reports only whether the CALL was well-formed. This mirrors
/// `props_si_many` in the JS bindings.
///
/// Passing `n = 0` is valid and does nothing.
///
/// # Safety
///
/// String arguments must be NUL-terminated; `vals1` and `vals2` must each
/// point to `n` readable doubles, and `out` to `n` writable doubles. `out` may
/// alias neither input.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rustprop_props_si_many(
    output: *const c_char,
    name1: *const c_char,
    vals1: *const f64,
    name2: *const c_char,
    vals2: *const f64,
    fluid: *const c_char,
    n: usize,
    out: *mut f64,
) -> c_int {
    guard(|| {
        clear_error();
        if n == 0 {
            return RUSTPROP_OK;
        }
        if vals1.is_null() || vals2.is_null() || out.is_null() {
            return set_error(
                RUSTPROP_NULL_ARGUMENT,
                "`vals1`, `vals2` and `out` must all be non-NULL when n > 0",
            );
        }
        // SAFETY: the caller guarantees these are NUL-terminated or NULL.
        let (output, name1, name2, fluid) = unsafe {
            str_args!(
                output => "output",
                name1 => "name1",
                name2 => "name2",
                fluid => "fluid",
            )
        };
        #[cfg(any(
            feature = "heos",
            feature = "if97",
            feature = "cubics",
            feature = "incompressible",
            feature = "pcsaft"
        ))]
        {
            // SAFETY: the caller guarantees n readable/writable elements.
            let (v1, v2, o) = unsafe {
                (
                    std::slice::from_raw_parts(vals1, n),
                    std::slice::from_raw_parts(vals2, n),
                    std::slice::from_raw_parts_mut(out, n),
                )
            };
            for i in 0..n {
                o[i] = rustprop::props_si(output, name1, v1[i], name2, v2[i], fluid)
                    .unwrap_or(f64::NAN);
            }
            RUSTPROP_OK
        }
        #[cfg(not(any(
            feature = "heos",
            feature = "if97",
            feature = "cubics",
            feature = "incompressible",
            feature = "pcsaft"
        )))]
        {
            let _ = (output, name1, name2, fluid);
            unavailable("heos")
        }
    })
}

// ---------------------------------------------------------------------------
// HAPropsSI
// ---------------------------------------------------------------------------

/// `HAPropsSI(output, name1, val1, name2, val2, name3, val3)` — humid air.
/// Same contract as [`rustprop_props_si`].
///
/// # Safety
///
/// Every string argument must be a NUL-terminated C string valid for the call;
/// `out` must point to a writable `double`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rustprop_ha_props_si(
    output: *const c_char,
    name1: *const c_char,
    val1: f64,
    name2: *const c_char,
    val2: f64,
    name3: *const c_char,
    val3: f64,
    out: *mut f64,
) -> c_int {
    guard(|| {
        clear_error();
        if out.is_null() {
            return set_error(RUSTPROP_NULL_ARGUMENT, "argument `out` was NULL");
        }
        // SAFETY: the caller guarantees these are NUL-terminated or NULL.
        let (output, name1, name2, name3) = unsafe {
            str_args!(
                output => "output",
                name1 => "name1",
                name2 => "name2",
                name3 => "name3",
            )
        };
        #[cfg(feature = "humid-air")]
        {
            match rustprop::ha_props_si(output, name1, val1, name2, val2, name3, val3) {
                Ok(v) => {
                    // SAFETY: checked non-NULL above.
                    unsafe { *out = v };
                    RUSTPROP_OK
                }
                Err(e) => set_error(status_of(&e), e.to_string()),
            }
        }
        #[cfg(not(feature = "humid-air"))]
        {
            let _ = (output, name1, val1, name2, val2, name3, val3);
            unavailable("humid-air")
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::CString;

    /// Call `rustprop_props_si` the way C does, from Rust.
    fn props(output: &str, n1: &str, v1: f64, n2: &str, v2: f64, fluid: &str) -> (c_int, f64) {
        let (o, a, b, f) = (
            CString::new(output).unwrap(),
            CString::new(n1).unwrap(),
            CString::new(n2).unwrap(),
            CString::new(fluid).unwrap(),
        );
        let mut out = f64::NAN;
        let status = unsafe {
            rustprop_props_si(
                o.as_ptr(),
                a.as_ptr(),
                v1,
                b.as_ptr(),
                v2,
                f.as_ptr(),
                &mut out,
            )
        };
        (status, out)
    }

    fn last_message() -> String {
        let n = unsafe { rustprop_last_error_message(std::ptr::null_mut(), 0) };
        let mut buf = vec![0u8; n + 1];
        unsafe { rustprop_last_error_message(buf.as_mut_ptr() as *mut c_char, buf.len()) };
        CStr::from_bytes_until_nul(&buf)
            .unwrap()
            .to_string_lossy()
            .into_owned()
    }

    /// The golden value the README and the Rust doctest both quote, reached
    /// through the C signature instead of the Rust one.
    #[test]
    #[cfg(feature = "heos")]
    fn props_si_matches_the_rust_api() {
        let (status, v) = props("Dmolar", "T", 300.0, "P", 101325.0, "Water");
        assert_eq!(status, RUSTPROP_OK);
        let expected = rustprop::props_si("Dmolar", "T", 300.0, "P", 101325.0, "Water").unwrap();
        assert_eq!(
            v.to_bits(),
            expected.to_bits(),
            "C ABI must not reshape the number"
        );
        assert!(((v - 55317.35277350119) / v).abs() < 1e-8);
    }

    /// An error must arrive as a status AND a message, with `*out` untouched.
    #[test]
    #[cfg(feature = "heos")]
    fn a_bad_fluid_reports_rather_than_writes() {
        let (o, a, b, f) = (
            CString::new("Dmolar").unwrap(),
            CString::new("T").unwrap(),
            CString::new("P").unwrap(),
            CString::new("NoSuchFluid").unwrap(),
        );
        let mut out = -1.0;
        let status = unsafe {
            rustprop_props_si(
                o.as_ptr(),
                a.as_ptr(),
                300.0,
                b.as_ptr(),
                101325.0,
                f.as_ptr(),
                &mut out,
            )
        };
        assert_ne!(status, RUSTPROP_OK);
        assert_eq!(out, -1.0, "`out` must not be written on failure");
        assert_eq!(rustprop_last_error_code(), status);
        assert!(
            last_message().contains("NoSuchFluid"),
            "message should name the offending key, got {:?}",
            last_message()
        );
    }

    /// NULL is a status, never a crash.
    #[test]
    fn null_arguments_are_refused() {
        let mut out = 0.0;
        let s = unsafe {
            rustprop_props_si(
                std::ptr::null(),
                std::ptr::null(),
                0.0,
                std::ptr::null(),
                0.0,
                std::ptr::null(),
                &mut out,
            )
        };
        assert_eq!(s, RUSTPROP_NULL_ARGUMENT);

        let ok = CString::new("T").unwrap();
        let s = unsafe {
            rustprop_props_si(
                ok.as_ptr(),
                ok.as_ptr(),
                0.0,
                ok.as_ptr(),
                0.0,
                ok.as_ptr(),
                std::ptr::null_mut(),
            )
        };
        assert_eq!(s, RUSTPROP_NULL_ARGUMENT);
    }

    /// Truncation must be reported, NUL-terminated, and never split a UTF-8
    /// sequence.
    #[test]
    #[cfg(feature = "heos")]
    fn error_message_truncates_safely() {
        let _ = props("Dmolar", "T", 300.0, "P", 101325.0, "NoSuchFluid");
        let needed = unsafe { rustprop_last_error_message(std::ptr::null_mut(), 0) };
        assert!(needed > 8, "expected a real message, got {needed} bytes");

        let mut small = [0xAAu8; 8];
        let reported =
            unsafe { rustprop_last_error_message(small.as_mut_ptr() as *mut c_char, small.len()) };
        assert_eq!(
            reported, needed,
            "must report the FULL length, not what fit"
        );
        assert_eq!(small[7], 0, "must NUL-terminate within the buffer");
        assert!(std::str::from_utf8(&small[..7]).is_ok());
    }

    /// `len = 0` must not write, and must still report the length.
    #[test]
    #[cfg(feature = "heos")]
    fn error_message_with_zero_length_writes_nothing() {
        let _ = props("Dmolar", "T", 300.0, "P", 101325.0, "NoSuchFluid");
        let mut canary = [0x5Au8; 4];
        let n = unsafe { rustprop_last_error_message(canary.as_mut_ptr() as *mut c_char, 0) };
        assert!(n > 0);
        assert_eq!(canary, [0x5A; 4], "len = 0 must not touch the buffer");
    }

    /// A success must clear a previous failure's status.
    #[test]
    #[cfg(feature = "heos")]
    fn success_clears_the_previous_error() {
        let _ = props("Dmolar", "T", 300.0, "P", 101325.0, "NoSuchFluid");
        assert_ne!(rustprop_last_error_code(), RUSTPROP_OK);
        let (status, _) = props("Dmolar", "T", 300.0, "P", 101325.0, "Water");
        assert_eq!(status, RUSTPROP_OK);
        assert_eq!(rustprop_last_error_code(), RUSTPROP_OK);
        assert_eq!(last_message(), "");
    }

    /// The batch path must agree with the scalar path element for element,
    /// bitwise, and must NaN-fill rather than abort on a bad cell.
    #[test]
    #[cfg(feature = "heos")]
    fn batch_agrees_with_scalar_and_survives_a_bad_cell() {
        let temps = [300.0, 400.0, 500.0, -1.0, 600.0];
        let pressures = [101325.0; 5];
        let mut out = [0.0f64; 5];
        let (o, a, b, f) = (
            CString::new("Dmolar").unwrap(),
            CString::new("T").unwrap(),
            CString::new("P").unwrap(),
            CString::new("Water").unwrap(),
        );
        let status = unsafe {
            rustprop_props_si_many(
                o.as_ptr(),
                a.as_ptr(),
                temps.as_ptr(),
                b.as_ptr(),
                pressures.as_ptr(),
                f.as_ptr(),
                temps.len(),
                out.as_mut_ptr(),
            )
        };
        assert_eq!(status, RUSTPROP_OK, "a bad cell must not fail the call");
        assert!(
            out[3].is_nan(),
            "T = -1 K should land as NaN, got {}",
            out[3]
        );
        for i in [0usize, 1, 2, 4] {
            let (s, want) = props("Dmolar", "T", temps[i], "P", pressures[i], "Water");
            assert_eq!(s, RUSTPROP_OK);
            assert_eq!(out[i].to_bits(), want.to_bits(), "batch cell {i} diverged");
        }
    }

    /// n = 0 must be a no-op, not a NULL-pointer complaint.
    #[test]
    fn empty_batch_is_ok() {
        let (o, a, b, f) = (
            CString::new("Dmolar").unwrap(),
            CString::new("T").unwrap(),
            CString::new("P").unwrap(),
            CString::new("Water").unwrap(),
        );
        let status = unsafe {
            rustprop_props_si_many(
                o.as_ptr(),
                a.as_ptr(),
                std::ptr::null(),
                b.as_ptr(),
                std::ptr::null(),
                f.as_ptr(),
                0,
                std::ptr::null_mut(),
            )
        };
        assert_eq!(status, RUSTPROP_OK);
    }

    /// The module docs claim every entry point is safe to call concurrently.
    /// This is the check behind that claim: the same states, computed on eight
    /// threads at once, must come back bit-identical to the serial answers.
    #[test]
    #[cfg(feature = "heos")]
    fn concurrent_props_si_agrees_with_serial() {
        let states: Vec<f64> = (0..200).map(|i| 280.0 + f64::from(i)).collect();
        let serial: Vec<u64> = states
            .iter()
            .map(|&t| props("Dmolar", "T", t, "P", 101325.0, "Water").1.to_bits())
            .collect();

        let threads: Vec<_> = (0..8)
            .map(|_| {
                let states = states.clone();
                std::thread::spawn(move || {
                    states
                        .iter()
                        .map(|&t| props("Dmolar", "T", t, "P", 101325.0, "Water").1.to_bits())
                        .collect::<Vec<u64>>()
                })
            })
            .collect();

        for t in threads {
            assert_eq!(t.join().unwrap(), serial, "concurrent answers diverged");
        }
    }

    /// The last-error slot is per-thread: one thread's failure must be
    /// invisible to another.
    #[test]
    #[cfg(feature = "heos")]
    fn last_error_does_not_leak_across_threads() {
        let _ = props("Dmolar", "T", 300.0, "P", 101325.0, "NoSuchFluid");
        assert_ne!(rustprop_last_error_code(), RUSTPROP_OK);
        let other = std::thread::spawn(|| rustprop_last_error_code())
            .join()
            .unwrap();
        assert_eq!(other, RUSTPROP_OK, "error leaked into another thread");
    }

    /// Introspection must describe THIS build, and `rustprop_backends()` must
    /// be a well-formed C string that agrees with `rustprop_has_backend`.
    #[test]
    fn introspection_describes_this_build() {
        let backends = unsafe { CStr::from_ptr(rustprop_backends()) }
            .to_str()
            .expect("backend list is valid UTF-8");
        assert!(
            !backends.ends_with(','),
            "the list must not carry a trailing separator, got {backends:?}"
        );
        for name in backends.split(',').filter(|s| !s.is_empty()) {
            let c = CString::new(name).unwrap();
            assert_eq!(
                unsafe { rustprop_has_backend(c.as_ptr()) },
                1,
                "`{name}` is listed but has_backend says no"
            );
        }
        let absent = CString::new("no-such-backend").unwrap();
        assert_eq!(unsafe { rustprop_has_backend(absent.as_ptr()) }, 0);
        assert_eq!(unsafe { rustprop_has_backend(std::ptr::null()) }, 0);

        let version = unsafe { CStr::from_ptr(rustprop_version()) }
            .to_str()
            .unwrap();
        assert_eq!(version, env!("CARGO_PKG_VERSION"));
        let upstream = unsafe { CStr::from_ptr(rustprop_upstream_version()) }
            .to_str()
            .unwrap();
        assert_eq!(upstream, rustprop::UPSTREAM_VERSION);
    }

    /// Fluid enumeration must be in range, NUL-terminated, and resolvable —
    /// a name this reports must be a name `props_si` accepts.
    #[test]
    #[cfg(feature = "heos")]
    fn fluid_enumeration_round_trips() {
        let n = rustprop_fluid_count();
        assert!(n > 0, "the heos build must carry at least one fluid");
        assert!(
            rustprop_fluid_name(n).is_null(),
            "index n must be out of range"
        );

        for i in 0..n {
            let p = rustprop_fluid_name(i);
            assert!(!p.is_null(), "index {i} < count must resolve");
            let name = unsafe { CStr::from_ptr(p) }.to_str().unwrap();
            assert!(!name.is_empty());
        }
        // Stability: the table is built once, so two reads give one pointer.
        assert_eq!(rustprop_fluid_name(0), rustprop_fluid_name(0));
    }

    /// The IF97 route through the C ABI, on the golden the README quotes.
    /// Also the reason `props` is exercised by an if97-only build, where the
    /// HEOS tests below are all compiled out.
    #[test]
    #[cfg(feature = "if97")]
    fn if97_route_matches_the_rust_api() {
        let (status, v) = props("H", "T", 300.0, "P", 101325.0, "IF97::Water");
        assert_eq!(status, RUSTPROP_OK, "{}", last_message());
        assert!(((v - 112665.04341853978) / v).abs() < 1e-11);
    }

    /// Call `rustprop_ha_props_si` the way C does, from Rust.
    fn ha(o: &str, n1: &str, v1: f64, n2: &str, v2: f64, n3: &str, v3: f64) -> (c_int, f64) {
        let (o, a, b, c) = (
            CString::new(o).unwrap(),
            CString::new(n1).unwrap(),
            CString::new(n2).unwrap(),
            CString::new(n3).unwrap(),
        );
        let mut out = f64::NAN;
        let status = unsafe {
            rustprop_ha_props_si(
                o.as_ptr(),
                a.as_ptr(),
                v1,
                b.as_ptr(),
                v2,
                c.as_ptr(),
                v3,
                &mut out,
            )
        };
        (status, out)
    }

    /// The point of exporting every symbol from every build: a call into an
    /// absent engine must be a STATUS, not a missing symbol. Only compiled
    /// where the engine really is absent, which is what makes it meaningful.
    #[test]
    #[cfg(not(any(
        feature = "heos",
        feature = "if97",
        feature = "cubics",
        feature = "incompressible",
        feature = "pcsaft"
    )))]
    fn props_si_without_an_engine_is_unavailable_not_absent() {
        let (status, out) = props("Dmolar", "T", 300.0, "P", 101325.0, "Water");
        assert_eq!(status, RUSTPROP_UNAVAILABLE);
        assert!(out.is_nan(), "`out` must not be written");
        let msg = last_message();
        assert!(
            msg.contains("--features"),
            "the message must tell the caller how to fix it, got {msg:?}"
        );
    }

    /// The humid-air half of the same contract.
    #[test]
    #[cfg(not(feature = "humid-air"))]
    fn ha_props_si_without_the_engine_is_unavailable_not_absent() {
        let (status, out) = ha("W", "T", 300.0, "P", 101325.0, "R", 0.5);
        assert_eq!(status, RUSTPROP_UNAVAILABLE);
        assert!(out.is_nan(), "`out` must not be written");
        assert!(last_message().contains("humid-air"));
    }

    /// And when it IS compiled in, it answers. The golden is upstream's
    /// `HAPropsSI("W", "T", 300, "P", 101325, "R", 0.5)`.
    #[test]
    #[cfg(feature = "humid-air")]
    fn ha_props_si_matches_the_rust_api() {
        let (status, v) = ha("W", "T", 300.0, "P", 101325.0, "R", 0.5);
        assert_eq!(status, RUSTPROP_OK, "{}", last_message());
        let expected = rustprop::ha_props_si("W", "T", 300.0, "P", 101325.0, "R", 0.5).unwrap();
        assert_eq!(
            v.to_bits(),
            expected.to_bits(),
            "C ABI must not reshape the number"
        );
    }

    /// The returned status and the recorded message must describe the SAME
    /// failure. They can disagree if the entry point evaluates every argument
    /// before deciding, because each failing `str_arg` overwrites the
    /// last-error slot while the returned code comes from a different one.
    #[test]
    fn status_and_message_describe_the_same_argument() {
        // `output` is invalid UTF-8 (INVALID_UTF8); `fluid` is NULL
        // (NULL_ARGUMENT). Two different codes, so a mismatch is visible.
        let bad_utf8: &[u8] = b"\xff\xfe\0";
        let ok = CString::new("T").unwrap();
        let mut out = 0.0;
        let status = unsafe {
            rustprop_props_si(
                bad_utf8.as_ptr() as *const c_char,
                ok.as_ptr(),
                300.0,
                ok.as_ptr(),
                101325.0,
                std::ptr::null(),
                &mut out,
            )
        };
        assert_eq!(
            status,
            rustprop_last_error_code(),
            "returned {status} but the last-error slot says {} — the caller \
             would read a message for a different argument",
            rustprop_last_error_code()
        );
        // And it should be the FIRST bad argument that is reported.
        assert_eq!(status, RUSTPROP_INVALID_UTF8);
        assert!(
            last_message().contains("output"),
            "got {:?}",
            last_message()
        );
    }

    /// Every status code must have a name, including one we never defined.
    #[test]
    fn every_status_has_a_string() {
        for s in [
            RUSTPROP_OK,
            RUSTPROP_SOLUTION,
            RUSTPROP_VALUE,
            RUSTPROP_UNAVAILABLE,
            RUSTPROP_PANIC,
            RUSTPROP_ERROR,
            -12345,
        ] {
            let p = rustprop_status_string(s);
            assert!(!p.is_null());
            assert!(!unsafe { CStr::from_ptr(p) }.to_bytes().is_empty());
        }
    }
}
