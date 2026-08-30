/*
 * rustprop — thermophysical properties, CoolProp 8.0.0 semantics, pure Rust.
 *
 * C ABI. Also the interface for C++, Python (ctypes/cffi), C#, Go (cgo),
 * Julia, MATLAB, Fortran (iso_c_binding) — anything that can call a C
 * function.
 *
 *   https://github.com/ernsoylu/RustProp
 *
 * SPDX-License-Identifier: MIT
 * Derivative work of CoolProp (MIT, (c) 2012-2018 Ian H. Bell and the
 * CoolProp developers).
 *
 * -------------------------------------------------------------------------
 * THIS HEADER IS HAND-WRITTEN, NOT GENERATED.
 *
 * It is the contract. `crates/rustprop-capi/src/lib.rs` is the other half of
 * it, and a change to a signature there is a change here in the same edit.
 * What keeps the two honest is not a diff: `examples/smoke.c` and
 * `examples/smoke.cc` compile against THIS file and link the built library,
 * so a mismatch is a link error in CI rather than a surprise in your program.
 * -------------------------------------------------------------------------
 *
 * WHICH ENGINES DOES MY COPY HAVE?
 *
 * rustprop selects calculation engines and fluid data at COMPILE time, so a
 * build carries only what was asked for — that is the point of the project,
 * and it is why a 128 KB build and a 4.2 MB build are both "rustprop".
 *
 * Every function below is exported by EVERY build regardless. A call into an
 * engine your copy does not carry returns RUSTPROP_UNAVAILABLE; it never
 * fails to link. Ask rustprop_backends() or rustprop_has_backend() what you
 * actually have. Release artifacts named `all-backends` carry everything.
 *
 * ERRORS
 *
 * Every calculation returns an int status; RUSTPROP_OK (0) means the `out`
 * argument was written, and nothing else writes it. The matching human-
 * readable message is rustprop_last_error_message(), which is per-thread.
 *
 * THREADS
 *
 * Every function here is safe to call concurrently from any number of
 * threads, with no initialisation call and no lock held by the caller. The
 * last-error slot is thread-local, so a failure on one thread is invisible to
 * another.
 *
 * MEMORY
 *
 * No function here ever returns memory the caller must free. `const char *`
 * returns are static or process-lifetime; everything else is written into
 * buffers you supply.
 */

#ifndef RUSTPROP_H
#define RUSTPROP_H

#include <stddef.h> /* size_t */

#ifdef __cplusplus
extern "C" {
#endif

/* Exported-symbol decoration. Consumers of the shared library on Windows
 * need __declspec(dllimport); the static library and every ELF/Mach-O target
 * need nothing. Define RUSTPROP_STATIC when linking rustprop.lib statically. */
#if defined(_WIN32) && !defined(RUSTPROP_STATIC)
#  define RUSTPROP_API __declspec(dllimport)
#else
#  define RUSTPROP_API
#endif

/* ------------------------------------------------------------------ */
/* Status codes                                                        */
/* ------------------------------------------------------------------ */

/* 1..11 correspond one-to-one with CoolProp's exception types, which this
 * port reproduces as error parity: the same refusal, in the same state, for
 * the same reason. 100+ exist only at this boundary. */
#define RUSTPROP_OK                   0
#define RUSTPROP_NOT_IMPLEMENTED      1  /* NotImplementedError              */
#define RUSTPROP_SOLUTION             2  /* SolutionError: solver diverged   */
#define RUSTPROP_ATTRIBUTE            3  /* AttributeError                   */
#define RUSTPROP_OUT_OF_RANGE         4  /* OutOfRangeError                  */
#define RUSTPROP_VALUE                5  /* ValueError: bad name or value    */
#define RUSTPROP_WRONG_FLUID          6  /* WrongFluidError                  */
#define RUSTPROP_COMPOSITION          7  /* CompositionError                 */
#define RUSTPROP_INPUT                8  /* InputError                       */
#define RUSTPROP_NOT_AVAILABLE        9  /* property absent for this fluid   */
#define RUSTPROP_KEY                 10  /* KeyError                         */
#define RUSTPROP_MULTIPLE_SOLUTIONS  11  /* input maps to >1 state           */
#define RUSTPROP_ERROR               99  /* an error this ABI predates       */

#define RUSTPROP_NULL_ARGUMENT      100  /* a required pointer was NULL      */
#define RUSTPROP_INVALID_UTF8       101  /* a string argument was not UTF-8  */
#define RUSTPROP_UNAVAILABLE        102  /* engine not in THIS build         */
#define RUSTPROP_PANIC              104  /* caught at the boundary; a bug    */

/* A short static name for a status ("SolutionError", "Unavailable", ...).
 * Never NULL, including for a code this header does not define. */
RUSTPROP_API const char *rustprop_status_string(int status);

/* ------------------------------------------------------------------ */
/* Errors (per-thread)                                                 */
/* ------------------------------------------------------------------ */

/* The status of the most recent call on THIS thread. */
RUSTPROP_API int rustprop_last_error_code(void);

/* Copy the most recent error message on this thread into `buf`.
 *
 * Returns the message length in bytes, NOT counting the NUL. So:
 *
 *   - a return value >= len means it was truncated; call again with
 *     (returned + 1) bytes to get all of it;
 *   - what is written is always NUL-terminated when len > 0, and is always
 *     valid UTF-8 (a truncation never splits a character);
 *   - buf = NULL / len = 0 asks for the length and writes nothing.
 *
 * Example:
 *
 *   char msg[256];
 *   if (rustprop_props_si("D", "T", 300, "P", 101325, "Water", &d)) {
 *       rustprop_last_error_message(msg, sizeof msg);
 *       fprintf(stderr, "rustprop: %s\n", msg);
 *   }
 */
RUSTPROP_API size_t rustprop_last_error_message(char *buf, size_t len);

/* ------------------------------------------------------------------ */
/* What is in this build                                               */
/* ------------------------------------------------------------------ */

/* This library's version, e.g. "0.1.0". Static, never NULL. */
RUSTPROP_API const char *rustprop_version(void);

/* The upstream CoolProp release these numbers reproduce, e.g. "8.0.0".
 * Static, never NULL. */
RUSTPROP_API const char *rustprop_upstream_version(void);

/* The backends compiled into THIS binary, comma-separated with no trailing
 * separator, e.g. "heos,if97,humid-air". Empty string if none. Never NULL,
 * and valid for the life of the process. */
RUSTPROP_API const char *rustprop_backends(void);

/* 1 if this binary carries `name`, else 0 (including for NULL or an unknown
 * name). Names are: "heos", "heos-mixtures", "if97", "cubics",
 * "incompressible", "pcsaft", "tabular", "svdsbtl", "humid-air". */
RUSTPROP_API int rustprop_has_backend(const char *name);

/* How many HEOS fluids this binary carries data for; 0 without "heos".
 * Fluid data is opt-in per fluid, so this is as build-dependent as the
 * backend list. */
RUSTPROP_API size_t rustprop_fluid_count(void);

/* The name of compiled-in HEOS fluid `index`, or NULL when out of range.
 * Sorted, stable, and valid for the life of the process. Iterate
 * 0 .. rustprop_fluid_count() to see everything this binary can answer for;
 * every name it reports is a name rustprop_props_si() accepts. */
RUSTPROP_API const char *rustprop_fluid_name(size_t index);

/* ------------------------------------------------------------------ */
/* Calculations                                                        */
/* ------------------------------------------------------------------ */

/* PropsSI(output, name1, val1, name2, val2, fluid).
 *
 * On RUSTPROP_OK, *out holds the result. On any other status *out is left
 * exactly as it was, and rustprop_last_error_message() says why.
 *
 *   double d;
 *   int rc = rustprop_props_si("Dmolar", "T", 300, "P", 101325, "Water", &d);
 *   // d == 55317.35277350119
 *
 * `fluid` takes the upstream spellings, including the backend prefix:
 * "Water", "R134a", "IF97::Water", "SRK::Propane", "HEOS::Methane&Ethane".
 */
RUSTPROP_API int rustprop_props_si(const char *output,
                                   const char *name1, double val1,
                                   const char *name2, double val2,
                                   const char *fluid,
                                   double *out);

/* PropsSI over `n` states at once: one output, one fluid, two input vectors,
 * `n` results written to `out`.
 *
 * Worth using when you are sweeping a grid: the names are parsed once
 * instead of n times, which for a large sweep dominates the thermodynamics.
 *
 * A state that fails yields NaN in its slot and does NOT stop the batch — a
 * sweep routinely clips the phase envelope, and one bad cell should not cost
 * you the other ten thousand. The returned status therefore reports only
 * whether the CALL was well-formed; check individual cells with isnan().
 *
 * `vals1`, `vals2` and `out` must each hold `n` doubles, and `out` must not
 * overlap the inputs. n = 0 is valid and does nothing.
 */
RUSTPROP_API int rustprop_props_si_many(const char *output,
                                        const char *name1, const double *vals1,
                                        const char *name2, const double *vals2,
                                        const char *fluid,
                                        size_t n,
                                        double *out);

/* HAPropsSI(output, name1, val1, name2, val2, name3, val3) — humid air /
 * psychrometrics. Same contract as rustprop_props_si.
 *
 *   double w;
 *   rustprop_ha_props_si("W", "T", 300, "P", 101325, "R", 0.5, &w);
 */
RUSTPROP_API int rustprop_ha_props_si(const char *output,
                                      const char *name1, double val1,
                                      const char *name2, double val2,
                                      const char *name3, double val3,
                                      double *out);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* RUSTPROP_H */
