/*
 * Compiles against include/rustprop.h and links the built library. This is
 * the check that the hand-written header and the Rust source agree: a
 * signature that drifted is a compile or link error here, not a surprise in
 * somebody's program.
 *
 * It is also the worked example. Every call a consumer needs is below, in the
 * order you would actually write them.
 *
 * Build (see ctest.sh, which is what CI runs):
 *   cc -I include smoke.c -L <libdir> -lrustprop -o smoke && ./smoke
 *
 * Exits 0 if every check passed, 1 otherwise.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "rustprop.h"

static int failures = 0;

static void check(int ok, const char *what) {
    printf("  %-58s %s\n", what, ok ? "ok" : "FAILED");
    if (!ok) failures++;
}

/* Relative agreement, the way the Rust suites state it. */
static int close_to(double got, double want, double rtol) {
    if (want == 0.0) return fabs(got) <= rtol;
    return fabs((got - want) / want) <= rtol;
}

static void print_last_error(const char *prefix) {
    /* The two-call idiom: ask for the length, then fetch. */
    size_t need = rustprop_last_error_message(NULL, 0);
    char buf[512];
    rustprop_last_error_message(buf, sizeof buf);
    printf("    %s [%s] %s%s\n", prefix,
           rustprop_status_string(rustprop_last_error_code()), buf,
           need >= sizeof buf ? " (truncated)" : "");
}

int main(void) {
    printf("rustprop %s (CoolProp %s)\n", rustprop_version(),
           rustprop_upstream_version());
    printf("backends: [%s]\n", rustprop_backends());
    printf("heos fluids compiled in: %zu\n\n", rustprop_fluid_count());

    /* --- introspection ------------------------------------------------ */
    printf("introspection\n");
    check(rustprop_version()[0] != '\0', "rustprop_version() is non-empty");
    check(strcmp(rustprop_upstream_version(), "8.0.0") == 0,
          "tracks CoolProp 8.0.0");
    check(rustprop_has_backend("no-such-backend") == 0,
          "an unknown backend is absent");
    check(rustprop_has_backend(NULL) == 0, "NULL backend name is absent");
    check(strcmp(rustprop_status_string(RUSTPROP_OK), "ok") == 0,
          "status_string(OK) == \"ok\"");
    check(rustprop_status_string(-424242) != NULL,
          "status_string of an unknown code is still a string");

    /* --- PropsSI ------------------------------------------------------ */
    printf("\nPropsSI\n");
    if (rustprop_has_backend("heos")) {
        double d = 0.0;
        int rc = rustprop_props_si("Dmolar", "T", 300.0, "P", 101325.0, "Water", &d);
        if (rc != RUSTPROP_OK) print_last_error("unexpected:");
        check(rc == RUSTPROP_OK, "PropsSI(Dmolar, T=300, P=101325, Water)");
        /* The golden the README and the Rust doctest both quote. */
        check(close_to(d, 55317.35277350119, 1e-8), "  == 55317.35277350119");

        /* Enumeration must round-trip: a name we are told about must work. */
        size_t n = rustprop_fluid_count();
        check(n > 0, "at least one fluid is compiled in");
        check(rustprop_fluid_name(n) == NULL, "index == count is out of range");
        int all_resolve = 1;
        for (size_t i = 0; i < n; i++) {
            const char *name = rustprop_fluid_name(i);
            double t = 0.0;
            if (name == NULL ||
                rustprop_props_si("Tcrit", "T", 300.0, "P", 101325.0, name, &t)
                    != RUSTPROP_OK) {
                printf("    fluid %zu (%s) did not resolve\n", i,
                       name ? name : "NULL");
                all_resolve = 0;
                break;
            }
        }
        check(all_resolve, "every enumerated fluid answers Tcrit");
    } else {
        /* The contract that matters for a prebuilt library: the symbol is
         * here even when the engine is not, so this LINKS and returns a
         * status rather than failing to load. */
        double d = -1.0;
        int rc = rustprop_props_si("Dmolar", "T", 300.0, "P", 101325.0, "Water", &d);
        check(rc == RUSTPROP_UNAVAILABLE, "absent engine reports UNAVAILABLE");
        check(d == -1.0, "  and leaves *out alone");
        print_last_error("message:");
    }

    /* --- IF97 --------------------------------------------------------- */
    if (rustprop_has_backend("if97")) {
        printf("\nIF97\n");
        double h = 0.0;
        int rc = rustprop_props_si("H", "T", 300.0, "P", 101325.0, "IF97::Water", &h);
        check(rc == RUSTPROP_OK, "PropsSI(H, T=300, P=101325, IF97::Water)");
        check(close_to(h, 112665.04341853978, 1e-11), "  == 112665.04341853978");
    }

    /* --- errors ------------------------------------------------------- */
    printf("\nerrors\n");
    {
        double d = -1.0;
        int rc = rustprop_props_si("Dmolar", "T", 300.0, "P", 101325.0,
                                   "NoSuchFluid", &d);
        check(rc != RUSTPROP_OK, "an unknown fluid fails");
        check(d == -1.0, "  and leaves *out alone");
        check(rustprop_last_error_code() == rc, "  last_error_code agrees");
        size_t need = rustprop_last_error_message(NULL, 0);
        check(need > 0, "  and carries a message");
        print_last_error("message:");

        /* Truncation must be reported and terminated, never overrun. */
        char tiny[8];
        memset(tiny, 0x7F, sizeof tiny);
        size_t reported = rustprop_last_error_message(tiny, sizeof tiny);
        check(reported == need, "  truncation reports the FULL length");
        check(tiny[sizeof tiny - 1] == '\0', "  and NUL-terminates in-buffer");

        /* A success must clear the slot. */
        if (rustprop_has_backend("heos")) {
            rustprop_props_si("Dmolar", "T", 300.0, "P", 101325.0, "Water", &d);
            check(rustprop_last_error_code() == RUSTPROP_OK,
                  "  a later success clears the error");
        }
    }
    {
        double d = 0.0;
        int rc = rustprop_props_si(NULL, "T", 300.0, "P", 101325.0, "Water", &d);
        check(rc == RUSTPROP_NULL_ARGUMENT, "NULL string is refused, not a crash");
        rc = rustprop_props_si("Dmolar", "T", 300.0, "P", 101325.0, "Water", NULL);
        check(rc == RUSTPROP_NULL_ARGUMENT, "NULL out is refused, not a crash");
    }

    /* --- batch -------------------------------------------------------- */
    if (rustprop_has_backend("heos")) {
        printf("\nbatch\n");
        double temps[5] = {300.0, 400.0, 500.0, -1.0 /* deliberately bad */, 600.0};
        double press[5] = {101325.0, 101325.0, 101325.0, 101325.0, 101325.0};
        double out[5] = {0};
        int rc = rustprop_props_si_many("Dmolar", "T", temps, "P", press,
                                        "Water", 5, out);
        check(rc == RUSTPROP_OK, "a bad cell does not fail the call");
        check(isnan(out[3]), "  the bad cell is NaN");

        int agree = 1;
        for (int i = 0; i < 5; i++) {
            if (i == 3) continue;
            double one = 0.0;
            rustprop_props_si("Dmolar", "T", temps[i], "P", press[i], "Water", &one);
            if (out[i] != one) { agree = 0; break; }
        }
        check(agree, "  every good cell equals the scalar answer exactly");

        rc = rustprop_props_si_many("Dmolar", "T", NULL, "P", NULL, "Water", 0, NULL);
        check(rc == RUSTPROP_OK, "n = 0 is a no-op, not a NULL complaint");
    }

    /* --- humid air ---------------------------------------------------- */
    printf("\nhumid air\n");
    {
        double w = 0.0;
        int rc = rustprop_ha_props_si("W", "T", 300.0, "P", 101325.0, "R", 0.5, &w);
        if (rustprop_has_backend("humid-air")) {
            check(rc == RUSTPROP_OK, "HAPropsSI(W, T=300, P=101325, R=0.5)");
            check(w > 0.0 && w < 1.0, "  humidity ratio is physical");
        } else {
            check(rc == RUSTPROP_UNAVAILABLE, "absent humid-air reports UNAVAILABLE");
        }
    }

    printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures,
           failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
