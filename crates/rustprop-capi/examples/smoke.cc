// C++ side of the header check: the same header, compiled by a C++ compiler.
//
// Two things it proves that smoke.c cannot:
//
//   1. include/rustprop.h is usable from C++ at all — the extern "C" guard is
//      right, and nothing in it needs a C-only construct.
//   2. The thread-safety the docs promise survives contact with real threads.
//      `rustprop` keeps per-fluid caches in `static`s, so the claim is worth
//      checking from the outside rather than only reasoning about: this runs
//      the same states on eight threads and demands bit-identical answers.
//
// It doubles as the idiomatic-C++ example — the small wrapper below is what
// a C++ consumer would actually write, and is the recommended shape:
// exceptions on the scalar path, status codes left to the batch path.
//
// Build (see ctest.sh, which is what CI runs):
//   c++ -std=c++17 -I include smoke.cc -L <libdir> -lrustprop -o smoke_cc

#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "rustprop.h"

namespace rustprop {

// The message for whatever just failed on this thread.
inline std::string last_error() {
    std::string s(rustprop_last_error_message(nullptr, 0), '\0');
    if (!s.empty()) rustprop_last_error_message(&s[0], s.size() + 1);
    return s;
}

struct Error : std::runtime_error {
    int status;
    explicit Error(int status_)
        : std::runtime_error(std::string(rustprop_status_string(status_)) + ": " +
                             last_error()),
          status(status_) {}
};

inline double props_si(const char* output, const char* name1, double val1,
                       const char* name2, double val2, const char* fluid) {
    double out = 0.0;
    int rc = rustprop_props_si(output, name1, val1, name2, val2, fluid, &out);
    if (rc != RUSTPROP_OK) throw Error(rc);
    return out;
}

inline double ha_props_si(const char* output, const char* n1, double v1,
                          const char* n2, double v2, const char* n3, double v3) {
    double out = 0.0;
    int rc = rustprop_ha_props_si(output, n1, v1, n2, v2, n3, v3, &out);
    if (rc != RUSTPROP_OK) throw Error(rc);
    return out;
}

// The batch path keeps its status-code shape: a failing cell is a NaN, not an
// exception, so a sweep that clips the phase envelope still returns.
inline std::vector<double> props_si_many(const char* output, const char* name1,
                                         const std::vector<double>& vals1,
                                         const char* name2,
                                         const std::vector<double>& vals2,
                                         const char* fluid) {
    if (vals1.size() != vals2.size())
        throw std::invalid_argument("props_si_many: input vectors differ in length");
    std::vector<double> out(vals1.size());
    int rc = rustprop_props_si_many(output, name1, vals1.data(), name2, vals2.data(),
                                    fluid, vals1.size(), out.data());
    if (rc != RUSTPROP_OK) throw Error(rc);
    return out;
}

inline bool has_backend(const char* name) { return rustprop_has_backend(name) != 0; }

}  // namespace rustprop

static int failures = 0;

static void check(bool ok, const char* what) {
    std::printf("  %-58s %s\n", what, ok ? "ok" : "FAILED");
    if (!ok) failures++;
}

static bool close_to(double got, double want, double rtol) {
    return want == 0.0 ? std::fabs(got) <= rtol
                       : std::fabs((got - want) / want) <= rtol;
}

int main() {
    std::printf("rustprop %s (CoolProp %s) via C++\n", rustprop_version(),
                rustprop_upstream_version());
    std::printf("backends: [%s]\n\n", rustprop_backends());

    std::printf("scalar\n");
    if (rustprop::has_backend("heos")) {
        double d = rustprop::props_si("Dmolar", "T", 300.0, "P", 101325.0, "Water");
        check(close_to(d, 55317.35277350119, 1e-8),
              "PropsSI(Dmolar, T=300, P=101325, Water)");
    }
    if (rustprop::has_backend("if97")) {
        double h = rustprop::props_si("H", "T", 300.0, "P", 101325.0, "IF97::Water");
        check(close_to(h, 112665.04341853978, 1e-11),
              "PropsSI(H, T=300, P=101325, IF97::Water)");
    }

    std::printf("\nexceptions\n");
    {
        bool threw = false;
        std::string what;
        try {
            rustprop::props_si("Dmolar", "T", 300.0, "P", 101325.0, "NoSuchFluid");
        } catch (const rustprop::Error& e) {
            threw = true;
            what = e.what();
        }
        check(threw, "a bad fluid throws rustprop::Error");
        check(what.find("NoSuchFluid") != std::string::npos,
              "  and the message names the offending key");
        std::printf("    %s\n", what.c_str());
    }

    if (rustprop::has_backend("heos")) {
        std::printf("\nbatch\n");
        std::vector<double> t, p;
        for (int i = 0; i < 64; i++) {
            t.push_back(280.0 + i);
            p.push_back(101325.0);
        }
        auto out = rustprop::props_si_many("Dmolar", "T", t, "P", p, "Water");
        check(out.size() == t.size(), "batch returns one result per state");
        bool agree = true;
        for (size_t i = 0; i < out.size(); i++) {
            if (out[i] != rustprop::props_si("Dmolar", "T", t[i], "P", p[i], "Water")) {
                agree = false;
                break;
            }
        }
        check(agree, "  every cell equals the scalar answer exactly");

        // The claim under test: concurrent calls are safe AND identical.
        std::printf("\nthreads\n");
        std::vector<std::vector<double>> results(8);
        std::vector<std::thread> pool;
        for (int k = 0; k < 8; k++) {
            pool.emplace_back([&results, &t, &p, k] {
                results[k] = rustprop::props_si_many("Dmolar", "T", t, "P", p, "Water");
            });
        }
        for (auto& th : pool) th.join();
        bool identical = true;
        for (auto& r : results)
            if (r != out) identical = false;
        check(identical, "8 threads agree bit-for-bit with the serial answer");

        // The last-error slot is thread-local: a failure over there must not
        // be visible over here.
        try {
            rustprop::props_si("Dmolar", "T", 300.0, "P", 101325.0, "NoSuchFluid");
        } catch (const rustprop::Error&) {
        }
        int seen_elsewhere = RUSTPROP_PANIC;
        std::thread([&seen_elsewhere] {
            seen_elsewhere = rustprop_last_error_code();
        }).join();
        check(seen_elsewhere == RUSTPROP_OK, "an error does not leak across threads");
    }

    if (rustprop::has_backend("humid-air")) {
        std::printf("\nhumid air\n");
        double w = rustprop::ha_props_si("W", "T", 300.0, "P", 101325.0, "R", 0.5);
        check(w > 0.0 && w < 1.0, "HAPropsSI(W, T=300, P=101325, R=0.5) is physical");
    }

    std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures,
                failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
