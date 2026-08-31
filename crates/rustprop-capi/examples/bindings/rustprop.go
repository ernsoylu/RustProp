// rustprop from Go, via cgo.
//
// Build with the SDK's include and lib directories on the cgo flags:
//
//	CGO_CFLAGS="-I<sdk>/include" \
//	CGO_LDFLAGS="-L<sdk>/lib -lrustprop" \
//	LD_LIBRARY_PATH=<sdk>/lib go run rustprop.go
//
// Or embed the paths in the #cgo lines below with ${SRCDIR}.
package main

/*
#include <stdlib.h>
#include "rustprop.h"
*/
import "C"

import (
	"errors"
	"fmt"
	"math"
	"os"
	"strings"
	"unsafe"
)

// Error carries the status code and the message rustprop recorded for it.
type Error struct {
	Status  int
	Kind    string
	Message string
}

func (e *Error) Error() string { return e.Kind + ": " + e.Message }

// lastError must be called on the same goroutine, and before any other
// rustprop call: the message slot is thread-local, and Go may move a
// goroutine between OS threads at a call boundary. Reading it immediately
// after the failing call, as every wrapper here does, is the safe pattern.
func lastError(status C.int) *Error {
	need := C.rustprop_last_error_message(nil, 0)
	buf := (*C.char)(C.malloc(C.size_t(need) + 1))
	defer C.free(unsafe.Pointer(buf))
	C.rustprop_last_error_message(buf, need+1)
	return &Error{
		Status:  int(status),
		Kind:    C.GoString(C.rustprop_status_string(status)),
		Message: C.GoString(buf),
	}
}

func cstr(s string) (*C.char, func()) {
	p := C.CString(s)
	return p, func() { C.free(unsafe.Pointer(p)) }
}

// PropsSI computes one property at one state.
func PropsSI(output, name1 string, val1 float64, name2 string, val2 float64, fluid string) (float64, error) {
	co, f1 := cstr(output)
	defer f1()
	c1, f2 := cstr(name1)
	defer f2()
	c2, f3 := cstr(name2)
	defer f3()
	cf, f4 := cstr(fluid)
	defer f4()

	var out C.double
	rc := C.rustprop_props_si(co, c1, C.double(val1), c2, C.double(val2), cf, &out)
	if rc != C.RUSTPROP_OK {
		return 0, lastError(rc)
	}
	return float64(out), nil
}

// HAPropsSI computes one humid-air property.
func HAPropsSI(output, n1 string, v1 float64, n2 string, v2 float64, n3 string, v3 float64) (float64, error) {
	co, f1 := cstr(output)
	defer f1()
	c1, f2 := cstr(n1)
	defer f2()
	c2, f3 := cstr(n2)
	defer f3()
	c3, f4 := cstr(n3)
	defer f4()

	var out C.double
	rc := C.rustprop_ha_props_si(co, c1, C.double(v1), c2, C.double(v2), c3, C.double(v3), &out)
	if rc != C.RUSTPROP_OK {
		return 0, lastError(rc)
	}
	return float64(out), nil
}

// PropsSIMany computes one output over many states. A state that fails yields
// NaN in its slot rather than failing the call.
func PropsSIMany(output, name1 string, vals1 []float64, name2 string, vals2 []float64, fluid string) ([]float64, error) {
	if len(vals1) != len(vals2) {
		return nil, errors.New("rustprop: input slices differ in length")
	}
	if len(vals1) == 0 {
		return nil, nil
	}
	co, f1 := cstr(output)
	defer f1()
	c1, f2 := cstr(name1)
	defer f2()
	c2, f3 := cstr(name2)
	defer f3()
	cf, f4 := cstr(fluid)
	defer f4()

	out := make([]float64, len(vals1))
	rc := C.rustprop_props_si_many(co, c1, (*C.double)(&vals1[0]), c2,
		(*C.double)(&vals2[0]), cf, C.size_t(len(vals1)), (*C.double)(&out[0]))
	if rc != C.RUSTPROP_OK {
		return nil, lastError(rc)
	}
	return out, nil
}

func Backends() []string {
	s := C.GoString(C.rustprop_backends())
	if s == "" {
		return nil
	}
	return strings.Split(s, ",")
}

func HasBackend(name string) bool {
	c, free := cstr(name)
	defer free()
	return C.rustprop_has_backend(c) == 1
}

func Fluids() []string {
	n := int(C.rustprop_fluid_count())
	out := make([]string, 0, n)
	for i := 0; i < n; i++ {
		out = append(out, C.GoString(C.rustprop_fluid_name(C.size_t(i))))
	}
	return out
}

func Version() string         { return C.GoString(C.rustprop_version()) }
func UpstreamVersion() string { return C.GoString(C.rustprop_upstream_version()) }

func main() {
	fmt.Printf("rustprop %s (CoolProp %s)\n", Version(), UpstreamVersion())
	fmt.Printf("backends: %s\n", strings.Join(Backends(), ","))
	fmt.Printf("fluids compiled in: %d\n", len(Fluids()))

	failures := 0
	check := func(ok bool, what string) {
		status := "ok"
		if !ok {
			status = "FAILED"
			failures++
		}
		fmt.Printf("  %-52s %s\n", what, status)
	}

	if HasBackend("heos") {
		d, err := PropsSI("Dmolar", "T", 300, "P", 101325, "Water")
		check(err == nil && math.Abs((d-55317.35277350119)/d) < 1e-8,
			fmt.Sprintf("PropsSI Dmolar Water = %.10g", d))

		temps := []float64{300, 400, 500}
		press := []float64{101325, 101325, 101325}
		many, err := PropsSIMany("Dmolar", "T", temps, "P", press, "Water")
		same := err == nil
		for i, t := range temps {
			one, _ := PropsSI("Dmolar", "T", t, "P", press[i], "Water")
			if many[i] != one {
				same = false
			}
		}
		check(same, "batch equals scalar exactly")
	}

	if HasBackend("if97") {
		h, err := PropsSI("H", "T", 300, "P", 101325, "IF97::Water")
		check(err == nil && math.Abs((h-112665.04341853978)/h) < 1e-11,
			fmt.Sprintf("PropsSI H IF97::Water = %.12g", h))
	}

	if HasBackend("humid-air") {
		w, err := HAPropsSI("W", "T", 300, "P", 101325, "R", 0.5)
		check(err == nil && w > 0 && w < 1, fmt.Sprintf("HAPropsSI W = %.6g", w))
	}

	_, err := PropsSI("Dmolar", "T", 300, "P", 101325, "NoSuchFluid")
	check(err != nil && strings.Contains(err.Error(), "NoSuchFluid"),
		"a bad fluid errors, naming the key")

	if failures > 0 {
		fmt.Printf("\nFAILED (%d failures)\n", failures)
		os.Exit(1)
	}
	fmt.Printf("\nPASSED (0 failures)\n")
}
