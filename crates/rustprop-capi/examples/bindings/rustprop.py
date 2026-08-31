"""rustprop from Python, via ctypes — no build step, no extension module.

Set RUSTPROP_LIB to the shared library, or drop it beside this file.
Run directly to check your copy:  python3 rustprop.py
"""

import ctypes
import os
import sys
from ctypes import POINTER, c_char_p, c_double, c_int, c_size_t, create_string_buffer

OK = 0
UNAVAILABLE = 102


def _find_library():
    if "RUSTPROP_LIB" in os.environ:
        return os.environ["RUSTPROP_LIB"]
    names = {"linux": "librustprop.so", "darwin": "librustprop.dylib", "win32": "rustprop.dll"}
    name = names.get(sys.platform, "librustprop.so")
    here = os.path.join(os.path.dirname(os.path.abspath(__file__)), name)
    return here if os.path.exists(here) else name


_lib = ctypes.CDLL(_find_library())

# Declaring argtypes/restype is not optional: without them ctypes assumes int
# arguments, and every double you pass arrives as garbage.
_lib.rustprop_props_si.argtypes = [c_char_p, c_char_p, c_double, c_char_p, c_double,
                                   c_char_p, POINTER(c_double)]
_lib.rustprop_props_si.restype = c_int
_lib.rustprop_ha_props_si.argtypes = [c_char_p, c_char_p, c_double, c_char_p, c_double,
                                      c_char_p, c_double, POINTER(c_double)]
_lib.rustprop_ha_props_si.restype = c_int
_lib.rustprop_props_si_many.argtypes = [c_char_p, c_char_p, POINTER(c_double), c_char_p,
                                        POINTER(c_double), c_char_p, c_size_t,
                                        POINTER(c_double)]
_lib.rustprop_props_si_many.restype = c_int
_lib.rustprop_last_error_message.argtypes = [c_char_p, c_size_t]
_lib.rustprop_last_error_message.restype = c_size_t
_lib.rustprop_status_string.argtypes = [c_int]
_lib.rustprop_status_string.restype = c_char_p
_lib.rustprop_backends.restype = c_char_p
_lib.rustprop_version.restype = c_char_p
_lib.rustprop_upstream_version.restype = c_char_p
_lib.rustprop_has_backend.argtypes = [c_char_p]
_lib.rustprop_has_backend.restype = c_int
_lib.rustprop_fluid_count.restype = c_size_t
_lib.rustprop_fluid_name.argtypes = [c_size_t]
_lib.rustprop_fluid_name.restype = c_char_p


class RustpropError(Exception):
    def __init__(self, status):
        self.status = status
        kind = _lib.rustprop_status_string(status).decode()
        need = _lib.rustprop_last_error_message(None, 0)
        buf = create_string_buffer(need + 1)
        _lib.rustprop_last_error_message(buf, need + 1)
        super().__init__(f"{kind}: {buf.value.decode()}")


def _enc(s):
    return s.encode("utf-8")


def props_si(output, name1, val1, name2, val2, fluid):
    """PropsSI(output, name1, val1, name2, val2, fluid)."""
    out = c_double()
    rc = _lib.rustprop_props_si(_enc(output), _enc(name1), val1, _enc(name2), val2,
                                _enc(fluid), ctypes.byref(out))
    if rc != OK:
        raise RustpropError(rc)
    return out.value


def ha_props_si(output, n1, v1, n2, v2, n3, v3):
    """HAPropsSI — humid air."""
    out = c_double()
    rc = _lib.rustprop_ha_props_si(_enc(output), _enc(n1), v1, _enc(n2), v2,
                                   _enc(n3), v3, ctypes.byref(out))
    if rc != OK:
        raise RustpropError(rc)
    return out.value


def props_si_many(output, name1, vals1, name2, vals2, fluid):
    """One output over many states. A failing state comes back as nan."""
    n = len(vals1)
    if n != len(vals2):
        raise ValueError("vals1 and vals2 must be the same length")
    arr = c_double * n
    a, b, out = arr(*vals1), arr(*vals2), arr()
    rc = _lib.rustprop_props_si_many(_enc(output), _enc(name1), a, _enc(name2), b,
                                     _enc(fluid), n, out)
    if rc != OK:
        raise RustpropError(rc)
    return list(out)


def backends():
    return _lib.rustprop_backends().decode().split(",") if _lib.rustprop_backends() else []


def has_backend(name):
    return _lib.rustprop_has_backend(_enc(name)) == 1


def fluids():
    return [_lib.rustprop_fluid_name(i).decode() for i in range(_lib.rustprop_fluid_count())]


def version():
    return _lib.rustprop_version().decode()


def upstream_version():
    return _lib.rustprop_upstream_version().decode()


if __name__ == "__main__":
    print(f"rustprop {version()} (CoolProp {upstream_version()})")
    print(f"backends: {','.join(backends())}")
    print(f"fluids compiled in: {len(fluids())}")

    failures = 0

    def check(ok, what):
        global failures
        print(f"  {what:<52} {'ok' if ok else 'FAILED'}")
        if not ok:
            failures += 1

    if has_backend("heos"):
        d = props_si("Dmolar", "T", 300.0, "P", 101325.0, "Water")
        check(abs((d - 55317.35277350119) / d) < 1e-8, f"PropsSI Dmolar Water = {d:.10g}")

        many = props_si_many("Dmolar", "T", [300.0, 400.0, 500.0], "P", [101325.0] * 3, "Water")
        one = [props_si("Dmolar", "T", t, "P", 101325.0, "Water") for t in (300.0, 400.0, 500.0)]
        check(many == one, "batch equals scalar exactly")

    if has_backend("if97"):
        h = props_si("H", "T", 300.0, "P", 101325.0, "IF97::Water")
        check(abs((h - 112665.04341853978) / h) < 1e-11, f"PropsSI H IF97::Water = {h:.12g}")

    if has_backend("humid-air"):
        w = ha_props_si("W", "T", 300.0, "P", 101325.0, "R", 0.5)
        check(0.0 < w < 1.0, f"HAPropsSI W = {w:.6g}")

    try:
        props_si("Dmolar", "T", 300.0, "P", 101325.0, "NoSuchFluid")
        check(False, "a bad fluid raises")
    except RustpropError as e:
        check("NoSuchFluid" in str(e), "a bad fluid raises, naming the key")

    print(f"\n{'FAILED' if failures else 'PASSED'} ({failures} failures)")
    sys.exit(1 if failures else 0)
