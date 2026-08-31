// rustprop from Java, via the Foreign Function & Memory API (stable in Java 22+).
//
// No JNI and no wrapper library to compile — the JVM calls the C ABI directly.
//
//   RUSTPROP_LIB=<sdk>/lib/librustprop.so \
//     java --enable-native-access=ALL-UNNAMED Rustprop.java
//
// Or pass -Drustprop.lib=<path>.

import java.lang.foreign.Arena;
import java.lang.foreign.FunctionDescriptor;
import java.lang.foreign.Linker;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.SymbolLookup;
import java.lang.foreign.ValueLayout;
import java.lang.invoke.MethodHandle;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public final class Rustprop implements AutoCloseable {

    public static final int OK = 0;
    public static final int UNAVAILABLE = 102;

    private final Arena arena = Arena.ofShared();
    private final MethodHandle propsSi, haPropsSi, propsSiMany, lastErrorMessage,
            statusString, backends, version, upstreamVersion, hasBackend,
            fluidCount, fluidName;

    /** Thrown for any non-OK status, carrying rustprop's own message. */
    public static final class RustpropException extends RuntimeException {
        public final int status;
        RustpropException(int status, String kind, String message) {
            super(kind + ": " + message);
            this.status = status;
        }
    }

    public Rustprop(Path library) {
        var linker = Linker.nativeLinker();
        var lookup = SymbolLookup.libraryLookup(library, arena);
        var A = ValueLayout.ADDRESS;
        var D = ValueLayout.JAVA_DOUBLE;
        var I = ValueLayout.JAVA_INT;
        var L = ValueLayout.JAVA_LONG; // size_t, on a 64-bit platform

        propsSi = down(linker, lookup, "rustprop_props_si",
                FunctionDescriptor.of(I, A, A, D, A, D, A, A));
        haPropsSi = down(linker, lookup, "rustprop_ha_props_si",
                FunctionDescriptor.of(I, A, A, D, A, D, A, D, A));
        propsSiMany = down(linker, lookup, "rustprop_props_si_many",
                FunctionDescriptor.of(I, A, A, A, A, A, A, L, A));
        lastErrorMessage = down(linker, lookup, "rustprop_last_error_message",
                FunctionDescriptor.of(L, A, L));
        statusString = down(linker, lookup, "rustprop_status_string",
                FunctionDescriptor.of(A, I));
        backends = down(linker, lookup, "rustprop_backends", FunctionDescriptor.of(A));
        version = down(linker, lookup, "rustprop_version", FunctionDescriptor.of(A));
        upstreamVersion = down(linker, lookup, "rustprop_upstream_version",
                FunctionDescriptor.of(A));
        hasBackend = down(linker, lookup, "rustprop_has_backend",
                FunctionDescriptor.of(I, A));
        fluidCount = down(linker, lookup, "rustprop_fluid_count", FunctionDescriptor.of(L));
        fluidName = down(linker, lookup, "rustprop_fluid_name", FunctionDescriptor.of(A, L));
    }

    private static MethodHandle down(Linker linker, SymbolLookup lookup, String name,
                                     FunctionDescriptor fd) {
        return linker.downcallHandle(
                lookup.find(name).orElseThrow(() ->
                        new UnsatisfiedLinkError("rustprop: missing symbol " + name)),
                fd);
    }

    /** A `const char *` return arrives with zero size; give it one before reading. */
    private static String str(MemorySegment p) {
        return p.equals(MemorySegment.NULL) ? null
                : p.reinterpret(Long.MAX_VALUE).getString(0);
    }

    private RustpropException error(int status) {
        try (var a = Arena.ofConfined()) {
            long need = (long) lastErrorMessage.invokeExact(MemorySegment.NULL, 0L);
            var buf = a.allocate(need + 1);
            long written = (long) lastErrorMessage.invokeExact(buf, need + 1);
            assert written == need;
            String kind = str((MemorySegment) statusString.invokeExact(status));
            return new RustpropException(status, kind, buf.getString(0));
        } catch (Throwable t) {
            throw new RuntimeException(t);
        }
    }

    public double propsSi(String output, String n1, double v1, String n2, double v2,
                          String fluid) {
        try (var a = Arena.ofConfined()) {
            var out = a.allocate(ValueLayout.JAVA_DOUBLE);
            int rc = (int) propsSi.invokeExact(a.allocateFrom(output), a.allocateFrom(n1), v1,
                    a.allocateFrom(n2), v2, a.allocateFrom(fluid), out);
            if (rc != OK) throw error(rc);
            return out.get(ValueLayout.JAVA_DOUBLE, 0);
        } catch (RustpropException e) {
            throw e;
        } catch (Throwable t) {
            throw new RuntimeException(t);
        }
    }

    public double haPropsSi(String output, String n1, double v1, String n2, double v2,
                            String n3, double v3) {
        try (var a = Arena.ofConfined()) {
            var out = a.allocate(ValueLayout.JAVA_DOUBLE);
            int rc = (int) haPropsSi.invokeExact(a.allocateFrom(output), a.allocateFrom(n1),
                    v1, a.allocateFrom(n2), v2, a.allocateFrom(n3), v3, out);
            if (rc != OK) throw error(rc);
            return out.get(ValueLayout.JAVA_DOUBLE, 0);
        } catch (RustpropException e) {
            throw e;
        } catch (Throwable t) {
            throw new RuntimeException(t);
        }
    }

    /** One output over many states; a state that fails comes back as NaN. */
    public double[] propsSiMany(String output, String n1, double[] v1, String n2,
                                double[] v2, String fluid) {
        if (v1.length != v2.length)
            throw new IllegalArgumentException("input arrays differ in length");
        if (v1.length == 0) return new double[0];
        try (var a = Arena.ofConfined()) {
            var arr1 = a.allocateFrom(ValueLayout.JAVA_DOUBLE, v1);
            var arr2 = a.allocateFrom(ValueLayout.JAVA_DOUBLE, v2);
            var out = a.allocate(ValueLayout.JAVA_DOUBLE, v1.length);
            int rc = (int) propsSiMany.invokeExact(a.allocateFrom(output), a.allocateFrom(n1),
                    arr1, a.allocateFrom(n2), arr2, a.allocateFrom(fluid),
                    (long) v1.length, out);
            if (rc != OK) throw error(rc);
            return out.toArray(ValueLayout.JAVA_DOUBLE);
        } catch (RustpropException e) {
            throw e;
        } catch (Throwable t) {
            throw new RuntimeException(t);
        }
    }

    public String version() {
        try { return str((MemorySegment) version.invokeExact()); }
        catch (Throwable t) { throw new RuntimeException(t); }
    }

    public String upstreamVersion() {
        try { return str((MemorySegment) upstreamVersion.invokeExact()); }
        catch (Throwable t) { throw new RuntimeException(t); }
    }

    public String backends() {
        try { return str((MemorySegment) backends.invokeExact()); }
        catch (Throwable t) { throw new RuntimeException(t); }
    }

    public boolean hasBackend(String name) {
        try (var a = Arena.ofConfined()) {
            return (int) hasBackend.invokeExact(a.allocateFrom(name)) == 1;
        } catch (Throwable t) { throw new RuntimeException(t); }
    }

    public List<String> fluids() {
        try {
            long n = (long) fluidCount.invokeExact();
            var out = new ArrayList<String>((int) n);
            for (long i = 0; i < n; i++)
                out.add(str((MemorySegment) fluidName.invokeExact(i)));
            return out;
        } catch (Throwable t) { throw new RuntimeException(t); }
    }

    @Override public void close() { arena.close(); }

    // ------------------------------------------------------------------ demo

    private static int failures = 0;

    private static void check(boolean ok, String what) {
        System.out.printf("  %-52s %s%n", what, ok ? "ok" : "FAILED");
        if (!ok) failures++;
    }

    public static void main(String[] args) {
        String path = System.getProperty("rustprop.lib",
                System.getenv().getOrDefault("RUSTPROP_LIB", "librustprop.so"));
        try (var rp = new Rustprop(Path.of(path))) {
            System.out.printf("rustprop %s (CoolProp %s)%n", rp.version(), rp.upstreamVersion());
            System.out.printf("backends: %s%n", rp.backends());
            System.out.printf("fluids compiled in: %d%n", rp.fluids().size());

            if (rp.hasBackend("heos")) {
                double d = rp.propsSi("Dmolar", "T", 300, "P", 101325, "Water");
                check(Math.abs((d - 55317.35277350119) / d) < 1e-8,
                        String.format("PropsSI Dmolar Water = %.10g", d));

                double[] t = {300, 400, 500};
                double[] p = {101325, 101325, 101325};
                double[] many = rp.propsSiMany("Dmolar", "T", t, "P", p, "Water");
                boolean same = true;
                for (int i = 0; i < t.length; i++)
                    if (many[i] != rp.propsSi("Dmolar", "T", t[i], "P", p[i], "Water"))
                        same = false;
                check(same, "batch equals scalar exactly");
            }

            if (rp.hasBackend("if97")) {
                double h = rp.propsSi("H", "T", 300, "P", 101325, "IF97::Water");
                check(Math.abs((h - 112665.04341853978) / h) < 1e-11,
                        String.format("PropsSI H IF97::Water = %.12g", h));
            }

            if (rp.hasBackend("humid-air")) {
                double w = rp.haPropsSi("W", "T", 300, "P", 101325, "R", 0.5);
                check(w > 0 && w < 1, String.format("HAPropsSI W = %.6g", w));
            }

            try {
                rp.propsSi("Dmolar", "T", 300, "P", 101325, "NoSuchFluid");
                check(false, "a bad fluid throws");
            } catch (RustpropException e) {
                check(e.getMessage().contains("NoSuchFluid"),
                        "a bad fluid throws, naming the key");
            }
        }
        System.out.printf("%n%s (%d failures)%n", failures == 0 ? "PASSED" : "FAILED", failures);
        System.exit(failures == 0 ? 0 : 1);
    }
}
