//! C `printf` conversions the ported error strings are built with.
//!
//! Upstream composes its messages through `format()` (`CoolPropTools.h`), a
//! `sprintf` wrapper, so a message the port must reproduce verbatim needs C's
//! conversion rules rather than Rust's `Display`.

/// C's `%g` at the default precision of 6 — the conversion upstream's error
/// strings are built with. Rust's `Display` agrees for many values but not
/// all (`%g` switches to exponential outside `[1e-4, 1e6)`, pads the exponent
/// to two digits, and drops trailing zeros), and these messages are asserted
/// against the wheel.
pub fn fmt_g(v: f64) -> String {
    if v.is_nan() {
        return "nan".into();
    }
    if v.is_infinite() {
        return if v < 0.0 { "-inf".into() } else { "inf".into() };
    }
    const P: i32 = 6;
    // The decimal exponent %e would use.
    let exp = if v == 0.0 {
        0
    } else {
        let e = format!("{:.*e}", (P - 1) as usize, v);
        e.split('e')
            .nth(1)
            .and_then(|s| s.parse::<i32>().ok())
            .unwrap_or(0)
    };
    let trim = |s: String| -> String {
        if !s.contains('.') {
            return s;
        }
        let s = s.trim_end_matches('0');
        s.trim_end_matches('.').to_string()
    };
    if !(-4..P).contains(&exp) {
        let mantissa = format!("{:.*e}", (P - 1) as usize, v);
        let (m, e) = mantissa.split_once('e').expect("scientific form");
        let e: i32 = e.parse().expect("exponent");
        format!(
            "{}e{}{:02}",
            trim(m.to_string()),
            if e < 0 { '-' } else { '+' },
            e.abs()
        )
    } else {
        trim(format!("{:.*}", (P - 1 - exp).max(0) as usize, v))
    }
}

#[cfg(test)]
mod tests {
    use super::fmt_g;

    #[test]
    fn g_format_matches_c() {
        // Values checked against printf("%g", x).
        assert_eq!(fmt_g(0.001), "0.001");
        assert_eq!(fmt_g(400.0), "400");
        assert_eq!(fmt_g(101325.0), "101325");
        assert_eq!(fmt_g(27965757.789598826), "2.79658e+07");
        assert_eq!(fmt_g(1e-5), "1e-05");
        assert_eq!(fmt_g(0.0001), "0.0001");
        assert_eq!(fmt_g(0.0), "0");
        assert_eq!(fmt_g(-3.5), "-3.5");
        assert_eq!(fmt_g(1234567.0), "1.23457e+06");
        // The alpha0 validity message quotes tau and delta this way.
        assert_eq!(fmt_g(6.470959999999999e-28), "6.47096e-28");
        assert_eq!(fmt_g(6.818239999999999e-31), "6.81824e-31");
    }
}
