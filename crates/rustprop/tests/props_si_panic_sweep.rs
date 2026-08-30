//! `props_si` must never panic, whatever it is handed.
//!
//! The release profile sets `panic = "abort"`, and the primary target is
//! wasm32, so a panic anywhere under this entry point does not unwind into a
//! catchable JS exception — it takes the whole module down for the host. That
//! makes "no panics" a hard property of the public surface rather than a
//! nicety, and CLAUDE.md already records that randomized coverage finds what
//! hand-chosen goldens do not.
//!
//! This is a seeded sweep, not a fuzzer: the same inputs run every time, so a
//! failure is reproducible from the printed case alone. It is deliberately
//! cheap enough to gate every push.

#![cfg(feature = "heos")]

use std::panic;

const OUTPUTS: &[&str] = &[
    "T",
    "P",
    "D",
    "Dmolar",
    "H",
    "Hmass",
    "S",
    "Smass",
    "U",
    "Q",
    "C",
    "Cvmass",
    "A",
    "V",
    "L",
    "surface_tension",
    "Prandtl",
    "Z",
    "Phase",
    "isothermal_compressibility",
];
const INPUTS: &[&str] = &[
    "T", "P", "D", "Dmolar", "H", "Hmass", "S", "Smass", "U", "Q",
];

/// Every boundary a property calculation has: non-finite, signed zero, the
/// subnormal floor, the overflow ceiling, and a few real physical values.
const VALUES: &[f64] = &[
    f64::NAN,
    f64::INFINITY,
    f64::NEG_INFINITY,
    -1.0e300,
    -1.0,
    -0.0,
    0.0,
    f64::MIN_POSITIVE,
    1e-300,
    1e-30,
    1.0,
    273.15,
    300.0,
    373.9,
    647.096,
    1000.0,
    1e5,
    101_325.0,
    1e8,
    22.064e6,
    1e30,
    1e300,
    f64::MAX,
    0.5,
    2.0,
    -0.5,
];

/// One fluid per engine, plus the pseudo-pure blends that ship with
/// `superancillary: None` — those exercise the `.expect(...)` sites that
/// assume a superancillary is present.
const FLUIDS: &[&str] = &[
    "Water",
    "R134a",
    "Nitrogen",
    "CO2",
    "Ammonia",
    "Air",
    "R404A",
    "R407C",
    "R410A",
    "R507A",
    "SES36",
    #[cfg(feature = "if97")]
    "IF97::Water",
    #[cfg(feature = "cubics")]
    "SRK::Propane",
    #[cfg(feature = "cubics")]
    "PR::Propane",
    #[cfg(feature = "pcsaft")]
    "PCSAFT::Propane",
    #[cfg(feature = "incompressible")]
    "INCOMP::MEG-20%",
    #[cfg(feature = "heos-mixtures")]
    "HEOS::Water&Ethanol",
];

/// xorshift64*, so the sweep is deterministic without a dev-dependency.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    fn pick<'a, T>(&mut self, s: &'a [T]) -> &'a T {
        &s[(self.next() % s.len() as u64) as usize]
    }
}

#[test]
fn props_si_never_panics() {
    // The default keeps the debug gate quick; CI's `safety` job re-runs this
    // in release with a much larger FUZZ_N, which is where the real coverage
    // comes from. Raise it locally to hunt.
    let n: u64 = std::env::var("FUZZ_N")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(2_000);

    let prev = panic::take_hook();
    panic::set_hook(Box::new(|_| {})); // the assertion below does the reporting

    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let mut found: Vec<String> = Vec::new();
    let mut seen = std::collections::HashSet::new();

    for _ in 0..n {
        let (out, n1, n2) = (*rng.pick(OUTPUTS), *rng.pick(INPUTS), *rng.pick(INPUTS));
        let (v1, v2, fluid) = (*rng.pick(VALUES), *rng.pick(VALUES), *rng.pick(FLUIDS));
        if n1 == n2 {
            continue;
        }
        let r = panic::catch_unwind(|| rustprop::props_si(out, n1, v1, n2, v2, fluid));
        if let Err(e) = r {
            let msg = e
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_else(|| "<non-string panic>".into());
            if seen.insert(msg.clone()) {
                found.push(format!(
                    "props_si({out:?}, {n1:?}, {v1:e}, {n2:?}, {v2:e}, {fluid:?}) panicked: {msg}"
                ));
            }
        }
    }

    panic::set_hook(prev);
    assert!(
        found.is_empty(),
        "{} distinct panic site(s) reachable from props_si across {n} seeded calls:\n  {}",
        found.len(),
        found.join("\n  ")
    );
}
