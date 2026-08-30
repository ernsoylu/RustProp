//! `ha_props_si` must never panic, for the same reason `props_si` must not:
//! `panic = "abort"` plus a wasm32 target means a panic is a dead module, not
//! a catchable exception.
//!
//! This is the second public entry point, and the original audit sweep only
//! covered the first — a gap in that audit's own coverage, closed here.
//! Psychrometrics is cheap per call (no Helmholtz cascade), so this reaches a
//! large number of real computations in a fraction of a second.

#![cfg(feature = "humid-air")]

use std::panic;

const OUTPUTS: &[&str] = &[
    "T", "P", "W", "Hda", "Sda", "B", "D", "R", "Vda", "Cda", "M", "K", "Tdp", "Twb", "Psi_w",
];
const INPUTS: &[&str] = &[
    "T", "P", "W", "R", "B", "D", "Hda", "Sda", "Tdp", "Twb", "Psi_w", "Vda",
];

/// Boundaries plus the psychrometric range that actually computes: sub-zero
/// through boiling, and pressures either side of one atmosphere.
const VALUES: &[f64] = &[
    f64::NAN,
    f64::INFINITY,
    f64::NEG_INFINITY,
    -1e300,
    -1.0,
    0.0,
    1e-300,
    0.5,
    1.0,
    253.15,
    273.15,
    300.0,
    320.0,
    373.15,
    500.0,
    101_325.0,
    1e5,
    1e8,
    1e300,
    f64::MAX,
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
fn ha_props_si_never_panics() {
    let n: u64 = std::env::var("FUZZ_N")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(50_000);

    let prev = panic::take_hook();
    panic::set_hook(Box::new(|_| {}));

    let mut rng = Rng(0x1234_5678_9ABC_DEF0);
    let mut found: Vec<String> = Vec::new();
    let mut seen = std::collections::HashSet::new();
    let mut solved = 0u64;

    for _ in 0..n {
        let (o, n1, n2, n3) = (
            *rng.pick(OUTPUTS),
            *rng.pick(INPUTS),
            *rng.pick(INPUTS),
            *rng.pick(INPUTS),
        );
        let (v1, v2, v3) = (*rng.pick(VALUES), *rng.pick(VALUES), *rng.pick(VALUES));
        if n1 == n2 || n2 == n3 || n1 == n3 {
            continue;
        }
        let r = panic::catch_unwind(|| rustprop::ha_props_si(o, n1, v1, n2, v2, n3, v3));
        match &r {
            Ok(Ok(_)) => solved += 1,
            Ok(Err(_)) => {}
            Err(_) => {}
        }
        if let Err(e) = r {
            let msg = e
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_else(|| "<non-string panic>".into());
            if seen.insert(msg.clone()) {
                found.push(format!(
                    "ha_props_si({o:?}, {n1:?}, {v1:e}, {n2:?}, {v2:e}, {n3:?}, {v3:e}) \
                     panicked: {msg}"
                ));
            }
        }
    }

    panic::set_hook(prev);
    assert!(
        found.is_empty(),
        "{} distinct panic site(s) reachable from ha_props_si across {n} seeded calls:\n  {}",
        found.len(),
        found.join("\n  ")
    );
    // A sweep that never reaches the solver proves nothing. Random name
    // triples are mostly rejected on sight, so this floor guards the sweep
    // itself: about 12% of draws resolve to a real psychrometric solve.
    let floor = n / 50;
    assert!(
        solved >= floor,
        "only {solved} of {n} draws reached a real computation (floor {floor}); the sweep \
         is bouncing off input validation and testing nothing"
    );
}
