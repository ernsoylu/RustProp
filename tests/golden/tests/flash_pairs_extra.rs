//! Tier-2 input-pair goldens (PLAN.md 4.6 deferrals): `(Hmolar, T)` and
//! `(T, Umolar)` — upstream `DHSU_T_flash` — and `(P, Umolar)` — upstream
//! `HSU_P_flash` — through `props_si`, including the mass-basis variants.
//! The T-based pairs solve density directly (1e-9); the (P, X) pair runs at
//! the established 1e-8 policy (upstream resolves T at ~30 bits).

use rustprop::props_si;
use rustprop_golden_tests::load_jsonl;
use std::path::Path;

#[test]
fn extra_flash_pairs_match_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/flash_pairs_extra.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(records.len(), 308);

    let mut failures = Vec::new();
    for rec in &records {
        let actual = match props_si(
            &rec.out,
            &rec.name1,
            rec.val1,
            &rec.name2,
            rec.val2,
            &format!("HEOS::{}", rec.fluid),
        ) {
            Ok(v) => v,
            Err(e) => {
                failures.push(format!("{} {}: error {e}", rec.fluid, rec.id()));
                continue;
            }
        };
        // (P, X): upstream resolves T at ~30 bits. P outputs of the T-based
        // caloric pairs: the density is Halley-solved to a 1e-8 residual and
        // the liquid's stiff dp/drho amplifies that into P (absolute error
        // stays ~1e-4 Pa) — the established solver-dependent 1e-8 tier.
        let rtol = if rec.name1 == "P" || rec.out == "P" {
            1e-8
        } else {
            1e-9
        };
        // Q = -1 sentinels compare exactly.
        if rec.expected == -1.0 && rec.out == "Q" {
            if actual != -1.0 {
                failures.push(format!(
                    "{} {}: got {actual}, expected -1",
                    rec.fluid,
                    rec.id()
                ));
            }
            continue;
        }
        let rel = ((actual - rec.expected) / rec.expected).abs();
        if rel > rtol || rel.is_nan() {
            failures.push(format!(
                "{} {}: actual {actual:e}, expected {:e}, rel {rel:e}",
                rec.fluid,
                rec.id(),
                rec.expected
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} failures:\n{}",
        failures.len(),
        records.len(),
        failures.join("\n")
    );
}
