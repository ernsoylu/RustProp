//! Regression tests for the `.svds` reader's handling of malformed blobs.
//!
//! The `.svds` container is this port's own invention — upstream caches its
//! surfaces as zlib-compressed msgpack and `artifact.rs` deliberately does not
//! parse that at runtime — so there is no upstream behaviour to be faithful to
//! here, and the reader is free to refuse bad input rather than trust it.
//!
//! Every case below was a live defect before the audit fixes landed. They are
//! cheap and they run in CI, so a future edit to the reader cannot quietly
//! reopen one.

use rustprop_core::params::Param;
use rustprop_svdsbtl::artifact;
use rustprop_svdsbtl::region::{AxisScale, AxisTransform, BoundaryCurve, Region, RegionAtlas};
use rustprop_svdsbtl::surface::SvdSurface;
use rustprop_svdsbtl::svd::{OutputTransform, SlopeSource, SvdDecomposition};

/// A bare 24-byte header with attacker-chosen counts.
fn header(n_props: u32, n_regions: u32) -> Vec<u8> {
    let mut b = Vec::new();
    b.extend_from_slice(b"RPSVDS01");
    b.extend_from_slice(&0i32.to_le_bytes()); // input_pair
    b.extend_from_slice(&n_props.to_le_bytes());
    b.extend_from_slice(&n_regions.to_le_bytes());
    b.extend_from_slice(&0u32.to_le_bytes()); // reserved
    b
}

/// Resident set, in KiB. Linux-only; the assertion below is skipped elsewhere.
#[cfg(target_os = "linux")]
fn vm_size_kb() -> Option<u64> {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmSize"))
                .and_then(|l| l.split_whitespace().nth(1))
                .and_then(|v| v.parse().ok())
        })
}

/// A count in the header must not size an allocation before the bytes behind
/// it have been checked.
///
/// Before the fix this reserved ~16 GB of address space from a 24-byte input
/// (measured via `/proc/self/status`) — survivable under Linux overcommit,
/// fatal on wasm32 where `usize` is 32 bits and linear memory is bounded, and
/// with `panic = "abort"` that takes the whole module down.
#[test]
fn header_counts_do_not_drive_allocation() {
    for n_props in [u32::MAX, u32::MAX / 2, 1 << 24] {
        let b = header(n_props, 0);
        assert_eq!(b.len(), 24);

        #[cfg(target_os = "linux")]
        let before = vm_size_kb();

        let r = artifact::load("x", &b);

        #[cfg(target_os = "linux")]
        if let (Some(before), Some(after)) = (before, vm_size_kb()) {
            let grew_mb = after.saturating_sub(before) / 1024;
            assert!(
                grew_mb < 64,
                "n_props = {n_props} grew the address space by {grew_mb} MB from a 24-byte \
                 input; the count must be bounded against the remaining bytes before any \
                 `with_capacity`"
            );
        }

        // `SvdSurface` is not `Debug`, so match rather than `expect_err`.
        let msg = match r {
            Ok(_) => panic!("a 24-byte blob cannot hold {n_props} properties"),
            Err(e) => e.to_string(),
        };
        assert!(
            msg.contains("entries") || msg.contains("truncated"),
            "expected a length-vs-input refusal, got: {msg}"
        );
    }
}

/// The same rule for the region count.
#[test]
fn region_count_does_not_drive_allocation() {
    let b = header(0, u32::MAX);
    let msg = match artifact::load("x", &b) {
        Ok(_) => panic!("a 24-byte blob cannot hold u32::MAX regions"),
        Err(e) => e.to_string(),
    };
    assert!(
        msg.contains("entries"),
        "expected a length-vs-input refusal, got: {msg}"
    );
}

/// A truncated blob is still a clean `Err`, not a panic.
#[test]
fn truncation_at_every_offset_is_an_error() {
    let full = header(2, 1);
    for cut in 0..full.len() {
        let r = artifact::load("x", &full[..cut]);
        assert!(r.is_err(), "a {cut}-byte blob must not load");
    }
}

/// Bad magic is refused before anything else is read.
#[test]
fn bad_magic_is_refused() {
    let mut b = header(0, 0);
    b[0] = b'X';
    let msg = match artifact::load("x", &b) {
        Ok(_) => panic!("a blob with bad magic must not load"),
        Err(e) => e.to_string(),
    };
    assert!(msg.contains("magic"), "got: {msg}");
}

/// A constant boundary curve: `kind = 0`, then `a_lo`, `a_hi`, `b`.
fn constant_curve(out: &mut Vec<u8>) {
    out.extend_from_slice(&0u32.to_le_bytes());
    for v in [0.0f64, 1.0, 0.5] {
        out.extend_from_slice(&v.to_le_bytes());
    }
}

/// One well-formed LINEAR/LINEAR region over a = [0, 1].
fn one_region(out: &mut Vec<u8>) {
    out.push(0); // primary scale: Linear
    out.push(0); // secondary scale: Linear
    out.extend_from_slice(&0u16.to_le_bytes()); // pad
    for v in [0.0f64, 1.0, 0.0, 1.0, 1.0] {
        // a_lo, a_hi, a_lo_t, a_hi_t, inv_span_t
        out.extend_from_slice(&v.to_le_bytes());
    }
    constant_curve(out); // b_lo
    constant_curve(out); // b_hi
}

/// `nx * rank` must not wrap into a length that then indexes out of bounds.
///
/// On a 32-bit target `usize` is 32 bits and the release profile has
/// `overflow-checks` off, so `nx * rank` with `nx = rank = 65536` wraps to 0.
/// A zero-length coefficient block then MATCHED the wrapped product in
/// `SvdDecomposition::validate`, so the decomposition passed its own
/// consistency check and the panic surfaced later, during evaluation, on a
/// surface that had been declared valid.
///
/// The header alone is enough to reach the check — the reader refuses before
/// it consumes any grid bytes — so this stays a small blob.
#[test]
fn grid_product_overflow_is_refused_not_wrapped() {
    let mut b = header(1, 1);
    b.extend_from_slice(&19i32.to_le_bytes()); // one property: upstream key 19 = T
    one_region(&mut b);
    // The decomposition header: nx = rank = 65536, so nx * rank = 2^32.
    b.extend_from_slice(&65536u32.to_le_bytes()); // nx
    b.extend_from_slice(&2u32.to_le_bytes()); // ny
    b.extend_from_slice(&65536u32.to_le_bytes()); // rank
    b.push(0); // out_transform: Identity
    b.push(1); // slope_source: HermiteFd
    b.extend_from_slice(&0u16.to_le_bytes()); // pad

    let msg = match artifact::load("x", &b) {
        Ok(_) => panic!(
            "a decomposition declaring a 65536 x 65536 coefficient block must not load from \
             a {}-byte blob",
            b.len()
        ),
        Err(e) => e.to_string(),
    };

    // On a 32-bit target the product overflows and must be named as such. On
    // a 64-bit target it is merely enormous, and the reader refuses on the
    // bytes instead — either way, an error rather than a wrapped length.
    if usize::BITS == 32 {
        assert!(
            msg.contains("overflows"),
            "on a {}-bit target the product must be refused as an overflow, got: {msg}",
            usize::BITS
        );
    } else {
        assert!(
            msg.contains("truncated") || msg.contains("overflows"),
            "expected a refusal, got: {msg}"
        );
    }
}

// ---------------------------------------------------------------------------
// Surface-level invariants
// ---------------------------------------------------------------------------

/// A rank-1 decomposition on an `n x n` grid over [0, 1]^2.
fn decomp(n: usize) -> SvdDecomposition {
    let grid: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
    SvdDecomposition {
        nx: n,
        ny: n,
        rank: 1,
        out_transform: OutputTransform::Identity,
        slope_source: SlopeSource::HermiteFd,
        x_grid: grid.clone(),
        y_grid: grid,
        u: vec![1.0; n],
        du_dx: vec![0.0; n],
        v_s: vec![1.0; n],
        dv_s_dy: vec![0.0; n],
        s: vec![1.0],
    }
}

fn atlas_with_one_region() -> RegionAtlas {
    let primary = AxisTransform::make(AxisScale::Linear, 0.0, 1.0).unwrap();
    let region = Region::new(
        primary,
        BoundaryCurve::Constant {
            a_lo: 0.0,
            a_hi: 1.0,
            b: 0.0,
        },
        BoundaryCurve::Constant {
            a_lo: 0.0,
            a_hi: 1.0,
            b: 1.0,
        },
        AxisScale::Linear,
    )
    .unwrap();
    let mut atlas = RegionAtlas::default();
    atlas.add(region);
    atlas
}

/// `eval_with_region*` take a caller-supplied index and return `Result`, so an
/// out-of-range index is an error, not an index panic.
#[test]
fn out_of_range_region_index_is_an_error() {
    let s = match SvdSurface::new(
        "poc".into(),
        0,
        vec![Param::T],
        atlas_with_one_region(),
        vec![vec![decomp(8)]],
    ) {
        Ok(s) => s,
        Err(e) => panic!("a uniform single-property region must seal: {e}"),
    };
    assert_eq!(s.region_count(), 1);

    let e = s
        .eval_with_region(Param::T, 99, 0.5, 0.5)
        .expect_err("region 99 does not exist");
    assert!(e.to_string().contains("out of range"), "got: {e}");

    let mut out = [0.0; 1];
    let e = s
        .eval_with_region_multi(99, 0.5, 0.5, &[Param::T], &mut out)
        .expect_err("region 99 does not exist");
    assert!(e.to_string().contains("out of range"), "got: {e}");
}

/// `eval_with_region_multi` builds one context from `row[0]` and reuses it, so
/// a region whose properties sit on different grids must be refused at seal
/// time — otherwise the two public entry points disagree and the multi path
/// indexes one decomposition's coefficients with another's cell index.
#[test]
fn mismatched_grids_within_a_region_are_refused() {
    // `SvdSurface` is not `Debug`, so match rather than `expect_err`.
    let e = match SvdSurface::new(
        "poc".into(),
        0,
        vec![Param::T, Param::P],
        atlas_with_one_region(),
        vec![vec![decomp(64), decomp(4)]],
    ) {
        Ok(_) => panic!("a region whose properties are on 64x64 and 4x4 grids is malformed"),
        Err(e) => e,
    };
    assert!(
        e.to_string().contains("share"),
        "expected a shared-grid refusal, got: {e}"
    );
}

/// Matching grids still seal, and both entry points agree on the result.
#[test]
fn matching_grids_seal_and_both_entry_points_agree() {
    let s = match SvdSurface::new(
        "ok".into(),
        0,
        vec![Param::T, Param::P],
        atlas_with_one_region(),
        vec![vec![decomp(16), decomp(16)]],
    ) {
        Ok(s) => s,
        Err(e) => panic!("uniform grids are well-formed: {e}"),
    };

    let single = s.eval_with_region(Param::P, 0, 0.5, 0.5).unwrap();
    let mut out = [0.0; 1];
    s.eval_with_region_multi(0, 0.5, 0.5, &[Param::P], &mut out)
        .unwrap();
    assert_eq!(
        single, out[0],
        "eval_with_region and eval_with_region_multi must agree"
    );
}

/// The guards `AxisTransform::make` enforces are reachable on the load path
/// too — the reader builds the struct literally to keep upstream's bits, so
/// `Region::new` is where they get applied.
#[test]
fn degenerate_axis_bounds_are_refused() {
    let bad = AxisTransform {
        scale: AxisScale::Linear,
        a_lo: 1.0,
        a_hi: 1.0, // not > a_lo
        a_lo_t: 1.0,
        a_hi_t: 1.0,
        inv_span_t: f64::INFINITY,
    };
    let e = Region::new(
        bad,
        BoundaryCurve::Constant {
            a_lo: 0.0,
            a_hi: 1.0,
            b: 0.0,
        },
        BoundaryCurve::Constant {
            a_lo: 0.0,
            a_hi: 1.0,
            b: 1.0,
        },
        AxisScale::Linear,
    )
    .expect_err("a_hi must exceed a_lo");
    assert!(e.to_string().contains("a_hi"), "got: {e}");
}
