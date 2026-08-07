//! Converts one upstream CoolProp `.svd.bin.z` surface into the flat
//! little-endian `.svds` blob `rustprop-svdsbtl` embeds.
//!
//! This is the ingestion half of the Phase 13.1 design decision: upstream's
//! coefficients come from Eigen's BDCSVD, which no independent
//! implementation reproduces bitwise, so they are copied rather than
//! recomputed. Every f64 is written as the exact bit pattern upstream
//! stored, and the conversion is checked by a golden test that evaluates
//! the result against the wheel's own `SVDSBTL&HEOS` backend.
//!
//! ```text
//! cargo run -p rustprop-svdgen -- <in.svd.bin.z> <out.svds>
//! ```
//!
//! Upstream's on-wire form (msgpack inside zlib), from
//! `include/CoolProp/sbtl/SVDSurfaceSerializer.h`:
//!
//! ```text
//! ["SVDS", revision, fluid_name, [surface, ...]]
//! surface = [input_pair, [[prop_key], ...], [region, ...], [decomp, ...]]
//! region  = [primary_scale, a_lo, a_hi, a_lo_t, a_hi_t, inv_span_t,
//!            curve_lo, curve_hi, secondary_scale]
//! curve   = [0, a_lo, a_hi, b]                        (constant)
//!         | [1, a[], b[], M[], b_min, b_max]          (cubic spline)
//! decomp  = [region_idx, prop_key, NX, NY, rank, out_transform,
//!            slope_source, x_grid[], y_grid[], U[], dU_dx[], V_S[],
//!            dV_S_dy[], S[]]
//! ```

use flate2::read::ZlibDecoder;
use rmpv::Value;
use rustprop_svdsbtl::artifact::MAGIC;
use std::io::Read;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: rustprop-svdgen <in.svd.bin.z> <out.svds>");
        std::process::exit(2);
    }
    match convert(&args[1], &args[2]) {
        Ok(n) => println!("wrote {} ({n} bytes)", args[2]),
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(1);
        }
    }
}

fn convert(src: &str, dst: &str) -> Result<usize, String> {
    let packed = std::fs::read(src).map_err(|e| format!("reading {src}: {e}"))?;
    let mut raw = Vec::new();
    ZlibDecoder::new(&packed[..])
        .read_to_end(&mut raw)
        .map_err(|e| format!("inflating {src}: {e}"))?;
    let root = rmpv::decode::read_value(&mut &raw[..]).map_err(|e| format!("msgpack: {e}"))?;

    let top = arr(&root, "root")?;
    if top.len() < 4 || top[0].as_str() != Some("SVDS") {
        return Err("not an SVDS artifact (bad magic)".into());
    }
    let surfaces = arr(&top[3], "surfaces")?;
    if surfaces.len() != 1 {
        return Err(format!(
            "expected exactly one surface per file, found {}",
            surfaces.len()
        ));
    }
    let s = arr(&surfaces[0], "surface")?;
    let input_pair = int(&s[0], "input_pair")? as i32;
    let props = arr(&s[1], "properties")?;
    let regions = arr(&s[2], "regions")?;
    let decomps = arr(&s[3], "decomps")?;

    let n_props = props.len();
    let n_regions = regions.len();
    if decomps.len() != n_props * n_regions {
        return Err(format!(
            "{} decompositions for {n_regions} regions x {n_props} properties",
            decomps.len()
        ));
    }

    let mut out = Vec::with_capacity(raw.len());
    out.extend_from_slice(MAGIC);
    out.extend_from_slice(&input_pair.to_le_bytes());
    out.extend_from_slice(&(n_props as u32).to_le_bytes());
    out.extend_from_slice(&(n_regions as u32).to_le_bytes());
    out.extend_from_slice(&0u32.to_le_bytes());

    // Property keys, in surface order. Each entry is a length-1 array.
    let mut prop_keys = Vec::with_capacity(n_props);
    for p in props {
        let key = int(&arr(p, "property")?[0], "prop key")? as i32;
        prop_keys.push(key);
        out.extend_from_slice(&key.to_le_bytes());
    }

    for r in regions {
        let r = arr(r, "region")?;
        let primary_scale = int(&r[0], "primary_scale")? as u8;
        let secondary_scale = int(&r[8], "secondary_scale")? as u8;
        out.push(primary_scale);
        out.push(secondary_scale);
        out.extend_from_slice(&0u16.to_le_bytes());
        for v in &r[1..=5] {
            out.extend_from_slice(&flt(v, "region bound")?.to_le_bytes());
        }
        write_curve(&mut out, &r[6])?;
        write_curve(&mut out, &r[7])?;
    }

    // Decompositions, region-major and in the surface's property order —
    // the artifact carries (region_idx, prop_key) on every entry, so index
    // rather than trust the file's ordering.
    for ri in 0..n_regions {
        for pk in &prop_keys {
            let d = decomps
                .iter()
                .map(|d| arr(d, "decomp"))
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .find(|d| {
                    int(&d[0], "region_idx").is_ok_and(|v| v as usize == ri)
                        && int(&d[1], "prop_key").is_ok_and(|v| v as i32 == *pk)
                })
                .ok_or_else(|| format!("no decomposition for region {ri}, property {pk}"))?;
            write_decomp(&mut out, d)?;
        }
    }

    std::fs::write(dst, &out).map_err(|e| format!("writing {dst}: {e}"))?;
    Ok(out.len())
}

fn write_curve(out: &mut Vec<u8>, v: &Value) -> Result<(), String> {
    let c = arr(v, "curve")?;
    let kind = int(&c[0], "curve kind")? as u32;
    out.extend_from_slice(&kind.to_le_bytes());
    match kind {
        0 => {
            for v in &c[1..=3] {
                out.extend_from_slice(&flt(v, "constant curve")?.to_le_bytes());
            }
        }
        1 => {
            let a = f64s(&c[1], "spline a")?;
            let b = f64s(&c[2], "spline b")?;
            let m = f64s(&c[3], "spline M")?;
            if b.len() != a.len() || m.len() != a.len() {
                return Err("spline a/b/M size mismatch".into());
            }
            out.extend_from_slice(&(a.len() as u32).to_le_bytes());
            out.extend_from_slice(&0u32.to_le_bytes());
            out.extend_from_slice(&flt(&c[4], "b_min")?.to_le_bytes());
            out.extend_from_slice(&flt(&c[5], "b_max")?.to_le_bytes());
            for v in a.iter().chain(b.iter()).chain(m.iter()) {
                out.extend_from_slice(&v.to_le_bytes());
            }
        }
        other => return Err(format!("boundary-curve kind {other} is not supported")),
    }
    Ok(())
}

fn write_decomp(out: &mut Vec<u8>, d: &[Value]) -> Result<(), String> {
    let nx = int(&d[2], "NX")? as u32;
    let ny = int(&d[3], "NY")? as u32;
    let rank = int(&d[4], "rank")? as u32;
    out.extend_from_slice(&nx.to_le_bytes());
    out.extend_from_slice(&ny.to_le_bytes());
    out.extend_from_slice(&rank.to_le_bytes());
    out.push(int(&d[5], "out_transform")? as u8);
    out.push(int(&d[6], "slope_source")? as u8);
    out.extend_from_slice(&0u16.to_le_bytes());
    let want = [
        (7usize, nx as usize),
        (8, ny as usize),
        (9, (nx * rank) as usize),
        (10, (nx * rank) as usize),
        (11, (ny * rank) as usize),
        (12, (ny * rank) as usize),
        (13, rank as usize),
    ];
    for (idx, n) in want {
        let v = f64s(&d[idx], "decomp array")?;
        if v.len() != n {
            return Err(format!(
                "decomp field {idx}: {} values, expected {n}",
                v.len()
            ));
        }
        for x in &v {
            out.extend_from_slice(&x.to_le_bytes());
        }
    }
    Ok(())
}

fn arr<'a>(v: &'a Value, what: &str) -> Result<&'a Vec<Value>, String> {
    v.as_array()
        .ok_or_else(|| format!("{what}: expected a msgpack array"))
}

fn int(v: &Value, what: &str) -> Result<i64, String> {
    v.as_i64()
        .ok_or_else(|| format!("{what}: expected an integer"))
}

fn flt(v: &Value, what: &str) -> Result<f64, String> {
    v.as_f64()
        .ok_or_else(|| format!("{what}: expected a float"))
}

fn f64s(v: &Value, what: &str) -> Result<Vec<f64>, String> {
    arr(v, what)?
        .iter()
        .map(|x| flt(x, what))
        .collect::<Result<Vec<_>, _>>()
}
