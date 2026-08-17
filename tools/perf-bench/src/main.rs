//! Timing harness for the hot HEOS paths (PLAN.md perf work). Standalone
//! package — NOT a workspace member — so it never affects the shipped crates.
//!
//! Usage:
//!   perf-bench                     # fast benches only
//!   perf-bench --grids             # + LogPT 200x200 and LogPH 60x60 builds
//!   perf-bench --logph200          # + LogPH 200x200 build (~100 s pre-memo)
//!   perf-bench --dump-grids DIR    # write bitwise grid dumps (LogPT 200x200,
//!                                  #   LogPH 60x60) as hex lines into DIR
//!
//! Machine note: this box runs the powersave governor — the first ~2 s of
//! work is 3-5x noisy, so a global warmup spin precedes every measurement.

use rustprop::props_si;
use rustprop_heos::flash_pt::PtFlash;
use rustprop_tabular::tables::{GridKind, GriddedTable};
use std::hint::black_box;
use std::io::Write;
use std::time::{Duration, Instant};

fn fluid(name: &str) -> &'static rustprop_core::fluid::FluidData {
    rustprop_data::fluids::all()
        .into_iter()
        .find(|(n, _)| *n == name)
        .expect("fluid in registry")
        .1
}

/// Run `f` for ~0.3 s of warmup, then time it over ~1 s; report mean ns/op.
fn bench<F: FnMut() -> f64>(name: &str, mut f: F) {
    let warm_until = Instant::now() + Duration::from_millis(300);
    while Instant::now() < warm_until {
        black_box(f());
    }
    // Calibrate an iteration count for ~1 s total.
    let t0 = Instant::now();
    let probe = 32;
    for _ in 0..probe {
        black_box(f());
    }
    let per = t0.elapsed().as_secs_f64() / probe as f64;
    let n = ((1.0 / per) as u64).clamp(16, 5_000_000);
    let t0 = Instant::now();
    for _ in 0..n {
        black_box(f());
    }
    let ns = t0.elapsed().as_secs_f64() * 1e9 / n as f64;
    println!("{name:<28} {ns:>14.1} ns/op   ({n} iters)");
}

fn bench_once<F: FnOnce() -> u64>(name: &str, f: F) {
    let t0 = Instant::now();
    let sink = f();
    let s = t0.elapsed().as_secs_f64();
    println!("{name:<28} {s:>14.3} s      (checksum {sink:016x})");
}

/// Fold every f64 of a grid into a checksum so builds cannot be optimized out.
fn grid_checksum(t: &GriddedTable) -> u64 {
    let mut acc = 0u64;
    let mut mix = |v: f64| acc = acc.rotate_left(7) ^ v.to_bits();
    for v in t.xvec.iter().chain(t.yvec.iter()) {
        mix(*v);
    }
    for g in [&t.t, &t.p, &t.rhomolar, &t.hmolar, &t.smolar, &t.umolar] {
        for m in [&g.val, &g.dx, &g.dy, &g.dxx, &g.dxy, &g.dyy] {
            for row in m {
                for v in row {
                    mix(*v);
                }
            }
        }
    }
    acc
}

/// Bitwise dump: every f64 as 16 hex digits, one per line, with section
/// headers — the byte-identity witness for the memoization change.
fn dump_grid(t: &GriddedTable, path: &std::path::Path) {
    let mut out = std::io::BufWriter::new(std::fs::File::create(path).expect("create dump"));
    let mut wf = |name: &str, m: &Vec<Vec<f64>>| {
        writeln!(out, "# {name}").unwrap();
        for row in m {
            for v in row {
                writeln!(out, "{:016x}", v.to_bits()).unwrap();
            }
        }
    };
    wf("xvec", &vec![t.xvec.clone()]);
    wf("yvec", &vec![t.yvec.clone()]);
    for (gname, g) in [
        ("t", &t.t),
        ("p", &t.p),
        ("rhomolar", &t.rhomolar),
        ("hmolar", &t.hmolar),
        ("smolar", &t.smolar),
        ("umolar", &t.umolar),
    ] {
        for (mname, m) in [
            ("val", &g.val),
            ("dx", &g.dx),
            ("dy", &g.dy),
            ("dxx", &g.dxx),
            ("dxy", &g.dxy),
            ("dyy", &g.dyy),
        ] {
            wf(&format!("{gname}.{mname}"), m);
        }
    }
    wf("visc", &t.visc);
    wf("cond", &t.cond);
    writeln!(out, "# nearest").unwrap();
    for i in 0..t.nx {
        for j in 0..t.ny {
            writeln!(out, "{} {}", t.nearest_i[i][j], t.nearest_j[i][j]).unwrap();
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let grids = args.iter().any(|a| a == "--grids");
    let logph200 = args.iter().any(|a| a == "--logph200");
    let dump_dir = args
        .iter()
        .position(|a| a == "--dump-grids")
        .map(|i| std::path::PathBuf::from(&args[i + 1]));

    // Global warmup: spin ~2 s so the powersave governor ramps up.
    let until = Instant::now() + Duration::from_secs(2);
    while Instant::now() < until {
        black_box(props_si("Dmolar", "T", 300.0, "P", 101325.0, "Water").unwrap());
    }

    // Self-consistent flash targets (derived, so they are always valid).
    let h_liq = props_si("Hmass", "T", 320.0, "P", 5e6, "Water").unwrap();
    let h_gas = props_si("Hmass", "T", 500.0, "P", 101325.0, "Water").unwrap();
    let s_gas = props_si("Smass", "T", 500.0, "P", 101325.0, "Water").unwrap();

    bench("warm PT Water (Dmolar)", || {
        props_si("Dmolar", "T", 300.0, "P", 101325.0, "Water").unwrap()
    });
    bench("HP liquid Water (T)", || {
        props_si("T", "P", 5e6, "Hmass", h_liq, "Water").unwrap()
    });
    bench("HP gas Water (T)", || {
        props_si("T", "P", 101325.0, "Hmass", h_gas, "Water").unwrap()
    });
    bench("HS gas Water (T)", || {
        props_si("T", "Hmass", h_gas, "Smass", s_gas, "Water").unwrap()
    });
    bench("QT Water (P)", || {
        props_si("P", "T", 300.0, "Q", 0.0, "Water").unwrap()
    });
    bench("mixture PT CH4/C2H6 (Dmolar)", || {
        props_si(
            "Dmolar",
            "T",
            250.0,
            "P",
            5e6,
            "HEOS::Methane[0.6]&Ethane[0.4]",
        )
        .unwrap()
    });

    // Raw kernels.
    let water = fluid("Water");
    let eos = rustprop_heos::alpha::HelmholtzEos::new(water);
    let (tau, delta) = (eos.t_reducing / 300.0, 55000.0 / eos.rhomolar_reducing);
    bench("alphar_all Water", || {
        eos.alphar_all(black_box(tau), black_box(delta)).d10
    });
    bench("alpha0_all Water", || {
        eos.alpha0_all(black_box(tau), black_box(delta)).d01
    });

    if grids || dump_dir.is_some() {
        let flash = PtFlash::new(water);
        let mut logpt = None;
        let mut logph60 = None;
        bench_once("LogPT 200x200 build", || {
            let t = GriddedTable::build(&flash, GridKind::LogPT, 200, 200, None).expect("LogPT");
            let c = grid_checksum(&t);
            logpt = Some(t);
            c
        });
        bench_once("LogPH 60x60 build", || {
            let t = GriddedTable::build(&flash, GridKind::LogPH, 60, 60, None).expect("LogPH");
            let c = grid_checksum(&t);
            logph60 = Some(t);
            c
        });
        if logph200 {
            bench_once("LogPH 200x200 build", || {
                let t =
                    GriddedTable::build(&flash, GridKind::LogPH, 200, 200, None).expect("LogPH");
                grid_checksum(&t)
            });
        }
        if let Some(dir) = dump_dir {
            std::fs::create_dir_all(&dir).expect("dump dir");
            dump_grid(logpt.as_ref().unwrap(), &dir.join("logpt-200x200.hex"));
            dump_grid(logph60.as_ref().unwrap(), &dir.join("logph-60x60.hex"));
            println!("dumps written to {}", dir.display());
        }
    }
}
