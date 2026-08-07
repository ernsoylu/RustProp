//! Example CLI that makes the rustprop libraries and calculations available
//! over stdout, in `PropsSI` argument order (PLAN.md 2.5/5.3).

use std::process::ExitCode;

const USAGE: &str = "usage: rustprop-cli props <OUT> <NAME1> <val1> <NAME2> <val2> <fluid>
examples: rustprop-cli props T P 101325 Q 0 Water
          rustprop-cli props H T 300 P 101325 IF97::Water
          rustprop-cli ha H T 298.15 P 101325 R 0.5";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match run(&args) {
        Ok(value) => {
            println!("{value}");
            ExitCode::SUCCESS
        }
        Err(msg) => {
            eprintln!("error: {msg}");
            eprintln!("{USAGE}");
            ExitCode::FAILURE
        }
    }
}

fn run(args: &[String]) -> Result<f64, String> {
    match args {
        [cmd, out, n1, v1, n2, v2, fluid] if cmd == "props" => {
            let v1 = parse_value(v1)?;
            let v2 = parse_value(v2)?;
            rustprop::props_si(out, n1, v1, n2, v2, fluid).map_err(|e| e.to_string())
        }
        [cmd, out, n1, v1, n2, v2, n3, v3] if cmd == "ha" => {
            let v1 = parse_value(v1)?;
            let v2 = parse_value(v2)?;
            let v3 = parse_value(v3)?;
            rustprop::ha_props_si(out, n1, v1, n2, v2, n3, v3).map_err(|e| e.to_string())
        }
        _ => Err("expected the props or ha command".into()),
    }
}

fn parse_value(s: &str) -> Result<f64, String> {
    s.parse().map_err(|_| format!("bad number {s:?}"))
}
