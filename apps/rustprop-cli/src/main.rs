//! Example CLI that makes the rustprop libraries and calculations available
//! over stdout, in `PropsSI` argument order (PLAN.md 2.5).

use rustprop::Param;
use std::process::ExitCode;

const USAGE: &str = "usage: rustprop-cli props <OUT> <NAME1> <val1> <NAME2> <val2> <fluid>
example: rustprop-cli props T P 101325 Q 0 IF97::Water";

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
            let out = parse_param(out)?;
            let n1 = parse_param(n1)?;
            let n2 = parse_param(n2)?;
            let v1 = parse_value(v1)?;
            let v2 = parse_value(v2)?;
            match fluid.as_str() {
                "IF97::Water" => {
                    rustprop::if97_api::props(out, n1, v1, n2, v2).map_err(|e| e.to_string())
                }
                _ => Err(format!(
                    "unknown fluid {fluid:?} (only IF97::Water is ported so far)"
                )),
            }
        }
        _ => Err("expected the props command".into()),
    }
}

fn parse_param(s: &str) -> Result<Param, String> {
    Param::parse(s).ok_or_else(|| format!("unknown parameter {s:?}"))
}

fn parse_value(s: &str) -> Result<f64, String> {
    s.parse().map_err(|_| format!("bad number {s:?}"))
}
