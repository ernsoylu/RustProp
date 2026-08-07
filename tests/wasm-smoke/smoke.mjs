// Node smoke test for the wasm-pack output (PLAN.md 14.1).
//
// Loads the built bundle and checks values against the SAME oracle records
// the Rust golden suite uses, so a JS-side regression cannot hide behind a
// green Rust suite. Run after:
//
//   wasm-pack build crates/rustprop-wasm --target nodejs --out-dir pkg-node \
//     --features if97,heos,water,humid-air
//
//   node tests/wasm-smoke/smoke.mjs
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";

const here = dirname(fileURLToPath(import.meta.url));
const root = join(here, "..", "..");
const wasm = await import(join(root, "crates/rustprop-wasm/pkg-node/rustprop_wasm.js"));

let checks = 0;
let failures = [];

function near(got, want, rtol, what) {
  checks += 1;
  const rel = Math.abs((got - want) / want);
  if (!(rel <= rtol)) {
    failures.push(`${what}: got ${got}, want ${want} (rel ${rel.toExponential(3)})`);
  }
}

// 1. The version the port tracks.
checks += 1;
if (wasm.upstream_version() !== "8.0.0") {
  failures.push(`upstream_version: got ${wasm.upstream_version()}`);
}

// 2. Golden records, read straight from the committed fixtures.
function jsonl(name) {
  return readFileSync(join(root, "tests/golden/fixtures", name), "utf8")
    .split("\n")
    .filter((l) => l.trim().length > 0)
    .map((l) => JSON.parse(l));
}

// HEOS water, through the string API.
const smoke = jsonl("water_propssi_smoke.jsonl");
for (const r of smoke) {
  const fluid = `${r.backend}::${r.fluid}`;
  const got = wasm.props_si(r.out, r.name1, r.val1, r.name2, r.val2, fluid);
  near(got, r.expected, r.rtol ?? 1e-9, `props_si ${r.out} ${fluid}`);
}

// IF97 steam, whose values the Rust suite pins separately.
near(
  wasm.props_si("D", "T", 400, "P", 101325, "IF97::Water"),
  wasm.props_si("D", "T", 400, "P", 101325, "IF97::Water"),
  0,
  "IF97 determinism",
);

// 3. The typed fast path agrees with the scalar one, element for element.
const temps = [300, 350, 400, 450, 500];
const press = temps.map(() => 101325);
const batch = wasm.props_si_many("D", "T", new Float64Array(temps), "P", new Float64Array(press), "IF97::Water");
if (batch.length !== temps.length) {
  failures.push(`props_si_many returned ${batch.length} values, expected ${temps.length}`);
} else {
  for (let i = 0; i < temps.length; i += 1) {
    const one = wasm.props_si("D", "T", temps[i], "P", press[i], "IF97::Water");
    checks += 1;
    if (batch[i] !== one) {
      failures.push(`props_si_many[${i}]: ${batch[i]} !== scalar ${one}`);
    }
  }
}

// 4. A failing state yields NaN in its slot without aborting the batch.
const withBad = wasm.props_si_many("D", "T", new Float64Array([400, -5]), "P", new Float64Array([101325, 101325]), "IF97::Water");
checks += 1;
if (!(withBad[0] > 0 && Number.isNaN(withBad[1]))) {
  failures.push(`props_si_many bad-state handling: ${Array.from(withBad)}`);
}

// 5. Errors arrive as thrown exceptions carrying rustprop's message.
checks += 1;
try {
  wasm.props_si("D", "T", 400, "P", 101325, "NoSuchFluid");
  failures.push("props_si on an unknown fluid did not throw");
} catch (e) {
  if (!String(e).length) failures.push("thrown error carried no message");
}

// 6. Humid air, when the bundle was built with it.
if (typeof wasm.ha_props_si === "function") {
  const ha = jsonl("humid_air.jsonl").slice(0, 40);
  for (const r of ha) {
    const got = wasm.ha_props_si(r.out, r.name1, r.val1, r.name2, r.val2, r.name3, r.val3);
    near(got, r.expected, r.rtol ?? 1e-8, `ha_props_si ${r.out}`);
  }
}

if (failures.length > 0) {
  console.error(`FAILED ${failures.length} of ${checks} checks:`);
  for (const f of failures.slice(0, 20)) console.error("  " + f);
  process.exit(1);
}
console.log(`wasm smoke: ${checks} checks passed`);
