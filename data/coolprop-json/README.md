# coolprop-json

Fluid runtime-JSON documents dumped **verbatim** from the CoolProp 8.0.0 oracle
wheel via `get_fluid_param_string(fluid, "JSON")` — see
`tools/golden-gen/dump_fluid_json.py`. This is the post-injection data the
reference implementation computes with (surface tension, environmental data,
and ancillaries already folded in by upstream's build pipeline; numeric values
exact, since the embedded CBOR carries binary IEEE doubles). Pinned data —
never edited by hand.

`tools/rustprop-datagen` consumes these documents to generate the
feature-gated modules in `crates/rustprop-data`; CI regenerates and diffs.

Derived from CoolProp (https://github.com/CoolProp/CoolProp), MIT License,
Copyright (c) 2012-2018 Ian H. Bell and other CoolProp developers — see the
repository `LICENSE` file.
