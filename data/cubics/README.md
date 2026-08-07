# cubics

`all_cubic_fluids.json` copied **verbatim** from the pinned upstream checkout
(`dev/cubics/all_cubic_fluids.json`, CoolProp tag v8.0.0) — the same document
upstream bakes into its binary via `all_cubics_JSON.h` and loads through
`CubicsLibrary`. Pinned data — never edited by hand.

`tools/rustprop-datagen` consumes it to generate
`crates/rustprop-data/src/cubics.rs` (feature `cubic-fluids`); CI regenerates
and diffs.

Derived from CoolProp (https://github.com/CoolProp/CoolProp), MIT License,
Copyright (c) 2012-2018 Ian H. Bell and other CoolProp developers — see the
repository `LICENSE` file.
