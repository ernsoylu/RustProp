# incompressible

The 126 incompressible-fluid documents copied **verbatim** from the pinned
upstream checkout (`dev/incompressible_liquids/json/*.json`, CoolProp tag
v8.0.0) — the same documents upstream bakes into its binary through
`all_incompressibles_JSON.h` and loads via `JSONIncompressibleLibrary`.
Pinned data — never edited by hand.

`tools/rustprop-datagen` consumes them to generate
`crates/rustprop-data/src/incompressible.rs` (feature
`incompressible-fluids`); CI regenerates and diffs.

Derived from CoolProp (https://github.com/CoolProp/CoolProp), MIT License,
Copyright (c) 2012-2018 Ian H. Bell and other CoolProp developers — see the
repository `LICENSE` file.
