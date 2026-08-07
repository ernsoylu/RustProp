# mixtures

`mixture_binary_pairs.json` (888 GERG-2008/Lemmon binary interaction records)
and `mixture_departure_functions.json` (28 departure functions) copied
**verbatim** from the pinned upstream checkout (`dev/mixtures/`, CoolProp tag
v8.0.0) — the same documents upstream bakes into its binary through the
generated `*_JSON.h` headers and loads via `MixtureParameters`. Pinned data —
never edited by hand.

`tools/rustprop-datagen` consumes them to generate
`crates/rustprop-data/src/mixtures.rs` (feature `mixture-data`); CI
regenerates and diffs. The six Lemmon `xi`/`zeta` records are converted to
GERG form at datagen time exactly as upstream's load-time
`LemmonAirHFCReducingFunction::convert_to_GERG` does (using the two fluids'
`EOS.reducing` states).

Derived from CoolProp (https://github.com/CoolProp/CoolProp), MIT License,
Copyright (c) 2012-2018 Ian H. Bell and other CoolProp developers — see the
repository `LICENSE` file.
