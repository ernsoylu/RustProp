# Archived oracle

One file: a byte-identical copy of the PyPI artifact that produced every golden
fixture in this repository.

```
coolprop-8.0.0-cp312-abi3-manylinux2014_x86_64.manylinux_2_17_x86_64.whl
sha256 8ca1aefd1873b14f3e7e9122a08a32fe8c62ac03a1fc4496512a31a3bba8f611
```

Verify it, install it, and understand why it is here: `../ORACLE.md`.

```bash
sha256sum -c <<<"8ca1aefd1873b14f3e7e9122a08a32fe8c62ac03a1fc4496512a31a3bba8f611  coolprop-8.0.0-cp312-abi3-manylinux2014_x86_64.manylinux_2_17_x86_64.whl"
```

CoolProp is MIT-licensed (Copyright 2012-2018 Ian H. Bell and other CoolProp
developers); the wheel carries its own licence at
`coolprop-8.0.0.dist-info/licenses/LICENSE`. This directory is dev-only —
`tools/golden-gen` is excluded from the Cargo workspace, so nothing here
reaches a published crate, a release tarball, or a wasm bundle.

Do not add other builds here. A second wheel would make "the oracle"
ambiguous, which is the exact problem this directory exists to end.
