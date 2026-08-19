# The oracle

Every number this project claims to be correct is correct *relative to one
compiled binary*: the CoolProp 8.0.0 Python extension module in the wheel
identified below. All 41,629 records in `tests/golden/fixtures/` came from it,
and the port deliberately follows **it** rather than the `v8.0.0` tag source
wherever the two disagree — see NEXT-STEPS.md, "Wheel-vs-tag discoveries", for
the two places they do (IF97 `set_phase`, HEOS DmolarT phase labels).

That is why `CoolProp==8.0.0` was never a sufficient pin. A version string
names a release; it does not name a build.

## Identity

| | |
|---|---|
| Distribution | `coolprop-8.0.0-cp312-abi3-manylinux2014_x86_64.manylinux_2_17_x86_64.whl` |
| Wheel sha256 | `8ca1aefd1873b14f3e7e9122a08a32fe8c62ac03a1fc4496512a31a3bba8f611` |
| Wheel size | 10,814,565 bytes |
| Wheel tags | `cp312-abi3-manylinux_2_17_x86_64`, `cp312-abi3-manylinux2014_x86_64` |
| Extension module | `CoolProp/CoolProp.abi3.so`, 9,050,856 bytes |
| **Extension sha256** | `05d85591871524e83bb23170e22ee149e1167fb3ef5deaf81b9d248283089be5` |
| `CoolProp.__gitrevision__` | `ae81610e7d23efc57f9d051c8e70a4d66e87537f` (the `v8.0.0` tag commit) |
| Python | 3.12.3 (CPython, Linux) |
| Host | `Linux-x86_64` |
| Verified | 2026-08-19 — the wheel above was downloaded fresh from PyPI and its `CoolProp.abi3.so` hashes to the digest in bold, matching the module the fixtures were generated with |

The extension digest is also written into
`tests/golden/fixtures/manifest.json` (`oracle_sha256`) on every generator run,
so a regeneration always records which binary answered.

An `abi3` wheel runs on CPython 3.12 **and later**. The bytes of the extension
do not change with the interpreter, so 3.13 or 3.14 will load the same oracle;
3.12.3 is simply what was used.

## Archived in-repo

`vendor/coolprop-8.0.0-cp312-abi3-manylinux2014_x86_64.manylinux_2_17_x86_64.whl`
is a byte-identical copy of the PyPI artifact.

**Why it is committed.** The project's entire correctness argument reduces to
this one file, and the fixtures cannot be reproduced without it — not even from
CoolProp's own source, because the shipped wheel demonstrably differs from the
`v8.0.0` tag in behaviour the port follows. PyPI does not normally delete
released files, but "normally" is not a guarantee worth resting a whole test
suite on. 10.3 MB against a 64 MB `.git` is a 16% cost, paid once, for the
difference between a reproducible project and one that silently becomes
unverifiable the day an upstream URL stops resolving.

**Licence.** CoolProp is MIT (Copyright 2012-2018 Ian H. Bell and other
CoolProp developers), and the wheel carries its own licence text at
`coolprop-8.0.0.dist-info/licenses/LICENSE`. Redistribution is explicitly
permitted; nothing in this project's own MIT licence conflicts.

**Scope.** The wheel is a dev-only artifact. `tools/golden-gen` is `exclude`d
from the Cargo workspace, so it is in no published crate, in no release
tarball, and in no wasm bundle — it costs `git clone` and nothing else.

## Cold start

From a machine with nothing on it, to regenerated fixtures.

```bash
# 1. A CPython 3.12+ interpreter. The wheel is abi3; 3.12 is what was used.
python3.12 --version

# 2. The venv the generator expects, at exactly this path.
python3.12 -m venv tools/golden-gen/.venv

# 3a. Preferred — install from PyPI with hashes enforced.
tools/golden-gen/.venv/bin/pip install --require-hashes \
    -r tools/golden-gen/requirements.txt

# 3b. Fallback — PyPI gone, offline, or the hash check failed: install the
#     archived copy. Identical bytes, so identical fixtures.
tools/golden-gen/.venv/bin/pip install --no-index \
    tools/golden-gen/vendor/coolprop-8.0.0-cp312-abi3-manylinux2014_x86_64.manylinux_2_17_x86_64.whl

# 4. Prove it is the right binary BEFORE generating anything.
tools/golden-gen/verify-oracle.sh

# 5. Regenerate. Deterministic by construction: on the recorded platform this
#    must leave `git status` clean.
tools/golden-gen/.venv/bin/python tools/golden-gen/gen_fixtures.py
git diff --stat tests/golden/fixtures/     # expected: empty
```

Step 5 rebuilds the tabular tables, which takes roughly 100 s each. To refresh
only the manifest, or one fixture, without paying that:

```bash
cd tools/golden-gen
./.venv/bin/python -c "import gen_fixtures as g; g.write_manifest()"
./.venv/bin/python -c "import gen_fixtures as g; g.write_jsonl('if97_water.jsonl', g.gen_if97_water())"
```

## If the digest does not match

`verify-oracle.sh` fails. Do not regenerate anything. Either:

- **You are on a different platform.** A macOS or Windows or aarch64 wheel is a
  different binary, and floating-point results are not guaranteed identical
  across them — that is precisely the risk `ci.yml`'s `platform` job exists to
  measure. Regenerate on `Linux-x86_64`, or be prepared to justify every diff
  the regeneration produces, one record at a time.
- **You installed a different build.** Reinstall from `vendor/`, which is the
  only copy whose provenance this repository can vouch for.

## If a future release replaces this oracle

Moving to CoolProp 8.1 (or any other build) is not a version bump. It changes
what "correct" means for every fixture in the tree, and the two wheel-vs-tag
divergences the port encodes are pinned to *this* binary's behaviour. Such a
move needs its own PLAN.md Decisions entry, a full regeneration, and a
record-by-record account of every fixture that moved — the same standard the
divergence table in NEXT-STEPS.md is held to.
