# golden-gen

Generates the committed golden fixtures in `tests/golden/fixtures/` by running the real
CoolProp 8.0.0 Python wheel — the project's primary verification oracle (see PLAN.md).

**`ORACLE.md` is the important file here**: which exact binary the fixtures came from, why
a version number does not identify it, the archived copy in `vendor/`, and the cold-start
walkthrough. Read it before touching a fixture.

Setup (once):

```bash
python3.12 -m venv tools/golden-gen/.venv
tools/golden-gen/.venv/bin/pip install --require-hashes -r tools/golden-gen/requirements.txt
tools/golden-gen/verify-oracle.sh     # confirm it is THE oracle, not just CoolProp 8.0.0
```

Regenerate fixtures (deterministic — rerunning must produce byte-identical files; run from
the repo root):

```bash
tools/golden-gen/.venv/bin/python tools/golden-gen/gen_fixtures.py
```

Fixtures are committed and never edited by hand.
