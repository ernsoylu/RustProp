# golden-gen

Generates the committed golden fixtures in `tests/golden/fixtures/` by running the real
CoolProp 8.0.0 Python wheel — the project's primary verification oracle (see PLAN.md).

Setup (once):

```bash
python3 -m venv tools/golden-gen/.venv
tools/golden-gen/.venv/bin/pip install -r tools/golden-gen/requirements.txt
```

Regenerate fixtures (deterministic — rerunning must produce byte-identical files; run from
the repo root):

```bash
tools/golden-gen/.venv/bin/python tools/golden-gen/gen_fixtures.py
```

Fixtures are committed and never edited by hand.
