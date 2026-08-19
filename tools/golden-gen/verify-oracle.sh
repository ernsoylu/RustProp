#!/bin/sh
# Confirm the installed CoolProp is the exact binary every committed golden
# came from. See tools/golden-gen/ORACLE.md for why that is not the same
# question as "is CoolProp 8.0.0 installed".
#
# The expected digest is read from tests/golden/fixtures/manifest.json, which
# the generator writes on every run, so this check cannot go stale.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
py="$here/.venv/bin/python"

if [ ! -x "$py" ]; then
    echo "no oracle venv at $here/.venv — see tools/golden-gen/ORACLE.md" >&2
    exit 2
fi

exec "$py" - "$here" <<'PY'
import hashlib
import json
import platform
import sys
from pathlib import Path

here = Path(sys.argv[1])
manifest = json.loads((here.parents[1] / "tests/golden/fixtures/manifest.json").read_text())

import CoolProp
import CoolProp.CoolProp as CP

got = hashlib.sha256(Path(CP.__file__).read_bytes()).hexdigest()
want = manifest["oracle_sha256"]

print(f"extension   {CP.__file__}")
print(f"version     {CoolProp.__version__}  (fixtures: {manifest['coolprop_version']})")
print(f"gitrevision {CoolProp.__gitrevision__}")
print(f"sha256      {got}")
print(f"expected    {want}")
print(f"python      {platform.python_version()}  (fixtures: {manifest['python']})")
print(f"platform    {platform.system()}-{platform.machine()}  (fixtures: {manifest['platform']})")

if got == want:
    print("\nOK — this is the binary the committed fixtures were generated from.")
    sys.exit(0)

this_platform = f"{platform.system()}-{platform.machine()}"
if this_platform != manifest["platform"]:
    print(
        f"\nMISMATCH, and this is a {this_platform} host while the fixtures were\n"
        f"generated on {manifest['platform']}. A different platform's wheel is a\n"
        "different binary by construction, so the digest CANNOT match here. Fixtures\n"
        "regenerated on this host are not guaranteed byte-identical to the committed\n"
        "ones; regenerate on the recorded platform, or expect to justify every diff.",
        file=sys.stderr,
    )
    sys.exit(1)

print(
    "\nMISMATCH on the recorded platform. The installed CoolProp is NOT the binary\n"
    "the fixtures came from — do not regenerate anything until this is resolved.\n"
    "See tools/golden-gen/ORACLE.md.",
    file=sys.stderr,
)
sys.exit(1)
PY
