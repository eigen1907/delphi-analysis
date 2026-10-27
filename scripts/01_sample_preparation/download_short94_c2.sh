#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

# Download the official short94_c2 files; retain originals and resume .part files.
python3 -B - <<'PY'
from pathlib import Path
import json
import subprocess
import urllib.request
import zlib

output = Path("data/20260828_florian/20261026_short94_c2")
with urllib.request.urlopen("https://opendata.cern.ch/api/records/81431", timeout=30) as response:
    files = json.load(response)["metadata"]["files"]


def valid(path, info):
    if not path.exists() or path.stat().st_size != info["size"]:
        return False
    checksum = 1
    with path.open("rb") as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
            checksum = zlib.adler32(block, checksum)
    return f"adler32:{checksum & 0xffffffff:08x}" == info["checksum"]


for index, info in enumerate(files, 1):
    name = Path(info["key"]).name
    target = output / Path(name).stem / name
    target.parent.mkdir(parents=True, exist_ok=True)
    if valid(target, info):
        print(f"[{index}/{len(files)}] Verified, skipping: {name}", flush=True)
        continue

    partial = target.with_name(name + ".part")
    if not valid(partial, info):
        # A complete but corrupt temporary file needs a fresh transfer.
        if partial.exists() and partial.stat().st_size >= info["size"]:
            partial.unlink()
        url = info["uri"].replace(
            "root://eospublic.cern.ch//", "https://eospublic.cern.ch/"
        )
        print(f"[{index}/{len(files)}] Downloading: {name}", flush=True)
        subprocess.run([
            "curl", "-fL", "--retry", "5", "--connect-timeout", "30",
            "-C", "-", "-o", str(partial), url,
        ], check=True)
    if not valid(partial, info):
        raise RuntimeError(f"Checksum mismatch: {partial}; rerun to download it again")
    partial.replace(target)
    print(f"[{index}/{len(files)}] Verified: {name}", flush=True)

print(f"All {len(files)} files downloaded and verified.")
PY
