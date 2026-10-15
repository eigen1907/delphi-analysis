#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.env"
cd "$PROJECT"
mkdir -p "$CACHE/build" "$OUTPUT"

docker build --platform linux/amd64 --build-context "generator=$SIMULATION/pythia8" \
    -t "$IMAGE" "$SAMPLE"
docker pull --platform linux/arm64 "$SIM_IMAGE"
# DELPHI conditions are stored as compressed CVMFS chunks.
uv run python - "$CACHE/condition-data" <<'PY'
"""Download the pinned DELPHI conditions without mounting CVMFS."""

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import hashlib
import os
import sys
import sqlite3
import tempfile
import time
import urllib.request
import zlib


SOURCE = "http://cvmfs-stratum-one.cern.ch:8000/cvmfs/delphi.cern.ch/data/"
CATALOG = "368a3fa512b044c9748aa12ec6fe40999ea262b3"
FILES = {
    ".cvmfscatalog": (0, "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"),
    "DBmisc.dat": (40960, "1afbce404793da90aa16f37fe764701b196081bf92b224990570a7fb7a9aaf8f"),
    "DBsysf.dat": (131072, "eabe267a5a9a6548681acfe25e9c09594cfdb0e11c502d115c3d1d010049c630"),
    "DBscon.dat": (450560, "7040b07fb11006e2e117b0737b6a638904659c5709e2b4d87b408716977e88e4"),
    "DBrunt.dat": (20307968, "67893c971ae1e27a885cfec71d5315a8433cc3f6f4005049923f788e5374f7d5"),
    "DBgeom.dat": (32432128, "67b7c8bbd666742fd1997b7303c339fd404aea36b45236d48249af10af250856"),
    "DBcalb.dat": (1497251840, "7530be62eda22296f15ce876a82b0d49149768288b5b7de20aec1592e76b5274"),
    "DBlepm.dat": (1655250944, "6f0f3a96324e1952e93037d19cf8cc7f39c6abfb895ef14e085569ce4ca26b70"),
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download_object(digest, suffix=""):
    url = SOURCE + digest[:2] + "/" + digest[2:] + suffix
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                compressed = response.read()
            if hashlib.sha1(compressed).hexdigest() != digest:
                raise ValueError(f"CVMFS object checksum mismatch: {url}")
            return zlib.decompress(compressed)
        except (OSError, ValueError, zlib.error):
            if attempt == 2:
                raise
            time.sleep(attempt + 1)


def fetch(destination):
    """Return a verified conditions directory; existing valid files are reused."""
    destination = Path(destination).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    missing = []
    for name, (size, checksum) in FILES.items():
        path = destination / name
        if (path.is_file() and path.stat().st_size == size
                and sha256(path) == checksum):
            print(f"Verified {name}", flush=True)
        else:
            missing.append(name)
    if not missing:
        return destination

    catalog = download_object(CATALOG, "C")
    with tempfile.TemporaryDirectory(prefix="delphi-conditions-") as temporary:
        catalog_path = Path(temporary) / "catalog.db"
        catalog_path.write_bytes(catalog)
        connection = sqlite3.connect(catalog_path)
        connection.row_factory = sqlite3.Row
        transfers = []
        handles = {}
        for name in missing:
            row = connection.execute("SELECT * FROM catalog WHERE name=?", (name,)).fetchone()
            if row is None or row["size"] != FILES[name][0]:
                raise ValueError(f"Unexpected condition catalog entry: {name}")
            handle = os.open(destination / (name + ".partial"),
                             os.O_CREAT | os.O_RDWR | os.O_TRUNC, 0o644)
            os.ftruncate(handle, row["size"])
            handles[name] = handle
            chunks = list(connection.execute(
                "SELECT offset, size, hash FROM chunks WHERE md5path_1=? AND md5path_2=? ORDER BY offset",
                (row["md5path_1"], row["md5path_2"]),
            ))
            offset = 0
            for chunk in chunks:
                if chunk["offset"] != offset:
                    raise ValueError(f"Non-contiguous condition chunks: {name}")
                transfers.append((name, offset, chunk["size"], chunk["hash"].hex(), "P"))
                offset += chunk["size"]
            if chunks and offset != row["size"]:
                raise ValueError(f"Incomplete condition chunks: {name}")
            if not chunks:
                transfers.append((name, 0, row["size"], row["hash"].hex(), ""))
        connection.close()

        def write_chunk(transfer):
            name, offset, size, digest, suffix = transfer
            content = download_object(digest, suffix)
            if len(content) != size or os.pwrite(handles[name], content, offset) != size:
                raise ValueError(f"Condition chunk size mismatch: {name} at {offset}")

        print(f"Downloading {len(missing)} condition files ({len(transfers)} objects)", flush=True)
        try:
            with ThreadPoolExecutor(max_workers=6) as workers:
                futures = [workers.submit(write_chunk, item) for item in transfers]
                for count, future in enumerate(as_completed(futures), 1):
                    future.result()
                    if count % 100 == 0 or count == len(transfers):
                        print(f"Condition objects: {count}/{len(transfers)}", flush=True)
        finally:
            for handle in handles.values():
                os.close(handle)
        for name in missing:
            partial = destination / (name + ".partial")
            if sha256(partial) != FILES[name][1]:
                raise ValueError(f"Condition file checksum mismatch: {name}")
            partial.replace(destination / name)
            print(f"Downloaded and verified {name}", flush=True)
    return destination


fetch(Path(sys.argv[1]))
PY

git -C "$SIMULATION" log -1 --oneline
git -C "$NANOAOD" log -1 --oneline
git -C "$NANOAOD" diff HEAD > "$OUTPUT/nanoaod.patch"
docker run --init --rm --platform linux/amd64 \
    -v "$NANOAOD:/src:ro" -v "$CACHE/build:/build" "$IMAGE" -c '
    set -e
    cmake -S /src -B /build -DFASTJETDIR=/usr/local
    cmake --build /build --parallel 2 --target delphi-nanoaod delphi-raw-nanoaod
    '
