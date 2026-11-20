"""Read-only SSH inventory and compact rsync staging for static remote owners."""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess

from .neb_umbrella_dispatch import ROOT, digest, read, write
from .neb_umbrella_sync import (
    CRITICAL,
    context,
    identity,
    publish_block,
    publish_latest,
    publish_view,
    published,
    verify_hashes,
)

# Executed through SSH stdin-free argv; this code only reads the remote campaign.
REMOTE_INVENTORY = r"""
import hashlib,json,pathlib,re,sys
root=pathlib.Path(sys.argv[1]); request=json.loads(sys.argv[2])
def digest(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(1048576),b""): h.update(block)
    return h.hexdigest()
def hashes(folder,names):
    result={}
    for name in names:
        p=folder/name
        if p.is_symlink() or not p.is_file(): raise ValueError("Missing/linked payload: "+str(p))
        result[name]=digest(p)
    return result
prepared=pathlib.Path(request["prepared"])
if digest(prepared)!=request["prepared_sha256"]: raise ValueError("Remote prepared identity changed")
if (root/"prepared.sha256").read_text().strip()!=request["prepared_sha256"]: raise ValueError("Remote prepared receipt changed")
result={"blocks":[],"views":[],"statuses":{}}
for window,stages in request["blocks"].items():
    for stage in stages:
        folder=root/window/stage
        if (folder/"analysis_request.json").is_file():
            result["blocks"].append({"window":window,"stage":stage,"sha256":hashes(folder,request["critical"])})
    for generation in sorted((root/window/"progress").glob("block_*")):
        manifest=generation/"MANIFEST.sha256"
        if not re.fullmatch(r"block_\d{2}",generation.name) or not manifest.is_file(): continue
        raw=manifest.read_text()
        if not raw or not raw.endswith("\n"): continue
        entries={}
        for line in raw.splitlines():
            expected,name=line.split("  ",1)
            if name not in request["view_files"] or name in entries: raise ValueError("Unexpected view payload")
            entries[name]=expected
        if not {"view.pdb","view_umbrella.vmd","manifest.json"} <= entries.keys(): raise ValueError("Incomplete view manifest")
        measured=hashes(generation,entries)
        if measured!=entries: raise ValueError("Remote view manifest hash mismatch")
        measured["MANIFEST.sha256"]=digest(manifest)
        result["views"].append({"window":window,"generation":generation.name,"sha256":measured})
    status=root/window/"window_status.json"
    if status.is_file():
        raw=status.read_bytes()
        if raw.endswith(b"\n"): result["statuses"][window+"/window_status.json"]=json.loads(raw)
for index in range(1,25):
    name=f"dispatch/{request['remote']}/worker{index:02d}/status.json"
    status=root/name
    if status.is_file():
        raw=status.read_bytes()
        if raw.endswith(b"\n"): result["statuses"][name]=json.loads(raw)
sys.stdout.write(json.dumps(result,allow_nan=False)+"\n")
"""
VIEW_FILES = {"view.pdb", "endpoints.xyz", "view_umbrella.vmd", "manifest.json"}
SOURCE_NAMES = [
    "scripts/sync_neb_umbrella_results.py",
    "src/nyl12_covalent/__init__.py",
    "src/nyl12_covalent/jobs.py",
    "src/nyl12_covalent/neb_umbrella_dispatch.py",
    "src/nyl12_covalent/neb_umbrella_sync.py",
    "src/nyl12_covalent/neb_umbrella_sync_transport.py",
]


def inventory(remote: str, root: Path, request: dict) -> dict:
    """Read remote hashes/statuses only, failing explicitly on SSH or identity errors."""
    command = shlex.join(
        ["python3", "-c", REMOTE_INVENTORY, str(root), json.dumps(request)]
    )
    result = subprocess.run(
        ["ssh", "-oBatchMode=yes", "-oConnectTimeout=20", remote, command],
        check=True,
        capture_output=True,
        text=True,
        timeout=120,
    )
    return json.loads(result.stdout)


def fetch(remote: str, root: Path, incoming: Path, names: list[str]) -> None:
    """Rsync only explicit verified relative names, preserving all remote raw data."""
    if not names:
        return
    incoming.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "rsync",
            "-rt",
            "--checksum",
            "--protect-args",
            "--from0",
            "--files-from=-",
            "-e",
            "ssh -oBatchMode=yes -oConnectTimeout=20",
            "--",
            f"{remote}:{root}/",
            str(incoming) + "/",
        ],
        input=("\0".join(names) + "\0").encode(),
        check=True,
        timeout=3600,
    )


def validate_inventory(data: dict, blocks: dict, remote: str) -> None:
    """Reject foreign ownership, traversal, duplicate blocks, and noncompact files."""
    seen = set()
    for row in data["blocks"]:
        key = (row["window"], row["stage"])
        if (
            key in seen
            or row["window"] not in blocks
            or row["stage"] not in blocks[row["window"]]
            or set(row["sha256"]) != set(CRITICAL)
        ):
            raise ValueError(
                "Inventory violates production-block ownership or payload contract"
            )
        seen.add(key)
    seen.clear()
    for row in data["views"]:
        key = (row["window"], row["generation"])
        if (
            key in seen
            or row["window"] not in blocks
            or not re.fullmatch(r"block_\d{2}", row["generation"])
            or not set(row["sha256"]) <= VIEW_FILES | {"MANIFEST.sha256"}
            or not {"MANIFEST.sha256", "manifest.json", "view.pdb", "view_umbrella.vmd"}
            <= set(row["sha256"])
        ):
            raise ValueError(
                "Inventory violates compact-view ownership or payload contract"
            )
        seen.add(key)
    for name in data["statuses"]:
        parts = Path(name).parts
        window = (
            len(parts) == 2 and parts[0] in blocks and parts[1] == "window_status.json"
        )
        worker = bool(
            re.fullmatch(
                rf"dispatch/{re.escape(remote)}/worker(?:0[1-9]|1\d|2[0-4])/status.json",
                name,
            )
        )
        if not window and not worker:
            raise ValueError("Inventory status violates remote ownership")


def sync_once(path: Path, remote: str) -> dict[str, int]:
    """Stage, verify and publish completed remote-owned compact artifacts once."""
    plan, cfg, prepared, blocks = context(path, remote)
    root = prepared.parent
    request = {
        "prepared": plan["hosts"][remote]["prepared"],
        "prepared_sha256": plan["hosts"][remote]["sha256"],
        "blocks": blocks,
        "critical": CRITICAL,
        "remote": remote,
        "view_files": sorted(VIEW_FILES),
    }
    data = inventory(remote, root, request)
    validate_inventory(data, blocks, remote)
    pending = []
    for kind, rows in (("block", data["blocks"]), ("view", data["views"])):
        for row in rows:
            relative = Path(row["window"]) / (
                row["stage"] if kind == "block" else "progress/" + row["generation"]
            )
            expected = identity(path, remote, str(root / relative), row["sha256"])
            gate = "acceptance.json" if kind == "block" else "MANIFEST.sha256"
            if not published(root / relative, row["sha256"], expected, gate):
                pending.append((kind, relative, row["sha256"], expected))
    incoming = path.parent / "incoming" / remote
    names = [
        str(relative / name) for _, relative, hashes, _ in pending for name in hashes
    ]
    fetch(remote, root, incoming, names)
    counts = {"blocks": 0, "views": 0, "statuses": len(data["statuses"])}
    for kind, relative, hashes, expected in pending:
        if kind == "block":
            changed = publish_block(
                incoming / relative, root / relative, hashes, cfg, expected
            )
        else:
            changed = publish_view(
                incoming / relative, root / relative, hashes, expected
            )
        counts["blocks" if kind == "block" else "views"] += int(changed)
    for window in sorted({row["window"] for row in data["views"]}):
        newest = max(
            row["generation"] for row in data["views"] if row["window"] == window
        )
        publish_latest(root / window / "progress" / newest)
    for name, status in data["statuses"].items():
        destination = path.parent / "remote_status" / remote / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        write(destination, status)
    write(
        path.parent / f"sync_{remote}_status.json",
        {
            "status": "completed_cycle",
            **counts,
            "plan_sha256": digest(path),
            "observed_utc": datetime.now(timezone.utc).isoformat(),
        },
    )
    return counts


def freeze_sources(path: Path, remote: str) -> Path:
    """Freeze the small synchronizer source tree without editing the dispatch plan."""
    context(path, remote)
    snapshot = path.parent / "sync_source_snapshot"
    snapshot.mkdir(exist_ok=False)
    for name in SOURCE_NAMES:
        destination = snapshot / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / name, destination)
    hashes = {str(snapshot / name): digest(snapshot / name) for name in SOURCE_NAMES}
    write(snapshot / "source_sha256.json", hashes)
    (snapshot / "source_sha256.sha256").write_text(
        digest(snapshot / "source_sha256.json") + "\n"
    )
    return snapshot / SOURCE_NAMES[0]


def verify_sources(script: Path) -> None:
    """Require execution from its retained source snapshot with pinned manifest."""
    snapshot = script.resolve().parents[1]
    manifest = snapshot / "source_sha256.json"
    if digest(manifest) != (snapshot / "source_sha256.sha256").read_text().strip():
        raise ValueError("Synchronizer source manifest changed")
    verify_hashes(
        Path("/"),
        {
            str(Path(name).relative_to("/")): expected
            for name, expected in read(manifest).items()
        },
    )
