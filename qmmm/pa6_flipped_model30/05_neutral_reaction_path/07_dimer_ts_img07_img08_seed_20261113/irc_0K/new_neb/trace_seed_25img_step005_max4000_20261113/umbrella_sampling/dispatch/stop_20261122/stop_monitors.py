"""Stop only verified completed-NEB-US monitoring processes; retain all data."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import signal
import time

ROOT = Path(__file__).resolve().parents[2]
TARGETS = [(4100242, "analyze_neb_umbrella.py"), (4098079, "sync_neb_umbrella_results.py")]


def main() -> None:
    """Validate command/cwd identities, signal monitors and retain closeout receipt."""
    rows = []
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT / "prepared.json", ROOT / "analysis/pmf_latest/provenance.json"]}
    for pid, name in TARGETS:
        proc = Path("/proc") / str(pid)
        args = (proc / "cmdline").read_bytes().replace(b"\0", b" ").decode()
        cwd = os.readlink(proc / "cwd")
        if cwd != str(ROOT) or name not in args or str(ROOT) not in args:
            raise RuntimeError(f"Monitor identity differs: {pid}")
        rows.append({"pid": pid, "command": args, "cwd": cwd, "signal": "SIGTERM"})
    receipt = {"observed_utc": datetime.now(timezone.utc).isoformat(), "reason": "PI-requested closeout after all 63 windows reached 30 ps; sampling convergence not established", "sampling_workers_interrupted": 0, "monitors": rows, "protected_sha256": hashes}
    output = Path(__file__).with_name("violet_monitor_stop_receipt.json")
    if output.exists():
        raise FileExistsError(output)
    for row in rows:
        os.kill(row["pid"], signal.SIGTERM)
    time.sleep(1)
    for row in rows:
        row["absent_after_signal"] = not Path(f"/proc/{row['pid']}").exists()
    for filename, digest in hashes.items():
        if hashlib.sha256(Path(filename).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"Protected artifact changed: {filename}")
    output.write_text(json.dumps(receipt, indent=2) + "\n")
    if not all(row["absent_after_signal"] for row in rows):
        raise RuntimeError("A signaled monitor is still present; inspect receipt")


if __name__ == "__main__":
    main()
