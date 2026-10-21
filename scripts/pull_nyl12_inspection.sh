#!/usr/bin/env bash
# Run on the laptop. This downloads only the curated inspection package.
set -euo pipefail

usage() {
    echo "Usage: bash $0 [--dry-run] USER@WORKSTATION [DESTINATION]"
    echo "Example: bash $0 jpc@violet ./nyl12-inspection-20261021"
}

dry_run=0
if [[ ${1:-} == --dry-run ]]; then
    dry_run=1
    shift
fi
if [[ ${1:-} == --help || ${1:-} == -h ]]; then
    usage
    exit 0
fi
if [[ $# -lt 1 || $# -gt 2 ]]; then
    usage >&2
    exit 2
fi
remote=$1
destination=${2:-./nyl12-inspection-20261021}
# Keep the SSH/rsync remote specification unambiguous; SSH aliases are supported.
if [[ ! $remote =~ ^[a-zA-Z0-9_][a-zA-Z0-9_.@-]*$ ]]; then
    echo "Use an SSH alias or USER@HOST (no spaces, colons, or options)." >&2
    exit 2
fi
if [[ -z $destination || $destination == -* ]]; then
    echo "Destination must be a directory path, not an option." >&2
    exit 2
fi
command -v rsync >/dev/null
command -v ssh >/dev/null
if command -v shasum >/dev/null; then
    checksum=(shasum -a 256 -c)
elif command -v sha256sum >/dev/null; then
    checksum=(sha256sum -c)
else
    echo "Install a SHA-256 checksum utility before downloading." >&2
    exit 2
fi

source_dir=/data/jpc/covalent-enzyme-design/output_laptop_inspection_20261021
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$destination"
fi
# No --delete: files already on the laptop are never deleted.
# These exclusions are defensive; the workstation package is already curated.
rsync_args=(-rtvz --safe-links --max-size=10m)
if [[ $dry_run -eq 1 ]]; then
    rsync_args+=(--dry-run)
fi
rsync "${rsync_args[@]}" \
    --exclude='*.chk' --exclude='*.fchk' --exclude='*.log' \
    --exclude='*.prmtop' --exclude='*.inpcrd' --exclude='*.rst7' \
    --exclude='*.nc' --exclude='*.dcd' --exclude='scratch/' \
    "$remote:$source_dir/" "$destination/"
if [[ $dry_run -eq 0 ]]; then
    (cd "$destination" && "${checksum[@]}" MANIFEST.sha256)
    echo "Download verified. Read $destination/README.md for PyMOL instructions."
fi
