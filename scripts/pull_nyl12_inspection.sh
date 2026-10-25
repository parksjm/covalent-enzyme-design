#!/usr/bin/env bash
# Run on the laptop. Download the curated inspection and VMD packages.
set -euo pipefail

usage() {
    echo "Usage: $0 [--dry-run] [USER@HOST [DESTINATION]]"
    echo "By default, pull from $(id -un)@violet into ./nyl12-inspection-20261021"
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
if [[ $# -gt 2 ]]; then
    usage >&2
    exit 2
fi
remote=${1:-$(id -un)@violet}
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
rsync_args=(-rtpvz --safe-links --max-size=10m)
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
fi

# This separately curated package contains only selected VMD inspection frames.
# Its DCD files are intentional; exclusions above apply to the older package.
vmd_source=/data/jpc/covalent-enzyme-design/output_vmd_inspection_20261024
vmd_destination="$destination/vmd_initial_final"
if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$vmd_destination"
fi
vmd_rsync_args=(-rtpvz --safe-links --max-size=100m)
if [[ $dry_run -eq 1 ]]; then
    vmd_rsync_args+=(--dry-run)
fi
rsync "${vmd_rsync_args[@]}" \
    --include='*.gz' --include='*.dcd' --include='*.vmd' \
    --include='*.tcl' --include='*.json' --include='*.md' \
    --include='MANIFEST.sha256' --exclude='*' \
    "$remote:$vmd_source/" "$vmd_destination/"
if [[ $dry_run -eq 0 ]]; then
    (cd "$vmd_destination" && "${checksum[@]}" MANIFEST.sha256)
    echo "Download verified. Read $vmd_destination/README.md for VMD instructions."
fi
