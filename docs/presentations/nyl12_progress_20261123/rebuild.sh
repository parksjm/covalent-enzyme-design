#!/usr/bin/env bash
# Rebuild from the source checkout; no project dependency changes.
set -euo pipefail
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
deck_dir="$repo_root/docs/presentations/nyl12_progress_20261123"
cd "$repo_root"
uv run --no-project --with python-pptx==1.0.2 \
    python scripts/build_nyl12_progress_deck.py \
    --source "$deck_dir/slides.json" \
    --output "$deck_dir/Nyl12_PA6_progress_20261123.pptx"
libreoffice -env:UserInstallation=file:///tmp/nyl12_progress_deck_rebuild \
    --headless --convert-to pdf --outdir "$deck_dir" \
    "$deck_dir/Nyl12_PA6_progress_20261123.pptx"
