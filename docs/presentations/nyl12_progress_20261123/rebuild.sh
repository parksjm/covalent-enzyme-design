#!/usr/bin/env bash
# Build the dated presentation and previews; no project dependencies are added.
set -euo pipefail
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
deck_dir="$repo_root/docs/presentations/nyl12_progress_20261123"
build_date=$(TZ=America/New_York date +%Y%m%d)
build_stamp=$(TZ=America/New_York date +%Y%m%dT%H%M%S)
output_base="$deck_dir/Nyl12_PA6_progress_${build_date}"
archive_dir="$repo_root/clutter/presentation_rebuilds/$build_stamp"
cd "$repo_root"

# Preserve older dated outputs and renders centrally before replacing anything.
for old_artifact in "$deck_dir"/Nyl12_PA6_progress_*.pptx \
    "$deck_dir"/Nyl12_PA6_progress_*.pdf \
    "$deck_dir"/Nyl12_PA6_progress_*.build.json; do
    [[ -f "$old_artifact" ]] || continue
    mkdir -p "$archive_dir"
    mv "$old_artifact" "$archive_dir/"
done
if [[ -d "$deck_dir/rendered" ]]; then
    mkdir -p "$archive_dir"
    mv "$deck_dir/rendered" "$archive_dir/rendered"
fi

uv run python scripts/plot_nyl12_presentation.py
uv run --offline --no-project --with python-pptx==1.0.2 \
    python -m scripts.build_nyl12_progress_deck \
    --source "$deck_dir/slides.json" \
    --output "$output_base.pptx"
libreoffice -env:UserInstallation=file:///tmp/nyl12_progress_deck_rebuild \
    --headless --convert-to pdf --outdir "$deck_dir" "$output_base.pptx"
mkdir -p "$deck_dir/rendered"
pdftoppm -png -r 140 "$output_base.pdf" "$deck_dir/rendered/slide"
uv run --offline --no-project --with python-pptx==1.0.2 \
    python scripts/make_presentation_contact_sheet.py \
    --source-dir "$deck_dir/rendered" \
    --output "$deck_dir/rendered/contact_sheet.png"
