#!/usr/bin/env bash
# Rebuild reviewed frozen evidence without rerunning scientific calculations.
set -euo pipefail
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
deck_dir="$repo_root/docs/presentations/nyl12_progress_20261127"
build_date=$(TZ=America/New_York date +%Y%m%d)
build_stamp=$(TZ=America/New_York date +%Y%m%dT%H%M%S)
output_base="$deck_dir/Nyl12_PA6_progress_${build_date}"
archive_dir="$repo_root/clutter/presentation_rebuilds/$build_stamp"
cd "$repo_root"
for artifact in "$deck_dir"/Nyl12_PA6_progress_*.pptx "$deck_dir"/Nyl12_PA6_progress_*.pdf "$deck_dir"/Nyl12_PA6_progress_*.build.json "$deck_dir/rendered"; do
    [[ -e "$artifact" ]] || continue
    mkdir -p "$archive_dir"
    mv "$artifact" "$archive_dir/"
done
uv run --offline python "$deck_dir/plot_evidence.py"
PYTHONPATH="$repo_root" uv run --offline --no-project --with python-pptx==1.0.2 python "$deck_dir/build_deck.py" --source "$deck_dir/slides.json" --output "$output_base.pptx"
libreoffice -env:UserInstallation=file:///tmp/nyl12_progress_20261127_rebuild --headless --convert-to pdf --outdir "$deck_dir" "$output_base.pptx"
mkdir -p "$deck_dir/rendered"
pdftoppm -png -r 110 "$output_base.pdf" "$deck_dir/rendered/slide"
uv run --offline --no-project --with python-pptx==1.0.2 python scripts/make_presentation_contact_sheet.py --source-dir "$deck_dir/rendered" --output "$deck_dir/rendered/contact_sheet.png"
