#!/usr/bin/env bash
# Refresh small azalea inspection artifacts on violet; never delete local files.
set -euo pipefail
root=/data/jpc/covalent-enzyme-design/qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/reaction_paths/03_acylenzyme_neb
rsync -a --safe-links --max-size=100m "$@" \
  --exclude='clutter/' --exclude='source_restarts/' --include='*/' --include='latest' \
  --include='*.json' --include='*.csv' --include='*.svg' --include='*.png' \
  --include='*.dcd' --include='*.parm7' --include='*.vmd' \
  --include='*.cpptraj' --include='*_export.log' --include='MANIFEST.sha256' \
  --include='*.py' --include='*.sh' --include='*.mdin' --include='*.RST' \
  --exclude='*' "azalea:$root/" "$root/"
