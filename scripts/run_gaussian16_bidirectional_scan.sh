#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
    echo "usage: $0 FORWARD_GJF FORWARD_LOG FORWARD_SCRATCH REVERSE_GJF REVERSE_LOG REVERSE_SCRATCH" >&2
    exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

"${script_dir}/run_gaussian16_file.sh" "$1" "$2" "$3"
"${script_dir}/run_gaussian16_file.sh" "$4" "$5" "$6"
