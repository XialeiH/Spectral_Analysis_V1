#!/bin/bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 STAGING_ROOT [--force]" >&2
    exit 2
fi

staging_root="$1"
force_export="${2:-}"
matlab_bin="/Applications/MATLAB_R2025b.app/bin/matlab"

for figure_file in "$staging_root"/*_16condition_full_eigenspectrum_maps.fig; do
    pdf_file="${figure_file%.fig}.pdf"
    pdf_size=0
    if [[ -f "$pdf_file" ]]; then
        pdf_size=$(stat -f '%z' "$pdf_file")
    fi
    if [[ "$force_export" != "--force" ]] && (( pdf_size >= 100000 )); then
        continue
    fi

    temporary_pdf="${pdf_file%.pdf}.fresh.pdf"
    escaped_figure=${figure_file//\'/\'\'}
    escaped_pdf=${temporary_pdf//\'/\'\'}
    "$matlab_bin" -batch \
        "f=openfig('$escaped_figure','invisible'); set(f,'Units','pixels','Position',[100 100 1800 2600]); exportgraphics(f,'$escaped_pdf','ContentType','vector','Resolution',600); close(f);"

    temporary_size=$(stat -f '%z' "$temporary_pdf")
    if (( temporary_size < 100000 )); then
        echo "Fresh export is still empty: $temporary_pdf" >&2
        exit 1
    fi
    mv "$temporary_pdf" "$pdf_file"
    echo "Re-exported $(basename "$pdf_file")"
done
