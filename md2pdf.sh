#!/usr/bin/env bash
# Convert one or more Markdown files to PDF via pandoc + xelatex.
# Uses DejaVu fonts (Unicode-compatible) and lang=cs for Czech typography
# (správné dělení slov, uvozovky, mezery).
#
# Použití:
#   md2pdf.sh file.md [file2.md ...]      # bez obsahu
#   md2pdf.sh -t file.md [file2.md ...]   # s obsahem (table of contents)

set -euo pipefail

toc_args=()
if [[ "${1:-}" == "-t" ]]; then
    toc_args=(--toc --toc-depth=2)
    shift
fi

if [[ $# -eq 0 ]]; then
    echo "Usage: $0 [-t] file1.md [file2.md ...]" 1>&2
    exit 2
fi

command -v pandoc   >/dev/null || { echo "ERROR: pandoc not installed" 1>&2; exit 1; }
command -v xelatex  >/dev/null || { echo "ERROR: xelatex not installed" 1>&2; exit 1; }

for name in "$@"; do
    if [[ "$name" != *.md ]]; then
        echo "ERROR: $name is not a .md file" 1>&2
        exit 1
    fi

    out="${name%.md}.pdf"
    echo "Convert ${name} -> ${out}"

    pandoc "$name" -o "$out" \
        --pdf-engine=xelatex \
        "${toc_args[@]}" \
        --variable lang=cs \
        --variable geometry:margin=2.4cm \
        --variable fontsize=11pt \
        --variable colorlinks=true \
        --variable linkcolor=blue \
        --variable urlcolor=blue \
        --variable mainfont="DejaVu Serif" \
        --variable monofont="DejaVu Sans Mono"
done
