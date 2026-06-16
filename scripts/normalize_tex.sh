#!/usr/bin/env bash
# normalize_tex.sh
# Automates the cleanup of squished/dense subscripts and un-nested math operators

set -euo pipefail

usage() {
    echo "Usage: $0 [path ...]" >&2
    echo "  Normalize LaTeX syntax in one or more .tex files or directories." >&2
    echo "  Defaults to: drafts" >&2
}

TARGETS=("$@")
if [ ${#TARGETS[@]} -eq 0 ]; then
    TARGETS=("drafts")
fi

for target in "${TARGETS[@]}"; do
    if [ ! -e "$target" ]; then
        echo "Error: Path '$target' not found." >&2
        usage
        exit 1
    fi
done

collect_tex_files() {
    local target
    for target in "${TARGETS[@]}"; do
        if [ -d "$target" ]; then
            find "$target" -type f -name "*.tex"
        elif [ -f "$target" ] && [[ "$target" == *.tex ]]; then
            printf '%s\n' "$target"
        fi
    done
}

TEX_FILES=()
while IFS= read -r tex_file; do
    TEX_FILES+=("$tex_file")
done < <(collect_tex_files | sort -u)

if [ ${#TEX_FILES[@]} -eq 0 ]; then
    echo "No .tex files found in target path(s): ${TARGETS[*]}"
    exit 0
fi

echo "==> Running LaTeX subscript and operator cleanup on ${#TEX_FILES[@]} file(s)..."

changed_files=0

# Find all .tex files and execute in-place substitutions
for tex_file in "${TEX_FILES[@]}"; do
    echo "  Processing: $tex_file"

    # macOS and GNU-compatible in-place edit via explicit backup suffix.
    cp "$tex_file" "${tex_file}.pre-normalize.bak"

    # 1. Compress long parameter subscripts to crisp abbreviations
    sed -i.bak 's/z_{\\mathrm{center}}/z_c/g' "$tex_file"
    sed -i.bak 's/\\alpha_{\\mathrm{stretch}}/\\alpha_s/g' "$tex_file"
    sed -i.bak 's/\\alpha_s\\text{stretch}/\\alpha_s/g' "$tex_file"
    sed -i.bak 's/L_{\\mathrm{domain}}/L_d/g' "$tex_file"
    sed -i.bak 's/\\epsilon_{\\text{mach}}/\\epsilon_{\\mathrm{mach}}/g' "$tex_file"

    # 2. Convert standard limits to cleaner syntax
    sed -i.bak 's/z_{\\mathrm{min}}/z_{\\min}/g' "$tex_file"
    sed -i.bak 's/z_{\\mathrm{max}}/z_{\\max}/g' "$tex_file"

    # 3. Inject proper operational spacing wrappers
    sed -i.bak 's/\\text{sech}/\\operatorname{sech}/g' "$tex_file"
    sed -i.bak 's/\\text{diag}/\\operatorname{diag}/g' "$tex_file"

    # 4. Harmonize vector profiles (A_obs to a_obs)
    sed -i.bak 's/\\mathbf{A}_{\\mathrm{obs}}/\\mathbf{a}_{\\mathrm{obs}}/g' "$tex_file"
    sed -i.bak 's/\\|\\mathbf{A}_{\\mathrm{obs}}\\|/\\|\\mathbf{a}_{\\mathrm{obs}}\\|/g' "$tex_file"

    if cmp -s "$tex_file" "${tex_file}.pre-normalize.bak"; then
        rm -f "${tex_file}.pre-normalize.bak"
    else
        changed_files=$((changed_files + 1))
        rm -f "${tex_file}.pre-normalize.bak"
    fi

    rm -f "${tex_file}.bak"
done

echo "==> Complete! Normalized ${#TEX_FILES[@]} TeX file(s); changed ${changed_files}."