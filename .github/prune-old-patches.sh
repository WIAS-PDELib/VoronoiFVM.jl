#!/usr/bin/env bash
#
# prune-old-patches.sh
#
# Purpose
# -------
# Documenter's deploydocs() adds a new vX.Y.Z/ directory to the gh-pages
# branch for every release and never deletes old ones. The whole branch is
# uploaded as ONE GitHub Pages artifact, which must stay below 1 GB.
# This script keeps only the NEWEST patch release of each minor version
# (e.g. of v1.2.0, v1.2.1 it keeps v1.2.1) and removes the rest.
#
# What it touches
# ---------------
#   - Real directories named vMAJOR.MINOR.PATCH   -> removed unless newest
#   - Symlinks vMAJOR.MINOR (e.g. v1.2 -> v1.2.1) -> repointed if they
#     referenced a removed directory
#   - Anything else (stable, latest, dev, v1, index.html, versions.js, ...)
#     is left alone. stable/v1 normally point at vX.Y symlinks, which stay valid.
#
# Usage (from the ROOT of a gh-pages checkout)
# --------------------------------------------
#   ./prune-old-patches.sh           # dry run, only prints what would happen
#   ./prune-old-patches.sh --apply   # actually remove / relink via git
#
# Requires: bash >= 4 (associative arrays), GNU find (for -printf), sort -V.
# On macOS use a GNU toolchain (brew install bash findutils coreutils).
########################################################################
# How to run it and push the result
# 1. Work in a separate clone or worktree so your normal checkout is untouched
#   git clone --branch gh-pages --single-branch https://github.com/WIAS-PDELib/VoronoiFVM.jl.git VoronoiFVM-gh-pages
#   cd VoronoiFVM-gh-pages
# 
# # 2. Dry run, then review the output
#   chmod +x ../prune-old-patches.sh
#   ../prune-old-patches.sh
# 
# # 3. Apply, check, and commit
#    ../prune-old-patches.sh --apply
#    git status --short | head
#    du -sh --exclude=.git .          # should now be well below 1 GB
#    git commit -m "Prune old patch release docs (keep newest patch per minor)"
# 
# # 4. Push
#   git push origin gh-pages
#########################################################################
set -euo pipefail

apply=0
[[ "${1:-}" == "--apply" ]] && apply=1

# Safety check: must be a git work tree containing versions.js (Documenter site)
if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "error: not inside a git work tree" >&2; exit 1
fi
if [[ ! -f versions.js ]]; then
    echo "error: versions.js not found; run from the root of the gh-pages checkout" >&2
    exit 1
fi

# 1. Collect real vX.Y.Z directories (symlinks are excluded by -type d),
#    sorted by version so the last entry per minor is the newest.
mapfile -t dirs < <(find . -maxdepth 1 -type d -regextype posix-extended \
    -regex '\./v[0-9]+\.[0-9]+\.[0-9]+' -printf '%f\n' | sort -V)

# 2. Remember the newest patch dir per minor (v1.2 -> v1.2.1).
declare -A newest
for d in "${dirs[@]}"; do
    newest[${d%.*}]="$d"      # input is sorted ascending: last one wins
done

# 3. Remove every vX.Y.Z that is not the newest of its minor.
removed=0
for d in "${dirs[@]}"; do
    if [[ "${newest[${d%.*}]}" != "$d" ]]; then
        echo "remove $d (keeping ${newest[${d%.*}]})"
        removed=$((removed + 1))
        if ((apply)); then git rm -rq -- "$d"; fi
    fi
done

# 4. Repoint vX.Y symlinks that do not reference the kept directory.
for minor in "${!newest[@]}"; do
    if [[ -L "$minor" && "$(readlink "$minor")" != "${newest[$minor]}" ]]; then
        echo "relink $minor -> ${newest[$minor]}"
        if ((apply)); then
            ln -sfn "${newest[$minor]}" "$minor"
            git add -- "$minor"
        fi
    fi
done

echo "$removed directories selected for removal."
((apply)) || echo "Dry run only; pass --apply to modify the working tree."


######
