#!/usr/bin/env bash
set -euo pipefail
# Uses YOUR git authentication. No token or password is stored in this script.
# Clone current main first; preserve history and avoid force pushes.
script_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
destination="${script_dir}-github"
git clone --branch main --single-branch \
  https://github.com/corne00/aagi-growth-tutorial.git "$destination"
if command -v rsync >/dev/null 2>&1; then
  rsync -a --exclude='.git/' --exclude='node_modules/' "$script_dir/" "$destination/"
else
  cp -R "$script_dir/." "$destination/"
fi
git -C "$destination" add -A
git -C "$destination" -c user.name="Corné Verburg" \
  -c user.email="c.verburg@tudelft.nl" commit \
  -m "Build reproducible APSIM wheat Bayesian R workshop"
git -C "$destination" push origin main
printf '\nPushed. GitHub Pages will publish https://corne00.github.io/aagi-growth-tutorial/\n'
