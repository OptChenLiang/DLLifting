#!/usr/bin/env bash
set -euo pipefail

REPO_URL="${1:-git@159.226.92.34:wangxintong/dllifting.git}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
WORK_DIR="${TMPDIR:-/tmp}/dllifting-publish-$$"

echo "Package: ${PKG_DIR}"
echo "Remote:  ${REPO_URL}"
echo "Work:    ${WORK_DIR}"

rm -rf "${WORK_DIR}"
mkdir -p "${WORK_DIR}"

rsync -a \
  --exclude 'build/' \
  --exclude 'test_dllifting' \
  --exclude 'test_dllifting_reduction' \
  --exclude 'test_isgeq' \
  --exclude 'test_mixed' \
  --exclude 'test_mixed_vars' \
  --exclude 'test_mixed_vars_r' \
  --exclude 'results_mixed_*.txt' \
  --exclude 'examples/example_basic' \
  --exclude 'examples/example_c_api' \
  "${PKG_DIR}/" "${WORK_DIR}/"

cd "${WORK_DIR}"
git init
git checkout -b main 2>/dev/null || git branch -M main

git add -A
git commit -m "Release DLLifting library (DL/DP hybrid knapsack lifting)"

if git remote | grep -q '^origin$'; then
  git remote set-url origin "${REPO_URL}"
else
  git remote add origin "${REPO_URL}"
fi

echo "Pushing to origin main ..."
git push -u origin main

echo "Done: ${REPO_URL}"
