#!/usr/bin/env bash
# Push the DLLifting package to GitHub (run on a machine with GitHub access).
set -euo pipefail

REPO_URL="${1:-https://github.com/wangxintong216-Cinty/Sequential-lifting-for-general-knapsack-set.git}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
WORK_DIR="${TMPDIR:-/tmp}/sequential-lifting-publish-$$"

echo "Package: ${PKG_DIR}"
echo "Remote:  ${REPO_URL}"
echo "Work:    ${WORK_DIR}"

rm -rf "${WORK_DIR}"
mkdir -p "${WORK_DIR}"

rsync -a --exclude 'build/' --exclude 'test_dllifting' --exclude 'test_dllifting_reduction' \
  --exclude 'test_isgeq' --exclude 'examples/example_basic' --exclude 'examples/example_c_api' \
  "${PKG_DIR}/" "${WORK_DIR}/"

cd "${WORK_DIR}"
if [[ ! -d .git ]]; then
  git init
  git checkout -b main 2>/dev/null || git branch -M main 2>/dev/null || true
  git add -A
  git commit -F - <<'ENDMSG' || true
Initial release: DLLifting library (DL/DP hybrid knapsack lifting)
ENDMSG
fi

if git remote | grep -q '^origin$'; then
  git remote set-url origin "${REPO_URL}"
else
  git remote add origin "${REPO_URL}"
fi

echo "Pushing to origin main ..."
git push -u origin main

echo "Done: ${REPO_URL}"
