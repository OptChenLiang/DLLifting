#!/usr/bin/env bash
# Publish/update the DLLifting package on GitLab (preserves remote history).
set -euo pipefail

REPO_URL="${1:-git@159.226.92.34:wangxintong/dllifting.git}"
COMMIT_MSG="${2:-Release 1.1.0: forced isdl_mode (DL/DP) for bounded and unbounded lifting}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
WORK_DIR="${TMPDIR:-/tmp}/dllifting-publish-$$"
MSG_FILE="${WORK_DIR}/commit_msg.txt"

echo "Package: ${PKG_DIR}"
echo "Remote:  ${REPO_URL}"
echo "Work:    ${WORK_DIR}"

rm -rf "${WORK_DIR}"
git clone "${REPO_URL}" "${WORK_DIR}"

rsync -a --delete \
  --exclude 'build/' \
  --exclude '.git/' \
  --exclude 'libdllifting.so' \
  --exclude 'test_dllifting' \
  --exclude 'test_isgeq' \
  --exclude 'test_mixed_vars' \
  --exclude 'test_mixed_vars_r' \
  --exclude 'example' \
  --exclude 'results_mixed_*.txt' \
  --exclude 'tests/results_mixed_*.txt' \
  "${PKG_DIR}/" "${WORK_DIR}/"

cd "${WORK_DIR}"
git add -A

if git diff --cached --quiet; then
  echo "No changes to push."
  rm -rf "${WORK_DIR}"
  exit 0
fi

cat > "${MSG_FILE}" <<EOF
${COMMIT_MSG}
EOF

# Use -F to avoid shell/git alias issues with -m on older git.
/usr/bin/git commit -F "${MSG_FILE}" || git commit -F "${MSG_FILE}"

git push origin main

echo "Done: ${REPO_URL}"
git log -1 --oneline
rm -rf "${WORK_DIR}"
