#!/bin/bash
# Compare DL, DP, DL+R, DP+R on mixed binary/bounded/unbounded instances.
set -euo pipefail
cd "$(dirname "$0")/.."

make test_mixed_vars test_mixed_vars_r

OUT_N=tests/results_mixed_N.txt
OUT_R=tests/results_mixed_R.txt

./test_mixed_vars  | tee "$OUT_N"
./test_mixed_vars_r | tee "$OUT_R"

echo ""
echo "=== Four-algorithm summary (time seconds, rhs) ==="
python3 - <<'PY'
import re
from collections import defaultdict

pat = re.compile(
    r"case=(\S+) build=(\S+) algo=(\S+) n=(\d+) ok=(\d+) time=([\d.e+-]+) "
    r"rhs=([\d.e+-]+) red=(\d+) ub=(\d+)")

def load(path):
    d = defaultdict(dict)
    for line in open(path):
        m = pat.search(line)
        if not m:
            continue
        case, build, algo, n, ok, t, rhs, red, ub = m.groups()
        d[case][(build, algo)] = dict(
            n=int(n), ok=int(ok), t=float(t), rhs=float(rhs), red=int(red), ub=int(ub))
    return d

n = load("tests/results_mixed_N.txt")
r = load("results_mixed_R.txt")
cases = sorted(set(n) | set(r))

print(f"{'case':<20} {'n':>3} {'DL':>9} {'DP':>9} {'DL+R':>9} {'DP+R':>9}  agree")
print("-" * 78)
for case in cases:
    nval = 0
    for src in (n, r):
        for key in src.get(case, {}):
            nval = max(nval, src[case][key].get("n", 0))
    def g(build, algo, key):
        src = n if build == "N" else r
        row = src.get(case, {}).get((build, algo))
        if not row or not row["ok"]:
            return "   fail "
        return f"{row[key]:8.4f}"

    dl = n.get(case, {}).get(("N", "DL"))
    dp = n.get(case, {}).get(("N", "DP"))
    dlr = r.get(case, {}).get(("R", "DL"))
    dpr = r.get(case, {}).get(("R", "DP"))
    agree = "yes" if dl and dp and dl["ok"] and dp["ok"] and abs(dl["rhs"] - dp["rhs"]) < 1e-4 else "no"
    print(f"{case:<20} {nval:3d} {g('N','DL','t')} {g('N','DP','t')} {g('R','DL','t')} {g('R','DP','t')}  {agree}")

ok_n = sum(1 for line in open("tests/results_mixed_N.txt") if line.startswith("OK:"))
warn_n = sum(1 for line in open("tests/results_mixed_N.txt") if "WARN:" in line)
print(f"\nStats (build N): OK={ok_n} WARN={warn_n}")
print("Notes: build N=DL/DP; R=DL+R/DP+R; ub=1 => red off; times in seconds.")
PY

echo ""
echo "Logs: $OUT_N , $OUT_R"
