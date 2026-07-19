"""Tiny CLI:  python3 -m unitri [POINTS_FILE]      (reads stdin if no file)

Reads a 2-D integer point set in essentially any bracket/comma/whitespace
format -- a pasted numpy array, a list of [x, y], or "x y" per line -- and
prints the number of fine triangulations to stdout.  Progress/diagnostics go to
stderr (so the count alone is still pipeable); pass -q to silence them.

    python3 -m unitri points.txt
    pbpaste | python3 -m unitri          # paste a numpy array straight in
    echo "[[0,0],[3,0],[0,3],[3,3]]" | python3 -m unitri
"""
from __future__ import annotations

import argparse
import re
import sys
import time

from .profiles import count_triangulations, points_to_profiles


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        prog="python3 -m unitri",
        description="Count fine triangulations of a 2-D lattice point set.")
    ap.add_argument("file", nargs="?",
                    help="points file (any bracket/comma/space format); "
                         "reads stdin if omitted")
    ap.add_argument("-q", "--quiet", action="store_true",
                    help="print only the count (suppress stderr diagnostics)")
    args = ap.parse_args(argv)

    def log(msg: str) -> None:
        if not args.quiet:
            print(msg, file=sys.stderr, flush=True)

    try:
        text = open(args.file).read() if args.file else sys.stdin.read()
    except OSError as e:
        sys.exit(f"error: {e}")

    nums = [int(t) for t in re.findall(r"-?\d+", text)]
    if not nums or len(nums) % 2:
        sys.exit("error: expected an even count of integers (x y pairs)")
    points = list(zip(nums[0::2], nums[1::2]))
    log(f"read {len(points)} points")

    try:
        m, n, upper, lower = points_to_profiles(points)
        log(f"minimal-width box: m={m}, n={n}")
        log(f"  upper = {upper}")
        log(f"  lower = {lower}")
        log("counting (exact; may take a few seconds for large regions)...")
        t = time.perf_counter()
        count = count_triangulations(points)
        log(f"done in {time.perf_counter() - t:.1f}s")
    except (ValueError, RuntimeError) as e:
        sys.exit(f"error: {e}")

    print(count)


if __name__ == "__main__":
    main()
