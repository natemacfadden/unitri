#!/bin/sh
# wall time and peak memory of a single run; auto-detects macOS (BSD time)
# vs GNU/Linux time.
# usage:  ./profile.sh <command> [args...]      # stdin is forwarded
# e.g.:   printf '' | ./profile.sh ./na-query 0

tmp=$(mktemp 2>/dev/null || mktemp -t profile)

if [ "$(uname)" = "Darwin" ]; then
    # BSD time: -l prints "<sec> real ..." and max RSS in bytes
    /usr/bin/time -l "$@" 2>"$tmp"
    secs=$(awk '/ real /{print $1}' "$tmp" | tail -1)
    mb=$(awk '/maximum resident set size/{printf "%.1f", $1/1048576}' "$tmp")
else
    # GNU time: -f gives elapsed seconds (%e) and max RSS in kbytes (%M)
    /usr/bin/time -f '%e %M' "$@" 2>"$tmp"
    secs=$(awk 'END{print $1}' "$tmp")
    mb=$(awk 'END{printf "%.1f", $2/1024}' "$tmp")
fi

echo "--- profile ---"
echo "wall time : $secs s"
echo "peak mem  : $mb MB"
rm -f "$tmp"
