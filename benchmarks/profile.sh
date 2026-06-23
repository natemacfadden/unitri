#!/bin/sh
# wall time (high-resolution; min/mean over ITERS runs) and peak memory of a
# command.  auto-detects macOS (BSD time) vs GNU/Linux time.
#   usage:  [ITERS=N] ./profile.sh <command> [args...]      # stdin is forwarded
#   e.g.:   echo "8 8 8 8" | ITERS=20 ./profile.sh ./na-query
#
# stdin is captured once and replayed to every run.  the command runs ITERS+1
# times total (one extra run measures peak memory).

iters=${ITERS:-1}
in=$(mktemp 2>/dev/null  || mktemp -t profin)
tmp=$(mktemp 2>/dev/null || mktemp -t proftmp)
trap 'rm -f "$in" "$tmp"' EXIT
cat > "$in"

# one run for peak memory.  the command's stdout and stderr pass through (so you
# see its output, including any error), while time's own report is captured to
# $tmp: fd 3 carries the command's stderr to the terminal; time writes to $tmp.
exec 3>&2
if [ "$(uname)" = "Darwin" ]; then
    /usr/bin/time -l sh -c 'exec "$@" 2>&3' _ "$@" <"$in" 2>"$tmp"; status=$?
    mem=$(awk '/maximum resident set size/{printf "%.1f MB", $1/1048576}' "$tmp")
else
    /usr/bin/time -f '%M' sh -c 'exec "$@" 2>&3' _ "$@" <"$in" 2>"$tmp"; status=$?
    mem=$(awk 'END{printf "%.1f MB", $1/1024}' "$tmp")
fi
exec 3>&-

# timing: ITERS runs under a high-resolution clock; report min and mean per run
PROFIN="$in" perl -MTime::HiRes=time -e '
    my ($iters, @cmd) = @ARGV;
    open(my $out, ">&", \*STDOUT);                  # keep the real stdout
    open(STDOUT, ">", "/dev/null");                 # silence the command
    open(STDERR, ">", "/dev/null");
    my ($min, $sum) = (1e30, 0);
    for (1 .. $iters) {
        open(STDIN, "<", $ENV{PROFIN}) or die;      # fresh stdin each run
        my $s = time; system { $cmd[0] } @cmd; my $d = time - $s;
        $sum += $d; $min = $d if $d < $min;
    }
    printf {$out} "--- profile ---\niterations: %d\nmin /iter : %.4f s\nmean/iter : %.4f s\n",
        $iters, $min, $sum / $iters;
' "$iters" "$@"
echo "peak mem  : $mem"
if [ "$status" -ne 0 ]; then
    echo "WARNING: command exited with status $status (see its output above); these numbers are for a failed/partial run -- 137 = killed (often out of memory)"
fi
exit "$status"
