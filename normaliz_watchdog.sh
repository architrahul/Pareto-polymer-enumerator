#!/bin/bash
# normaliz_watchdog.sh
# Kills any `normaliz` process that has been running longer than $THRESHOLD
# seconds. Designed to be launched alongside a long-running benchmark sweep
# so a single long Hilbert basis computation can't hold up the whole
# run.
#
# When a normaliz child is killed, the parent (hilbert_pipeline.py) sees the
# child exit, calls cleanup_normaliz_files() before the next block (so no
# stale eqs.out is mis-read), and continues. The block contributes 0 vectors
# to the union.

THRESHOLD=${1:-10800}   # default 3 hours in seconds
LOG=${LOG:-/tmp/normaliz_watchdog.log}
INTERVAL=${INTERVAL:-60}

echo "[$(date)] watchdog started: threshold=${THRESHOLD}s, interval=${INTERVAL}s, log=$LOG" | tee -a "$LOG"

while true; do
    # ps -eo pid=,etime=,comm= prints PID, elapsed, command name.
    # etime format: [[DD-]HH:]MM:SS
    ps -eo pid=,etime=,comm= 2>/dev/null \
      | awk -v thr="$THRESHOLD" '
          $NF == "normaliz" {
              # Parse etime in $2 — split on - and :
              n = split($2, a, /[-:]/)
              if      (n == 2) secs = a[1]*60 + a[2]
              else if (n == 3) secs = a[1]*3600 + a[2]*60 + a[3]
              else if (n == 4) secs = a[1]*86400 + a[2]*3600 + a[3]*60 + a[4]
              else             secs = 0
              if (secs > thr) print $1 " " secs
          }' \
      | while read -r pid secs; do
            echo "[$(date)] killing normaliz PID=$pid elapsed=${secs}s (> ${THRESHOLD}s)" | tee -a "$LOG"
            kill -TERM "$pid" 2>/dev/null
            sleep 5
            # SIGKILL if still alive
            if kill -0 "$pid" 2>/dev/null; then
                echo "[$(date)] PID=$pid did not exit on SIGTERM, sending SIGKILL" | tee -a "$LOG"
                kill -KILL "$pid" 2>/dev/null
            fi
        done

    sleep "$INTERVAL"
done
