#!/usr/bin/env bash
# stop a running DRAM panel job: the wrapper, DRAM.py and their children
for pat in 'run_dram_panel' 'DRAM.py annotate'; do
  for pid in $(pgrep -f "$pat"); do
    [ "$pid" = "$$" ] && continue
    pkill -P "$pid" 2>/dev/null; kill "$pid" 2>/dev/null
  done
done
sleep 2
pgrep -fa 'dram_out' | grep -v _stop_dram | awk '{print $1}' | xargs -r kill 2>/dev/null
echo stopped
