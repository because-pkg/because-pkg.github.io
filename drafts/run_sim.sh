#!/bin/bash
# run_sim.sh: Robust runner for the coverage simulation

cd /Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because

while true; do
  echo "--- Starting Simulation @ $(date) ---"
  /usr/local/bin/Rscript drafts/sim_coverage_butterfly.R
  
  # Check exit code
  if [ $? -eq 0 ]; then
    echo "--- Simulation Completed Successfully @ $(date) ---"
    break
  else
    echo "--- Simulation Crashed / Stopped @ $(date). Restarting in 5s... ---"
    sleep 5
  fi
done
