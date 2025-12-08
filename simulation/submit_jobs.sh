#!/bin/bash

# Define arrays for options and methods
smoothings=("clusterbased" "adjacencybased" "distancebased") 
outcomemods=("linear" "linearinteraction" "nonlinearinteraction") 

# Loop through each combination and submit jobs
for smoothing in "${smoothings[@]}"; do
  for outcomemod in "${outcomemods[@]}"; do
  sbatch --job-name=${smoothing}_${outcomemod}\
    --output=output/${smoothing}_${outcomemod}.out \
    --error=error/${smoothing}_${outcomemod}.err \
    run_job.sh 1000 "$smoothing" "$outcomemod"  
  sleep 1 # pause to be kind to the scheduler
  done
done

