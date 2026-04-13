#!/bin/bash
#SBATCH --job-name=because_benchmark
#SBATCH --output=logs/bench_%a.out
#SBATCH --error=logs/bench_%a.err
#SBATCH --array=1-5
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.com

# --- Slurm Array Benchmarking Suite ---
# Task 1: because (JAGS)
# Task 2: because (NIMBLE)
# Task 3: brms (Spatial + Phylo PGLMM)
# Task 4: MCMCglmm (Sequential PGLMM)
# Task 5: Classical Baselines (phylopath + GLMs)

echo "Starting Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Running engine $(hostname)"

# Create logs directory if it doesn't exist
mkdir -p logs

# Load R module (uncomment/edit based on your cluster)
# module load R/4.2.0-foss-2022a

# Execute the task-aware R script
Rscript drafts/comparison_server.R ${SLURM_ARRAY_TASK_ID}

echo "Task ${SLURM_ARRAY_TASK_ID} complete."
