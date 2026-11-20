#!/bin/bash
#SBATCH --job-name=pcwde-paper-build
#SBATCH --account=2016456
#SBATCH --partition=general
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
set -euo pipefail
module load gcc/14.2.0-a75l texlive/20260301-lc7o
exec python3 -B "$(dirname "$0")/run_job.py" "$@"
