#!/usr/bin/env bash
#SBATCH --job-name=synteny-coord-batched
#SBATCH --account=2016456
#SBATCH --partition=general
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=recover_modelless_batched.%j.out
set -euo pipefail
cd "${MODELLESS_RUN_ROOT:?}"
test "$(git --git-dir="$HOME/repos/cynipoid-pcwde-evolution.git" rev-parse main)" = 6d7ae3d394e5bc313ab1ac45d7248766734ac348
export LD_LIBRARY_PATH=${MODELLESS_TOOL_BIN:?}/../lib:${LD_LIBRARY_PATH:-}
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
${MODELLESS_TOOL_BIN:?}/../bin/python -B "${MODELLESS_CODE_ROOT:?}/recover_modelless_batched.py"
