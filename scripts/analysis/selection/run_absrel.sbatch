#!/usr/bin/env bash
#SBATCH --job-name=pcwde_absrel
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-2%3
#SBATCH --output=logs/absrel_parallel_%A_%a.out
#SBATCH --error=logs/absrel_parallel_%A_%a.err

# Recovery-only route. It writes a separately named candidate and never
# overwrites the combined producer's accepted <family>_absrel.json. Promotion
# requires an external byte/hash receipt.
# Required environment: PCWDE_SELECTION_PROJECT_DIR, PCWDE_SELECTION_WORK,
# and PCWDE_HYPHY.
set -euo pipefail

PROJECT_DIR=${PCWDE_SELECTION_PROJECT_DIR:?PCWDE_SELECTION_PROJECT_DIR is required}
WORK=${PCWDE_SELECTION_WORK:?PCWDE_SELECTION_WORK is required}
HYPHY=${PCWDE_HYPHY:?PCWDE_HYPHY is required}
TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
CPUS=${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK is required}
test "$CPUS" -eq 16
test -d "$PROJECT_DIR"
test -d "$WORK"
test -x "$HYPHY"
cd "$PROJECT_DIR"

FAMILIES=(GH5_2 PL1_8 PL4_4)
case "$TASK_ID" in
  0|1|2) FAMILY=${FAMILIES[$TASK_ID]} ;;
  *) printf 'invalid SLURM_ARRAY_TASK_ID: %s\n' "$TASK_ID" >&2; exit 64 ;;
esac
ALN="$WORK/${FAMILY}_current_codon_aligned.fna"
TREE="$WORK/outputs/${FAMILY}/${FAMILY}_iqtree.treefile"
OUT="$WORK/outputs/${FAMILY}/${FAMILY}_absrel_parallel.json"
LOG="$WORK/outputs/${FAMILY}/${FAMILY}_absrel_parallel.log"

test -s "$ALN"
test -s "$TREE"
test ! -e "$OUT"
test ! -e "$LOG"
"$HYPHY" --version 2>&1 | grep -Eq '(^|[^0-9])2\.5\.85([^0-9]|$)'

"$HYPHY" absrel \
  --alignment "$ALN" \
  --tree "$TREE" \
  --branches All \
  --output "$OUT" \
  CPU="$CPUS" \
  2>&1 | tee "$LOG"

test -s "$OUT"
grep -q '"runtime":"2.5.85"' "$OUT"
