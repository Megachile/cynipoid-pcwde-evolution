#!/usr/bin/env bash
#SBATCH --job-name=pcwde_sel
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-2%3
#SBATCH --output=logs/selection_%A_%a.out
#SBATCH --error=logs/selection_%A_%a.err

# Submit from the explicit project directory after creating its logs directory:
#   cd "$PCWDE_SELECTION_PROJECT_DIR" && mkdir -p logs && sbatch run_selection.sbatch
# Required environment: PCWDE_SELECTION_PROJECT_DIR, PCWDE_SELECTION_WORK,
# PCWDE_IQTREE, and PCWDE_HYPHY.
set -euo pipefail
PROJECT_DIR=${PCWDE_SELECTION_PROJECT_DIR:?PCWDE_SELECTION_PROJECT_DIR is required}
WORK=${PCWDE_SELECTION_WORK:?PCWDE_SELECTION_WORK is required}
IQTREE=${PCWDE_IQTREE:?PCWDE_IQTREE is required}
HYPHY=${PCWDE_HYPHY:?PCWDE_HYPHY is required}
TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
CPUS=${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK is required}
test "$CPUS" -eq 16
test -d "$PROJECT_DIR"
test -d "$WORK"
test -x "$IQTREE"
test -x "$HYPHY"
cd "$PROJECT_DIR"

FAMILIES=(GH5_2 PL1_8 PL4_4)
case "$TASK_ID" in
  0|1|2) FAMILY=${FAMILIES[$TASK_ID]} ;;
  *) printf 'invalid SLURM_ARRAY_TASK_ID: %s\n' "$TASK_ID" >&2; exit 64 ;;
esac
ALN="$WORK/${FAMILY}_current_codon_aligned.fna"
OUT="$WORK/outputs/$FAMILY"
mkdir -p "$OUT" "$WORK/logs"

test -s "$ALN"
grep -q '^>' "$ALN"
"$IQTREE" --version 2>&1 | grep -Eq '(^|[^0-9])3\.0\.1([^0-9]|$)'
"$HYPHY" --version 2>&1 | grep -Eq '(^|[^0-9])2\.5\.85([^0-9]|$)'
for result in "$OUT/${FAMILY}_fubar.json" "$OUT/${FAMILY}_absrel.json"; do
  test ! -e "$result"
done

"$IQTREE" \
  -s "$ALN" \
  --prefix "$OUT/${FAMILY}_iqtree" \
  -m GTR+G4 \
  -B 1000 \
  -T "$CPUS" \
  --safe

TREE="$OUT/${FAMILY}_iqtree.treefile"
test -s "$TREE"

"$HYPHY" fubar \
  --alignment "$ALN" \
  --tree "$TREE" \
  --branches All \
  --output "$OUT/${FAMILY}_fubar.json" \
  CPU="$CPUS" \
  2>&1 | tee "$OUT/${FAMILY}_fubar.log"

"$HYPHY" absrel \
  --alignment "$ALN" \
  --tree "$TREE" \
  --branches All \
  --output "$OUT/${FAMILY}_absrel.json" \
  CPU="$CPUS" \
  2>&1 | tee "$OUT/${FAMILY}_absrel.log"

for result in "$OUT/${FAMILY}_fubar.json" "$OUT/${FAMILY}_absrel.json"; do
  test -s "$result"
  grep -q '"runtime":"2.5.85"' "$result"
done
