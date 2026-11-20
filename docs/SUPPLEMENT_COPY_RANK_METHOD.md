# Condensed copy-rank table

`scripts/tables/build_cazy_rank_table.py` reads the same-build `paired_ranks.tsv` and `primary_ecami_ranks.tsv` plus the committed 64-genome `config/cazy_lineage_groups.tsv`. The table uses only the `largest_group` primary eCAMI route. The producer checks the 192 genome-family pairs, all 576 route rows and the accepted-versus-scan sensitivity counts before summarizing.

Each of the six lineage groups has one row per GH5_2, PL1_8 and PL4_4. The four named main-text groups come directly from the lineage authority. Remaining Figitidae form the fifth group; other source families form the non-cynipoid Hymenoptera/Diptera sixth group. For each family, the table reports the number of genomes in the group, the number with a nonempty primary rank, the median integer rank, median midrank percentile and min–max integer rank among present genomes. An absent family has zero present genomes and blank rank statistics. No new cohort, threshold or rank reference is introduced.

The full 192-row per-genome table, including pooled, accepted-member and mixed-unit ranks, is written to `cazy_rank_table/per_genome_cazy_ranks.tsv` and copied to the data deposit. `build_numbers.py` reads the condensed reader table for main-text lineage values. The exact build records input and output SHA256.
