# Full competitive dbCAN-sub inputs

The 64 Figure 1 genomes use the scans pinned in `config/upstream_locked_inputs.tsv` and the source/proteome/scan hash mapping in the producer's `coverage.tsv`. The 48 frozen cynipoid scans were produced in `source_inputs/gh1e72_cynipoid_fill/scan_outputs`; the 15 external and one internal reference scan were produced in `analysis/figure6_uniform_gh1e72_scan/full_dbcan_outgroups/`. All use dbCAN 5.2.4, dbCAN-sub HMM SHA256 `fc431183fff93bd6a9e0a346dff1118f6bb9ecba02fb881b8c03a57d08ca5fe8`, with the competitive command:

```
run_dbcan CAZyme_annotation --input_raw_data INPUT.faa --mode protein --db_dir ${DBCAN_DB_DIR} --output_dir OUTPUT --methods dbCANsub --threads 8 --e_value_threshold_dbsub 1e-15 --coverage_threshold_dbsub 0.35 --log-file OUTPUT/run.log --log-level INFO
```

The producer's `inventory.py` verifies each current proteome SHA256 against the recorded scan input and checks the scan output SHA256. The initial inventory covered 64/64 inherited scans. Subsequent source reconciliation showed that AE, AQP and DL required their paper proteomes rather than the inherited scan inputs. These three were rescanned with the same command, database and thresholds; their paper source, processed input, scan result and validation hashes are pinned in `config/cazy_paper_proteomes.tsv`, `config/cazy_processed_proteomes.tsv` and `config/cazy_scans.tsv`.

The paper deposit now carries the 64 scan tables, 64 required proteomes, 64 processed scan proteomes and 17 mapping sources under `repository_inputs/cazy_inputs/`. Shared paths are stored once, for 148 checksum-pinned files. `config/logical_inputs.tsv` registers every file. The CAZy producers resolve these paths below `PCWDE_DATA_ROOT` and reject paths that escape it, including symlinks. The accepted-member table and sequence FASTAs remain in their existing deposit locations.

The three saved GH1_e72 tables in `repository_inputs/figures/presence_counts/` serve only as frozen comparison controls. The scan-count producer checks all 64 counts against them and writes `cazy_counts/derived_gh1e72_gene_counts.tsv`. Figure 1 and Figure S2 read that derived table after the producer has run in the same manuscript build.
