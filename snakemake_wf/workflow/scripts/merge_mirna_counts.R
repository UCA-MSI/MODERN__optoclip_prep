#!/usr/bin/env Rscript
# =====================================================================
# merge_mirna_counts.R — cross-sample miRNA RPM matrix (draft l.186-200)
# ---------------------------------------------------------------------
# Reproduces the draft's R block: read every per-group normalized count
# table and outer-join them on "Geneid" into one matrix, written as CSV.
# Driven by the Snakemake `script:` directive (snakemake@ S4 object).
# =====================================================================

counts_files <- snakemake@input[["counts"]]
csv_file <- snakemake@output[[1]]

tables <- lapply(counts_files, function(f) read.delim(f, check.names = FALSE))

merged_counts_SPORTS_miRNAs <- Reduce(
  function(...) merge(..., all = TRUE, by = "Geneid"),
  tables
)

write.csv(merged_counts_SPORTS_miRNAs, csv_file)
