#!/usr/bin/env bash
# =====================================================================
# restore_headers.sh — encode ' ' and '/' as '#' in FASTQ headers
# ---------------------------------------------------------------------
# Shared by trim.smk (rule restore_headers) and mirna.smk (rule
# mirna_restore_headers): both branches replay the exact same draft
# block (1_script l.345 / 2_script_miRNA l.34) on different inputs.
# Factored out once so the sed logic isn't duplicated (see
# docs/02_miRNA_script_uml_audit.md, point R3).
#
# Usage:
#   restore_headers.sh <in.fastq.gz> <out.fastq.gz>
# =====================================================================
set -euo pipefail

in=$1
out=$2

zcat "$in" | sed 's/ /#/g; s/\//#/g' | gzip > "$out"
