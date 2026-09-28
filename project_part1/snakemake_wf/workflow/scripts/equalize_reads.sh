#!/usr/bin/env bash
# =====================================================================
# equalize_reads.sh — P1 read equalization between R1 and R2
# ---------------------------------------------------------------------
# Reproduces the draft block verbatim: extract read IDs from the
# quality-filtered R1 and the header-encoded R2, keep the intersection
# (`comm -12`), re-attach the mate suffix, and subset both FASTQs with
# seqtk so R1 and R2 contain exactly the same reads in the same order.
#
# Usage:
#   equalize_reads.sh <R1_filtered.fq.gz> <R2_encoded.fq.gz> \
#                     <out_R1_subset.fq.gz> <out_R2_subset.fq.gz> <tmp_dir>
# =====================================================================
set -euo pipefail
export LC_ALL=C

r1_in=$1
r2_in=$2
r1_out=$3
r2_out=$4
tmp=$5

mkdir -p "$tmp"

# Read IDs (strip leading '@' and everything after the first '#')
zcat "$r1_in" | awk 'NR%4==1{split($0,a,"#"); print substr(a[1],2)}' | sort > "$tmp/R1_ids.txt"
zcat "$r2_in" | awk 'NR%4==1{split($0,a,"#"); print substr(a[1],2)}' | sort > "$tmp/R2_ids.txt"

comm -12 "$tmp/R1_ids.txt" "$tmp/R2_ids.txt" > "$tmp/common_ids.txt"

sed 's/$/#1:N:0:1/' "$tmp/common_ids.txt" > "$tmp/R1_final_ids.txt"
sed 's/$/#2:N:0:1/' "$tmp/common_ids.txt" > "$tmp/R2_final_ids.txt"

seqtk subseq "$r1_in" "$tmp/R1_final_ids.txt" | gzip > "$r1_out"
seqtk subseq "$r2_in" "$tmp/R2_final_ids.txt" | gzip > "$r2_out"
