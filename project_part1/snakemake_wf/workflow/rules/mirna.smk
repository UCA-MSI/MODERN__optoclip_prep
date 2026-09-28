# =====================================================================
# mirna.smk — miRNA profiling with SPORTS1.1 (migration of
#             draft_code/raw_bash/2_script_miRNA.txt)
# =====================================================================
# Independent post-processing branch of the same iCLIP preprocessing:
# reuses the P2b length-filtered intermediate (results/trim/lenfilter/),
# re-selects miRNA-sized reads (18-25 nt), re-annotates them against
# miRBase with SPORTS1.1 (Shi et al. 2018), and produces an
# RPM-normalized, multi-sample expression matrix.
#
# Full audit (steps, draft line numbers, open questions) is in
# docs/02_miRNA_script_uml_audit.md §2. Three deliberate deviations from
# a literal replay of the draft, all documented at the point they apply:
#   R2  size-selection now depends explicitly on the trim/lenfilter
#       temp() files (declared input, not a detached post-hoc script).
#   R3  header restoration reuses scripts/restore_headers.sh instead of
#       re-implementing the same sed pipeline (shared with trim.smk).
#   R4  the R2 merge branch is dropped: sports.pl runs single-end (-s)
#       and the draft never reads the merged R2 files back.
#   R1 (reference chain) and R5 (gene-count extraction) required an
#   explicit modeling choice where the draft was silent/ambiguous — see
#   the rule docstrings below for the reasoning; both are flagged
#   "ASSUMPTION" and should be checked against a real SPORTS run before
#   the numbers are trusted.
#
# One rule per step, wildcarded by {sample} (barcode unit, reused from
# common.smk) or {group} (biological sample) — every step can be run in
# isolation, e.g. `snakemake results/mirna/sports/1_R1`.


# ---- Reference preparation (independent of any sample) --------------


rule hsa_mature_species_filter:
    """Keep only the target species' entries from miRBase's mature.fa (draft l.4-5)."""
    input:
        config["mirna"]["mirbase_mature_fa"],
    output:
        "results/mirna/ref/hsa_mature.fa",
    params:
        species=config["mirna"]["species_tag"],
    log:
        "results/logs/mirna/ref/species_filter.log",
    conda:
        "../envs/base.yaml"
    shell:
        """grep -A 1 "{params.species}" {input} | grep -v "^--$" > {output} 2> {log}"""


rule hsa_mature_u2t:
    """RNA -> DNA alphabet (U->T) so the sequences can seed a Bowtie1 index (draft l.7)."""
    input:
        "results/mirna/ref/hsa_mature.fa",
    output:
        "results/mirna/ref/hsa_mature_final.fa",
    log:
        "results/logs/mirna/ref/u2t.log",
    conda:
        "../envs/base.yaml"
    shell:
        "sed 's/U/T/g' {input} > {output} 2> {log}"


rule mirna_bowtie_index:
    """Build the Bowtie1 miRBase index consumed by sports.pl -m (draft l.86-88).

    ASSUMPTION (audit point R1): the draft builds `miRBase_21-hsa` from a
    `miRBase_21-hsa.fa` obtained "separately", while `hsa_mature_final.fa`
    (this rule's real input) is produced just above and never referenced
    again. The U->T conversion only makes sense as prep for a Bowtie
    index, so the two are almost certainly the same file under two names
    (renamed/copied before `bowtie-build` on the original host). This
    rule makes that link explicit instead of leaving it implicit.
    """
    input:
        "results/mirna/ref/hsa_mature_final.fa",
    output:
        directory(config["reference"]["mirbase_bowtie_index"]),
    params:
        prefix=lambda w, output: f"{output[0]}/miRBase_hsa",
    threads: threads_for("bowtie_index_threads", 4)
    log:
        "results/logs/mirna/ref/mirna_bowtie_index.log",
    conda:
        "../envs/bowtie1.yaml"
    shell:
        "mkdir -p {output} && bowtie-build {input} {params.prefix} > {log} 2>&1"


rule mirna_genome_bowtie_index:
    """Build the Bowtie1 genome index consumed by sports.pl -g (draft l.135, `-g .../hg38/genome`).

    Separate from `star_index` (mapping.smk) and `novoalign_index`: SPORTS
    hard-requires a Bowtie1 index, distinct from STAR's and Novoalign's.
    """
    input:
        config["reference"]["genome_fasta"],
    output:
        directory(config["reference"]["mirna_genome_bowtie_index"]),
    params:
        prefix=lambda w, output: f"{output[0]}/genome",
    threads: threads_for("bowtie_index_threads", 4)
    log:
        "results/logs/mirna/ref/genome_bowtie_index.log",
    conda:
        "../envs/bowtie1.yaml"
    shell:
        "mkdir -p {output} && bowtie-build {input} {params.prefix} > {log} 2>&1"


# ---- Per barcode-unit: miRNA size selection + header restore --------


rule mirna_size_select:
    """Keep reads of miRNA length (18-25 nt) from the length-filtered pair (draft l.13).

    Consumes results/trim/lenfilter/{sample}_{read}.fastq.gz, the same
    temp() intermediate that trim.smk's restore_headers reads (audit
    point R2): declaring it here as an explicit rule input, rather than
    running this as a detached post-hoc script, tells Snakemake there
    are two consumers, so the temp file survives until both have run.
    """
    input:
        r1="results/trim/lenfilter/{sample}_1.fastq.gz",
        r2="results/trim/lenfilter/{sample}_2.fastq.gz",
    output:
        r1=temp("results/mirna/sizesel/{sample}_1.fastq.gz"),
        r2=temp("results/mirna/sizesel/{sample}_2.fastq.gz"),
    params:
        min_len=config["mirna"]["min_length"],
        max_len=config["mirna"]["max_length"],
    threads: threads_for("cutadapt_threads", 4)
    log:
        "results/logs/mirna/size_select/{sample}.log",
    conda:
        "../envs/cutadapt.yaml"
    shell:
        "cutadapt -j {threads} --pair-filter=any "
        "--minimum-length {params.min_len} --maximum-length {params.max_len} "
        "-o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log} 2>&1"


rule mirna_restore_headers:
    """Header ' '/'/'-> '#' on the size-selected FASTQ (draft l.34).

    Only R1 is restored: sports.pl runs single-end (-s) and the draft's
    merged R2 (audit point R4) is never read again downstream, so the R2
    branch is dropped here to save compute/disk.
    """
    input:
        fastq="results/mirna/sizesel/{sample}_1.fastq.gz",
        script="workflow/scripts/restore_headers.sh",
    output:
        temp("results/mirna/final/{sample}_1.fastq.gz"),
    log:
        "results/logs/mirna/restore_headers/{sample}.log",
    conda:
        "../envs/base.yaml"
    shell:
        "bash {input.script} {input.fastq} {output} 2> {log}"


# ---- Per biological sample: merge barcodes, annotate, count ---------


rule mirna_merge_group:
    """Merge a biological sample's barcode-unit R1 files (draft l.95-104).

    The draft `cat`s barcode1+barcode2 for every group except sample 11,
    which is `cp`'d from a single file. samples.tsv shows sample 11 has
    two demultiplexed barcode units (11_1, 11_2) exactly like every other
    group, and 1_script merges both of them before mapping — so the
    draft's single-file shortcut for sample 11 in this script looks like
    an inconsistency rather than a deliberate choice (see
    docs/02_miRNA_script_uml_audit.md §2.4). This rule merges every
    group uniformly via `cat` (a no-op concatenation when there is only
    one input), which both matches how sample 11 is treated everywhere
    else in the workflow and removes the need for a special case.
    """
    input:
        lambda w: [f"results/mirna/final/{u}_1.fastq.gz" for u in units_of_group(w.group)],
    output:
        "results/mirna/merged/{group}_R1.fastq.gz",
    log:
        "results/logs/mirna/merge_group/{group}.log",
    conda:
        "../envs/base.yaml"
    shell:
        "cat {input} > {output} 2> {log}"


rule mirna_gunzip:
    """Decompress the merged FASTQ: sports.pl expects plain FASTQ input (draft l.119-130)."""
    input:
        "results/mirna/merged/{group}_R1.fastq.gz",
    output:
        temp("results/mirna/merged/{group}_R1.fastq"),
    log:
        "results/logs/mirna/gunzip/{group}.log",
    conda:
        "../envs/base.yaml"
    shell:
        "gunzip -c {input} > {output} 2> {log}"


rule sports_annotate:
    """Annotate a sample's miRNA-sized reads with SPORTS1.1 (draft l.135-144).

    `-s` = single-end + antisense counts, `-k` = keep intermediates,
    `-z` = gzip SPORTS' own outputs (draft comment: "KEEP THIS ONE
    BECAUSE GIVES ALSO ANTISENSE"). SPORTS itself must be on PATH (see
    envs/sports.yaml) — it is not conda-packaged.
    """
    input:
        fastq="results/mirna/merged/{group}_R1.fastq",
        genome_index=config["reference"]["mirna_genome_bowtie_index"],
        mirna_index=config["reference"]["mirbase_bowtie_index"],
    output:
        directory("results/mirna/sports/{group}_R1"),
    params:
        genome_prefix=f"{config['reference']['mirna_genome_bowtie_index']}/genome",
        mirna_prefix=f"{config['reference']['mirbase_bowtie_index']}/miRBase_hsa",
    threads: threads_for("sports_threads", 4)
    log:
        "results/logs/mirna/sports/{group}.log",
    conda:
        "../envs/sports.yaml"
    shell:
        "mkdir -p {output} && "
        "sports.pl -i {input.fastq} -s -p {threads} "
        "-g {params.genome_prefix} -m {params.mirna_prefix} "
        "-o {output} -k -z > {log} 2>&1"


rule sports_extract_counts:
    """Extract per-miRNA read counts from SPORTS' native output (draft: implicit before l.154).

    ASSUMPTION (audit point R5): `gene_counts_sample{N}.txt` is consumed
    by the draft's normalization step but never produced by it — SPORTS
    writes several files per RNA class under tool-specific names, and the
    extraction/rename was done by hand outside the script. This rule
    reconstructs it explicitly and traceably from SPORTS' own documented
    6-column `<seqfile>_output.txt` table (columns: ID, Sequence, Length,
    Reads, Match_Genome, Annotation — see the SPORTS1.1 README,
    "annotation.pl" section): it sums `Reads` per `Annotation` for every
    row whose Annotation matches `mirna.annotation_regex` (default
    `^hsa-`, i.e. a mature miRBase hsa- name). Validate this against a
    real SPORTS run before trusting downstream counts.
    """
    input:
        "results/mirna/sports/{group}_R1",
    output:
        "results/mirna/counts/gene_counts_{group}.txt",
    params:
        table=lambda w, input: f"{input[0]}/1_{w.group}_R1/{w.group}_R1_result/{w.group}_R1_output.txt",
        regex=config["mirna"]["annotation_regex"],
    log:
        "results/logs/mirna/extract_counts/{group}.log",
    conda:
        "../envs/base.yaml"
    shell:
        r"""
        awk -F'\t' -v OFS='\t' -v re="{params.regex}" '
          $6 ~ re {{ sum[$6] += $4 }}
          END {{ print "Geneid", "Count"; for (g in sum) print g, sum[g] }}
        ' {params.table} > {output} 2> {log}
        """


rule mirna_normalize_rpm:
    """Reads-per-million normalization of a sample's miRNA counts (draft l.150-182)."""
    input:
        "results/mirna/counts/gene_counts_{group}.txt",
    output:
        "results/mirna/counts/normalized_gene_counts_{group}.txt",
    log:
        "results/logs/mirna/normalize/{group}.log",
    conda:
        "../envs/base.yaml"
    shell:
        r"""
        ( total=$(awk -F'\t' '{{ if (NR>1) sum+=$2 }} END {{ print sum }}' {input})
          awk -F'\t' -v total="$total" 'BEGIN {{ OFS="\t" }}
            NR==1 {{ print $1, $2, "Normalized_Count_{wildcards.group}" }}
            NR>1  {{ printf "%s\t%s\t%.6f\n", $1, $2, ($2/total)*1000000 }}' {input} \
          | sed 's/,/./g' > {output} ) 2> {log}
        """


rule merge_mirna_counts:
    """Outer-join every sample's normalized counts into one matrix (draft l.186-200)."""
    input:
        counts=expand("results/mirna/counts/normalized_gene_counts_{group}.txt", group=GROUPS),
    output:
        "results/mirna/merged_counts_SPORTS_miRNAs.csv",
    log:
        "results/logs/mirna/merge_counts.log",
    conda:
        "../envs/r-base.yaml"
    script:
        "../scripts/merge_mirna_counts.R"
