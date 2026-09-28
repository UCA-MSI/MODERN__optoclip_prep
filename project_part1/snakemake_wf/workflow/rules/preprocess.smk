# =====================================================================
# preprocess.smk — P1 : QC-oriented preprocessing on the whole library
# =====================================================================
# Operates on the single multiplexed library (no {sample} wildcard yet):
# header encoding, barcode-frequency diagnostics, barcode-region quality
# filtering, and R1/R2 read equalization. Byte-exact reproduction of the
# draft (sed/awk/comm kept verbatim to preserve results strictly).


rule encode_headers_r2:
    """Encode spaces and slashes in R2 headers as '#' (draft l.18)."""
    input:
        config["raw_reads"]["r2"],
    output:
        "results/preprocess/R2_encoded.fastq.gz",
    log:
        "results/logs/preprocess/encode_headers_r2.log",
    conda:
        "../envs/base.yaml"
    shell:
        r"zcat {input} | sed 's/ /#/g; s/\//#/g' | gzip > {output} 2> {log}"


rule barcode_frequencies:
    """Detected experimental barcodes and their frequencies (diagnostic)."""
    input:
        config["raw_reads"]["r1"],
    output:
        "results/preprocess/exp_barcodes_R1.detected",
    params:
        start=config["read_structure"]["umi1_len"] + 1,
        blen=config["read_structure"]["barcode_len"],
    log:
        "results/logs/preprocess/barcode_frequencies.log",
    conda:
        "../envs/base.yaml"
    shell:
        "zcat {input} | "
        "awk -v s={params.start} -v l={params.blen} "
        "'{{ if (FNR%4==2) print substr($1,s,l) }}' | "
        "sort | uniq -c | sort -k1,1rn > {output} 2> {log}"


rule qualfilter_barcode_region:
    """Quality-filter the R1 5' barcode region and subset R1 (Option A/B).

    barcode_qc_length = 15 (Option B, default) or 6 (Option A). Both draft
    variants wrote the same file; this single parametrized rule replaces
    them and makes the choice explicit via config.
    """
    input:
        config["raw_reads"]["r1"],
    output:
        fastq="results/preprocess/R1_filtered.fastq.gz",
        ids=temp("results/preprocess/R1_qualFilteredIDs.list"),
    params:
        qc_len=config["preprocess"]["barcode_qc_length"],
        q=config["preprocess"]["quality_threshold"],
        p=config["preprocess"]["quality_percent"],
    log:
        "results/logs/preprocess/qualfilter_barcode_region.log",
    conda:
        "../envs/fastx.yaml"
    shell:
        r"""
        ( zcat {input} \
            | fastx_trimmer -f 1 -l {params.qc_len} \
            | fastq_quality_filter -q {params.q} -p {params.p} \
            | awk 'FNR%4==1 {{ print substr($1,2) }}' > {output.ids}
          seqtk subseq {input} {output.ids} \
            | sed 's/ /#/g; s/\//#/g' | gzip > {output.fastq} ) 2> {log}
        """


rule equalize_reads:
    """Keep reads common to filtered R1 and encoded R2 (draft comm block)."""
    input:
        r1="results/preprocess/R1_filtered.fastq.gz",
        r2="results/preprocess/R2_encoded.fastq.gz",
        script="workflow/scripts/equalize_reads.sh",
    output:
        r1="results/preprocess/R1_subset.fastq.gz",
        r2="results/preprocess/R2_subset.fastq.gz",
    params:
        tmp="results/preprocess/tmp_equalize",
    log:
        "results/logs/preprocess/equalize_reads.log",
    conda:
        "../envs/fastx.yaml"
    shell:
        "bash {input.script} {input.r1} {input.r2} "
        "{output.r1} {output.r2} {params.tmp} > {log} 2>&1"
