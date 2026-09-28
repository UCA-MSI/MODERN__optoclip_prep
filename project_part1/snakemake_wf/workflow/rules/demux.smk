# =====================================================================
# demux.smk — P2a : demultiplexing with Flexbar
# =====================================================================
# The barcode FASTA is generated from the sample sheet so demux stays
# fully config-driven. Flexbar produces all per-unit files in one run;
# because the sample set is known a priori we enumerate the outputs
# (no checkpoint needed) and downstream rules consume them per {sample}.


rule barcodes_fasta:
    """Build the flexbar barcode FASTA from the sample sheet."""
    input:
        config["samples"],
    output:
        "resources/barcodes_antisense.fasta",
    run:
        with open(output[0], "w") as fh:
            for sid in UNITS:
                fh.write(f">{sid}\n{samples.loc[sid, 'flexbar_barcode']}\n")


rule flexbar_demux:
    """Demultiplex the equalized library into per-barcode paired FASTQs."""
    input:
        r1="results/preprocess/R1_subset.fastq.gz",
        r2="results/preprocess/R2_subset.fastq.gz",
        barcodes="resources/barcodes_antisense.fasta",
    output:
        r1=expand("results/demux/flexbarOut_barcode_{sample}_1.fastq.gz", sample=UNITS),
        r2=expand("results/demux/flexbarOut_barcode_{sample}_2.fastq.gz", sample=UNITS),
        unassigned=expand(
            "results/demux/flexbarOut_barcode_unassigned_{read}.fastq.gz",
            read=["1", "2"],
        ),
    params:
        target="results/demux/flexbarOut",
        trim_end=config["demux"]["barcode_trim_end"],
        error_rate=config["demux"]["barcode_error_rate"],
        min_len=config["demux"]["min_read_length"],
        umi=lambda w: "--umi-tags" if config["demux"].get("umi_tags", True) else "",
    threads: threads_for("align_threads", 8)
    log:
        "results/logs/demux/flexbar.log",
    conda:
        "../envs/flexbar.yaml"
    shell:
        "flexbar -r {input.r1} -p {input.r2} --target {params.target} "
        "--zip-output GZ --barcodes {input.barcodes} --barcode-unassigned "
        "--barcode-trim-end {params.trim_end} "
        "--barcode-error-rate {params.error_rate} "
        "--min-read-length {params.min_len} {params.umi} "
        "--threads {threads} > {log} 2>&1"
