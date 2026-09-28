# =====================================================================
# qc.smk — quality control (FastQC per checkpoint, MultiQC aggregation)
# =====================================================================
# A single generic FastQC rule serves every checkpoint. Each checkpoint
# has a short label ({qcid}); the actual FASTQ is looked up in
# FASTQC_INPUTS, and FastQC's output (named after the input) is renamed
# to the label so outputs never depend on input path depth. This removes
# the ~120 duplicated fastqc calls of the draft.

# ---- QC checkpoints: label -> FASTQ path ----------------------------
FASTQC_INPUTS = {
    "raw_R1": config["raw_reads"]["r1"],
    "raw_R2": config["raw_reads"]["r2"],
    "subset_R1": "results/preprocess/R1_subset.fastq.gz",
    "subset_R2": "results/preprocess/R2_subset.fastq.gz",
}
for _u in UNITS:
    for _r in ("1", "2"):
        FASTQC_INPUTS[f"final_{_u}_{_r}"] = f"results/trim/final/{_u}_{_r}.fastq.gz"


# FastQC on one checkpoint FASTQ; output renamed to the {qcid} label.
rule fastqc:
    wildcard_constraints:
        qcid=r"[A-Za-z0-9_]+",
    input:
        lambda w: FASTQC_INPUTS[w.qcid],
    output:
        html="results/qc/fastqc/{qcid}_fastqc.html",
        zip="results/qc/fastqc/{qcid}_fastqc.zip",
    params:
        outdir="results/qc/fastqc",
        stem=lambda w, input: os.path.basename(input[0])[: -len(".fastq.gz")],
    threads: threads_for("fastqc_threads", 1)
    log:
        "results/logs/fastqc/{qcid}.log",
    conda:
        "../envs/fastqc.yaml"
    shell:
        "mkdir -p {params.outdir} && "
        "fastqc --nogroup --threads {threads} --outdir {params.outdir} {input} > {log} 2>&1 && "
        "mv {params.outdir}/{params.stem}_fastqc.html {output.html} && "
        "mv {params.outdir}/{params.stem}_fastqc.zip {output.zip}"


def fastqc_targets():
    """All FastQC checkpoints (feed MultiQC)."""
    return [f"results/qc/fastqc/{qcid}_fastqc.zip" for qcid in FASTQC_INPUTS]


rule multiqc:
    """Aggregate all FastQC checkpoints into one report."""
    input:
        fastqc_targets(),
    output:
        "results/qc/multiqc/multiqc_report.html",
    params:
        indir="results/qc/fastqc",
        outdir="results/qc/multiqc",
    log:
        "results/logs/multiqc/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc --force --outdir {params.outdir} {params.indir} > {log} 2>&1"
