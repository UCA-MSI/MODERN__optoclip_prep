# =====================================================================
# fasta.smk — FASTA conversion and read-length histograms (diagnostic)
# =====================================================================
# Per {sample}/{read}: convert the clean FASTQ to FASTA and draw the
# read-length distribution (FASTX-Toolkit), as in the draft.


rule fastq_to_fasta:
    """Convert a clean FASTQ to FASTA (draft l.435)."""
    input:
        "results/trim/final/{sample}_{read}.fastq.gz",
    output:
        "results/fasta/{sample}_{read}.fasta.gz",
    log:
        "results/logs/fasta/fastq_to_fasta/{sample}_{read}.log",
    conda:
        "../envs/fastx.yaml"
    shell:
        "zcat {input} | fastq_to_fasta -n -r | gzip > {output} 2> {log}"


rule readlength_histogram:
    """Read-length histogram PNG from a FASTA (draft l.478)."""
    input:
        "results/fasta/{sample}_{read}.fasta.gz",
    output:
        "results/fasta/{sample}_{read}.readlength.png",
    log:
        "results/logs/fasta/readlength/{sample}_{read}.log",
    conda:
        "../envs/fastx.yaml"
    shell:
        "fasta_clipping_histogram.pl {input} {output} > {log} 2>&1"
