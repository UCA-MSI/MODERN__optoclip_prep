# =====================================================================
# peakcalling.smk — P6 : pyicoclip peak calling + R visualization
# =====================================================================
# pyicos / pyicoclip only run under Python 2 (see envs/pyicoclip.yaml).
# The annotation is converted once; peaks are called per biological
# sample (group) and summarized as a peak-length histogram.


rule gtf_to_bed:
    """Convert the GENCODE GTF to BED with gffread (draft l.816)."""
    input:
        config["reference"]["gtf"],
    output:
        temp("results/peaks/annotation.bed"),
    log:
        "results/logs/peaks/gtf_to_bed.log",
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread -E {input} -T -o {output} 2> {log}"


rule sort_annotation_bed:
    input:
        "results/peaks/annotation.bed",
    output:
        temp("results/peaks/annotation.sorted.bed"),
    log:
        "results/logs/peaks/sort_annotation.log",
    conda:
        "../envs/base.yaml"
    shell:
        "sort -k1,1 -k2,2n {input} > {output} 2> {log}"


rule annotation_to_bedpk:
    """Convert the sorted annotation BED to pyicos region format (draft l.818)."""
    input:
        "results/peaks/annotation.sorted.bed",
    output:
        "results/peaks/annotation.bedpk",
    log:
        "results/logs/peaks/annotation_to_bedpk.log",
    conda:
        "../envs/pyicoclip.yaml"
    shell:
        "pyicos convert {input} {output} -f bed -F bed_pk > {log} 2>&1"


rule pyicos_extend:
    """Extend reads by N nt before peak calling (draft l.801)."""
    input:
        "results/dedup/group/{group}.dedup.bed",
    output:
        temp("results/peaks/{group}.ext.bed"),
    params:
        n=config["peaks"]["extend"],
    log:
        "results/logs/peaks/extend/{group}.log",
    conda:
        "../envs/pyicoclip.yaml"
    shell:
        "pyicos extend {input} {output} {params.n} -f bed -F bed > {log} 2>&1"


rule sort_extended_bed:
    input:
        "results/peaks/{group}.ext.bed",
    output:
        temp("results/peaks/{group}.ext.sorted.bed"),
    log:
        "results/logs/peaks/sort_extended/{group}.log",
    conda:
        "../envs/base.yaml"
    shell:
        "sort -k1,1 -k2,2n {input} > {output} 2> {log}"


rule pyicoclip:
    """Call CLIP peaks per biological sample (draft l.834)."""
    input:
        bed="results/peaks/{group}.ext.sorted.bed",
        region="results/peaks/annotation.bedpk",
    output:
        "results/peaks/{group}.genome_sorted.pk",
    params:
        pval=config["peaks"]["p_value"],
    log:
        "results/logs/peaks/pyicoclip/{group}.log",
    conda:
        "../envs/pyicoclip.yaml"
    shell:
        "pyicoclip {input.bed} {output} -f bed --region {input.region} "
        "--stranded --p-value {params.pval} > {log} 2>&1"


rule peak_length_histogram:
    """Peak-length histogram per biological sample (draft l.849)."""
    input:
        pk="results/peaks/{group}.genome_sorted.pk",
    output:
        pdf="results/peaks/{group}.peaklength.pdf",
    params:
        plot=config["plot"],
    log:
        "results/logs/peaks/hist/{group}.log",
    conda:
        "../envs/r-viz.yaml"
    script:
        "../scripts/peak_length_hist.R"
