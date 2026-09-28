# =====================================================================
# dedup.smk — P5 : UMI deduplication and hierarchical merges
# =====================================================================
# unit  -> umi_tools dedup                    (results/dedup/unit)
# group -> merge a sample's barcode units     (results/dedup/group)
# cond. -> merge groups of a condition        (results/dedup/condition)


rule dedup_unit:
    """Remove PCR duplicates per mapping unit with UMI-tools (draft l.677)."""
    input:
        "results/mapping/merged/{munit}.sorted.bam",
    output:
        bam="results/dedup/unit/{munit}.dedup.sorted.bam",
    params:
        tmp="results/dedup/unit/{munit}.dedup.bam",
        method=config["dedup"]["method"],
        umi_method=config["dedup"]["extract_umi_method"],
        stats="results/dedup/unit/{munit}.duprm.log",
    threads: threads_for("samtools_threads", 4)
    log:
        "results/logs/dedup/unit/{munit}.log",
    conda:
        "../envs/umitools.yaml"
    shell:
        """
        ( samtools index {input}
          umi_tools dedup -I {input} -L {params.stats} --paired \
            -S {params.tmp} --extract-umi-method {params.umi_method} \
            --method {params.method}
          samtools sort -@ {threads} -o {output.bam} {params.tmp}
          rm -f {params.tmp} ) > {log} 2>&1
        """


rule merge_group_dedup:
    """Merge a biological sample's dedup'd barcode units (draft l.739)."""
    input:
        lambda w: group_dedup_bams(w.group),
    output:
        "results/dedup/group/{group}.dedup.sorted.bam",
    params:
        tmp="results/dedup/group/{group}.merged.bam",
    threads: threads_for("samtools_threads", 4)
    log:
        "results/logs/dedup/merge_group/{group}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        ( samtools merge -@ {threads} -f {params.tmp} {input}
          samtools sort -@ {threads} -o {output} {params.tmp}
          rm -f {params.tmp} ) > {log} 2>&1
        """


rule merge_condition:
    """Merge all groups of a condition into one BAM (draft l.765)."""
    input:
        lambda w: condition_group_bams(w.condition),
    output:
        bam="results/dedup/condition/{condition}.sorted.bam",
        bai="results/dedup/condition/{condition}.sorted.bam.bai",
        sam="results/dedup/condition/{condition}.sorted.sam",
    params:
        tmp="results/dedup/condition/{condition}.merged.bam",
    threads: threads_for("samtools_threads", 4)
    log:
        "results/logs/dedup/merge_condition/{condition}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        ( samtools merge -@ {threads} -f {params.tmp} {input}
          samtools sort -@ {threads} -o {output.bam} {params.tmp}
          samtools index {output.bam}
          samtools view -h -o {output.sam} {output.bam}
          rm -f {params.tmp} ) > {log} 2>&1
        """


rule bam_to_bed:
    """Convert a biological sample's dedup'd BAM to BED (draft l.781)."""
    input:
        "results/dedup/group/{group}.dedup.sorted.bam",
    output:
        "results/dedup/group/{group}.dedup.bed",
    log:
        "results/logs/dedup/bam_to_bed/{group}.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools bamtobed -i {input} > {output} 2> {log}"
