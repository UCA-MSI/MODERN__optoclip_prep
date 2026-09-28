# =====================================================================
# mapping.smk — P3/P4 : genomic mapping (STAR + Novoalign) and merge
# =====================================================================
# Per mapping unit {munit} (a barcode unit, or a whole group when listed
# in merge_before_mapping). STAR maps first; Novoalign rescues the reads
# STAR left unmapped; the two BAMs are sorted and merged.


rule star_index:
    """Build the STAR genome index (draft l.537)."""
    input:
        fasta=config["reference"]["genome_fasta"],
        gtf=config["reference"]["gtf"],
    output:
        directory(config["reference"]["star_index"]),
    params:
        overhang=config["star"]["sjdb_overhang"],
    threads: threads_for("star_index_threads", 8)
    log:
        "results/logs/mapping/star_index.log",
    conda:
        "../envs/star.yaml"
    shell:
        "mkdir -p {output} && "
        "STAR --runThreadN {threads} --runMode genomeGenerate "
        "--genomeDir {output} --genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} --sjdbOverhang {params.overhang} > {log} 2>&1"


if MERGE_BEFORE_MAPPING:

    # Concatenate a group's barcode units before mapping (draft l.523).
    rule premerge_fastq:
        wildcard_constraints:
            group="|".join(MERGE_BEFORE_MAPPING),
        input:
            lambda w: mapping_unit_fastq(w.group, w.read),
        output:
            temp("results/mapping/premerge/{group}_{read}.fastq.gz"),
        log:
            "results/logs/mapping/premerge/{group}_{read}.log",
        shell:
            "cat {input} > {output} 2> {log}"


def star_reads(wildcards):
    """R1/R2 feeding STAR for a mapping unit (plain unit vs merged group)."""
    m = wildcards.munit
    if "_" in m:
        base = "results/trim/final/{m}_{r}.fastq.gz"
    else:
        base = "results/mapping/premerge/{m}_{r}.fastq.gz"
    return {"r1": base.format(m=m, r="1"), "r2": base.format(m=m, r="2")}


rule star_align:
    """Map a unit with STAR; emit sorted BAM + unmapped mates (draft l.543)."""
    input:
        unpack(star_reads),
        index=config["reference"]["star_index"],
        gtf=config["reference"]["gtf"],
    output:
        bam="results/mapping/star/{munit}Aligned.sortedByCoord.out.bam",
        mate1="results/mapping/star/{munit}Unmapped.out.mate1",
        mate2="results/mapping/star/{munit}Unmapped.out.mate2",
    params:
        prefix="results/mapping/star/{munit}",
        overhang=config["star"]["sjdb_overhang"],
        extra=config["star"]["extra"],
    threads: threads_for("star_align_threads", 8)
    log:
        "results/logs/mapping/star_align/{munit}.log",
    conda:
        "../envs/star.yaml"
    shell:
        "STAR --runThreadN {threads} --runMode alignReads "
        "--genomeDir {input.index} --sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.overhang} --readFilesCommand zcat "
        "{params.extra} --outFileNamePrefix {params.prefix} "
        "--readFilesIn {input.r1} {input.r2} > {log} 2>&1"


rule novoindex:
    """Build the Novoalign index (draft l.567). Novoalign is proprietary."""
    input:
        config["reference"]["genome_fasta"],
    output:
        config["reference"]["novoalign_index"],
    params:
        k=config["novoalign"]["index_kmer"],
        s=config["novoalign"]["index_step"],
    log:
        "results/logs/mapping/novoindex.log",
    conda:
        "../envs/novoalign.yaml"
    shell:
        "novoindex -k {params.k} -s {params.s} {output} {input} > {log} 2>&1"


rule novoalign:
    """Rescue STAR-unmapped reads with Novoalign (draft l.569)."""
    input:
        mate1="results/mapping/star/{munit}Unmapped.out.mate1",
        mate2="results/mapping/star/{munit}Unmapped.out.mate2",
        index=config["reference"]["novoalign_index"],
    output:
        sam=temp("results/mapping/novoalign/{munit}.sam"),
        report="results/mapping/novoalign/{munit}.report.txt",
    params:
        t=config["novoalign"]["threshold"],
        l=config["novoalign"]["min_len"],
        r=config["novoalign"]["repeat"],
    log:
        "results/logs/mapping/novoalign/{munit}.log",
    conda:
        "../envs/novoalign.yaml"
    shell:
        "novoalign -t {params.t} -l {params.l} -F STDFQ -r {params.r} "
        "-o sam -d {input.index} -f {input.mate1} {input.mate2} "
        "> {output.sam} 2> {output.report}"


rule sort_novoalign:
    """Coordinate-sort the Novoalign SAM into BAM."""
    input:
        "results/mapping/novoalign/{munit}.sam",
    output:
        temp("results/mapping/novoalign/{munit}.sorted.bam"),
    threads: threads_for("samtools_threads", 4)
    log:
        "results/logs/mapping/sort_novoalign/{munit}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input} 2> {log}"


rule sort_star:
    """Re-sort the STAR BAM (draft re-sorts before merging)."""
    input:
        "results/mapping/star/{munit}Aligned.sortedByCoord.out.bam",
    output:
        temp("results/mapping/star/{munit}.sorted.bam"),
    threads: threads_for("samtools_threads", 4)
    log:
        "results/logs/mapping/sort_star/{munit}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input} 2> {log}"


rule merge_aligners:
    """Merge STAR + Novoalign alignments for a unit (draft l.632)."""
    input:
        star="results/mapping/star/{munit}.sorted.bam",
        novo="results/mapping/novoalign/{munit}.sorted.bam",
    output:
        temp("results/mapping/merged/{munit}.merged.bam"),
    threads: threads_for("samtools_threads", 4)
    log:
        "results/logs/mapping/merge_aligners/{munit}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools merge -@ {threads} -f {output} {input.star} {input.novo} 2> {log}"


rule sort_merged_unit:
    """Coordinate-sort the merged per-unit BAM (feeds dedup, draft l.652)."""
    input:
        "results/mapping/merged/{munit}.merged.bam",
    output:
        "results/mapping/merged/{munit}.sorted.bam",
    threads: threads_for("samtools_threads", 4)
    log:
        "results/logs/mapping/sort_merged_unit/{munit}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input} 2> {log}"
