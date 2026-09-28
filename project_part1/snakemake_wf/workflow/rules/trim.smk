# =====================================================================
# trim.smk — P2b : adapter / barcode trimming and length filtering
# =====================================================================
# Per {sample} (barcode-level unit). Reproduces the draft chain:
#   cutadapt(R1 adapter) -> header #->space -> cutadapt(R2 barcode)
#   -> paired length filter -> header space->#  (= clean `final_` FASTQ)
# All intermediates are temp(); only results/trim/final/ is kept.


rule cutadapt_r1_adapter:
    """Remove the 3' sequencing adapter from R1 (draft l.115)."""
    input:
        "results/demux/flexbarOut_barcode_{sample}_1.fastq.gz",
    output:
        temp("results/trim/cutadapt_r1/{sample}_1.fastq.gz"),
    params:
        adapter=config["trim"]["r1_adapter"],
    threads: threads_for("cutadapt_threads", 4)
    log:
        "results/logs/trim/cutadapt_r1/{sample}.log",
    conda:
        "../envs/cutadapt.yaml"
    shell:
        "cutadapt -j {threads} -a {params.adapter} -o {output} {input} > {log} 2>&1"


rule fix_header_space_r1:
    """Header '#1' -> ' 1' so cutadapt sees a standard read id (draft l.183)."""
    input:
        "results/trim/cutadapt_r1/{sample}_1.fastq.gz",
    output:
        temp("results/trim/hdr_space_r1/{sample}_1.fastq.gz"),
    log:
        "results/logs/trim/fix_header_space_r1/{sample}.log",
    conda:
        "../envs/base.yaml"
    shell:
        "zcat {input} | "
        "awk '{{if(NR%4==1){{sub(\"#1\",\" 1\",$1)}} print}}' | "
        "gzip > {output} 2> {log}"


rule fix_header_space_r2:
    """Header '#2' -> ' 2' on the demuxed R2 (draft l.204)."""
    input:
        "results/demux/flexbarOut_barcode_{sample}_2.fastq.gz",
    output:
        temp("results/trim/hdr_space_r2/{sample}_2.fastq.gz"),
    log:
        "results/logs/trim/fix_header_space_r2/{sample}.log",
    conda:
        "../envs/base.yaml"
    shell:
        "zcat {input} | "
        "awk '{{if(NR%4==1){{sub(\"#2\",\" 2\",$1)}} print}}' | "
        "gzip > {output} 2> {log}"


rule cutadapt_r2_barcode:
    """Remove the 3' UMI+barcode+UMI block from R2 (draft l.228).

    The adapter is per-sample: NNNN<sense barcode>NNNNN (from the sheet).
    """
    input:
        "results/trim/hdr_space_r2/{sample}_2.fastq.gz",
    output:
        temp("results/trim/cutadapt_r2/{sample}_2.fastq.gz"),
    params:
        adapter=r2_adapter,
    threads: threads_for("cutadapt_threads", 4)
    log:
        "results/logs/trim/cutadapt_r2/{sample}.log",
    conda:
        "../envs/cutadapt.yaml"
    shell:
        "cutadapt -j {threads} -a {params.adapter} -o {output} {input} > {log} 2>&1"


rule cutadapt_length_filter:
    """Paired length filter, drop pairs with a read < min_length (draft l.277)."""
    input:
        r1="results/trim/hdr_space_r1/{sample}_1.fastq.gz",
        r2="results/trim/cutadapt_r2/{sample}_2.fastq.gz",
    output:
        r1=temp("results/trim/lenfilter/{sample}_1.fastq.gz"),
        r2=temp("results/trim/lenfilter/{sample}_2.fastq.gz"),
    params:
        min_len=config["trim"]["min_length"],
    threads: threads_for("cutadapt_threads", 4)
    log:
        "results/logs/trim/length_filter/{sample}.log",
    conda:
        "../envs/cutadapt.yaml"
    shell:
        "cutadapt -j {threads} --pair-filter=any --minimum-length={params.min_len} "
        "-o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log} 2>&1"


rule restore_headers:
    """Header ' ' and '/' -> '#' on the clean FASTQ (draft l.345).

    Shared logic lives in scripts/restore_headers.sh (also used by the
    miRNA branch, see mirna.smk::mirna_restore_headers).
    """
    input:
        fastq="results/trim/lenfilter/{sample}_{read}.fastq.gz",
        script="workflow/scripts/restore_headers.sh",
    output:
        "results/trim/final/{sample}_{read}.fastq.gz",
    log:
        "results/logs/trim/restore_headers/{sample}_{read}.log",
    conda:
        "../envs/base.yaml"
    shell:
        "bash {input.script} {input.fastq} {output} 2> {log}"
