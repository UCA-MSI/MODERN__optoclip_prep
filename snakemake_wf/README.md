# optiCLIP — Snakemake workflow

Paired-end iCLIP-seq pipeline, migrated from `../draft_code/raw_bash/script.txt`.
It goes from a single multiplexed library to per-sample CLIP peaks:

```text
QC → preprocessing → demultiplexing → trimming → FASTA/QC
   → mapping (STAR + Novoalign) → dedup (UMI-tools) → merges
   → peak calling (pyicoclip) → peak-length plots
```

A second, independent branch profiles small RNAs from the same
preprocessing (migrated from `2_script_miRNA.txt`):

```text
trim/lenfilter → miRNA size selection (18-25nt) → header restore
   → per-sample barcode merge → SPORTS1.1 annotation (Bowtie1 + miRBase)
   → gene-count extraction → RPM normalization → merged expression matrix
```

The full design rationale (audit, architecture, DAG, dev plan) is in
[`docs/01_audit_architecture_plan.md`](docs/01_audit_architecture_plan.md)
(main branch) and
[`docs/02_miRNA_script_uml_audit.md`](docs/02_miRNA_script_uml_audit.md)
(miRNA branch).

## Layout

```text
Snakefile                 entrypoint: config, includes, targets
config/
  config.yaml             all scientific parameters
  samples.tsv             sample sheet (units, barcodes, conditions)
workflow/
  rules/*.smk             rules grouped by pipeline phase
  scripts/                equalize_reads.sh, restore_headers.sh,
                           peak_length_hist.R, merge_mirna_counts.R
  envs/*.yaml             one conda env per tool
  schemas/*.yaml          config / sample-sheet validation
profiles/{local,slurm}/   execution profiles
install/                  setup.sh, environment.yaml, Singularity.def
resources/                inputs (reads, genome) — not versioned
results/                  outputs (+ results/logs/)
```

## Quick start

```bash
# 1. install the controller env and per-rule envs (or the container)
bash install/setup.sh              # conda envs
bash install/setup.sh --container  # Singularity image instead

conda activate opticlip

# 2. stage inputs
#    resources/reads/   iClip_S1_R1_001.fastq.gz, iClip_S1_R2_001.fastq.gz
#    resources/genome/  GRCh38.p14.genome.fa, gencode.v44...gtf

# 3. verify config/samples.tsv  (see WARNING below)

# 4. dry-run and DAG
snakemake --workflow-profile profiles/local -n
snakemake --workflow-profile profiles/local --dag | dot -Tpng > dag.png

# 5. run
snakemake --workflow-profile profiles/local        # local
snakemake --workflow-profile profiles/slurm        # cluster
```

Handy sub-targets: `snakemake … preprocess | trimming | mapping | dedup | peaks`.
`mirna` (SPORTS1.1 branch) is separate and not part of `all` — see below.

## Configuration

- **`config/config.yaml`** — every parameter of the original script (read
  structure, quality thresholds, adapters, references, STAR/Novoalign options,
  peak parameters, resources). Notably `preprocess.barcode_qc_length` selects
  Option B = 15 nt (default) or Option A = 6 nt for the barcode-region QC.
- **`config/samples.tsv`** — one row per barcode unit (`<group>_<n>`) with its
  condition and barcodes. Groups are merged into biological samples, then into
  conditions (`control`, `KD`).

> ⚠️ **`flexbar_barcode` is a placeholder** (reverse-complement of the sense
> barcode). The original `barcodes_antisense.fasta` was not in the repo — verify
> or replace these before running the demux step. `r2_barcode` values are taken
> verbatim from the draft and are trusted.

## miRNA branch (SPORTS1.1)

`workflow/rules/mirna.smk` re-profiles the same iCLIP reads for small
RNAs. It is a separate target, not part of `all`:

```bash
# stage the miRNA reference (download from https://www.mirbase.org/)
#    resources/mirna/mature.fa

snakemake --workflow-profile profiles/local mirna
```

Each step is its own wildcarded rule (e.g. `snakemake results/mirna/sports/1_R1`
runs SPORTS for group 1 only), so any part of the branch can be run in
isolation the same way as the rest of the workflow. Configuration lives
under `mirna:` and `reference.mirna_genome_bowtie_index` /
`reference.mirbase_bowtie_index` in `config/config.yaml`.

Two points in the original script were ambiguous or left implicit and
required a modeling choice during migration (both documented in the
rule docstrings and in
[`docs/02_miRNA_script_uml_audit.md`](docs/02_miRNA_script_uml_audit.md)
§2.4, points R1 and R5) — **validate them against a real SPORTS1.1 run
before trusting the output counts**:

- which file actually feeds `bowtie-build` for the miRBase index
  (`mirna_bowtie_index` rule assumes it's `hsa_mature_final.fa`, the
  species-filtered, U→T-converted `mature.fa`);
- how `gene_counts_sample{N}.txt` is derived from SPORTS' output
  (`sports_extract_counts` rule reconstructs it from SPORTS' documented
  `*_output.txt` table, summing `Reads` per `Annotation` matching
  `mirna.annotation_regex`).

SPORTS1.1 itself (Perl/Python2, not conda-packaged) must be installed
manually and put on `PATH` — see `workflow/envs/sports.yaml`.

## Notes on reproducibility

- Software is deployed per rule via `workflow/envs/*.yaml` (conda/mamba) or via
  the Singularity image. **Versions are intentionally unpinned** — pin them
  before production.
- **Novoalign** is proprietary (licence required); **pyicoclip** is Python 2
  only (installed via pip in a dedicated env / container layer).
- Intermediate files are `temp()`; only deliverables are kept. Every rule writes
  to `results/logs/`.

## Fidelity to the original

The scientific logic is preserved exactly: header (`sed`/`awk`) and read-id
(`comm`) operations are reproduced byte-for-byte, and sample 11's pre-mapping
FASTQ merge is kept as a documented special case (`merge_before_mapping`). The
draft's ~120 duplicated commands become parametrized, wildcard-driven rules.
