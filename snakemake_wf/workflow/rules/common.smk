# =====================================================================
# common.smk — shared setup: sample sheet, wildcard constraints, helpers
# =====================================================================
# This module is included first by the Snakefile. It exposes the sample
# model (units / groups / conditions / mapping units) and small pure
# helper functions used across all rule modules. No rule lives here.

import os
import pandas as pd
from snakemake.utils import validate

# ---- Configuration --------------------------------------------------
validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#")
    .set_index("sample_id", drop=False)
    .sort_index()
)
validate(samples, schema="../schemas/samples.schema.yaml")

# ---- Derived sample model -------------------------------------------
# UNITS   : barcode-level units (e.g. "1_1", "1_2", ...) — the {sample} wildcard
# GROUPS  : biological samples (e.g. "1", "2", ...)
# MAPPING_UNITS : what actually goes into the aligner. A group listed in
#           `merge_before_mapping` is aligned as a single merged unit named
#           after the group; every other group aligns its barcode units.
UNITS = list(samples["sample_id"])
GROUPS = list(dict.fromkeys(samples["sample_group"]))
MERGE_BEFORE_MAPPING = set(config.get("merge_before_mapping", []))
CONDITIONS = config["conditions"]


def units_of_group(group):
    """Barcode-level units belonging to a biological sample."""
    return list(samples.loc[samples["sample_group"] == group, "sample_id"])


MAPPING_UNITS = []
for _g in GROUPS:
    if _g in MERGE_BEFORE_MAPPING:
        MAPPING_UNITS.append(_g)
    else:
        MAPPING_UNITS.extend(units_of_group(_g))


# ---- Wildcard constraints (disambiguate the DAG) --------------------
wildcard_constraints:
    sample=r"\d+_\d+",  # barcode-level unit
    group=r"\d+",  # biological sample
    munit=r"\d+(_\d+)?",  # mapping unit (unit or merged group)
    read=r"[12]",
    condition="|".join(CONDITIONS),


# ---- Per-sample parameter helpers -----------------------------------
def r2_adapter(wildcards):
    """cutadapt R2 3' adapter for a unit: <flank5><sense barcode><flank3>."""
    bc = samples.loc[wildcards.sample, "r2_barcode"]
    return (
        config["trim"]["r2_adapter_flank5"] + bc + config["trim"]["r2_adapter_flank3"]
    )


# ---- Mapping-unit input resolution ----------------------------------
def mapping_unit_fastq(munit, read):
    """final FASTQ(s) feeding a mapping unit for a given read (1|2).

    A plain unit maps its own final FASTQ; a merged group concatenates
    the final FASTQs of all its barcode units (draft: `cat` of sample 11).
    """
    if "_" in munit:
        return [f"results/trim/final/{munit}_{read}.fastq.gz"]
    return [f"results/trim/final/{u}_{read}.fastq.gz" for u in units_of_group(munit)]


def is_merged_group(group):
    return group in MERGE_BEFORE_MAPPING


def group_dedup_bams(group):
    """Deduplicated, sorted BAM(s) that make up a biological sample."""
    if is_merged_group(group):
        return [f"results/dedup/unit/{group}.dedup.sorted.bam"]
    return [f"results/dedup/unit/{u}.dedup.sorted.bam" for u in units_of_group(group)]


def condition_group_bams(condition):
    """Per-group dedup'd BAMs that make up a condition."""
    return [f"results/dedup/group/{g}.dedup.sorted.bam" for g in CONDITIONS[condition]]


# ---- Resource helper ------------------------------------------------
def threads_for(key, default=1):
    return int(config.get("resources", {}).get(key, default))
