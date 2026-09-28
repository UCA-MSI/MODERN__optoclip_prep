#!/usr/bin/env bash
# =====================================================================
# setup.sh — bootstrap the optiCLIP workflow on a Linux host
# ---------------------------------------------------------------------
# 1. create the controller conda env (Snakemake + Slurm plugin)
# 2. either pre-build every per-rule conda env, or build the container
#
# Usage:
#   bash install/setup.sh                 # controller env + per-rule conda envs
#   bash install/setup.sh --container     # controller env + Singularity image
# Run from the snakemake_wf/ project root.
# =====================================================================
set -euo pipefail

MODE="${1:-conda}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"

# ---- pick a conda frontend -----------------------------------------
if command -v mamba >/dev/null 2>&1; then CONDA=mamba; else CONDA=conda; fi
echo ">> using '$CONDA' as conda frontend"

# ---- 1. controller environment -------------------------------------
if ! conda env list | grep -qE '^opticlip\s'; then
  echo ">> creating controller env 'opticlip'"
  "$CONDA" env create -f install/environment.yaml
else
  echo ">> controller env 'opticlip' already exists — skipping"
fi

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate opticlip

# ---- 2. software deployment ----------------------------------------
case "$MODE" in
  --container|container)
    echo ">> building Singularity image install/opticlip.sif"
    apptainer build install/opticlip.sif install/Singularity.def
    echo ">> run with: snakemake --workflow-profile profiles/slurm \\"
    echo "             --software-deployment-method apptainer \\"
    echo "             --apptainer-args '-B /path/to/data'"
    ;;
  *)
    echo ">> pre-creating all per-rule conda envs (no jobs run)"
    snakemake --workflow-profile profiles/local --conda-create-envs-only
    echo ">> done. Dry-run with: snakemake --workflow-profile profiles/local -n"
    ;;
esac

cat <<'EOF'

Next steps
----------
1. Put the raw reads in           resources/reads/
2. Put genome + GTF in            resources/genome/
3. Verify config/samples.tsv      (flexbar_barcode column is a placeholder!)
4. Dry-run:   snakemake --workflow-profile profiles/local -n
5. DAG:       snakemake --workflow-profile profiles/local --dag | dot -Tpng > dag.png
6. Run:       snakemake --workflow-profile profiles/local          # local
        or:   snakemake --workflow-profile profiles/slurm          # cluster
EOF
