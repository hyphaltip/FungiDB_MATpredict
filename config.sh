#!/usr/bin/env bash
# Shared settings for every script in scripts/. Source it; do not run it.
#
# Every script is run from the repository root (the folder holding this file),
# or with PROJECT_DIR set to that folder. SLURM jobs get it from
# SLURM_SUBMIT_DIR, which sbatch sets to the folder sbatch was run in.
# BASH_SOURCE is not used: it resolves to the wrong path inside SLURM jobs.

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
if [ ! -f "$PROJECT_DIR/config.sh" ]; then
    echo "config.sh: PROJECT_DIR=$PROJECT_DIR does not hold config.sh. Run from the repository root." >&2
    exit 1
fi

# --- software versions (pinned) ---------------------------------------------
MATPREDICT_VERSION="v0.8.1"            # https://github.com/stajichlab/MATPredict/releases
GMG_VERSION="v0.3.2"                   # https://github.com/stajichlab/GmG/releases
GMG_INDEX_NAME="bfd-reps-2026.10.5-s20000.gmgi"   # the index MATPredict $MATPREDICT_VERSION is validated with
NCBI_TAXDUMP_SNAPSHOT="2026-10-01"     # https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/

# --- data source -------------------------------------------------------------
FUNGIDB_RELEASE="68"
# A local mirror of the FungiDB download tree, if one exists. When the files
# are not there, scripts/02_fetch_fungidb_data.sh downloads them from FungiDB.
FUNGIDB_LOCAL_MIRROR="${FUNGIDB_LOCAL_MIRROR:-/bigdata/stajichlab/shared/FungiDB/release-${FUNGIDB_RELEASE}}"
FUNGIDB_URL="https://fungidb.org/common/downloads/release-${FUNGIDB_RELEASE}"

# --- folders (all inside the repository) -------------------------------------
TOOLS_DIR="$PROJECT_DIR/tools"
MATPREDICT_DIR="$TOOLS_DIR/MATPredict"
GMG_DIR="$TOOLS_DIR/gmg"
DATA_DIR="$PROJECT_DIR/data"
RESULTS_DIR="$PROJECT_DIR/results"
RUNS_DIR="$RESULTS_DIR/matpredict"       # one detect --out-dir per sample
TABLES_DIR="$RESULTS_DIR/tables"
COMPARE_DIR="$RESULTS_DIR/compare"
LOGS_DIR="$PROJECT_DIR/logs"
SAMPLES_FILE="$PROJECT_DIR/samples.txt"  # one FungiDB organism abbreviation per line

# --- environment MATPredict reads ---------------------------------------------
export MATPREDICT_GMG_CACHE="$GMG_DIR/index"
export MATPREDICT_GMG_INDEX="$GMG_INDEX_NAME"
export MATPREDICT_TAXONOMY="$TOOLS_DIR/taxonomy/ncbi_taxonomy.tsv.zst"
export MATPREDICT_OFFLINE=1              # compute nodes have no internet; use the local taxonomy table
export MATPREDICT_CACHE_DIR="$TOOLS_DIR/cache"
export PATH="$GMG_DIR/gmg-${GMG_VERSION}-linux-64:$PATH"

# --- SLURM -------------------------------------------------------------------
SLURM_PARTITION="${SLURM_PARTITION:-epyc}"
SLURM_CPUS="${SLURM_CPUS:-8}"
SLURM_MEM="${SLURM_MEM:-16G}"
SLURM_TIME="${SLURM_TIME:-4:00:00}"

# Run a command inside the MATPredict pixi environment.
matpredict_run() {
    (cd "$MATPREDICT_DIR" && pixi run --manifest-path "$MATPREDICT_DIR/pixi.toml" "$@")
}

# Use the pixi module on UCR HPCC when pixi is not already on PATH.
if ! command -v pixi >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then module load pixi >/dev/null 2>&1 || true; fi
fi
