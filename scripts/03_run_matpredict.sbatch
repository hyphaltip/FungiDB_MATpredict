#!/usr/bin/env bash
# Run `matpredict detect --phylum auto` on every sample in samples.txt as a
# SLURM array. Each array task handles SAMPLES_PER_TASK consecutive samples
# (config.sh, default 25). Submit from the repository root with
#     bash scripts/03_submit.sh
# which computes the array size. Output: results/matpredict/<sample>/.
# A genome that already has detection_report.yaml is skipped, so re-submitting
# after a failure only runs what is missing.
#SBATCH --job-name=matpredict
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
set -euo pipefail
source "${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}/config.sh"

N="${SLURM_ARRAY_TASK_ID:-1}"
K="${SAMPLES_PER_TASK:-25}"
CPUS="${SLURM_CPUS_PER_TASK:-4}"
export MATPREDICT_POLISH_WORKERS="$CPUS"   # MATPredict has no --threads flag; this sets its worker pool
mkdir -p "$LOGS_DIR" "$RUNS_DIR"

run_one() {   # run_one <abbrev>
    local abbrev=$1
    R="FungiDB-${FUNGIDB_RELEASE}_${abbrev}"
    GENOME_GZ="$(genome_file "$abbrev")"   # .fasta.gz or .fasta
    PROTEIN_GZ="$DATA_DIR/protein/${R}_AnnotatedProteins.fasta.gz"
    OUT="$RUNS_DIR/$abbrev"

    [ -n "$GENOME_GZ" ] || { echo "$abbrev: no genome in $DATA_DIR/genome (run scripts/02_fetch_fungidb_data.sh)" >&2; return 1; }
    if [ -s "$OUT/detection_report.yaml" ]; then
        echo "$abbrev: $OUT/detection_report.yaml exists, skipping (delete it to re-run)"; return 0
    fi

    # Node-local scratch for the uncompressed FASTA and the work files.
    # $SCRATCH is set by SLURM on UCR HPCC; fall back to a temp dir for local runs.
    WORK="${SCRATCH:-$(mktemp -d)}/matpredict_$abbrev"
    mkdir -p "$WORK" "$OUT" "$LOGS_DIR"
    rm -f "$OUT/failed_reason.txt"
    gzip -cdf "$GENOME_GZ" > "$WORK/$abbrev.genome.fasta" || return 1   # -f passes a plain file through
    PROTEIN_ARG=()
    if [ -s "$PROTEIN_GZ" ]; then
        # MATPredict needs each protein defline to carry its genomic location as
        # 'contig:start-end:strand'. FungiDB deflines carry
        # '... | location=contig:start-end(-) | ...', so rewrite them to
        # '>protein_id contig:start-end:-'.
        zcat "$PROTEIN_GZ" | awk '
            /^>/ { id=substr($1,2); loc="";
                   if (match($0, /location=[^ |]+/)) { loc=substr($0, RSTART+9, RLENGTH-9);
                       sub(/\(/, ":", loc); sub(/\)/, "", loc) }
                   print ">" id (loc=="" ? "" : " " loc); next }
            { print }' > "$WORK/$abbrev.proteins.fasta"
        head -1 "$WORK/$abbrev.proteins.fasta"
        PROTEIN_ARG=(--proteins "$WORK/$abbrev.proteins.fasta")
    fi

    echo "$abbrev: start $(date)  cpus=$CPUS  work=$WORK"
    # errexit is off inside a function called from `if`, so check the status explicitly
    if ! matpredict_run matpredict detect \
            --genome "$WORK/$abbrev.genome.fasta" \
            "${PROTEIN_ARG[@]}" \
            --phylum auto \
            --sample "$abbrev" \
            --out-dir "$OUT" \
            2>&1 | tee "$OUT/detect.log"; then
        rm -rf "$WORK"; return 1
    fi
    [ -s "$OUT/detection_report.yaml" ] || { rm -rf "$WORK"; return 1; }
    echo "$abbrev: done $(date)"
    rm -rf "$WORK"
}

# samples (N-1)*K+1 .. N*K of samples.txt, ignoring comments and blank lines
first=$(( (N - 1) * K + 1 )); last=$(( N * K ))
mapfile -t CHUNK < <(grep -v '^\s*#' "$SAMPLES_FILE" | grep -v '^\s*$' | sed -n "${first},${last}p")
[ "${#CHUNK[@]}" -gt 0 ] || { echo "task $N: no samples in lines $first-$last of $SAMPLES_FILE" >&2; exit 1; }
echo "task $N: ${#CHUNK[@]} genomes (lines $first-$last)"
failed=0
for abbrev in "${CHUNK[@]}"; do
    if ! run_one "$abbrev"; then
        echo "$abbrev: FAILED" >&2
        mkdir -p "$RUNS_DIR/$abbrev"; echo "detect exited with an error; see logs" > "$RUNS_DIR/$abbrev/failed_reason.txt"
        failed=$((failed + 1))
    fi
done
echo "task $N: done, $failed of ${#CHUNK[@]} failed"
exit $(( failed > 0 ))
