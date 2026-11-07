#!/usr/bin/env bash
# For every sample in samples.txt, put the FungiDB genome, GFF and protein files
# under data/. A local mirror is symlinked when it has the file; otherwise the
# file is downloaded from fungidb.org (login node, needs internet).
set -euo pipefail
source "${PROJECT_DIR:-$PWD}/config.sh"

mkdir -p "$DATA_DIR/genome" "$DATA_DIR/gff" "$DATA_DIR/protein"

fetch() {   # fetch <abbrev> <subdir> <filename> <dest dir>
    local abbrev=$1 sub=$2 name=$3 dest=$4
    local target="$dest/$name"
    [ -s "$target" ] && { echo "have  $target"; return; }
    local local_src="$FUNGIDB_LOCAL_MIRROR/$abbrev/$sub/data/$name"
    if [ -s "$local_src" ]; then
        ln -sfn "$local_src" "$target"; echo "link  $target -> $local_src"
    else
        echo "fetch $FUNGIDB_URL/$abbrev/$sub/data/$name"
        curl -fsSL --retry 3 -o "$target.part" "$FUNGIDB_URL/$abbrev/$sub/data/$name"
        mv "$target.part" "$target"
    fi
}

grep -v '^\s*#' "$SAMPLES_FILE" | grep -v '^\s*$' | while read -r abbrev; do
    R="FungiDB-${FUNGIDB_RELEASE}_${abbrev}"
    fetch "$abbrev" fasta "${R}_Genome.fasta.gz"            "$DATA_DIR/genome"
    fetch "$abbrev" gff   "${R}.gff.gz"                     "$DATA_DIR/gff"
    fetch "$abbrev" fasta "${R}_AnnotatedProteins.fasta.gz" "$DATA_DIR/protein"
done
echo "--- data/:"; ls -lL "$DATA_DIR"/*/
