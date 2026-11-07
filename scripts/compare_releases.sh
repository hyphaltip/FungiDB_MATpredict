#!/usr/bin/env bash
# Which genomes differ between two FungiDB releases? Writes results/release_changes.tsv
# with one row per sample: status new_in_<B> / dropped_after_<A> / same_genome / changed_genome.
# Genome identity is the md5 of the uncompressed sequence, computed only for
# samples present in both releases.
#   bash scripts/compare_releases.sh [OLD_RELEASE NEW_RELEASE]   (default: 68 71)
set -euo pipefail
source "${PROJECT_DIR:-$PWD}/config.sh"
A="${1:-68}"; B="${2:-71}"
OUT="$PROJECT_DIR/results/release_changes.tsv"; mkdir -p "$(dirname "$OUT")"
list() { ls "$PROJECT_DIR/data/release-$1/genome" | sed -n "s/^FungiDB-$1_\(.*\)_Genome\.fasta\(\.gz\)\?$/\1/p" | sort -u; }
gfile() { local b="$PROJECT_DIR/data/release-$1/genome/FungiDB-$1_$2_Genome.fasta"; [ -s "$b.gz" ] && echo "$b.gz" || echo "$b"; }
md5() { gzip -cdf "$1" | md5sum | cut -c1-32; }
printf "sample\tin_release_%s\tin_release_%s\tstatus\tmd5_%s\tmd5_%s\n" "$A" "$B" "$A" "$B" > "$OUT"
sort -u <(list "$A") <(list "$B") | while read -r s; do
    inA=0; inB=0; [ -n "$(list "$A" | grep -x "$s")" ] && inA=1; [ -n "$(list "$B" | grep -x "$s")" ] && inB=1
    ma=""; mb=""
    if [ $inA = 1 ] && [ $inB = 1 ]; then
        ma=$(md5 "$(gfile "$A" "$s")"); mb=$(md5 "$(gfile "$B" "$s")")
        [ "$ma" = "$mb" ] && st=same_genome || st=changed_genome
    elif [ $inB = 1 ]; then st="new_in_$B"; else st="dropped_after_$A"; fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$s" $inA $inB "$st" "$ma" "$mb"
done >> "$OUT"
echo "wrote $OUT"; cut -f4 "$OUT" | tail -n +2 | sort | uniq -c
