#!/usr/bin/env bash
# Download (1) the GmG "Guess My Genome" binary, (2) the GmG reference index,
# and (3) the NCBI taxonomy dump, then build MATPredict's offline taxonomy table.
# Run on a login node (needs internet). Safe to re-run; existing files are kept.
set -euo pipefail
source "${PROJECT_DIR:-$PWD}/config.sh"

# 1. GmG binary, checked against its published SHA-256.
mkdir -p "$GMG_DIR"
P="linux-64"
if [ ! -x "$GMG_DIR/gmg-$GMG_VERSION-$P/gmg" ]; then
    ( cd "$GMG_DIR"
      curl -fsSLO "https://github.com/stajichlab/GmG/releases/download/$GMG_VERSION/gmg-$GMG_VERSION-$P.tar.gz"
      curl -fsSLO "https://github.com/stajichlab/GmG/releases/download/$GMG_VERSION/gmg-$GMG_VERSION-$P.tar.gz.sha256"
      sha256sum -c "gmg-$GMG_VERSION-$P.tar.gz.sha256"
      tar xzf "gmg-$GMG_VERSION-$P.tar.gz"
      rm "gmg-$GMG_VERSION-$P.tar.gz" "gmg-$GMG_VERSION-$P.tar.gz.sha256" )
fi
gmg --version

# 2. GmG reference index into the cache folder MATPredict reads ($MATPREDICT_GMG_CACHE),
#    checked against the SHA256SUMS of the GmG index release (index-<version>).
mkdir -p "$MATPREDICT_GMG_CACHE"
IDX_VERSION=$(echo "$GMG_INDEX_NAME" | sed -E 's/^bfd(-reps)?-([0-9.]+)-s[0-9]+\.gmgi$/\2/')
if [ ! -s "$MATPREDICT_GMG_CACHE/$GMG_INDEX_NAME" ]; then
    curl -fsSL -o "$MATPREDICT_GMG_CACHE/$GMG_INDEX_NAME.part" "https://gmg.fungalgenomes.org/index/$GMG_INDEX_NAME"
    mv "$MATPREDICT_GMG_CACHE/$GMG_INDEX_NAME.part" "$MATPREDICT_GMG_CACHE/$GMG_INDEX_NAME"
fi
curl -fsSL -o "$MATPREDICT_GMG_CACHE/SHA256SUMS.index-$IDX_VERSION" \
    "https://github.com/stajichlab/GmG/releases/download/index-$IDX_VERSION/SHA256SUMS"
( cd "$MATPREDICT_GMG_CACHE" && grep " \*\?$GMG_INDEX_NAME\$" "SHA256SUMS.index-$IDX_VERSION" | sed 's/ \*/  /' | sha256sum -c - )
# the site manifest, so MATPredict can say whether a newer index exists without network on a compute node
curl -fsSL -o "$MATPREDICT_GMG_CACHE/manifest.json" "https://gmg.fungalgenomes.org/index/manifest.json" || true
ls -l "$MATPREDICT_GMG_CACHE"

# 3. NCBI taxonomy snapshot -> MATPredict offline table ($MATPREDICT_TAXONOMY).
TAX_DIR="$(dirname "$MATPREDICT_TAXONOMY")"
mkdir -p "$TAX_DIR"
ZIP="$TAX_DIR/taxdmp_${NCBI_TAXDUMP_SNAPSHOT}.zip"
if [ ! -s "$MATPREDICT_TAXONOMY" ]; then
    [ -s "$ZIP" ] || curl -fsSL -o "$ZIP" \
        "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_${NCBI_TAXDUMP_SNAPSHOT}.zip"
    matpredict_run matpredict curate-db build-taxonomy --taxdump "$ZIP" \
        --snapshot "$NCBI_TAXDUMP_SNAPSHOT" --out "$MATPREDICT_TAXONOMY"
    rm -f "$ZIP"
fi
ls -l "$MATPREDICT_TAXONOMY"
echo "gmg $GMG_VERSION, index $GMG_INDEX_NAME, taxonomy $NCBI_TAXDUMP_SNAPSHOT" > "$TOOLS_DIR/gmg.version"
