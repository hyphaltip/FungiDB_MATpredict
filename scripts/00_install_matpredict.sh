#!/usr/bin/env bash
# Install MATPredict at the pinned release into tools/MATPredict with pixi.
# Run on a login node (needs internet). Safe to re-run.
set -euo pipefail
source "${PROJECT_DIR:-$PWD}/config.sh"

mkdir -p "$TOOLS_DIR"
if [ ! -d "$MATPREDICT_DIR/.git" ]; then
    git clone --branch "$MATPREDICT_VERSION" --depth 1 \
        https://github.com/stajichlab/MATPredict "$MATPREDICT_DIR"
else
    echo "MATPredict already cloned at $MATPREDICT_DIR ($(git -C "$MATPREDICT_DIR" describe --tags))"
fi

command -v pixi >/dev/null || { echo "pixi not found. Install: https://pixi.sh or 'module load pixi'" >&2; exit 1; }
# The pixi cache on a network filesystem is slow; keep it on local disk when a job scratch exists.
export PIXI_CACHE_DIR="${PIXI_CACHE_DIR:-${SCRATCH:-$HOME/.cache}/pixi}"
pixi install --manifest-path "$MATPREDICT_DIR/pixi.toml"

echo "--- installed:"
matpredict_run matpredict --version
git -C "$MATPREDICT_DIR" describe --tags > "$TOOLS_DIR/MATPredict.version"
cat "$TOOLS_DIR/MATPredict.version"
