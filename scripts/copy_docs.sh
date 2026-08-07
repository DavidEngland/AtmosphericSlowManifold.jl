#!/usr/bin/env bash
set -euo pipefail

DEST_DIR="docs/src"
mkdir -p "$DEST_DIR"

# Copy all .md files from src/ into docs/src/ (flattened)
find src -type f -name "*.md" -exec cp -v {} "$DEST_DIR/" \;

echo "Markdown files successfully copied to $DEST_DIR."