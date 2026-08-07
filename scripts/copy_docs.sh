#!/usr/bin/env bash
set -euo pipefail

DEST_DIR="docs/src"
mkdir -p "$DEST_DIR"

# Process all .md files inside src/
find src -type f -name "*.md" | sort | while IFS= read -r filepath; do
    # Extract relative path inside src/ (e.g., "Discovery/WSINDyEngine.md")
    rel_path="${filepath#src/}"

    # Replace slashes with underscores to guarantee unique, collision-free target filenames
    flat_name="${rel_path//\//_}"
    target_file="$DEST_DIR/$flat_name"

    echo "Processing $filepath -> $target_file"

    # 1. Inject relative origin path callout at the top
    cat <<EOF > "$target_file"
<!-- Auto-generated from package source -->
> **Source:** \`$filepath\`

EOF

    # 2. Append original file contents
    cat "$filepath" >> "$target_file"
done

echo "Successfully flattened and annotated Markdown files into $DEST_DIR."