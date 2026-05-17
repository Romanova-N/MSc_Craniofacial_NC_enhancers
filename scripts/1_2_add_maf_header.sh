#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <maf_file_or_dir>" >&2
  exit 1
fi

TARGET="$1"

add_header() {
  local f="$1"
  # skip empty files
  [[ ! -s "$f" ]] && return 0

  # if already has header, do nothing
  if head -n 1 "$f" | grep -q '^##maf'; then
    return 0
  fi

  # prepend header safely
  tmp="$(mktemp)"
  printf "##maf version=1\n" > "$tmp"
  cat "$f" >> "$tmp"
  mv "$tmp" "$f"
}

export -f add_header

if [[ -d "$TARGET" ]]; then
  # all .maf files in directory (adjust pattern if needed)
  find "$TARGET" -maxdepth 1 -type f -name "*.maf" -print0 \
    | while IFS= read -r -d '' f; do
        add_header "$f"
      done
elif [[ -f "$TARGET" ]]; then
  add_header "$TARGET"
else
  echo "Error: not a file or directory: $TARGET" >&2
  exit 1
fi

echo "Done."