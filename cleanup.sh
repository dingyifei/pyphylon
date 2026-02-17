#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"

phase=""
dir=""
force=false
dry_run=false

usage() {
    cat <<'EOF'
Usage: cleanup.sh [-p <phase>] [-d <dir>] [--force] [--dry-run]

Remove intermediate/output files from temp/ and output/ directories.

Flags:
  -p <phase>   Remove files matching a phase prefix (e.g. -p 1 for 1[a-z]_*, -p 1a for 1a_*)
  -d <dir>     Limit cleanup to 'temp' or 'output' only
  --force      Also remove files/folders containing "protected" in name
  --dry-run    Print what would be deleted without deleting

Files named readme.md are always preserved.
Items with "protected" in their name are skipped unless --force is used.
EOF
    exit "${1:-0}"
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -p)
            [[ $# -lt 2 ]] && { echo "Error: -p requires an argument"; usage 1; }
            phase="$2"; shift 2 ;;
        -d)
            [[ $# -lt 2 ]] && { echo "Error: -d requires an argument"; usage 1; }
            dir="$2"; shift 2 ;;
        --force)
            force=true; shift ;;
        --dry-run)
            dry_run=true; shift ;;
        -h|--help)
            usage 0 ;;
        *)
            echo "Error: unknown option '$1'"
            usage 1 ;;
    esac
done

# Validate -d value
if [[ -n "$dir" && "$dir" != "temp" && "$dir" != "output" ]]; then
    echo "Error: -d must be 'temp' or 'output', got '$dir'"
    exit 1
fi

# Validate -p value (single digit or digit+letters)
if [[ -n "$phase" && ! "$phase" =~ ^[0-9][a-zA-Z]*$ ]]; then
    echo "Error: -p must be a digit optionally followed by letters (e.g. 1, 1a, 2b), got '$phase'"
    exit 1
fi

# Build target directories
if [[ -n "$dir" ]]; then
    targets=("$SCRIPT_DIR/$dir")
else
    targets=("$SCRIPT_DIR/temp" "$SCRIPT_DIR/output")
fi

deleted=0
skipped_protected=0

for target in "${targets[@]}"; do
    [[ -d "$target" ]] || continue

    # Build find arguments
    find_args=("$target" -maxdepth 1 -mindepth 1)

    # Always exclude readme.md
    find_args+=(-not -iname "readme.md")

    # Phase filter
    if [[ -n "$phase" ]]; then
        if [[ ${#phase} -eq 1 ]]; then
            # Single digit: match <digit>[a-zA-Z]_*
            find_args+=(-name "${phase}[a-zA-Z]_*")
        else
            # Digit+letters: match <phase>_*
            find_args+=(-name "${phase}_*")
        fi
    fi

    while IFS= read -r item; do
        [[ -z "$item" ]] && continue
        basename="$(basename "$item")"

        # Check protected (case-insensitive)
        if [[ "${basename,,}" == *protected* && "$force" == false ]]; then
            skipped_protected=$((skipped_protected + 1))
            if [[ "$dry_run" == true ]]; then
                echo "[skip protected] $item"
            fi
            continue
        fi

        if [[ "$dry_run" == true ]]; then
            echo "[would delete] $item"
        else
            if [[ -d "$item" ]]; then
                rm -rf "$item"
            else
                rm -f "$item"
            fi
        fi
        deleted=$((deleted + 1))
    done < <(find "${find_args[@]}" 2>/dev/null || true)
done

# Summary
if [[ "$dry_run" == true ]]; then
    echo ""
    echo "Dry run: $deleted item(s) would be deleted, $skipped_protected protected item(s) skipped"
else
    echo "$deleted item(s) deleted, $skipped_protected protected item(s) skipped"
fi
