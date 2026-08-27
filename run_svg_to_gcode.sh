#!/usr/bin/env bash
#
# Convert an SVG file into G-code fit inside a fixed-size target box, without
# hand-computing --scale.
#
# Usage:
#   ./run_svg_to_gcode.sh drawing.svg
#
# Edit FIT_XMIN/FIT_YMIN/FIT_XMAX/FIT_YMAX and OUTPUT_PATH below to change the
# target box and output location. Set CLOSED_LOOP=true to append a return
# pass (--closed-loop); RETURN_FEED optionally sets its feed rate.
#
set -euo pipefail

# Resolve the project root (directory of this script) so it runs from anywhere.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PYTHON="${PYTHON:-python}"

FIT_XMIN=50
FIT_YMIN=50
FIT_XMAX=150
FIT_YMAX=150
OUTPUT_PATH="output/example_gcode_svg.gcode"
CLOSED_LOOP=false
RETURN_FEED=""

SVG_PATH="${1:-}"
if [[ -z "$SVG_PATH" ]]; then
    echo "error: missing required argument: path to source SVG file" >&2
    echo "usage: $0 <drawing.svg>" >&2
    exit 1
fi
if [[ ! -f "$SVG_PATH" ]]; then
    echo "error: no such file: $SVG_PATH" >&2
    exit 1
fi

EXTRA_ARGS=()
if [[ "$CLOSED_LOOP" == "true" ]]; then
    EXTRA_ARGS+=(--closed-loop)
    if [[ -n "$RETURN_FEED" ]]; then
        EXTRA_ARGS+=(--return-feed "$RETURN_FEED")
    fi
fi

# No --home or --z-offset-mm: this script never emits a homing or Z move.
"$PYTHON" svg_to_gcode.py "$SVG_PATH" -o "$OUTPUT_PATH" \
    --fit-box-mm "$FIT_XMIN" "$FIT_YMIN" "$FIT_XMAX" "$FIT_YMAX" \
    ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}
