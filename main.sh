
#!/usr/bin/env bash
#
# SpraySim launcher.
#
# Loads a config file from config/ and runs the simulation with those values,
# so you never have to pass parameters on the command line.
#
# Usage:
#   ./main.sh                       # uses config/default.conf
#   ./main.sh fine_mist             # uses config/fine_mist.conf
#   ./main.sh config/big_drops.conf # explicit path also works
#   ./main.sh --list                # list available configs
#
set -euo pipefail

# Resolve the project root (directory of this script) so it runs from anywhere.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CONFIG_DIR="config"
PYTHON="${PYTHON:-python3}"

if [[ "${1:-}" == "--list" || "${1:-}" == "-l" ]]; then
    echo "Available configs in ${CONFIG_DIR}/:"
    for f in "${CONFIG_DIR}"/*.conf; do
        [[ -e "$f" ]] && echo "  - $(basename "${f%.conf}")"
    done
    exit 0
fi

# Pick the config: default, a bare name (fine_mist), or an explicit path.
ARG="${1:-default}"
if [[ -f "$ARG" ]]; then
    CONFIG="$ARG"
elif [[ -f "${CONFIG_DIR}/${ARG}.conf" ]]; then
    CONFIG="${CONFIG_DIR}/${ARG}.conf"
elif [[ -f "${CONFIG_DIR}/${ARG}" ]]; then
    CONFIG="${CONFIG_DIR}/${ARG}"
else
    echo "Error: config '${ARG}' not found." >&2
    echo "Try: ./main.sh --list" >&2
    exit 1
fi

echo "Loading config: ${CONFIG}"
# shellcheck disable=SC1090
source "$CONFIG"

# Fall back to sensible defaults for anything the config didn't set.
DROPLETS="${DROPLETS:-4000}"
SPEED="${SPEED:-9.0}"
CONE="${CONE:-25.0}"
HEIGHT="${HEIGHT:-1.5}"
RADIUS_MM="${RADIUS_MM:-0.4}"
DT="${DT:-0.001}"
SEED="${SEED:-42}"
OUT="${OUT:-output/spray_summary.png}"
NO_PLOT="${NO_PLOT:-false}"

# Assemble the run.py invocation.
CMD=("$PYTHON" run.py
    --droplets "$DROPLETS"
    --speed "$SPEED"
    --cone "$CONE"
    --height "$HEIGHT"
    --radius-mm "$RADIUS_MM"
    --dt "$DT"
    --seed "$SEED"
    --out "$OUT")

if [[ "$NO_PLOT" == "true" ]]; then
    CMD+=(--no-plot)
fi

# Pass through any extra CLI args (e.g. ./main.sh default --no-plot).
if [[ $# -gt 1 ]]; then
    CMD+=("${@:2}")
fi

echo "Running: ${CMD[*]}"
echo
exec "${CMD[@]}"
