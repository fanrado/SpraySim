
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
MATERIAL="${MATERIAL:-water}"
DENSITY="${DENSITY:-}"
VISCOSITY="${VISCOSITY:-}"
SOLIDS_FRACTION="${SOLIDS_FRACTION:-1.0}"
DRAG_MODEL="${DRAG_MODEL:-clift_gauvin}"
PRESSURE_BAR="${PRESSURE_BAR:-3.0}"
ORIFICE_MM="${ORIFICE_MM:-0.8}"
NOZZLE_SHAPE="${NOZZLE_SHAPE:-full_cone}"
SPRAY_DURATION="${SPRAY_DURATION:-0.15}"
DISTRIBUTION="${DISTRIBUTION:-lognormal}"
MEAN_RADIUS_MM="${MEAN_RADIUS_MM:-0.4}"
RADIUS_STD_MM="${RADIUS_STD_MM:-0.12}"
CONE="${CONE:-25.0}"
HEIGHT="${HEIGHT:-1.5}"
SPEED_SPREAD="${SPEED_SPREAD:-0.15}"
DT="${DT:-0.001}"
SEED="${SEED:-42}"
OUT="${OUT:-output/spray_summary.png}"
NO_PLOT="${NO_PLOT:-false}"
DATA="${DATA:-output/spray_data.npz}"
NO_DATA="${NO_DATA:-false}"
DROPLETS="${DROPLETS:-}"

# Assemble the run.py invocation.
CMD=("$PYTHON" run.py
    --material "$MATERIAL"
    --solids-fraction "$SOLIDS_FRACTION"
    --pressure-bar "$PRESSURE_BAR"
    --orifice-mm "$ORIFICE_MM"
    --shape "$NOZZLE_SHAPE"
    --spray-duration "$SPRAY_DURATION"
    --distribution "$DISTRIBUTION"
    --mean-radius-mm "$MEAN_RADIUS_MM"
    --radius-std-mm "$RADIUS_STD_MM"
    --drag-model "$DRAG_MODEL"
    --cone "$CONE"
    --height "$HEIGHT"
    --speed-spread "$SPEED_SPREAD"
    --dt "$DT"
    --seed "$SEED"
    --out "$OUT"
    --data "$DATA")

# Only override the material's default density if the config set one.
if [[ -n "$DENSITY" ]]; then
    CMD+=(--density "$DENSITY")
fi

# Only override the material's default viscosity if the config set one.
if [[ -n "$VISCOSITY" ]]; then
    CMD+=(--viscosity "$VISCOSITY")
fi

# Only pin an explicit droplet count if the config set one.
if [[ -n "$DROPLETS" ]]; then
    CMD+=(--droplets "$DROPLETS")
fi

if [[ "$NO_PLOT" == "true" ]]; then
    CMD+=(--no-plot)
fi

if [[ "$NO_DATA" == "true" ]]; then
    CMD+=(--no-data)
fi

# Pass through any extra CLI args (e.g. ./main.sh default --no-plot).
if [[ $# -gt 1 ]]; then
    CMD+=("${@:2}")
fi

echo "Running: ${CMD[*]}"
echo
exec "${CMD[@]}"
