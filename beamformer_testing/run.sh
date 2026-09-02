#!/bin/bash
set -euo pipefail

# ---- Defaults ----
THREADS=4
BINS=64

usage() {
    cat <<EOF
Usage: $0 --data-dir DIR --output-dir DIR --par PSR --start-obsid ID --end-obsid ID \\
          --beam BEAM --start-chan N --end-chan N [--bins N] [--threads N]

  --data-dir     Path to the MWAX Beamformer VDIF files
  --output-dir   Path to output directory for .ar and .png files
  --par          Full path to the pulsar parameter file, e.g. "/path/to/J1752-2806.par"
  --start-obsid  Observation ID to process (start)
  --end-obsid    Observation ID to process (end)
  --beam         Beam number to process - normally 01,02,03... etc
  --start-chan   First receiver channel number to process, e.g. 109
  --end-chan     Last receiver channel number to process, e.g. 132
  --threads      Number of threads to use (default: $THREADS)
  --bins         Number of phase bins to use (default: $BINS)
  
  -h, --help     Show this help message
EOF
    exit 1
}

# ---- Parse arguments ----
while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir)    DATA="$2";           shift 2 ;;
        --output-dir)  OUTPUT="$2";         shift 2 ;;
        --par)         PSR="$2";            shift 2 ;;
        --start-obsid) OBSID_START="$2";    shift 2 ;;
        --end-obsid)   OBSID_END="$2";      shift 2 ;;
        --beam)        BEAM="$2";           shift 2 ;;
        --start-chan)  CHAN_START="$2";     shift 2 ;;
        --end-chan)    CHAN_END="$2";       shift 2 ;;
        --bins)        BINS="$2";           shift 2 ;;
        --threads)     THREADS="$2";        shift 2 ;;
        -h|--help)     usage ;;
        *)
            echo "Error: unknown argument '$1'" >&2
            usage
            ;;
    esac
done

# ---- Validation helpers ----
is_int() {
    [[ "$1" =~ ^[0-9]+$ ]]
}

fail() {
    echo "Error: $1" >&2
    exit 1
}

# ---- Required argument checks ----
: "${DATA:?Missing required argument: --data-dir}"
: "${OUTPUT:?Missing required argument: --output-dir}"
: "${PSR:?Missing required argument: --par}"
: "${OBSID_START:?Missing required argument: --start-obsid}"
: "${OBSID_END:?Missing required argument: --end-obsid}"
: "${BEAM:?Missing required argument: --beam}"
: "${CHAN_START:?Missing required argument: --start-chan}"
: "${CHAN_END:?Missing required argument: --end-chan}"

# ---- Type / value validation ----
[[ -d "$DATA" ]] || fail "--data-dir '$DATA' does not exist or is not a directory"
[[ -d "$OUTPUT" ]] || fail "--output-dir '$OUTPUT' does not exist or is not a directory"
[[ -f "$PSR" ]] || fail "--par '$PSR' does not exist or is not a file"

is_int "$OBSID_START" || fail "--start-obsid '$OBSID_START' must be an integer"
is_int "$OBSID_END"   || fail "--end-obsid '$OBSID_END' must be an integer"
(( OBSID_START <= OBSID_END )) || fail "--start-obsid ($OBSID_START) must be <= --end-obsid ($OBSID_END)"

[[ "$BEAM" =~ ^[0-9]{1,3}$ ]] || fail "--beam '$BEAM' must be numeric (e.g. 01, 02, 12)"

is_int "$CHAN_START" || fail "--start-chan '$CHAN_START' must be an integer"
is_int "$CHAN_END"   || fail "--end-chan '$CHAN_END' must be an integer"
(( CHAN_START <= CHAN_END )) || fail "--start-chan ($CHAN_START) must be <= --end-chan ($CHAN_END)"

is_int "$BINS"    || fail "--bins '$BINS' must be a positive integer"
is_int "$THREADS" || fail "--threads '$THREADS' must be a positive integer"
(( BINS > 0 ))    || fail "--bins must be greater than 0"
(( THREADS > 0 )) || fail "--threads must be greater than 0"

[[ -x ./fold_and_plot.sh ]] || fail "./fold_and_plot.sh not found or not executable in the current directory"

# ---- vdif_cat setup / validation ----
VDIF_CAT_DIR="$HOME/mwax_mover"
VDIF_CAT_PY="$VDIF_CAT_DIR/.venv/bin/python3"
VDIF_CAT_SCRIPT="$VDIF_CAT_DIR/src/mwax_mover/cli/vdif_cat.py"

[[ -x "$VDIF_CAT_PY" ]] || fail "venv python not found or not executable: $VDIF_CAT_PY"
[[ -f "$VDIF_CAT_SCRIPT" ]] || fail "vdif_cat script not found: $VDIF_CAT_SCRIPT"

# ---- Download the metafits files and run vdif_cat ----
METADATA_URL="http://ws.mwatelescope.org/metadata/fits?obs_id="

echo "Scanning $DATA for .vdif obsids between $OBSID_START and $OBSID_END..."

# Extract the first 10 digits of every .vdif filename in DATA, keep unique numeric values
FOUND_OBSIDS=()
while IFS= read -r obsid; do
    FOUND_OBSIDS+=("$obsid")
done < <(
    find "$DATA" -maxdepth 1 -type f -name '*.vdif' -printf '%f\n' \
        | grep -oE '^[0-9]{10}' \
        | sort -un
)

if (( ${#FOUND_OBSIDS[@]} == 0 )); then
    fail "No .vdif filenames starting with a 10-digit obsid were found in $DATA"
fi

for obsid in "${FOUND_OBSIDS[@]:-}"; do
    [[ -z "$obsid" ]] && continue
    if (( obsid < OBSID_START || obsid > OBSID_END )); then
        continue
    fi

    METAFITS="$DATA/${obsid}_metafits.fits"
    if [[ -f "$METAFITS" ]]; then
        echo "Metafits already exists for $obsid, skipping download and vdif_cat."
        continue
    fi

    echo "Downloading metafits for obsid $obsid..."
    if ! curl -sSL --fail -o "$METAFITS" "${METADATA_URL}${obsid}"; then
        rm -f "$METAFITS"
        fail "Failed to download metafits for obsid $obsid from ${METADATA_URL}${obsid}"
    fi

    echo "Running vdif_cat for obsid $obsid..."
    if ! "$VDIF_CAT_PY" "$VDIF_CAT_SCRIPT" -m "$METAFITS" -i "$DATA" -o "$DATA"; then
        fail "vdif_cat failed for obsid $obsid"
    fi
done

# ---- Main ----
for chan in $(seq "$CHAN_START" "$CHAN_END"); do
    ./fold_and_plot.sh \
        --data-dir "$DATA" \
        --output-dir "$OUTPUT" \
        --start-obsid "$OBSID_START" \
        --end-obsid "$OBSID_END" \
        --par "$PSR" \
        --beam "$BEAM" \
        --chan "$chan" \
        --threads "$THREADS" \
        --bins "$BINS"
done

echo "Done"