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

[[ "$PSR" =~ ^J[0-9]{4}[+-][0-9]{2,4}$ ]] || fail "--par '$PSR' doesn't look like a valid pulsar name (expected format e.g. J1752-2806)"

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

PARFILE="/voltdata/test_data/par/$PSR.par"
[[ -f "$PARFILE" ]] || fail "par file not found: $PARFILE"

[[ -x ./fold_and_plot.sh ]] || fail "./fold_and_plot.sh not found or not executable in the current directory"

# ---- Main ----
for chan in $(seq "$CHAN_START" "$CHAN_END"); do
    ./fold_and_plot.sh \
        --data-dir "$DATA" \
        --output-dir "$OUTPUT" \
        --start-obsid "$OBSID_START" \
        --end-obsid "$OBSID_END" \
        --par "$PARFILE" \
        --beam "$BEAM" \
        --chan "$chan" \
        --threads "$THREADS" \
        --bins "$BINS"
done

echo "Done"
