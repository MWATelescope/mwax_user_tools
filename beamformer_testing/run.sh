#!/bin/bash
set -euo pipefail

# ---- Defaults ----
BINS=64
THREADS=4

usage() {
    cat <<EOF
Usage: $0 --vdif-dir DIR --psr-name PSR --obsid-start ID --obsid-end ID \\
          --beam BEAM --chan-start N --chan-end N [--bins N] [--threads N]

  --vdif-dir     Path to the MWAX Beamformer VDIF files
  --psr-name     Name of the pulsar to be looked up, e.g. "J1752-2806"
  --obsid-start  Observation ID to process (start)
  --obsid-end    Observation ID to process (end)
  --beam         Beam number to process - normally 01,02,03... etc
  --chan-start   First receiver channel number to process, e.g. 109
  --chan-end     Last receiver channel number to process, e.g. 132
  --bins         Number of phase bins to use (default: $BINS)
  --threads      Number of threads to use (default: $THREADS)
  -h, --help     Show this help message
EOF
    exit 1
}

# ---- Parse arguments ----
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vdif-dir)    DATA="$2";           shift 2 ;;
        --psr-name)    PSR="$2";            shift 2 ;;
        --obsid-start) OBSID_START="$2";    shift 2 ;;
        --obsid-end)   OBSID_END="$2";      shift 2 ;;
        --beam)        BEAM="$2";           shift 2 ;;
        --chan-start)  CHAN_START="$2";     shift 2 ;;
        --chan-end)    CHAN_END="$2";       shift 2 ;;
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
: "${DATA:?Missing required argument: --vdif-dir}"
: "${PSR:?Missing required argument: --psr-name}"
: "${OBSID_START:?Missing required argument: --obsid-start}"
: "${OBSID_END:?Missing required argument: --obsid-end}"
: "${BEAM:?Missing required argument: --beam}"
: "${CHAN_START:?Missing required argument: --chan-start}"
: "${CHAN_END:?Missing required argument: --chan-end}"

# ---- Type / value validation ----
[[ -d "$DATA" ]] || fail "--vdif-dir '$DATA' does not exist or is not a directory"

[[ "$PSR" =~ ^J[0-9]{4}[+-][0-9]{2,4}$ ]] || fail "--psr-name '$PSR' doesn't look like a valid pulsar name (expected format e.g. J1752-2806)"

is_int "$OBSID_START" || fail "--obsid-start '$OBSID_START' must be an integer"
is_int "$OBSID_END"   || fail "--obsid-end '$OBSID_END' must be an integer"
(( OBSID_START <= OBSID_END )) || fail "--obsid-start ($OBSID_START) must be <= --obsid-end ($OBSID_END)"

[[ "$BEAM" =~ ^[0-9]{1,3}$ ]] || fail "--beam '$BEAM' must be numeric (e.g. 01, 02, 12)"

is_int "$CHAN_START" || fail "--chan-start '$CHAN_START' must be an integer"
is_int "$CHAN_END"   || fail "--chan-end '$CHAN_END' must be an integer"
(( CHAN_START <= CHAN_END )) || fail "--chan-start ($CHAN_START) must be <= --chan-end ($CHAN_END)"

is_int "$BINS"    || fail "--bins '$BINS' must be a positive integer"
is_int "$THREADS" || fail "--threads '$THREADS' must be a positive integer"
(( BINS > 0 ))    || fail "--bins must be greater than 0"
(( THREADS > 0 )) || fail "--threads must be greater than 0"

PARFILE="/voltdata/test_data/par/$PSR.par"
[[ -f "$PARFILE" ]] || fail "par file not found: $PARFILE"

[[ -x ./fold_and_plot.sh ]] || fail "./fold_and_plot.sh not found or not executable in the current directory"

# ---- Main ----
OUT="$DATA/$PSR"
mkdir -p "$OUT"

for chan in $(seq "$CHAN_START" "$CHAN_END"); do
    ./fold_and_plot.sh \
        --data-dir "$DATA" \
        --output-dir "$OUT" \
        --start-obsid "$OBSID_START" \
        --end-obsid "$OBSID_END" \
        --par "$PARFILE" \
        --beam "$BEAM" \
        --chan "$chan" \
        --threads "$THREADS" \
        --bins "$BINS"
done

echo "Done"
