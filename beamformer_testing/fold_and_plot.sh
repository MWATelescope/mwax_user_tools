#!/bin/bash

# Usage: ./fold_and_plot.sh --data-dir <dir> --output-dir <dir> --start-obsid <id> --end-obsid <id> --par <file> --beam <num> --chan <num> [--threads <num>] [--bins <num>]
# Example: ./fold_and_plot.sh --data-dir /path/to/data --output-dir /path/to/output --start-obsid 1234567890 --end-obsid 1234567899 --par pulsar.par --beam 01 --chan 109 --threads 4 --bins 256

# --- Defaults ---
THREADS=2
BINS=128

# --- Parse named arguments ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir)    HOST_DIR="$2";    shift 2 ;;
        --output-dir)  OUTPUT_DIR="$2";  shift 2 ;;
        --start-obsid) START_OBSID="$2"; shift 2 ;;
        --end-obsid)   END_OBSID="$2";   shift 2 ;;
        --par)         PAR_PATH="$2";    shift 2 ;;
        --beam)        BEAM="$2";        shift 2 ;;
        --chan)         CHANNEL="ch$2";   shift 2 ;;
        --threads)     THREADS="$2";     shift 2 ;;
        --bins)        BINS="$2";        shift 2 ;;
        --help)
            echo "Usage: $0 --data-dir <dir> --output-dir <dir> --start-obsid <id> --end-obsid <id> --par <file> --beam <num> --chan <num> [--threads <num>] [--bins <num>]"
            echo ""
            echo "  --data-dir      Directory containing VDIF, HDR and PAR files"
            echo "  --output-dir    Directory to write output .ar and .png files"
            echo "  --start-obsid   10-digit start observation ID (inclusive)"
            echo "  --end-obsid     10-digit end observation ID (inclusive)"
            echo "  --par           PAR filename (no path)"
            echo "  --beam          Zero-padded beam number (e.g. 01, 02)"
            echo "  --chan          3-digit receiver channel number (e.g. 109)"
            echo "  --threads       Number of dspsr threads (default: 2)"
            echo "  --bins          Number of profile bins for dspsr (default: 128)"
            exit 0
            ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

PAR_DIR="$(dirname "$PAR_PATH")"
PAR_FILE="$(basename "$PAR_PATH")"
PAR_BASE="${PAR_FILE%.*}"

# --- Validate inputs ---
MISSING=()
[[ -z "$HOST_DIR" ]]    && MISSING+=("--data-dir")
[[ -z "$OUTPUT_DIR" ]]  && MISSING+=("--output-dir")
[[ -z "$START_OBSID" ]] && MISSING+=("--start-obsid")
[[ -z "$END_OBSID" ]]   && MISSING+=("--end-obsid")
[[ -z "$PAR_PATH" ]]    && MISSING+=("--par")
[[ -z "$BEAM" ]]        && MISSING+=("--beam")
[[ -z "$CHANNEL" ]]     && MISSING+=("--chan")

if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "Error: Missing required arguments: ${MISSING[*]}"
    echo "Run '$0 --help' for usage."
    exit 1
fi

if [[ ! -f "$PAR_PATH" ]]; then
    echo "Error: PAR file '$PAR_PATH' does not exist."
    exit 1
fi

if [[ ! -d "$HOST_DIR" ]]; then
    echo "Error: Data directory '$HOST_DIR' does not exist."
    exit 1
fi

if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "Error: Output directory '$OUTPUT_DIR' does not exist."
    exit 1
fi

if [[ ! "$START_OBSID" =~ ^[0-9]{10}$ ]]; then
    echo "Error: --start-obsid must be a 10-digit integer (got '$START_OBSID')."
    exit 1
fi

if [[ ! "$END_OBSID" =~ ^[0-9]{10}$ ]]; then
    echo "Error: --end-obsid must be a 10-digit integer (got '$END_OBSID')."
    exit 1
fi

if (( 10#$END_OBSID < 10#$START_OBSID )); then
    echo "Error: --end-obsid ($END_OBSID) must be greater than or equal to --start-obsid ($START_OBSID)."
    exit 1
fi

# --- Find matching HDR files ---
HDR_FILES=()
for f in "$HOST_DIR"/*"${CHANNEL}"*beam"${BEAM}"*.hdr; do
    [[ -e "$f" ]] || continue
    fname="$(basename "$f")"
    obsid="${fname:0:10}"
    if [[ "$obsid" =~ ^[0-9]{10}$ ]] && (( 10#$obsid >= 10#$START_OBSID )) && (( 10#$obsid <= 10#$END_OBSID )); then
        HDR_FILES+=("$f")
    fi
done

NUM_HDR=${#HDR_FILES[@]}

if [[ $NUM_HDR -eq 0 ]]; then
    echo "Error: No .hdr files found in '$HOST_DIR' matching channel '$CHANNEL', beam '$BEAM', and obsid range $START_OBSID-$END_OBSID."
    exit 1
fi

echo "Found $NUM_HDR HDR file(s) to process."

# --- Process each HDR file ---
FIRST=true
AR_FILES=()

for HDR_PATH in "${HDR_FILES[@]}"; do
    HDR_FILENAME="$(basename "$HDR_PATH")"
    OBSID="${HDR_FILENAME:0:10}"
    AR_FILES+=("/output/${OBSID}_${CHANNEL}_beam${BEAM}.ar")

    echo "Processing: $HDR_FILENAME (OBSID: $OBSID)"

    if [[ "$FIRST" == true ]]; then
        docker run -it --rm --user "$(id -u):$(id -g)" --entrypoint dspsr \
            -v "$HOST_DIR":/data \
            -v "$OUTPUT_DIR":/output \
            -v "$PAR_DIR":/par \
            cirapulsarsandtransients/psr-analysis:latest \
            -t "$THREADS" -A -L 10 -F 32:D -b "$BINS" \
            -S 4 \
            -E /par/"$PAR_FILE" \
            -O /output/"${OBSID}_${CHANNEL}_beam${BEAM}" \
            /data/"$HDR_FILENAME"
        FIRST=false
    else
        docker run -it --rm --user "$(id -u):$(id -g)" --entrypoint dspsr \
            -v "$HOST_DIR":/data \
            -v "$OUTPUT_DIR":/output \
            -v "$PAR_DIR":/par \
            cirapulsarsandtransients/psr-analysis:latest \
            -t "$THREADS" -A -L 10 -F 32:D -b "$BINS" \
            -E /par/"$PAR_FILE" \
            -O /output/"${OBSID}_${CHANNEL}_beam${BEAM}" \
            /data/"$HDR_FILENAME"
    fi

    if [[ $? -ne 0 ]]; then
        echo "Error: dspsr failed for $HDR_FILENAME. Aborting."
        exit 1
    fi
done

# --- Combine .ar files (skip if only one HDR was processed) ---
if [[ $NUM_HDR -gt 1 ]]; then
    echo "Combining .ar files for ch${CHANNEL} beam $BEAM..."

    docker run -it --rm --user "$(id -u):$(id -g)" --entrypoint psradd \
        -v "$OUTPUT_DIR":/output \
        cirapulsarsandtransients/psr-analysis:latest \
        -o /output/${PAR_BASE}_${CHANNEL}_beam${BEAM}_combined.ar "${AR_FILES[@]}"

    if [[ $? -ne 0 ]]; then
        echo "Error: psradd failed. Aborting."
        exit 1
    fi
else
    echo "Only one HDR file processed, skipping psradd."
    SINGLE_OBSID="$(basename "${HDR_FILES[0]}")"
    SINGLE_OBSID="${SINGLE_OBSID:0:10}"
    cp "$OUTPUT_DIR/${SINGLE_OBSID}_${CHANNEL}_beam${BEAM}.ar" \
       "$OUTPUT_DIR/${PAR_BASE}_${CHANNEL}_beam${BEAM}_combined.ar"
fi

# --- Generate plots ---
echo "Generating profile plot for ch${CHANNEL} beam $BEAM..."

docker run -it --rm --user "$(id -u):$(id -g)" --entrypoint pav \
    -v "$OUTPUT_DIR":/output \
    cirapulsarsandtransients/psr-analysis:latest \
    -DFTp /output/${PAR_BASE}_${CHANNEL}_beam${BEAM}_combined.ar -g /output/${PAR_BASE}_${CHANNEL}_beam${BEAM}_combined_profile.png/png

if [[ $? -ne 0 ]]; then
    echo "Error: pav failed."
    exit 1
fi

echo "Generating waterfall plot for beam ch${CHANNEL} $BEAM..."

docker run -it --rm --user "$(id -u):$(id -g)" --entrypoint pav \
    -v "$OUTPUT_DIR":/output \
    cirapulsarsandtransients/psr-analysis:latest \
    -GTp /output/${PAR_BASE}_${CHANNEL}_beam${BEAM}_combined.ar -g /output/${PAR_BASE}_${CHANNEL}_beam${BEAM}_combined_waterfall.png/png

if [[ $? -ne 0 ]]; then
    echo "Error: pav failed."
    exit 1
fi

echo "Done! Output files:"
echo "  Combined archive : $OUTPUT_DIR/beam${BEAM}_${PAR_BASE}_combined.ar"
echo "  Profile plot     : $OUTPUT_DIR/beam${BEAM}_${PAR_BASE}_profile.png"
echo "  Waterfall plot   : $OUTPUT_DIR/beam${BEAM}_${PAR_BASE}_waterfall.png"