#!/bin/bash

NUM_TILES=2
TONE1=4800
TONE2=5311
TONE3=1458
NOISE_LEVEL=0.1
NBITS=8
FFT_SIZE=6400
SUBOB_LENGTH=8.0
HEADER_SIZE=0
HEADER_FILENAME="min_v2_header.txt"
SUBFILE_NAME="test_subfile_2T.sub"

# generate all the delayed tones files, both pols the same
# - at the same time construct a filenames file and a delays file
# - every tile/pol has a different delay, but the same tones and noise level (different noise seed for each)
echo
echo "Calling make_delayed_tones_plus_noise to generate test signal files for $NUM_TILES tiles..."
echo
NUM_SOURCES=$((NUM_TILES * 2))   # 2 pols
DELAY_STRIDE=$(echo "1.0 / $NUM_SOURCES" | bc -l)  # fractions of a sample - this allows all sources to have a different delay within [0.0..1.0]

BASENAME="${HEADER_SIZE}_$((2*FFT_SIZE))_${NBITS}_tones_${TONE1}_${TONE2}_${TONE3}_${NOISE_LEVEL}_length_${SUBOB_LENGTH}_D"
EXTENSION=".dat"

rm -f subfile_source_files.txt
rm -f subfile_source_delays.txt
rm -f "${BASENAME}"*"$EXTENSION"

for ((i=0; i<${NUM_SOURCES}; i++)); do
    DELAY=$(echo "$i * ${DELAY_STRIDE}" | bc -l)
    echo "$DELAY" >> subfile_source_delays.txt
    OUTFILE="${BASENAME}_${DELAY}${EXTENSION}"
    echo "$OUTFILE" >> subfile_source_files.txt
    ./make_delayed_tones_plus_noise -H ${HEADER_SIZE} -d $((2*FFT_SIZE)) -b ${NBITS} -D ${DELAY} -1 ${TONE1} -2 ${TONE2} -3 ${TONE3} -n ${NOISE_LEVEL} -s $(( i + 1 )) -l ${SUBOB_LENGTH} -o ${OUTFILE}
done

# build the test subfile, reading in filenames and delays from file, and the PSRDADA header
#./make_test_subfile -H 0 -d 12800 -b 8 -f filenames.txt -D delays.txt -o test_subfile.sub
echo
echo "Calling make_test_subfile to assemble a test subfile for $NUM_TILES tiles..."
echo
./make_test_subfile -T ${NUM_TILES} -H ${HEADER_FILENAME} -s subfile_source_files.txt -D subfile_source_delays.txt -o ${SUBFILE_NAME}
echo
