#!/bin/bash

DIR="$1"
ENV_DIR="$2"
ARRAY_SIZE=10000

# Count total fasta files
TOTAL=$(ls "${DIR}"/*.fasta | wc -l)

# Number of full sbatch submissions
BATCHES=$(( (TOTAL + ARRAY_SIZE - 1) / ARRAY_SIZE ))

for (( i=0; i<$BATCHES; i++ )); do
    OFFSET=$(( i * ARRAY_SIZE ))
    ARRAY_END=$(( TOTAL - OFFSET - 1 ))
    if [ $ARRAY_END -gt 9999 ]; then
        ARRAY_END=9999
    fi
    sbatch --array=0-${ARRAY_END} scripts/_decompress_array.sh "$DIR" "$ENV_DIR" "$OFFSET"
done
