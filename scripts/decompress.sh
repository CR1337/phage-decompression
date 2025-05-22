#!/bin/bash

# Directory containing *.fasta files
INPUT_DIR="$1"

NUM_FILES=$(ls "$INPUT_DIR"/*.fasta | wc -l)

sbatch --export=DIR="$INPUT_DIR",ENV_PATH=".venv" --array=0-$(($NUM_FILES - 1)) _decompress_array.sh
