#!/bin/bash

# Directory containing *.fasta files
INPUT_DIR="$1"

sbatch --export=DIR="$INPUT_DIR",ENV_PATH=".venv" scripts/_decompress_array.sh
