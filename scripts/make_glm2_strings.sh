#!/bin/bash

INPUT_DIR="$1"
NUM_FILES=$(ls "$INPUT_DIR"/*.fasta | wc -l)

sbatch --export=DIR="$INPUT_DIR",ENV_PATH=".Svenv" --array=0-$(($NUM_FILES - 1)) scripts/_make_glm2_strings_array.sh
