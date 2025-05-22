#!/bin/bash

INPUT_DIR="$1"

sbatch --export=DIR="$INPUT_DIR",ENV_PATH=".venv" scripts/_make_glm2_strings_array.sh
