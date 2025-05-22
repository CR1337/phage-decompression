#!/bin/bash
#SBATCH --job-name=decompress_genomes
#SBATCH --output=logs/decompress_%A_%a.out
#SBATCH --error=logs/decompress_%A_%a.err
#SBATCH --array=0-0  # overridden by caller
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=512M
#SBATCH --partition=cpu
#SBATCH --nice=5000

# Load input directory from the environment
DIR="$DIR"
ENV_PATH="$ENV_PATH"

# Get the list of FASTA files and select the one corresponding to the array task ID
FASTA_FILE=$(ls "$DIR"/*.fasta | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")

# Extract base name (e.g., "abc" from "abc.fasta")
BASENAME=$(basename "$FASTA_FILE" .fasta)

# Construct path to corresponding GFF file (e.g., DIR/abc/abc.gff)
SUBDIR="$DIR/$BASENAME"
GFF_FILE="$SUBDIR/$BASENAME.gff"

# Optional: Check if files exist
if [[ ! -f "$FASTA_FILE" ]]; then
    echo "FASTA file not found: $FASTA_FILE"
    exit 1
fi

if [[ ! -f "$GFF_FILE" ]]; then
    echo "GFF file not found: $GFF_FILE"
    exit 1
fi

# Activate virtual environment
source "$ENV_PATH/bin/activate"

# Run the Python script
python3 phagetools.py decompress -i "$FASTA_FILE" -g "$GFF_FILE" -l processed_files.txt -o "$SUBDIR"
