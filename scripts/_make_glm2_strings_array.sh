#!/bin/bash -ux
#SBATCH --job-name=glm2_strings
#SBATCH --output=logs/glm2_strings_%A_%a.out
#SBATCH --error=logs/glm2_strings_%A_%a.err
#SBATCH --array=0-30000
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=512M
#SBATCH --partition=cpu
#SBATCH --nice=5000

# Load input directory and virtualenv path
DIR="$DIR"
ENV_PATH="$ENV_PATH"

# Get the list of original FASTA files
FASTA_FILE=$(ls "$DIR"/*.fasta | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")

# Get the base name (e.g., "abc" from "abc.fasta")
BASENAME=$(basename "$FASTA_FILE" .fasta)
SUBDIR="$DIR/$BASENAME"
GFF_FILE="$SUBDIR/$BASENAME.gff"

# Check required files
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

# Output path for first run
OUT_FILE1="$SUBDIR/${BASENAME}.seq"

# Run first script
python3 phagetools.py glm2 -i "$FASTA_FILE" -g "$GFF_FILE" -l processed_files.txt -o "$OUT_FILE1"

# Construct decompressed filename
DECOMPRESSED_FASTA="${FASTA_FILE%.fasta}_decompressed.fasta"
DECOMPRESSED_BASENAME=$(basename "$DECOMPRESSED_FASTA" .fasta)

# Update GFF file path if needed — assuming GFF is still the same
OUT_FILE2="$SUBDIR/${DECOMPRESSED_BASENAME}.seq"

# Run second script with decompressed FASTA
if [[ -f "$DECOMPRESSED_FASTA" ]]; then
    python3 phagetools.py glm2 -i "$DECOMPRESSED_FASTA" -g "$GFF_FILE" -l processed_files.txt -o "$OUT_FILE2"
else
    echo "Decompressed FASTA not found: $DECOMPRESSED_FASTA"
    exit 1
fi
