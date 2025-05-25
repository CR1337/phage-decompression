#!/bin/bash -ux
#SBATCH --job-name=decompress_genomes
#SBATCH --output=logs/decompress_%A_%a.out
#SBATCH --error=logs/decompress_%A_%a.err
#SBATCH --array=0-9999
#SBATCH --ntasks=1
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=512M
#SBATCH --partition=cpu
#SBATCH --nice=5000 

DIR="$1"
ENV_DIR="$2"
OFFSET="$3"

# Activate the virtual environment
source "${ENV_DIR}/bin/activate"

# Get all fasta files and select the one for this task
FILES=($(ls "${DIR}"/*.fasta))
INDEX=$(($SLURM_ARRAY_TASK_ID + $OFFSET))
FASTA_FILE="${FILES[$INDEX]}"

# Derive NAME and GFF_FILE from FASTA file
NAME=$(basename "$FASTA_FILE" .fasta)
GFF_FILE="${NAME}.gff"

# Create output directory if needed
mkdir -p "${DIR}/${NAME}"

# Run the Python script
python3 phagetools.py decompress -i "${DIR}/${NAME}.fasta" -g "${DIR}/${NAME}/${GFF_FILE}" -l processed_files.txt -o "${DIR}/${NAME}"
