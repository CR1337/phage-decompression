#!/bin/bash -ux
#SBATCH --job-name=decompress_genomes
#SBATCH --output=logs/decompress_%A_%a.out
#SBATCH --error=logs/decompress_%A_%a.err
#SBATCH --account=sci-renard-student
#SBATCH --array=0-9999
#SBATCH --ntasks=1
#SBATCH --time=00:50:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --partition=cpu
#SBATCH --nice=5000 

DIR="$1"
OFFSET="$2"

# Activate the virtual environment
eval "$(~/miniforge3/bin/conda shell.bash hook)"

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
python3 phagetools.py glm2 -i "${DIR}/${NAME}.fasta" -g "${DIR}/${NAME}/${GFF_FILE}" -o "${DIR}/${NAME}/${NAME}_original.seq"
python3 phagetools.py glm2 -i "${DIR}/${NAME}.fasta" -g "${DIR}/${NAME}/${GFF_FILE}" -r -o "${DIR}/${NAME}/${NAME}_decompressed.seq"
