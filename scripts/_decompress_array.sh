#!/bin/bash -ux
#SBATCH --job-name=decompress_genomes
#SBATCH --account=sci-renard-student
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --array=0-9999
#SBATCH --ntasks=1
#SBATCH --time=00:40:00
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

mkdir -p logs
exec > "logs/decompress_${NAME}.out" 2> "logs/decompress_${NAME}.err"

# Create output directory if needed
mkdir -p "${DIR}/${NAME}"

# Run the Python script
python3 phagetools.py decompress -i "${DIR}/${NAME}.fasta" -g "${DIR}/${NAME}/${GFF_FILE}" -l processed_files.txt -o "${DIR}/${NAME}"
