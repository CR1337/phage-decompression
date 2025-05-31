
# Phage Genome Toolset

A command-line utility for decompressing, analyzing, transforming, and comparing phage genome sequences. Supports GenBank and FASTA formats with optional GFF annotation. This tool is designed to streamline genomic workflows, generate glm2-compatible output, visualize sequence statistics, and compute sequence similarities.

---

## 📦 Installation

Ensure Python 3.8+ is installed, and install required dependencies:

```bash
pip install -r requirements.txt
```

---

## 🚀 Usage Overview

```bash
python main.py <command> [options]
```

Available subcommands:

- [`decompress`](#decompress)
- [`summary`](#summary)
- [`glm2`](#glm2)
- [`similarity`](#similarity)

Use `-h` or `--help` with any subcommand for full usage instructions.

---

## 🔍 Subcommands

### 🧬 `decompress`

#### 1. What it's for
This command expands a compressed and feature-encoded genome into its full nucleotide sequence.

#### 2. What it does
- Loads a genome from either GenBank or FASTA + GFF.
- Decompresses annotated regions by deoverlapping them.
- Stores output as a FASTA and GFF file.
- Computes statistics on compression ratio, gene count, and overlaps.
- Outputs a JSON file with these metrics.
- Optionally appends the result file path to a list file for batch processing.

#### 3. How to use it

```bash
python main.py decompress \
    -i path/to/genome.fasta \
    -g path/to/genome.gff \
    -o output_directory/ \
    --input-format fasta \
    -l list_of_results.txt
```

If you're using a GenBank file, GFF input is not required:

```bash
python main.py decompress -i genome.gb -o output_directory/ --input-format genbank
```

---

### 📊 `summary`

#### 1. What it's for
This command summarizes a collection of decompressed genome results, producing descriptive statistics and visualization plots.

#### 2. What it does
- Reads a list of JSON files produced by the `decompress` command.
- Aggregates values like compression ratios, sequence lengths, gene counts, etc.
- Computes statistics: mean, median, standard deviation, percentiles, etc.
- Generates histograms (both linear and log-scale) for each metric.
- Outputs a CSV of raw data and a CSV summary.

#### 3. How to use it

```bash
python main.py summary \
    -l path/to/list_of_json_files.txt \
    -o output_directory/
```

The `list_of_json_files.txt` file should contain one JSON path per line.

---

### 🔁 `glm2`

#### 1. What it's for
Generates GLM2-compatible sequence output.

#### 2. What it does
- If `--raw` is used: encodes the raw DNA sequence.
- If `--raw` is not used: extracts features (like CDS regions) and encodes them as Aminoacids.
- Works with both GenBank files or FASTA + GFF combinations.

#### 3. How to use it

From GenBank:

```bash
python main.py glm2 \
    -i path/to/genome.gb \
    -o output.glm2 \
    --input-format genbank
```

From FASTA + GFF:

```bash
python main.py glm2 \
    -i path/to/genome.fasta \
    -g path/to/genome.gff \
    -o output.glm2 \
    --input-format fasta
```

Use raw mode to ignore annotations and convert the whole DNA (eg. for compressed sequences):

```bash
python main.py glm2 -i genome.fasta -o raw_output.glm2 --raw
```

---

### 🧪 `similarity`

#### 1. What it's for
Computes or visualizes similarity between DNA sequences, useful for comparative genomics or clustering related genomes.

#### 2. What it does
Two operation modes:
- **Similarity score mode** (direct comparison): loads two sequences and computes a numerical similarity score.
- **Plot mode**: loads a CSV file with pairwise similarity results and plots a heatmap, optionally with clustering.

#### 3. How to use it

**Compare two sequences:**

```bash
python main.py similarity \
    -i genome1.fasta genome2.gb \
    -o similarity_scores.csv
```

**Plot a similarity heatmap from a CSV file:**

```bash
python main.py similarity \
    -i pairwise_similarities.csv \
    -o similarity_plot.png \
    --plot --cluster
```

---

## 📂 File Format Support

| Format     | Description                        |
|------------|------------------------------------|
| **FASTA**  | DNA sequences; requires GFF        |
| **GenBank**| Contains both sequence and annotations |

---

## 📑 Example Workflow

```bash
# Step 1: Decompress
python main.py decompress -i genome.gb -o decompressed/

# Step 2: Convert to GLM2 format
python main.py glm2 -i genome.gb -o genome.glm2

# Step 3: Compare genomes
python main.py similarity -i genome1.gb genome2.gb -o results.csv

# Step 4: Summarize and visualize results
python main.py summary -l decompressed/list.txt -o analysis/
```

---

## 🧩 Dependencies

- Python 3.8+
- numpy
- pandas
- matplotlib
- Biopython

Install them using:

```bash
pip install -r requirements.txt
```

---
