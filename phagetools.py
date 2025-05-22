import os
import json
from typing import List

from modules.arguments import parse_arguments
from modules.file_io import load_annotated_sequence_from_fasta, load_annotated_sequence_from_genbank, store_annotated_sequence_as_fasta, load_sequence_from_fasta, load_sequence_from_genbank
from modules.decompress import decompress as decompress_genome, count_overlapping_cds
from modules.glm2 import create_glm2_string_from_dna, create_glm2_string_from_features
from modules.similarity import compute_similarity, plot_similarities

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def decompress(
    input_format: str, 
    input_filename: str,
    input_gff_filename: str | None,
    file_list_filename: str | None,
    output_directory: str
):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    if input_format == 'fasta':
        sequence = load_annotated_sequence_from_fasta(input_filename, input_gff_filename)
    elif input_format == 'genbank':
        sequence = load_annotated_sequence_from_genbank(input_filename)

    decompressed_sequence = decompress_genome(sequence)
    name = os.path.basename(input_filename).split(".")[0]
    
    fasta_filename = os.path.join(output_directory, f"{name}_decompressed.fasta")
    gff_filename = os.path.join(output_directory, f"{name}_decompressed.gff")
    store_annotated_sequence_as_fasta(decompressed_sequence, fasta_filename, gff_filename)


    compressed_length = len(sequence)
    decompressed_length = len(decompressed_sequence)
    compression_ratio = decompressed_length / compressed_length
    difference = decompressed_length - compressed_length
    n_genes = sum((1 if f.type == 'CDS' else 0) for f in sequence.features)
    n_overlapping_genes = count_overlapping_cds(sequence)

    json_data = {
        'compressed_length': compressed_length,
        'decompressed_length': decompressed_length,
        'compression_ratio': compression_ratio,
        'difference': difference,
        'n_genes': n_genes,
        'n_overlapping_genes': n_overlapping_genes
    }

    json_filename = os.path.join(output_directory, f"{name}_decompression.json")

    with open(json_filename, 'w') as file:
        json.dump(json_data, file)

    if file_list_filename:
        with open(file_list_filename, 'a') as file:
            file.write(f"{json_filename}\n")


def summary(file_list_filename: str, output_directory: str) :
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        
    with open(file_list_filename, 'r') as file:
        file_list = [line.strip() for line in file.readlines()]
    clean_file_list = [line for line in file_list if line]

    columns = [
        'compressed_length', 'decompressed_length', 'compression_ratio', 
        'difference', 'n_genes', 'n_overlapping_genes'
    ]

    data_df = pd.DataFrame(columns=columns)
    for filename in clean_file_list:
        with open(filename, 'r') as file:
            data = json.load(file)
        data_df.loc[len(data_df)] = data

    summary_funcs = {
        'count': data_df.count(numeric_only=True),
        'mean': data_df.mean(numeric_only=True),
        'median': data_df.median(numeric_only=True),
        'stddev': data_df.std(numeric_only=True),
        'min': data_df.min(numeric_only=True),
        'max': data_df.max(numeric_only=True),
        '25 percentile': data_df.quantile(0.25, numeric_only=True),
        '75 percentile': data_df.quantile(0.75, numeric_only=True)
    }

    summary_df = pd.DataFrame(summary_funcs).T
    summary_df = summary_df.reset_index().rename(columns={'index': 'stat'})

    data_filename = os.path.join(output_directory, "data.csv")
    summary_filename = os.path.join(output_directory, "summary.csv")

    data_df.to_csv(data_filename)
    summary_df.to_csv(summary_filename)

    for column in columns:
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        ax.hist(data_df[column])
        ax.set_title(column.replace("_", " ").capitalize())
        ax.set_xlabel(column.replace("_", " "))
        ax.set_ylabel("Count")
        filename = os.path.join(output_directory, f"{column}.png")
        fig.tight_layout()
        fig.savefig(filename)

        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        ax.hist(np.log(list(filter(lambda x: x > 0, data_df[column]))))
        ax.set_title("Log " + column.replace("_", " ").capitalize())
        ax.set_xlabel("log " + column.replace("_", " "))
        ax.set_ylabel("Count")
        filename = os.path.join(output_directory, f"log_{column}.png")
        fig.tight_layout()
        fig.savefig(filename)


def glm2(
    input_format: str, 
    input_filename: str, 
    input_gff_filename: str | None, 
    raw: bool, 
    output_filename: str
):
    if input_format == 'fasta':
        sequence = load_annotated_sequence_from_fasta(input_filename, input_gff_filename)
    elif input_format == 'genbank':
        sequence = load_annotated_sequence_from_genbank(input_filename)

    if raw:
        result = create_glm2_string_from_dna(sequence)
    else:
        result = create_glm2_string_from_features(sequence)

    with open(output_filename, 'w') as file:
        file.write(result)


def similarity(input_filenames: List[str], plot: bool, cluster: bool, output_filename: str):
    if plot:
        df = pd.read_csv(
            input_filenames[0], 
            header=0, 
            names=["filename_1", "filename_2", "similarity"], 
            index_col=False
        )
        plot_similarities(df, output_filename, cluster)

    else:
        sequences = []
        for filename in input_filenames[:2]:
            if filename.endswith(".fasta"):
                sequences.append(load_sequence_from_fasta(filename))
            else:
                sequences.append(load_sequence_from_genbank(filename))
        similarity = compute_similarity(sequences[0].seq, sequences[1].seq)

        with open(output_filename, 'a') as file:
            file.write(f"{input_filenames[0]},{input_filenames[1]},{similarity}\n")


def main():
    arguments = parse_arguments()

    if arguments.command == 'decompress':
        decompress(
            arguments.input_format,
            arguments.input_filename,
            (
                arguments.input_gff_filename 
                if hasattr(arguments, 'input_gff_filename') 
                else None
            ),
            (
                arguments.file_list_filename 
                if hasattr(arguments, 'file_list_filename') 
                else None
            ),
            arguments.output_directory
        )

    elif arguments.command == 'summary':
        summary(
            arguments.file_list_filename,
            arguments.output_directory
        )

    elif arguments.command == 'glm2':
        glm2(
            arguments.input_format,
            arguments.input_filename,
            (
                arguments.input_gff_filename 
                if hasattr(arguments, 'input_gff_filename') 
                else None
            ),
            arguments.raw,
            arguments.output_filename
        )

    elif arguments.command == 'similarity':
        similarity(
            arguments.input_filenames,
            arguments.plot,
            arguments.cluster,
            arguments.output_filename
        )


if __name__ == "__main__":
    main()
