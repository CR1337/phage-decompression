import os
import re
from pathlib import Path
from typing import List

from modules.arguments import parse_arguments
from modules.file_io import SequenceFile, load_annotated_sequence, load_sequence, store_decompressed_sequence
from modules.decompression import decompress as decompress_genome
from modules.plotting import plot_decompressed_plasmid, plot_analysis
from modules.analysis import analyze as analyze_genomes
from modules.analysis import summary, make_dataframe

import warnings
warnings.filterwarnings("ignore", message="Biotite does not support FASTA data mixed into GFF files")


def get_sequence_files(names: str, input_directory: str) -> List[SequenceFile]:
    pattern = f"^{names}$"
    regex = re.compile(pattern)
    genome_names = [
        ".".join(f.name.split(".")[:-1]) 
        for f in Path(input_directory).glob("*.fasta")
    ]
    return [
        SequenceFile(name, input_directory) for name in genome_names
        if regex.match(name)
    ]


def decompress(names: str, input_directory: str):
    for file in get_sequence_files(names, input_directory):
        genome, fasta_parameters, gff_parameters = load_annotated_sequence(file)
        decompressed_genome = decompress_genome(genome)
        store_decompressed_sequence(decompressed_genome, file, fasta_parameters, gff_parameters)
        print(file.name)


def plot(names: str, input_directory: str):
    for file in get_sequence_files(names, input_directory):
        genome, *_ = load_annotated_sequence(file)
        decompressed_genome, *_ = load_annotated_sequence(file, decompressed=True)
        plot_filename = file.make_path(f"{file.name}_plasmids.png")
        plot_decompressed_plasmid(
            genome, 
            decompressed_genome, 
            plot_filename,
            feature_format={
                'CDS': (True, 'orange', 'black', lambda f: f.qual.get('function')),
                'gap': (False, 'gray', 'black', lambda _: None)
            }
        )
        print(file.name)


def analyze(names: str, input_directory: str, output_directory: str):
    files = get_sequence_files(names, input_directory)
    names = (file.name for file in files)
    genomes = (load_sequence(file) for file in files)
    decompressed_genomes = (load_sequence(file, decompressed=True) for file in files)
    
    summary_filename = os.path.join(output_directory, "summary.txt")
    plot_filename = os.path.join(output_directory, "analysis.png")
    log_plot_filename = os.path.join(output_directory, "log_analysis.png")

    analyses = [
        analyze_genomes(name, genome, decompressed_genome) 
        for name, genome, decompressed_genome 
        in zip(names, genomes, decompressed_genomes)
    ]

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    plot_analysis(analyses, plot_filename)
    plot_analysis(analyses, log_plot_filename, log_scale=True)

    genomes_summary = summary(analyses)
    with open(summary_filename, 'w') as summary_file:
        summary_file.write(genomes_summary)
    print(genomes_summary)

    for analysis, file in zip(analyses, files):
        filename = file.make_path(f"{file.name}_analysis.json")
        analysis.save_as_json(filename)
        print(file.name)

    df = make_dataframe(analyses)
    filename = os.path.join(output_directory, "analysis.csv")
    df.to_csv(filename, index = False)


def main():
    arguments = parse_arguments()

    if 'decompress' in arguments.operations:
        decompress(arguments.names, arguments.input_directory)

    if 'plot' in arguments.operations:
        plot(arguments.names, arguments.input_directory)
    
    if 'analyze' in arguments.operations:
        analyze(arguments.names, arguments.input_directory, arguments.output_directory)


if __name__ == "__main__":
    main()
