import os
from biotite.sequence import AnnotatedSequence, Feature
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import biotite.sequence.graphics as graphics
from biotite.sequence.io import fasta, gff


GENOME_DIRECTORY: str = "GenomesDB"


def plot_plasmids(genome: AnnotatedSequence, decompressed_genome: AnnotatedSequence, filename: str):
    """
    Plot the plasmid maps of compressed and decompressed genomes.
    """
    def feature_formatter(feature: Feature):
        if feature.key == "CDS":
            return True, "orange", 'black', feature.qual.get('function')
        elif feature.key == "RBS":
            return False, "blue", 'black', None
        elif feature.key == "start_codon":
            return False, "green", 'black', None
        elif feature.key == "stop_codon":
            return False, "red", 'black', None
        elif feature.key == "gap":
            return False, "gray", 'black', None
        else:
            raise ValueError("Invalid Feature")

    fig, axs = plt.subplots(1, 2, subplot_kw={'projection': 'polar'}, figsize=(16, 8))
    
    graphics.plot_plasmid_map(
        axs[0], 
        genome.annotation, 
        len(genome.sequence), 
        label="compressed",
        feature_formatter=feature_formatter
    )

    graphics.plot_plasmid_map(
        axs[1], 
        decompressed_genome.annotation, 
        len(decompressed_genome.sequence), 
        label="decompressed",
        feature_formatter=feature_formatter
    )

    legend_elements = [
        Patch(facecolor='orange', edgecolor='black', label='CDS'),
        Patch(facecolor='blue', edgecolor='black', label='RBS'),
        Patch(facecolor='green', edgecolor='black', label='Start Codon'),
        Patch(facecolor='red', edgecolor='black', label='Stop Codon'),
        Patch(facecolor='gray', edgecolor='black', label='Gap'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', title='Feature Types')
    fig.savefig(filename)


if __name__ == "__main__":
    name = input("Enter name>").strip()

    original_fasta_filename = os.path.join(GENOME_DIRECTORY, f"{name}.fasta")
    original_gff_filename = os.path.join(GENOME_DIRECTORY, name, f"{name}.gff")
    decompressed_fasta_filename = os.path.join(GENOME_DIRECTORY, name, f"{name}_decompressed.fasta")
    decompressed_gff_filename = os.path.join(GENOME_DIRECTORY, name, f"{name}_decompressed.gff")

    for filename in [
        original_fasta_filename,
        original_gff_filename,
        decompressed_fasta_filename,
        decompressed_gff_filename
    ]:
        if not os.path.exists(filename):
            print(f"No such file: {filename}")
            exit(1)

    original_fasta_file = fasta.FastaFile.read(original_fasta_filename)
    original_gff_file = gff.GFFFile.read(original_gff_filename)
    original_sequence = AnnotatedSequence(
        gff.get_annotation(original_gff_file), 
        fasta.get_sequence(original_fasta_file)
    )

    decompressed_fasta_file = fasta.FastaFile.read(decompressed_fasta_filename)
    decompressed_gff_file = gff.GFFFile.read(decompressed_gff_filename)
    decompressed_sequence = AnnotatedSequence(
        gff.get_annotation(decompressed_gff_file), 
        fasta.get_sequence(decompressed_fasta_file)
    )

    plot_filename = os.path.join(GENOME_DIRECTORY, name, f"{name}_decompressed_plasmids.png")

    plot_plasmids(original_sequence, decompressed_sequence, plot_filename)
