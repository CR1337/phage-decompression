import argparse


def _parse_file_format(parser: argparse.ArgumentParser):
    """
    Adds the input file format argument to the parser.

    Allows selection between 'fasta' and 'genbank' as input formats.
    """
    ALLOWED_FILE_FORMATS = ["fasta", "genbank"]

    parser.add_argument(
        "-if",
        "--input-format",
        type=str,
        choices=ALLOWED_FILE_FORMATS,
        default="fasta",
        help="Input file format. Must be either 'fasta' or 'genbank'. Default is 'fasta'.",
    )


def _parse_output_directory(parser: argparse.ArgumentParser):
    """
    Adds the output directory argument to the parser.

    This is the location where output files (FASTA, GFF, plots, summaries, etc.) will be stored.
    """
    parser.add_argument(
        "-o",
        "--output-directory",
        type=str,
        required=True,
        help="Directory where output files will be saved.",
    )


def _parse_decompress_arguments(parser: argparse.ArgumentParser):
    """
    Adds arguments for the 'decompress' command.

    These arguments configure decompression of an annotated genome file into its full nucleotide representation.
    """
    _parse_file_format(parser)

    parser.add_argument(
        "-i",
        "--input-filename",
        type=str,
        required=True,
        help="Path to the input genome file (FASTA or GenBank).",
    )

    parser.add_argument(
        "-g",
        "--input-gff-filename",
        type=str,
        help="Optional GFF file with annotations (required if input is FASTA).",
    )

    parser.add_argument(
        "-l",
        "--file-list-filename",
        type=str,
        help="Optional path to a file where a list of JSON output files will be appended.",
    )

    _parse_output_directory(parser)


def _parse_summary_arguments(parser: argparse.ArgumentParser):
    """
    Adds arguments for the 'summary' command.

    Used to compute statistics and generate visualizations for a batch of decompression results.
    """
    parser.add_argument(
        "-l",
        "--file-list-filename",
        type=str,
        required=True,
        help="Path to a file containing a list of JSON files with decompression statistics.",
    )

    _parse_output_directory(parser)


def _parse_glm2_arguments(parser: argparse.ArgumentParser):
    """
    Adds arguments for the 'glm2' command.

    Converts annotated genome sequences into glm2-compatible format for downstream analysis.
    """
    _parse_file_format(parser)

    parser.add_argument(
        "-i",
        "--input-filename",
        type=str,
        required=True,
        help="Path to the input genome file (FASTA or GenBank).",
    )

    parser.add_argument(
        "-g",
        "--input-gff-filename",
        type=str,
        help="Optional GFF file with annotations (required if input is FASTA and --raw is not set).",
    )

    parser.add_argument(
        "-r",
        "--raw",
        action="store_true",
        help="Output the raw DNA sequence in glm2 format, ignoring annotations.",
    )

    parser.add_argument(
        "-o",
        "--output-filename",
        type=str,
        required=True,
        help="Path to the output file where glm2-formatted sequence will be saved.",
    )


def _parse_similarity_arguments(parser: argparse.ArgumentParser):
    """
    Adds arguments for the 'similarity' command.

    Compares sequences and computes similarity metrics. Can also visualize pairwise similarities.
    """
    parser.add_argument(
        "-i",
        "--input-filenames",
        type=str,
        nargs="+",
        required=True,
        help="List of two input genome files to compare or one csv file containing similarity values.",
    )

    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="Plot a similarity heatmap from a precomputed CSV file.",
    )

    parser.add_argument(
        "-c",
        "--cluster",
        action="store_true",
        help="Apply hierarchical clustering when plotting similarities.",
    )

    parser.add_argument(
        "-o",
        "--output-filename",
        type=str,
        required=True,
        help="Path to the output file (either for similarity results or plot image).",
    )


def _validate_arguments(arguments: argparse.Namespace, parser: argparse.ArgumentParser):
    """
    Validates logical constraints between arguments depending on the command.

    Ensures required input files are provided based on input format and selected operation.
    """
    if (
        arguments.command == "decompress"
        and arguments.input_format == "fasta"
        and not hasattr(arguments, "input_gff_filename")
    ):
        parser.error("You need to specify a gff file when decompressing a fasta file.")

    if (
        arguments.command == "glm2"
        and arguments.input_format == "fasta"
        and not hasattr(arguments, "input_gff_filename")
        and not arguments.raw
    ):
        parser.error("You need to specify a gff file when translating a fasta file.")


def parse_arguments():
    """
    Parses command-line arguments and validates constraints.

    Returns:
        argparse.Namespace: Parsed arguments object.
    """
    parser = argparse.ArgumentParser(
        description="Tool for decompressing, plotting and analyzing phage genomes."
    )

    subparsers = parser.add_subparsers(dest="command")

    decompress_parser = subparsers.add_parser(
        "decompress",
        required=True,
        help="Decompress a genome into full nucleotide sequence and save as FASTA + GFF.",
    )
    _parse_decompress_arguments(decompress_parser)

    summary_parser = subparsers.add_parser(
        "summary",
        required=True,
        help="Summarize and visualize decompression results from multiple JSON files.",
    )
    _parse_summary_arguments(summary_parser)

    glm2_parser = subparsers.add_parser(
        "glm2",
        required=True,
        help="Generate a glm2-compatible sequence file from an annotated genome.",
    )
    _parse_glm2_arguments(glm2_parser)

    similarity_parser = subparsers.add_parser(
        "similarity",
        required=True,
        help="Compute or plot similarity between genome sequences."
    )
    _parse_similarity_arguments(similarity_parser)

    arguments = parser.parse_args()
    _validate_arguments(arguments, parser)
    return arguments
