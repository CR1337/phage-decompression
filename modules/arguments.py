import argparse


def _parse_file_format(parser: argparse.ArgumentParser):
    ALLOWED_FILE_FORMATS = ['fasta', 'genbank']

    parser.add_argument(
        '-if', '--input-format',
        type=str,
        choices=ALLOWED_FILE_FORMATS,
        default='fasta',
        help=""
    )

def _parse_output_directory(parser: argparse.ArgumentParser):
    parser.add_argument(
        '-o', '--output-directory',
        type=str,
        required=True,
        help=""
    )
    

def _parse_decompress_arguments(parser: argparse.ArgumentParser):
    _parse_file_format(parser)

    parser.add_argument(
        '-i', '--input-filename',
        type=str,
        required=True,
        help=""
    )

    parser.add_argument(
        '-g', '--input-gff-filename',
        type=str,
        help=""
    )

    parser.add_argument(
        '-l', '--file-list-filename',
        type=str,
        help=""
    )

    _parse_output_directory(parser)


def _parse_summary_arguments(parser: argparse.ArgumentParser):
    parser.add_argument(
        '-l', '--file-list-filename',
        type=str,
        required=True,
        help=""
    )

    _parse_output_directory(parser)


def _parse_glm2_arguments(parser: argparse.ArgumentParser):
    _parse_file_format(parser)

    parser.add_argument(
        '-i', '--input-filename',
        type=str,
        required=True,
        help=""
    )

    parser.add_argument(
        '-g', '--input-gff-filename',
        type=str,
        help=""
    )

    parser.add_argument(
        '-r', '--raw',
        action="store_true",
        help=""
    )

    parser.add_argument(
        '-o', '--output-filename',
        type=str,
        required=True,
        help=""
    )


def _parse_similarity_arguments(parser: argparse.ArgumentParser):
    parser.add_argument(
        '-i', '--input-filenames',
        type=str,
        nargs='+',
        required=True,
        help=""
    )

    parser.add_argument(
        '-p', '--plot',
        action='store_true',
        help=""
    )

    parser.add_argument(
        '-c', '--cluster',
        action='store_true',
        help=""
    )

    parser.add_argument(
        '-o', '--output-filename',
        type=str,
        required=True,
        help=""
    )


def _validate_arguments(arguments: argparse.Namespace, parser: argparse.ArgumentParser):
    if (
        arguments.command == 'decompress'
        and arguments.input_format == 'fasta'
        and not hasattr(arguments, 'input_gff_filename')
    ):
        parser.error("You need to specify a gff file when decompressing a fasta file.")

    if (
        arguments.command == 'glm2'
        and arguments.input_format == 'fasta'
        and not hasattr(arguments, 'input_gff_filename')
        and not arguments.raw
    ):
        parser.error("You need to specify a gff file when translating a fasta file.")


def parse_arguments():
    parser = argparse.ArgumentParser(description="Tool for decompressing, plotting and analyzing phage genomes.")

    subparsers = parser.add_subparsers(dest='command')

    decompress_parser = subparsers.add_parser('decompress', help="")
    _parse_decompress_arguments(decompress_parser)

    summary_parser = subparsers.add_parser('summary', help="")
    _parse_summary_arguments(summary_parser)

    glm2_parser = subparsers.add_parser('glm2', help="")
    _parse_glm2_arguments(glm2_parser)

    similarity_parser = subparsers.add_parser('similarity', help="")
    _parse_similarity_arguments(similarity_parser)

    arguments = parser.parse_args()
    _validate_arguments(arguments, parser)
    return arguments
