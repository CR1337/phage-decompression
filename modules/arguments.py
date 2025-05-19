import argparse
from typing import List


ALLOWED_OPERATIONS: List[str] = ['decompress', 'plot', 'analyze']


def parse_arguments():
    parser = argparse.ArgumentParser(description="Tool for decompressing, plotting and analyzing phage genomes.")

    parser.add_argument(
        '-f', '--input-fasta-file',
        type=str,
        required=True,
        help=""
    )

    parser.add_argument(
        '-g', '--input-gff-file',
        type=str,
        required=True,
        help=""    
    )

    parser.add_argument(
        '-o', '--output-directory',
        type=str,
        help="Output directory for the analysis"
    )

    parser.add_argument(
        '-op', '--operations',
        nargs='+',
        choices=ALLOWED_OPERATIONS,
        required=True,
        help="List of operations to perform (decompress, plot, analyze)"
    )

    arguments = parser.parse_args()

    if 'analyze' in arguments.operations:
        if not getattr(arguments, 'output_directory'):
            parser.error("Analysis requires an output directory.")

    return arguments
