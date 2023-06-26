from typing import Iterator
import numpy as np
import argparse
import sys
import os

from classes.txgroup import Transcriptome, Gene, Bundle
from classes.transcript import Transcript

def orfclust(args):
    transcriptome = Transcriptome()
    transcriptome.build_from_file(args.gtf)
    transcriptome.load_expression(args.exp)
    transcriptome.gid_sort()

    for gene in transcriptome.gene_it():
        print(gene)

def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')
    parser.add_argument("-g",
                        "--gtf",
                        required=True,
                        type=str,
                        help="Annotation in a GFF/GTF format.")
    parser.add_argument("-e",
                        "--exp",
                        required=False,
                        type=str,
                        help="Optional file with expression data to be aggregated for transcripts that are groupped together. \
                            The file should be a tab delimited file with the first column containing transcript IDs \
                            and the second column containing expression values. The file should not contain a header.")
    parser.add_argument("--use_geneid",
                        action="store_true",
                        required=False,
                        help="If selected will use gene_id attribute to group transcripts. \
                            Without this flag transcripts are combined based on overlap instead.")
    parser.add_argument("-m",
                        "--mode",
                        required=False,
                        type=str,
                        default="ORF",
                        const='ORF',
                        nargs='?',
                        choices=['ORF', 'FUNC'],
                        help="Mode to group transcripts. Options are: ORF and FUNC. \
                            If ORF is selected - transcripts are groupped if they share the same ORF. \
                            If FUNC is selected - transcripts are split into two categories: with and without ORF. \
                            Default is (default: %(default)s).")
    parser.add_argument("-r",
                        "--rep",
                        required=False,
                        type=str,
                        default="EXP",
                        const='EXP',
                        nargs='?',
                        choices=['EXP', 'SHORT', "LONG", "ORF"],
                        help="Representative transcript. \
                              If EXP is selected - the top expressed transcript in each group will be chosen (requires expression data via '-e'). \
                              If SHORT is selected - the shorted transcript will be selected. \
                              IF ORF is selected - transcripts will be represented by the coordinates of their ORFs instead. Default is (default: %(default)s).")

    parser.set_defaults(func=orfclust)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])