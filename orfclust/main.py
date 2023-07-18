from typing import Iterator
import numpy as np
import argparse
import copy
import sys
import os

from classes.txgroup import Transcriptome, Gene, Bundle
from classes.transcript import Transcript

# TODO: if txcounts/expression is not provided - can not report highest expression - report mode should be set to something else

def orfclust(args):

    out_fname = args.output.rstrip(".gtf").rstrip(".GTF").rstrip(".gff").rstrip(".GFF")
    out_gtf_fp = open(out_fname+".gtf","w+") # clustered version of the annotation file
    out_exp_fp = open(out_fname+".txcounts","w+") # clustered version of the expression file
    out_grp_fp = open(out_fname+".groups","w+") # two columns: representative transcript_id, comma-separated list of all transcript_ids in the group


    # load header line from the expression file and output to the new file
    if args.exp is not None:
        with open(args.exp,"r") as exp_fp:
            # if the first line is a header - output it to the new file
            # header is defined as a line that starts with "tx_id" or line where all values other than the first one are numeric
            first_line = exp_fp.readline().strip()
            has_tx_id_col = first_line.startswith("tx_id")
            numeric_values = [x.replace(".","",1).isnumeric() for x in first_line.split("\t")[1:]]

            if has_tx_id_col or not all(numeric_values):
                out_exp_fp.write(first_line+"\n")

    transcriptome = Transcriptome()
    transcriptome.build_from_file(args.gtf)
    if args.exp is not None:
        transcriptome.load_expression(args.exp)
    transcriptome.gid_sort()

    representative_func = None # this function is applied to each group of transcripts to select a single representative transcript
    if args.rep == "SHORT":
        representative_func = lambda tx: tx.orf_len()
    elif args.rep == "LONG":
        representative_func = lambda tx: tx.orf_len()
    elif args.rep == "EXP":
        representative_func = lambda tx: max(tx.expression)
    elif args.rep == "ORF":
        representative_func = lambda tx: tx.orf_coords()
    else:
        raise Exception("Unknown representative transcript option: "+args.rep)

    for gene in transcriptome.gene_it():
        groupper = ["cds"]
        if args.mode == "FUNC":
            groupper = "cds"
        for grp_id,grp in gene.group_by(groupper):
            # select representative transcript for each group
            rep_tx = None
            if args.rep == "EXP":
                rep_tx = copy.deepcopy(max(grp.object_it(),key=lambda tx: tx.get_expression(sum)))
            elif args.rep == "SHORT":
                rep_tx = copy.deepcopy(min(grp.object_it(),key=lambda tx: tx.elen()))
            elif args.rep == "LONG":
                rep_tx = copy.deepcopy(max(grp.object_it(),key=lambda tx: tx.elen()))
            elif args.rep == "ORF":
                rep_tx = copy.deepcopy(next(grp.object_it())) # get first transcript as a model to modify
                rep_tx.clear_exons()
                rep_tx.set_exons(rep_tx.get_cds())
                rep_tx.finalize()
            else:
                raise Exception("Unknown representative transcript option: "+args.rep)

            # accumulate expression data
            total_exp = sum([tx.get_expression(sum) for tx in grp.object_it()])
            rep_tx.set_expression([total_exp])

            out_gtf_fp.write(rep_tx.to_gtf()+"\n")
            out_grp_fp.write(rep_tx.get_tid()+"\t"+",".join([tx.get_tid() for tx in grp.object_it()])+"\n")
            # lastly, we need to figure out how to handle the expression data
            # we need to sum up expression values for each sample in the txcounts file separately
            # the output txcounts should have the same samples preserved as in the input txcounts file
            out_exp_fp.write(rep_tx.get_tid()+"\t")
            cumulative_expressions = []
            for tx in grp.object_it():
                if not cumulative_expressions:
                    cumulative_expressions = tx.get_expression()
                else:
                    cumulative_expressions = [x + y for x, y in zip(cumulative_expressions, tx.get_expression())]

            for exp in cumulative_expressions:
                out_exp_fp.write(str(round(exp,2))+"\t")
            out_exp_fp.write("\n")

    out_gtf_fp.close()
    out_exp_fp.close()
    out_grp_fp.close()

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
    parser.add_argument("-o",
                        "--output",
                        required=True,
                        type=str,
                        help="Output base file name. Anythign after the last dot will be ignored.")
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