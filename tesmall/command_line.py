import argparse
import logging
import os
import sys

from abundance import *
from alignment import *
from annotation import *
from settings import *
from summary import *

def main():
    parser = argparse.ArgumentParser(prog="TEsmall")
    parser.add_argument("-a", "--adapter", metavar="STR",
        default='TGGAATTCTCGGGTGCCAAGG', help="Sequence of an adapter that was "
        "ligated to the 3' end. The adapter itself and anything that follows "
        "is trimmed. (default: %(default)s)")
    parser.add_argument("-m", "--minlen", metavar="INT", type=int,
        default=16, help="Discard trimmed reads that are shorter than INT. "
        "Reads that are too short even before adapter removal are also "
        "discarded. (default: %(default)d)")
    parser.add_argument("-M", "--maxlen", metavar="INT", type=int,
        default=36, help="Discard trimmed reads that are longer than INT. "
        "Reads that are too long even before adapter removal are also "
        "discarded. (default: %(default)d)")
    parser.add_argument("-g", "--genome", metavar="STR", default="hg19",
        choices=["dm3", "mm9", "hg19"], help="Version of reference genome "
        "(default: %(default)s)")
    parser.add_argument("--maxaln", metavar="INT", type=int, default=100,
        help="Suppress all alignments for a particular read if more than INT "
        "reportable alignments exist for it. (default: %(default)s)")
    parser.add_argument("--mismatch", metavar="INT", type=int, default=0,
        choices=[0, 1, 2, 3], help="Report alignments with at most INT "
        "mismatches. (default: %(default)s)")
    parser.add_argument("-o", "--order", metavar="STR", nargs="+",
    choices=["structural_RNA", "miRNA", "hairpin", "exon", "TE", "intron",
        "piRNA_cluster"], default=["structural_RNA", "miRNA", "hairpin", "exon",
        "TE", "intron", "piRNA_cluster"], help="Annotation priority. (default: structural_RNA miRNA "
        "hairpin exon TE intron piRNA_cluster)")
    parser.add_argument("-p", "--parallel", metavar="INT", type=int, default=1,
        help="Parallel execute by INT CPUs. (default: %(default)s)")
    parser.add_argument("-f", "--fastq", metavar="STR", nargs="+", help="Input in "
        "FASTQ format. Compressed input is supported and auto-detected from "
        "the filename extension (.gz).")
    parser.add_argument("-l", "--label", metavar="STR", nargs="+",
        help="Unique label for each sample.")
    parser.add_argument('--verbose', metavar="INT", type=int, nargs='?',
        default=2, help="Set verbose level. 0: only show critical message, 1: show additional "
        "warning message, 2: show process information, 3: show debug messages. DEFAULT:2")
    args = parser.parse_args()
    
    logging.basicConfig(level=(4 - args.verbose) * 10,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stderr, filemode="w")

    get_requirements(args.genome)
    annofiles = []
    if not args.label:
        args.label = map(lambda s: re.sub(r".f(ast)?q(.gz)?$", "", os.path.basename(s)), args.fastq)
    else:
        assert len(set(args.label)) == len(args.fastq)
    for label, fastq in zip(args.label, args.fastq):
        trimmed_fastq = trim_3p_adapters(fastq, label, args.adapter, args.minlen, args.maxlen)
        trimmed_fastq = trim_5p_adapters(trimmed_fastq, label, "GTTCAGAGTTCTACAGTCCGACGATC", args.minlen, args.maxlen)
        btidx = os.path.join(BOWTIE_INDEX.format(args.genome), "genome")
        rbtidx = os.path.join(BOWTIE_INDEX.format(args.genome), "rDNA")
        bestone_bam, unfile = map_bestone_reads(trimmed_fastq, rbtidx, 2)
        #bestone_bam, unfile = map_bestone_reads(fastq, rbtidx, args.mismatch)
        multi_bam = map_multi_reads(unfile, btidx, args.maxaln, args.mismatch)
        readinfo = get_read_info(multi_bam)
        annofiles.append(annotate_reads(args.genome, args.order, multi_bam))

    calc_composition(annofiles)
    gen_summary(args.label, args.order, args.maxaln)
    calc_abundance(annofiles)
