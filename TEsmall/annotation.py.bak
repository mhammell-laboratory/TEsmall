import argparse
import glob
import logging
import os

import pybedtools

from settings import *

def annotate_reads(genome, order, multi):
    root = os.path.splitext(os.path.basename(multi))[0]
    root = os.path.splitext(root)[0]
    logging.info("Assigning reads to genomic features...")
    bedfiles = map(lambda s: os.path.join(ANNOTATION.format(genome), "{0}.bed".format(s)), order)
    #bedfiles = sorted(glob.glob(os.path.join(ANNOTATION.format(genome), "*.bed")))
    outfname = "{0}.anno".format(root)
    with open(outfname, "w") as outfile:
        columns = ["rid", "rchr", "rstart", "rend", "rstrand", "rlen", "ftype", "fid", "fchr", "fstart", "fend", "fstrand", "overlap"]
        outfile.write("\t".join(columns) + "\n")
        anno_reads = set()
        bamfile = multi
        for bedfile in bedfiles:
            current_reads = set()
            f_type = os.path.splitext(os.path.basename(bedfile))[0]
            bam = pybedtools.BedTool(bamfile)
            bed = pybedtools.BedTool(bedfile)
            for line in bam.intersect(bed, wo=True, bed=True, f=0.9, sorted=False, stream=True):
                r_chrom = line[0]
                r_start = line[1]
                r_end = line[2]
                r_id = line[3]
                r_strand = line[5]
                r_len = line[10][:-1]
                f_chrom = line[12]
                f_start = line[13]
                f_end = line[14]
                f_id = line[15]
                f_strand = line[17]
                f_len = line[18]
                if f_type in ["miRNA", "hairpin"]:
                    if r_strand == f_strand:
                        current_reads.add(r_id)
                        if r_id not in anno_reads:
                            outfile.write("\t".join([r_id, r_chrom, r_start, r_end, r_strand, r_len, f_type, f_id, f_chrom, f_start, f_end, f_strand, f_len]) + "\n")
                else:
                    current_reads.add(r_id)
                    if r_id not in anno_reads:
                        if f_type == "TE":
                            if r_strand == f_strand:
                                outfile.write("\t".join([r_id, r_chrom, r_start, r_end, r_strand, r_len, "sense_TE", f_id, f_chrom, f_start, f_end, f_strand, f_len]) + "\n")
                            else:
                                outfile.write("\t".join([r_id, r_chrom, r_start, r_end, r_strand, r_len, "anti_TE", f_id, f_chrom, f_start, f_end, f_strand, f_len]) + "\n")
                        else:
                            outfile.write("\t".join([r_id, r_chrom, r_start, r_end, r_strand, r_len, f_type, f_id, f_chrom, f_start, f_end, f_strand, f_len]) + "\n")
            anno_reads = anno_reads.union(current_reads)
    return outfname
