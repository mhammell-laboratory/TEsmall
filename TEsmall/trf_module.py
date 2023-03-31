import argparse
import glob
import logging
import os

import subprocess
import pybedtools

import pandas as pd
import numpy as np

from .settings import *

import sys


def get_CCA(bam):
    """Makes a fasta file of CCA terminated mapped reads from multibam file"""
    root = os.path.splitext(os.path.basename(bam))[0]
    root = os.path.splitext(root)[0]
    cca_fasta = "{0}_cca.fa".format(root)
    grep_input1 = "samtools view {0} | awk -v OFS='\t' '{{print $1, $10}}' | uniq -c | grep 'CCA$' | awk -v OFS=',' '{{print $2, substr($3, 1, length($3)-3)}}'".format(bam)
    grep_process1 = subprocess.Popen(grep_input1, stdout=subprocess.PIPE, shell=True)
    grep_output1 = grep_process1.stdout.read()
    grep_input2 = "samtools view {0} | awk -v OFS='\t' '{{print $1, $10}}' | uniq -c | grep 'CC$' | awk -v OFS=',' '{{print $2, substr($3, 1, length($3)-2)}}'".format(bam)
    grep_process2 = subprocess.Popen(grep_input2, stdout=subprocess.PIPE, shell=True)
    grep_output2 = grep_process2.stdout.read()
    grep_input3 = "samtools view {0} | awk -v OFS='\t' '{{print $1, $10}}' | uniq -c | grep -E 'AC$|GC$|TC$' | awk -v OFS=',' '{{print $2, substr($3, 1, length($3)-1)}}'".format(bam)
    grep_process3 = subprocess.Popen(grep_input3, stdout=subprocess.PIPE, shell=True)
    grep_output3 = grep_process3.stdout.read()
    grep_lines1 = grep_output1.split()
    grep_lines2 = grep_output2.split()
    grep_lines3 = grep_output3.split()
    grep_lines = grep_lines1 + grep_lines2 + grep_lines3

    with open(cca_fasta, "w") as cca_out:
        for l in grep_lines:
            l = l.strip()
            l = l.split(b',')
            read_id = l[0].decode("utf-8")
            read_seq = l[1].decode("utf-8")
            cca_out.write(">{0}\n{1}\n".format(read_id, read_seq))
    return cca_fasta


def map_CCA(cca_fasta, ebwt, multibam, maxaln=50, mm=0):
    prefix = os.path.splitext(os.path.basename(multibam))[0]
    prefix = os.path.splitext(prefix)[0]
    prefix = os.path.splitext(prefix)[0]
    logfile = "{0}.tRNA.log".format(prefix)
    unfile = "{0}.unaligned.cca.fa".format(prefix)
    command = 'bowtie -f --chunkmbs 1024 -a -m {0} --best --strata -v {1} -S '\
              '--un "{2}" "{3}" "{4}" 2> "{5}" | samtools view -F 4 -buS - '\
              '2>> "{5}" | samtools sort -n -o "{6}.tRNA.bam" - 2>> "{5}"'
    command = command.format(maxaln, mm, unfile, ebwt, cca_fasta, logfile, prefix)

    logging.info("Aligning CCA trimmed reads to tRNA sequences...")
    subprocess.call(command, shell=True)

    return "{0}.tRNA.bam".format(prefix)


def extract_CCA_map(multibam, cca_trimmed_bam):
    """Makes file of join output of original CCA(+) mapping and CCA(-) mapping reads. Extracts these from multi bam file piped to downstream bed intersect"""
    prefix = os.path.splitext(os.path.basename(multibam))[0]
    prefix = os.path.splitext(prefix)[0]
    prefix = os.path.splitext(prefix)[0]

    command_single = "samtools view {0}.tRNA.bam | cut -f 1 | sort -u > {0}_3trf_rid.txt"
    command_single = command_single.format(prefix)
    command_sam_head = "samtools view -H {0} > {1}.header.txt"
    command_sam_head = command_sam_head.format(multibam, prefix)
    command_clean_multi_bam = "samtools view {0} | sort -k1,1 | join -t $'\t' -v 1 -1 1 -2 1 -  {1}_3trf_rid.txt | cat {1}.header.txt - | samtools view -bo {1}.3trf_free.bam -"
    command_clean_multi_bam = command_clean_multi_bam.format(multibam, prefix)
    command_3trf_mappers_bam = "samtools view {0} | sort -k1,1 | join -t $'\t' -1 1 -2 1 - {1}_3trf_rid.txt | cat {1}.header.txt - | samtools view -bo {1}.3trf.bam -"
    command_3trf_mappers_bam = command_3trf_mappers_bam.format(multibam, prefix)

    command_clean = "rm -f {0}.3trf_body.sam {0}.3trf.sam {0}_3trf_rid.txt"
    command_clean = command_clean.format(prefix)

    logging.info("Finding 3' tRF mappers...")
    subprocess.call(command_single, shell=True)
    subprocess.check_call(command_sam_head, shell=True)
    subprocess.Popen(command_clean_multi_bam, shell=True, executable='/bin/bash')
    subprocess.check_call(command_3trf_mappers_bam, shell=True, executable='/bin/bash')
    subprocess.call(command_clean, shell=True)
    out_trf_align = '{0}.3trf.bam'.format(prefix)
    out_residual_align = '{0}.3trf_free.bam'.format(prefix)
    return out_trf_align, out_residual_align

def fix_tRNA_mapped_coor(trnabam):
    root = os.path.splitext(os.path.basename(trnabam))[0]
    root = os.path.splitext(root)[0]
    command_bam_to_sam = 'samtools view -h {0} > {1}.change_coor.sam'
    command_bam_to_sam = command_bam_to_sam.format(trnabam, root)
    subprocess.call(command_bam_to_sam, shell=True)

    genome_head = '{0}.header.txt'.format(root)
    sam = '{0}.change_coor.sam'.format(root)
    files = [genome_head, sam]
    sam_out = '{0}.trna_for_intersect.sam'.format(root)
    with open(sam_out, 'w') as outfile:
        with open(genome_head) as header:
            for lines in header:
                outfile.write(lines)
        with open(sam) as infile:
            for line in infile:
                if line[0] == '@':
                    continue
                else:
                    read_line = line.split()
                    location = read_line[2].split(':')
                    name = location[0]
                    chrom = location[1]
                    coor = location[2].split('-')
                    strand = location[3]
                    if(strand == "-"):
                        if(int(read_line[1]) & 16):
                            read_line[1] = str(int(read_line[1]) - 16)
                        else:
                            read_line[1] = str(int(read_line[1]) + 16)
                    newcoor = int(coor[0]) + int(read_line[3])  # seems like pybedtools takes 1 indexed bam and turns it to 0 index upon intersect check this
                    end_read = "\t".join(read_line[4:])
                    outfile.write("\t".join([read_line[0], read_line[1], chrom, str(newcoor), end_read]) + "\n")

    command_sam_bam = "samtools view -b {0}.trna_for_intersect.sam > {0}.trna_for_intersect.bam".format(root)
    subprocess.call(command_sam_bam, shell=True)
    bam_out = "{0}.trna_for_intersect.bam".format(root)
    command_clean = "rm -f {0}.trna_for_intersect.sam {0}.change_coor.sam {0}.header.txt".format(root)
    subprocess.call(command_clean, shell=True)
    return bam_out


def intersect_CCA_map(trnaorgbam, trnabam, annot_dir, bed):
    root = os.path.splitext(os.path.basename(trnaorgbam))[0]
    root = os.path.splitext(root)[0]
    logging.info("Assigning 3' tRFs to transposable elements...")
    if bed == "TE":
        bedfile =  os.path.join(annot_dir, "TE.bed")
        outfname = "{0}.3trf.TE.mapper.anno".format(root)
        bamfile = trnaorgbam
        f_type = "TE"
    else:
        bedfile = os.path.join(annot_dir, "structural_RNA.bed")
        outfname = "{0}.3trf.struc.mapper.anno".format(root)
        bamfile = trnabam
        f_type = "structural_RNA"

    with open(outfname, "w") as outfile:
        columns = ["rid", "rchr", "rstart", "rend", "rstrand", "rlen", "ftype", "fid", "fchr", "fstart", "fend", "fstrand", "overlap"]
        outfile.write("\t".join(columns) + "\n")
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
            if f_type == "TE":
                if r_strand == f_strand:
                    outfile.write("\t".join([r_id, r_chrom, r_start, r_end, r_strand, r_len, "sense_TE", f_id, f_chrom, f_start, f_end, f_strand, f_len]) + "\n")
                else:
                    outfile.write("\t".join([r_id, r_chrom, r_start, r_end, r_strand, r_len, "anti_TE", f_id, f_chrom, f_start, f_end, f_strand, f_len]) + "\n")
            else:
                outfile.write("\t".join([r_id, r_chrom, r_start, r_end, r_strand, r_len, f_type, f_id, f_chrom, f_start, f_end, f_strand, f_len]) + "\n")
    return outfname


def handle_cca(multibam, ebwt, annot_dir):
    # filters CCA terminated reads and cleaves CCA tail makes fasta for align
    trna_fa = get_CCA(multibam)
    # aligns fasta
    cca_bam = map_CCA(trna_fa, ebwt, multibam)
    # cleans original map file of mapped CCA trfs
    out_trf_align, out_residual_align = extract_CCA_map(multibam, cca_bam)
    trnabam = fix_tRNA_mapped_coor(cca_bam)
    trf_map = intersect_CCA_map(out_trf_align, trnabam, annot_dir, "structural_RNA") # change to out cca_bam after you fix where it maps
    trf_te_map = intersect_CCA_map(out_trf_align, trnabam, annot_dir, "TE")
    return trf_map, out_residual_align


# def main():
#     if len(sys.argv) < 3:
#         print("Usage: %s [multi BAM] [tRNA BAM]" % sys.argv[0])
#         sys.exit(1)
#     extract_CCA_map(sys.argv[1], sys.argv[2])

# if __name__ == "__main__":
#     main()
