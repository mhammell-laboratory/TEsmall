import argparse
import glob
import logging
import os

import subprocess
import pybedtools

import pandas as pd
import numpy as np
from .EMAlgorithm import *
from .abundance import *
from .settings import *


def flatten(multi_reads, container=None):
    if container is None:
        container = []
    for s in multi_reads:
        container.extend(s)
    return container


def annotate_reads(multi, annot_dir, order):
    root = os.path.splitext(os.path.basename(multi))[0]
    root = os.path.splitext(root)[0]
    logging.info("Assigning reads to genomic features...")
    bedfiles = [os.path.join(annot_dir, "{0}.bed".format(s)) for s in order]
    #bedfiles = ['/sonas-hs/mhammell/hpc/home/koneill/.tesmall/genomes/hg19/annotation/structural_RNA.bed', '/sonas-hs/mhammell/hpc/home/koneill/.tesmall/genomes/hg19/annotation/miRNA.bed', '/sonas-hs/mhammell/hpc/home/koneill/.tesmall/genomes/hg19/annotation/hairpin.bed', '/sonas-hs/mhammell/hpc/data/kat/tesmall_em_test/tes_em_test_2/TE.bed', '/sonas-hs/mhammell/hpc/home/koneill/.tesmall/genomes/hg19/annotation/exon.bed', '/sonas-hs/mhammell/hpc/home/koneill/.tesmall/genomes/hg19/annotation/intron.bed', '/sonas-hs/mhammell/hpc/home/koneill/.tesmall/genomes/hg19/annotation/piRNA_cluster.bed']
    outfname = "{0}.anno".format(root)
    with open(outfname, "w") as outfile:
        columns = ["rid", "rchr", "rstart", "rend", "rstrand", "rlen", "ftype", "fid", "fchr", "fstart", "fend", "fstrand", "overlap"]
        outfile.write("\t".join(columns) + "\n")
        anno_reads = set()
        bamfile = multi
        te_s_multi = []
        te_s_index = {}
        te_s_len = {}
        te_s_idx_ct = 0
        te_as_multi = []
        te_as_index = {}
        te_as_len = {}
        te_as_idx_ct = 0
        mir_multi = []
        mir_index = {}
        mir_len = {}
        mir_idx_ct = 0
        te_s_read_ref = {}
        te_as_read_ref = {}
        mir_read_ref = {}
        # deciding size of numpy array for bedgraph file
        array_count = 0
        for bedfile in bedfiles:
            current_reads = set()
            f_type = os.path.splitext(os.path.basename(bedfile))[0]
            bam = pybedtools.BedTool(bamfile)
            bed = pybedtools.BedTool(bedfile)
            te_s_int = []
            te_as_int = []
            mir_int = []
            previousread = ""
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
                currentread = r_id
                if f_type in ["miRNA", "hairpin"]:
                    if r_strand == f_strand:
                        current_reads.add(r_id)
                        if r_id not in anno_reads:
                            x = 0
                            if f_id in mir_read_ref:
                                x = 1
                                mir_read_ref[f_id].append(r_id)
                                array_count += 1
                            if x != 1:
                                mir_read_ref[f_id] = [r_id]
                                array_count += 1
                            if f_id not in list(mir_index.values()):
                                mir_index[mir_idx_ct] = str(f_id)
                                mir_len[mir_idx_ct] = int(f_len)
                                mir_idx_ct += 1
                            # read.append(r_id)
                            if currentread == previousread or previousread == "":
                                mir_int.append(f_id)
                            if currentread != previousread:
                                mir_multi.append(mir_int)
                                mir_int = []
                                mir_int.append(f_id)
                            outfile.write("\t".join([r_id, r_chrom, r_start, r_end, r_strand, r_len, f_type, f_id, f_chrom, f_start, f_end, f_strand, f_len]) + "\n")
                else:
                    current_reads.add(r_id)
                    if r_id not in anno_reads:
                        if f_type == "TE":
                            if r_strand == f_strand:
                                y = 0
                                if f_id in te_s_read_ref:
                                    y = 1
                                    te_s_read_ref[f_id].append(r_id)
                                    array_count += 1
                                if y != 1:
                                    te_s_read_ref[f_id] = [r_id]
                                    array_count += 1
                                if f_id not in list(te_s_index.values()):
                                    te_s_index[te_s_idx_ct] = str(f_id)
                                    te_s_len[te_s_idx_ct] = int(f_len)
                                    te_s_idx_ct += 1
                                # read.append(r_id)  # test for order, is ordered
                                if currentread == previousread or previousread == "":
                                    te_s_int.append(f_id)
                                if currentread != previousread:
                                    te_s_multi.append(te_s_int)
                                    te_s_int = []
                                    te_s_int.append(f_id)
                                outfile.write("\t".join([r_id, r_chrom, r_start, r_end, r_strand, r_len, "sense_TE", f_id, f_chrom, f_start, f_end, f_strand, f_len]) + "\n")

                            else:
                                z = 0
                                if f_id in te_as_read_ref:
                                    z = 1
                                    te_as_read_ref[f_id].append(r_id)
                                    array_count += 1
                                if z != 1:
                                    te_as_read_ref[f_id] = [r_id]
                                    array_count += 1
                                if f_id not in list(te_as_index.values()):
                                    te_as_index[te_as_idx_ct] = str(f_id)
                                    te_as_len[te_as_idx_ct] = int(f_len)
                                    te_as_idx_ct += 1
                                # read.append(r_id)  # test for order, is ordered
                                if currentread == previousread or previousread == "":
                                    te_as_int.append(f_id)
                                if currentread != previousread:
                                    te_as_multi.append(te_as_int)
                                    te_as_int = []
                                    te_as_int.append(f_id)
                                    outfile.write("\t".join([r_id, r_chrom, r_start, r_end, r_strand, r_len, "anti_TE", f_id, f_chrom, f_start, f_end, f_strand, f_len]) + "\n")
                        else:
                            outfile.write("\t".join([r_id, r_chrom, r_start, r_end, r_strand, r_len, f_type, f_id, f_chrom, f_start, f_end, f_strand, f_len]) + "\n")
                        previousread = currentread
            anno_reads = anno_reads.union(current_reads)
    return outfname, te_s_multi, te_s_index, te_s_len, te_as_multi, te_as_index, te_as_len, mir_multi, mir_index, mir_len, te_s_read_ref, te_as_read_ref, mir_read_ref, array_count


def getcountsEM(te_multi, te_idx_dict):
    # make separate lists of lists from unique and multimapper reads
    # returns list ordered by te_idx of counts for unique and multi mappers
    unique_count = []
    multi_count = []
    te_names = list(te_idx_dict.values())  # get list of all te int names
    multi_te_weights = {key: 0 for key in te_names}  # make a dic of all TEs assigned a value of 0

    for i in te_multi:
        if len(i) > 1:
            weight = 1.0 / len(i)
            multi_count.append(i)
            for te in i:
                multi_te_weights[te] += weight
        else:
            unique_count.append(i)
    # go from dictionary of multimapper derived weights to list of counts with order consistent with te_idx key
    # invert dict to make te_names keys
    # inv_te_idx = {v: k for k, v in te_idx_dict.iteritems()}
    mcounts = []
    for i in range(len(te_idx_dict)):
        teoi = te_idx_dict[i]
        mcounts.append(multi_te_weights[teoi])

    # flatten unique list of lists to count number of instances
    flat_unique = flatten(unique_count)

    # list of counts whose order represents index of associated tid for EM feed in
    ucounts = []
    for i in range(len(te_idx_dict)):
        uq = flat_unique.count(te_idx_dict[i])
        ucounts.append(uq)

    return ucounts, mcounts, multi_count


def merge_weight_dict(EM_output, em_index):
    count_dic = {}
    for i in range(len(EM_output)):
        k = em_index[i]
        v = EM_output[i]
        count_dic[k] = v
    return count_dic


def make_list_for_bed_grph(te_s_weights, te_s_read_ref, te_as_weights, te_as_read_ref, array_count):
    idx = array_count - 1
    result_array = np.empty((array_count, 4), dtype='object')
    x = 'sense_TE'
    y = 'anti_TE'
    z = 'miRNA'
    # turn these into a numpy array indexed
    for k, v in list(te_s_read_ref.items()):
        read_num = len(v)
        weight_feature = te_s_weights[k]
        read_weight = weight_feature / read_num
        for i in te_s_read_ref[k]:
            result_array[idx] = [x, k, i, read_weight]
            idx -= 1
    for k, v in list(te_as_read_ref.items()):
        read_num = len(v)
        weight_feature = te_as_weights[k]
        read_weight = weight_feature / read_num
        for i in te_as_read_ref[k]:
            result_array[idx] = [y, k, i, read_weight]
            idx -= 1
            # result_array = np.append(result_array, [q], axis=0)
#    for k, v in list(mir_read_ref.items()):
#        read_num = len(v)
#        weight_feature = mir_weights[k]
#        read_weight = weight_feature / read_num
#        for i in mir_read_ref[k]:
#            result_array[idx] = [z, k, i, read_weight]
#            idx -= 1
            # result_array = np.append(result_array, [q], axis=0)
    weight_df = pd.DataFrame(result_array, columns=['ftype', 'fid', 'rid', 'weight'])
    return weight_df


def handle_annotation_EM(trf_free_multi, annot_dir, order):
    anno_fi, te_s_read_list, te_s_index, te_s_len, te_as_read_list, te_as_index, te_as_len, mir_read_list, mir_index, mir_len, te_s_ref_dic, te_as_ref_dic, mir_ref_dic, array_count = annotate_reads(trf_free_multi, annot_dir, order)
    sys.stderr.write("Annotation done. Fixing counts for EM\n")
#    mir_uniquecs, mir_multics, mir_multi_list = getcountsEM(mir_read_list, mir_index)
#    te_as_uniquecs, te_as_multics, te_as_multi_list = getcountsEM(te_as_read_list, te_as_index)
#    te_s_uniquecs, te_s_multics, te_s_multi_list = getcountsEM(te_s_read_list, te_s_index)
    # make function to not hardcode average read length to 20
#    EM_te_as_multi = EMestimate(te_as_len, te_as_multi_list, te_as_index, te_as_uniquecs, te_as_multics, 100, 20)
#    EM_te_s_multi = EMestimate(te_s_len, te_s_multi_list, te_s_index, te_s_uniquecs, te_s_multics, 100, 20)
#    EM_mir_multi = EMestimate(mir_len, mir_multi_list, mir_index, mir_uniquecs, mir_multics, 100, 20)
#    EM_te_as = list(map(operator.add, te_as_uniquecs, EM_te_as_multi))
#    EM_te_s = list(map(operator.add, te_s_uniquecs, EM_te_s_multi))
#    EM_mir = list(map(operator.add, mir_uniquecs, EM_mir_multi))
#    sys.stderr.write("EM done. Matching weights to features\n")
#    mir_weights = merge_weight_dict(EM_mir, mir_index)
#    te_s_weights = merge_weight_dict(EM_te_s, te_s_index)
#    te_as_weights = merge_weight_dict(EM_te_as, te_as_index)
#    em_weight_df_for_bed = make_list_for_bed_grph(te_s_weights, te_s_ref_dic, te_as_weights, te_as_ref_dic, mir_weights, mir_ref_dic, array_count)
#    em_weight_df_for_bed = make_list_for_bed_grph(te_s_weights, te_s_ref_dic, te_as_weights, te_as_ref_dic, array_count)
#    return anno_fi, em_weight_df_for_bed, mir_weights, te_s_weights, te_as_weights
#    return anno_fi, em_weight_df_for_bed, te_s_weights, te_as_weights
    return anno_fi
