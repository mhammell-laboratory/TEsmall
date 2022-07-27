import logging
import re
from os.path import splitext
import pandas as pd
from functools import reduce

# for library compostion not count table, dont need EM to calculate becasue reads aren't assigned to switch categories


def calc_composition(anno_list, cca_anno_list):
    for anno in anno_list:
        logging.info("Calculating read composition...")
        for cca_anno in cca_anno_list:
            if cca_anno_list.index(cca_anno) == anno_list.index(anno):
                df = pd.read_csv(anno, sep="\t", usecols=["rid", "rlen", "ftype"])
                df_cca = pd.read_csv(cca_anno, sep="\t", usecols=["rid", "rlen", "ftype"])
                df = pd.concat([df, df_cca], sort=True)
                df = df.groupby(["rlen", "ftype"]).rid.nunique()
                df = df.unstack("rlen")
        root = splitext(anno)[0]
        outfile = "{0}.comp".format(root)
        df.T.to_csv(outfile, sep="\t", na_rep=0, float_format="%.0f")  # df.T is a transpose

# for count table


def em_dic_to_df(em_dict, srna_class, file_root):
    if em_dict:
        df = pd.DataFrame.from_dict(em_dict, orient='index')
        df = pd.DataFrame.reset_index(df)
        te_name_df = pd.DataFrame(columns=['ftype'], index=list(range(len(em_dict))))
        te_name_df = te_name_df.fillna(srna_class)
        te_table = pd.concat([te_name_df, df], sort=True, axis=1)
        te_table.columns = ['ftype', 'fid', file_root]
        te_table.iloc[:, 2] = te_table.iloc[:, 2].round(10)
        #te_table.to_csv('/Users/koneill/Desktop/fix_tesmall/df_{0}_{1}.txt'.format(file_root, srna_class), sep=' ', mode='a', index=False)
        return te_table
    else:
        empty_df = pd.DataFrame(columns=['ftype', 'fid', file_root])  # in case you get an empty dictionary from the EM, there are no reads in that class
        return empty_df


# join EM derived read weight and abundance.py derived weights with anno file to make bedgraph file


def calc_abundance(em_in):
    dfs = []  # holds dfs from each sample to merge for final count table
    for anno, cca_anno, em_weight_df, mirna_dic, te_s_dic, te_as_dic in em_in:
        logging.info("Calculating feature abundances...")
        root = splitext(anno)[0]
        # create data frames from annotation files (1/n) and EM weighted dictionaries
        df = pd.read_csv(anno, sep="\t", usecols=["rid", "ftype", "fid"])  # get necessary columns
        df = df[(df["ftype"] != "anti_TE") & (df["ftype"] != "sense_TE") & (df["ftype"] != "miRNA")]  # filter out EM treated categories
        rweight = 1 / df.groupby("rid").fid.nunique()  # calculate 1/n weight per read
        ftable = df.groupby(["ftype", "fid"]).rid.unique()  # group on sRNA id in class
        count = ftable.apply(lambda l: round(sum([rweight[s] for s in l])))  # hierarchecal df with ftype and fid (fid contains fid and weight tab delim)
        temp_dic = {root: count}  # makes new key with sample name
        count = pd.DataFrame(temp_dic).reset_index()  # resets weight column name to sample root name
        te_s_table = em_dic_to_df(te_s_dic, 'sense_TE', root)  # get matching EM dfs
        te_as_table = em_dic_to_df(te_as_dic, 'anti_TE', root)
        mir_table = em_dic_to_df(mirna_dic, 'miRNA', root)

        em_count = pd.concat([te_s_table, te_as_table, mir_table], sort=True)  # put EM dfs together
        #em_count.to_csv('/Users/koneill/Desktop/fix_tesmall/em_merge_df_{0}.txt'.format(root), sep=' ', mode='a')
        em_count = em_count[["ftype", "fid", root]]  # order columns properly
        #em_count.to_csv('/Users/koneill/Desktop/fix_tesmall/em_merge_df_named_{0}.txt'.format(root), sep=' ', mode='a')

        cca_df = pd.read_csv(cca_anno, sep="\t", usecols=["rid", "ftype", "fid"])  # do everything from before on cca annotation
        cca_rweight = 1 / cca_df.groupby("rid").fid.nunique()
        cca_ftable = cca_df.groupby(["ftype", "fid"]).rid.unique()
        cca_count = cca_ftable.apply(lambda l: round(sum([cca_rweight[s] for s in l])))
        temp_cca_dic = {root: cca_count}  # fix weight column
        cca_count = pd.DataFrame(temp_cca_dic).reset_index()
        # put main (1/n) anno counts, EM weighted counts, and (1/n) 3' trf counts together
        count_out = pd.concat([count, em_count, cca_count], sort=True, ignore_index=True)
        count_out = count_out.loc[count_out[root] > 0.1, :]  # filter species with 0 counts
        count_out = count_out[['fid', 'ftype', root]]
        count_out = count_out.groupby(['fid','ftype'], as_index = False)[root].sum()
        count_out = count_out[['fid', 'ftype', root]]
        #count_out.to_csv("test_sample_df_final_{0}.txt".format(root), sep="\t", na_rep=0, float_format="%.0f", index=False)
        dfs.append(count_out)

        # make a bedgraph file too
        coor = pd.read_csv(anno, sep="\t", usecols=["rid", "ftype", "fid", 'rchr', 'rstart', 'rend'])  # get the right columns from anno file
        cca_coor = pd.read_csv(cca_anno, sep="\t", usecols=["rid", "ftype", "fid", 'rchr', 'rstart', 'rend'])
        bed = coor.merge(rweight.to_frame(), on='rid')  # get their already calculated associated weights
        bed.columns = ['rid', 'rchr', 'rstart', 'rend', 'ftype', 'fid', 'weight']  # name columns
        bed = bed[(bed["ftype"] != "anti_TE") & (bed["ftype"] != "sense_TE") & (bed["ftype"] != "miRNA")]  # filter out those that arent TEs or miRNAs since these weights are from EM
        em_bed = pd.merge(coor, em_weight_df, left_on=['ftype', 'fid', 'rid'], right_on=['ftype', 'fid', 'rid'])  # merge in EM weights
        cca_bed = cca_coor.merge(cca_rweight.to_frame(), on='rid')  # get 3' trf weights from parallel annotation counting
        cca_bed.columns = ['rid', 'rchr', 'rstart', 'rend', 'ftype', 'fid', 'weight']  # rename columns

        # merge bed dataframes to one big one
        for_bed_graph = pd.concat([bed, em_bed, cca_bed])
        for_bed_graph = for_bed_graph[['rchr', 'rstart', 'rend', 'weight']]  # pick the columns an actual .bed file uses
        for_bed_graph = for_bed_graph.sort_values(by=['rchr', 'rstart'])  # sort 'em for your viewer
        for_bed_graph = for_bed_graph.groupby(['rchr', 'rstart', 'rend'], as_index=False)['weight'].sum()  # sum within bins of read length denoted by annotation start stop
        for_bed_graph.to_csv("{0}.bedgraph".format(root), sep='\t', na_rep=0, index=False)  # export to file

    df_final = reduce(lambda left, right: pd.merge(left, right, how='outer', on=['ftype', 'fid']), dfs)  # merge dfs from each sample
    df_final = df_final.sort_values(by=['ftype'])
    df_final.to_csv("count_summary.txt", sep="\t", na_rep=0, float_format="%.0000f", index=False, index_label=None)

    return df_final, for_bed_graph
