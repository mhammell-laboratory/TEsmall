import logging
import re
from os.path import splitext
import pandas as pd

def calc_composition(anno_list):
    for anno in anno_list:
        logging.info("Calculating read composition...")
        df = pd.read_table(anno, usecols=["rid", "rlen", "ftype"])
        df = df.groupby(["rlen", "ftype"]).rid.nunique()
        df = df.unstack("rlen")

        root = splitext(anno)[0]
        outfile = "{0}.comp".format(root)
        df.T.to_csv(outfile, sep="\t", na_rep=0, float_format="%.0f")

def calc_abundance(anno_list):
    count_dict = {}
    for anno in anno_list:
        logging.info("Calculating feature abundances...")
        df = pd.read_table(anno, usecols=["rid", "ftype", "fid"])
        rweight = 1/df.groupby("rid").fid.nunique()
        ftable = df.groupby(["ftype", "fid"]).rid.unique()
        count_dict[splitext(anno)[0]] = ftable.apply(lambda l: round(sum([rweight[s] for s in l])))

        count = pd.DataFrame(count_dict).reset_index()
        #strcount = count[count["ftype"] == "structural_RNA"]
        #count = count[count["ftype"] != "structural_RNA"]
        #strcount.to_csv("structural_RNA.count.txt", sep="\t", na_rep=0, float_format="%.0f", index=False)
        tecount = count[(count["ftype"] == "anti_TE") | (count["ftype"] == "sense_TE")]
        count = count[(count["ftype"] != "anti_TE") & (count["ftype"] != "sense_TE")]
        tecount["fid"] = [re.sub(r"_dup\d+$", "", s) for s in tecount["fid"]]
        tecount = tecount.groupby(["ftype", "fid"]).sum().reset_index()

        count = pd.concat([count, tecount])
        root = splitext(anno)[0]
        count.to_csv("{0}.count".format(root), sep="\t", na_rep=0, float_format="%.0f", index=False)
