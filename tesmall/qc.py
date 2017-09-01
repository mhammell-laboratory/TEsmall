import argparse
from collections import Counter
from distutils.spawn import find_executable
import os
import subprocess

from bokeh.charts import Bar, output_file, save
from bokeh.charts.attributes import cat
import pandas as pd
import pysam

def get_read_info(multibam):
    root = os.path.splitext(os.path.basename(multibam))[0]
    root = os.path.splitext(root)[0]
    output = "{0}.rinfo".format(root)
    with open(output, "w") as outfile:
        outfile.write("id\tlength\tcount\tmapper\n")
        fname = multibam
        rcounter = Counter()
        bamfile = pysam.AlignmentFile(fname, "rb", check_sq=False)
        for read in bamfile:
            rid = read.query_name
            rlen = read.query_length
            rcounter.update([(rid, rlen)])

        for rid, rlen in rcounter:
            num = rcounter[(rid, rlen)]
            outfile.write("{0}\t{1}\t{2}\t{3}\n".format(rid, rlen, num, "unique" if num == 1 else "multi"))
    return output

def plot_read_dist(rinfo):
    df = pd.read_table(rinfo)
    data = df[['length', 'mapper', 'count']].groupby(['length', 'mapper']).count()
    data = data.apply(lambda s: s/data.sum()*100, axis=1).reset_index()
    root = os.path.splitext(os.path.basename(rinfo))[0]
    outfile = "{0}.rdist.html".format(root)
    output_file(outfile)
    p = Bar(data, plot_width=600, plot_height=400,
        label='length', values='count', agg='sum',
        stack=cat(columns='mapper', sort=False), legend='top_right',
        xlabel='read length (nt)', ylabel='percent (%)',
        tooltips=[('length', '@length'), ('mapper', '@mapper'), ('percent', '@height')],
        logo=None)
    p.outline_line_color = None
    save(p)
    return outfile
