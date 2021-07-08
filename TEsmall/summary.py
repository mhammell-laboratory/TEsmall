import datetime
import argparse
import logging
import os
import string
from math import pi, cos, sin
from collections import defaultdict, Counter
import pandas as pd
import bokeh
from bokeh.layouts import column, row, gridplot
from bokeh.plotting import figure, output_file, save
from bokeh.models import HoverTool, ColumnDataSource, Title
from bokeh.models.widgets import DataTable, TableColumn, NumberFormatter
from bokeh.embed import components
import seaborn as sns
import pysam

bokeh_version = bokeh.__version__

header1 = """
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <!-- The above 3 meta tags *must* come first in the head; any other head content must come *after* these tags -->
    <meta name="description" content="">
    <meta name="author" content="">
    <link rel="icon" href="../../favicon.ico">

    <title>TEsmall Dashboard</title>

    <!-- Bootstrap core CSS -->
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css" integrity="sha384-1q8mTJOASx8j1Au+a5WDVnPi2lkFfwwEAa8hDDdjZlpLegxhjVME1fgjWPGmkzs7" crossorigin="anonymous">
"""

bokeh_import = ["<link rel=\"stylesheet\" href=\"https://cdn.pydata.org/bokeh/release/bokeh-{0}.min.css\" type=\"text/css\" />".format(bokeh_version)]
bokeh_import.append("<link rel=\"stylesheet\" href=\"https://cdn.pydata.org/bokeh/release/bokeh-widgets-{0}.min.css\" type=\"text/css\" />".format(bokeh_version))
bokeh_import.append("<script type=\"text/javascript\" src=\"https://cdn.pydata.org/bokeh/release/bokeh-{0}.min.js\"></script>".format(bokeh_version))
bokeh_import.append("<script type=\"text/javascript\" src=\"https://cdn.pydata.org/bokeh/release/bokeh-widgets-{0}.min.js\"></script>".format(bokeh_version))
bokeh_import.append("<script type=\"text/javascript\" src=\"https://cdn.pydata.org/bokeh/release/bokeh-tables-{0}.min.js\"></script>".format(bokeh_version))
bokeh_importLines = "\n".join(bokeh_import)

header2 = """
    <style>
      /*
       * Base structure
       */

      /* Move down content because we have a fixed navbar that is 50px tall */
      body {
        padding-top: 50px;
      }


      /*
       * Global add-ons
       */

      .sub-header {
        padding-bottom: 10px;
        border-bottom: 1px solid #eee;
      }

      /*
       * Top navigation
       * Hide default border to remove 1px line.
       */
      .navbar-fixed-top {
        border: 0;
      }

      /*
       * Sidebar
       */

      /* Hide for mobile, show later */
      .sidebar {
        display: none;
      }
      @media (min-width: 768px) {
        .sidebar {
          position: fixed;
          top: 51px;
          bottom: 0;
          left: 0;
          z-index: 1000;
          display: block;
          padding: 20px;
          overflow-x: hidden;
          overflow-y: auto; /* Scrollable contents if viewport is shorter than content. */
          background-color: #f5f5f5;
          border-right: 1px solid #eee;
        }
      }

      /* Sidebar navigation */
      .nav-sidebar {
        margin-right: -21px; /* 20px padding + 1px border */
        margin-bottom: 20px;
        margin-left: -20px;
      }
      .nav-sidebar > li > a {
        padding-right: 20px;
        padding-left: 20px;
      }
      .nav-sidebar > .active > a,
      .nav-sidebar > .active > a:hover,
      .nav-sidebar > .active > a:focus {
        color: #fff;
        background-color: #428bca;
      }


      /*
       * Main content
       */

      .main {
        padding: 20px;
      }
      @media (min-width: 768px) {
        .main {
          padding-right: 40px;
          padding-left: 40px;
        }
      }
      .main .page-header {
        margin-top: 0;
      }


      /*
       * Placeholder dashboard ideas
       */

      .placeholders {
        margin-bottom: 30px;
        text-align: center;
      }
      .placeholders h4 {
        margin-bottom: 0;
      }
      .placeholder {
        margin-bottom: 20px;
      }
      .placeholder img {
        display: inline-block;
        border-radius: 50%;
      }
    </style>

  </head>

  <body>

    <nav class="navbar navbar-inverse navbar-fixed-top">
      <div class="container-fluid">
        <div class="navbar-header">
          <a class="navbar-brand" href="#">TEsmall</a>
        </div>
        <div id="navbar" class="navbar-collapse collapse">
          <ul class="nav navbar-nav navbar-right">
"""
navh = """
          </ul>
        </div>
      </div>
    </nav>

    <div class="container-fluid">
      <div class="row">
        <div class="col-sm-3 col-md-2 sidebar">
          <ul class="nav nav-sidebar">
"""
navt = """
          </ul>
        </div>
"""

footer = """
      </div>
    </div>

    <!-- Bootstrap core JavaScript
    ================================================== -->
    <!-- Placed at the end of the document so the pages load faster -->
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>
    <script src="../../dist/js/bootstrap.min.js"></script>
    <script>
    window.addEventListener("hashchange", function() { scrollBy(0, -50) })
    $(".nav li").on("click", function() {
      $(".nav li").removeClass("active");
      $(this).addClass("active");});
    window.jQuery || document.write('<script src="../../assets/js/vendor/jquery.min.js"><\/script>')
    </script>
  </body>
</html>
"""

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

def get_stat(prefix, maxaln):
    stat = {"Statistics": [], "Number of reads": [], "Proportion": []}
    with open("{0}.cutadapt1.log".format(prefix)) as infile:
        for line in infile:
            if "Total reads processed:" in line:
                raw_reads = int(line.strip().split()[-1].replace(",", ""))
                stat["Statistics"].append("Raw reads")
                stat["Number of reads"].append(raw_reads)
                stat["Proportion"].append(None)
                break
    with open("{0}.rRNA.log".format(prefix)) as infile:
        for line in infile:
            if "# reads processed:" in line:
                trimmed_reads = int(line.strip().split()[-1])
                stat["Statistics"].append("After trimming adapters")
                stat["Number of reads"].append(trimmed_reads)
                stat["Proportion"].append(None)
    with open("{0}.log".format(prefix)) as infile:
        for line in infile:
            if "# reads processed:" in line:
                rm_reads = int(line.strip().split()[-1])
                stat["Statistics"].append("After removing rRNAs")
                stat["Number of reads"].append(rm_reads)
                stat["Proportion"].append(rm_reads/rm_reads)
            elif "# reads with at least one" in line:
                up_reads = int(line.strip().split()[-2])
                stat["Statistics"].append("Aligned reads (up to {0} alignments)".format(maxaln))
                stat["Number of reads"].append(up_reads)
                stat["Proportion"].append(up_reads/rm_reads)
            elif "# reads with alignments suppressed due to -m:" in line:
                over_reads = int(line.strip().split()[-2])
                stat["Statistics"].append("Over {0} alignments".format(maxaln))
                stat["Number of reads"].append(over_reads)
                stat["Proportion"].append(over_reads/rm_reads)
            elif "# reads that failed to align:" in line:
                un_reads = int(line.strip().split()[-2])
                stat["Statistics"].append("Unaligned reads".format(maxaln))
                stat["Number of reads"].append(un_reads)
                stat["Proportion"].append(un_reads/rm_reads)

    df = pd.read_table("{0}.anno".format(prefix), usecols=["rid"])
    anno_reads = len(df.rid.unique())
    stat["Statistics"].append("Annotated reads of aligned reads")
    stat["Number of reads"].append(anno_reads)
    stat["Proportion"].append(anno_reads/up_reads)

    return stat

def output_components(prefix, order, maxaln):
    rinfo = "{0}.rinfo".format(prefix)
    comp = "{0}.comp".format(prefix)

    def plot_read_dist(rinfo):
        df = pd.read_table(rinfo)
        data = df[['length', 'mapper', 'count']].groupby(['length', 'mapper']).count()
        data = data.apply(lambda s: s/data.sum()*100, axis=1).reset_index()
        unique = data.loc[data['mapper'] == "unique"]
        multi = data.loc[data['mapper'] == "multi"]
        data = unique.set_index('length').join(multi.set_index('length'),how='outer',rsuffix="_multi")
        data = data[['count','count_multi']]
        p = figure(x_range=(min(list(data.index)),max(list(data.index))),
                   plot_width=500,
                   plot_height=400)
        p.vbar_stack(['count','count_multi'],
                     source = data,
                     x = 'length',
                     width = 0.9,
                     color = ["#f98283","#a4a4a4"],
                     legend_label = ["Unique","Multi"],
#                     tooltips=[('length', '@length'), (['unique','multi'], ['@unique','@multi'])]
                     )
        p.xaxis.axis_label='read length (nt)'
        p.yaxis.axis_label='proportion (%)'
#        p.x_range.range_padding = 0.1
        p.ygrid.grid_line_color = None
        p.legend.location = "top_right"
        p.toolbar.logo = None
        p.outline_line_color = None
        return p

    rdist = plot_read_dist(rinfo)
    
    df = pd.read_table(comp, index_col=0)
    total = df.sum()
    total = total*100/total.sum()
    df = df.apply(lambda s: s*100/s.sum(), axis=1)
    ftypes = list(df.columns)
    df = df.reset_index()
    colors = sns.color_palette("hls", len(ftypes)).as_hex()
    bar = figure(x_range=rdist.x_range,
                 y_range = (0, 100),
                 plot_width=500,
                 plot_height=400)
    bar.vbar_stack(ftypes,
                   x='rlen',
                   color=colors,
                   width = 0.9,
                   source=df,
                   legend_label=ftypes,
#                   tooltips=[('length', '@rlen'), ('feature', '@ftypes'), ('percent', '@height')]
                   )
    bar.xaxis.axis_label='read length (nt)'
    bar.yaxis.axis_label='proportion (%)'
    bar.ygrid.grid_line_color = None
    bar.legend.location = "top_right"
#    bar.x_range.range_padding = 0.1
    bar.toolbar.logo = None
    bar.outline_line_color = None

    start_angles = {}
    end_angles = {}
    start = 0
    for ftype in ftypes:
        end = 2*pi*total[ftype]/100
        start_angles[ftype] = start
        end_angles[ftype] = start + end
        start += end

    colors = dict(list(zip(ftypes, colors)))
    df = pd.DataFrame(total).reset_index()
    df.columns = ["ftype", "percent"]
    df["start"] = df.apply(lambda s: start_angles[s["ftype"]], axis=1)
    df["end"] = df.apply(lambda s: end_angles[s["ftype"]], axis=1)
    df["color"] = df.apply(lambda s: colors[s["ftype"]], axis=1)
    df["x"] = df.apply(lambda s: 1.2*cos((start_angles[s["ftype"]] + end_angles[s["ftype"]])/2), axis=1)
    df["y"] = df.apply(lambda s: 1.2*sin((start_angles[s["ftype"]] + end_angles[s["ftype"]])/2), axis=1)
    df["text"] = df.apply(lambda s: "{0:.3f}%".format(s["percent"]), axis=1)
    source = ColumnDataSource(data=df)

    pie = figure(width=400, height=400, x_range=(-1.4, 1.4), y_range=(-1.4, 1.4))
    pie.toolbar.logo = None
    wr = pie.annular_wedge(x=0, y=0, inner_radius=0.5, outer_radius=1,
        start_angle="start", end_angle="end", fill_color="color",
        line_color="#ffffff", line_width=0.5, source=source,legend_field = 'ftype')


    pie.axis.visible = False
    pie.grid.grid_line_color = None
    pie.outline_line_color = None

    hover = HoverTool(tooltips=[("feature", "@ftype"), ("percent", "@percent")], renderers=[wr])
    pie.add_tools(hover)

    text_props = {
        "source": source,
        "angle": 0,
        "color": "black",
        "text_align": "center",
        "text_baseline": "middle",
        "text_font_size": "8pt",
        "text_font_style": "normal"
    }
    pie.text(x="x", y="y", text="text", **text_props)

    empty = figure(width=400, height=400, x_range=(-1.1, 1.1), y_range=(-1.1, 1.1))
    empty.toolbar.logo = None
    empty.axis.visible = False
    empty.grid.grid_line_color = None
    empty.outline_line_color = None
    empty.toolbar_location = None

    stat = get_stat(prefix, maxaln)
    source = ColumnDataSource(data=stat)
    columns = [
        TableColumn(field="Statistics", title="Statistics", width=200),
        TableColumn(field="Number of reads", title="Number of reads", formatter=NumberFormatter(format="0,0"), width=150),
        TableColumn(field="Proportion", title="Proportion", formatter=NumberFormatter(format="0.000%"), width=100),
    ]
    data_table = DataTable(source=source, columns=columns, width=450)

    script, div = components(gridplot([[data_table, rdist], [pie, bar]], sizing_mode="scale_width"))
    return script, div

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--prefix", nargs="+")
    parser.add_argument("-o", "--order", nargs="+")
    parser.add_argument("-m", "--maxaln", type=int)
    return parser

def gen_summary(prefix, order, maxaln):
    logging.info("Generating summary report...")
    idx = order.index("TE")
    order[idx] = "sense_TE"
    order.insert(idx, "anti_TE")
    date = datetime.datetime.now().strftime("%B %d, %Y")
    with open("report.html", "w") as outfile:
        outfile.write(header1)
        outfile.write(bokeh_importLines)
        outfile.write(header2)
        outfile.write('            <li><a href="#">{0}</a></li>'.format(date))
        outfile.write(navh)
        for i, p in enumerate(prefix):
            if i == 0:
                outfile.write('            <li class="active"><a href="#{0}">{0}</a></li>\n'.format(p))
            else:
                outfile.write('            <li><a href="#{0}">{0}</a></li>\n'.format(p))
        outfile.write(navt)
        for p in prefix:
            script, div = output_components(p, order, maxaln)
            outfile.write('    <div class="col-sm-9 col-sm-offset-3 col-md-10 col-md-offset-2 main">\n')
            outfile.write('      <h2 class="sub-header" id="{0}">{0}</h2>\n'.format(p))
            outfile.write(div + "\n")
            outfile.write(script + "\n")
            outfile.write("    </div>\n")
        outfile.write(footer)

if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    gen_summary(args.prefix, args.order, args.maxaln)
