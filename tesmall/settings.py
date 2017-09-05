"""
Check the requirements
"""

# from contextlib import closing
import os
from os.path import basename, expanduser, isdir, isfile, join
# import shutil
import subprocess
import tarfile
#from distutils.spawn import find_executable
# import urllib2

HOME = expanduser("~")
TESMALL = join(HOME, ".tesmall")
WHOLE_GENOME = join(TESMALL, "genomes/{0}/sequence")
ANNOTATION = join(TESMALL, "genomes/{0}/annotation")
BOWTIE_INDEX = join(TESMALL, "genomes/{0}/sequence/bowtie_index")


def check_executions():
    try:
        subprocess.check_call(["cutadapt", "-h"])
    except OSError:
        print "Could not find cutadapt, please reinstall."

    try:
        subprocess.check_call(["bowtie", "-h"])
    except OSError:
        print "Could not find bowtie, please reinstall."

    try:
        subprocess.check_call(["samtools"])
    except OSError:
        print "Could not find samtools, please reinstall."

    try:
        subprocess.check_call(["bedtools"])
    except OSError:
        print "Could not find bedtools, please reinstall."


def check_requirements(genome):
    if (isdir(TESMALL) and
       isdir(WHOLE_GENOME.format(genome)) and
       isfile(join(WHOLE_GENOME.format(genome), "genome.fa")) and
       isfile(join(WHOLE_GENOME.format(genome), "genome.fa.fai")) and
       isdir(ANNOTATION.format(genome)) and
       isfile(join(ANNOTATION.format(genome), "structural_RNA.bed")) and
       isfile(join(ANNOTATION.format(genome), "miRNA.bed")) and
       isfile(join(ANNOTATION.format(genome), "hairpin.bed")) and
       isfile(join(ANNOTATION.format(genome), "exon.bed")) and
       isfile(join(ANNOTATION.format(genome), "TE.bed")) and
       isfile(join(ANNOTATION.format(genome), "intron.bed")) and
       isfile(join(ANNOTATION.format(genome), "piRNA_cluster.bed")) and
       isdir(BOWTIE_INDEX.format(genome)) and
       isfile(join(BOWTIE_INDEX.format(genome), "genome.1.ebwt")) and
       isfile(join(BOWTIE_INDEX.format(genome), "genome.2.ebwt")) and
       isfile(join(BOWTIE_INDEX.format(genome), "genome.3.ebwt")) and
       isfile(join(BOWTIE_INDEX.format(genome), "genome.4.ebwt")) and
       isfile(join(BOWTIE_INDEX.format(genome), "genome.rev.1.ebwt")) and
       isfile(join(BOWTIE_INDEX.format(genome), "genome.rev.2.ebwt"))):
        return True
    return False


def get_requirements(genome):
    if not check_requirements(genome):
        url = "https://www.dropbox.com/s/2mclhrejogig07u/{0}.tar.gz"
        url = url.format(genome)
        path = join(TESMALL, "genomes")
        filepath = join(path, basename(url))
        """
        with closing(urllib2.urlopen(url)) as infile:
            try:
                os.makedirs(path)
            except OSError:
                if not os.path.isdir(path):
                    raise
            with open(filepath, "wb") as outfile:
                shutil.copyfileobj(infile, outfile)
        """

        cmd = "wget --no-check-certificate {0} -O {1}".format(url, filepath)
        subprocess.call(cmd, shell=True)

        tar = tarfile.open(filepath, "r:gz")
        tar.extractall(path)
        tar.close()
        os.remove(filepath)
