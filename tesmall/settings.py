"""
Check the requirements
"""

# from contextlib import closing
import logging
import os
from os.path import basename, expanduser, isdir, isfile, join
# import shutil
import subprocess
import tarfile
#from distutils.spawn import find_executable
# import urllib2

try:
    HOME = os.environ["TESMALLROOT"]
except KeyError:
    HOME = expanduser("~")
TESMALL = join(HOME, ".tesmall")
WHOLE_GENOME = join(TESMALL, "genomes/{0}/sequence")
ANNOTATION = join(TESMALL, "genomes/{0}/annotation")
BOWTIE_INDEX = join(TESMALL, "genomes/{0}/sequence/bowtie_index")

def check_requirements(genome):
    logging.info("Checking if reference genome and annotation files exist...")
    if (isdir(TESMALL) and
       isdir(WHOLE_GENOME.format(genome)) and
       isfile(join(WHOLE_GENOME.format(genome), "genome.fa")) and
       isfile(join(WHOLE_GENOME.format(genome), "genome.fa.fai")) and
       isfile(join(WHOLE_GENOME.format(genome), "rDNA.fa")) and
       isfile(join(WHOLE_GENOME.format(genome), "rDNA.fa.fai")) and
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
       isfile(join(BOWTIE_INDEX.format(genome), "genome.rev.2.ebwt")) and
       isfile(join(BOWTIE_INDEX.format(genome), "rDNA.1.ebwt")) and
       isfile(join(BOWTIE_INDEX.format(genome), "rDNA.2.ebwt")) and
       isfile(join(BOWTIE_INDEX.format(genome), "rDNA.3.ebwt")) and
       isfile(join(BOWTIE_INDEX.format(genome), "rDNA.4.ebwt")) and
       isfile(join(BOWTIE_INDEX.format(genome), "rDNA.rev.1.ebwt")) and
       isfile(join(BOWTIE_INDEX.format(genome), "rDNA.rev.2.ebwt"))):
        return True
    return False

def get_requirements(genome):
    if not check_requirements(genome):
        if not isdir(join(TESMALL, "genomes")):
            os.makedirs(join(TESMALL, "genomes"))
        logging.info("Downloading reference genome and annotation files...")
        url = {"hg19": "http://labshare.cshl.edu/shares/mhammelllab/www-data/TEsmall/hg19.tar.gz",
               "dm3": "http://labshare.cshl.edu/shares/mhammelllab/www-data/TEsmall/dm3.tar.gz",
               "mm9": "http://labshare.cshl.edu/shares/mhammelllab/www-data/TEsmall/mm9.tar.gz",
               "hg38": "http://labshare.cshl.edu/shares/mhammelllab/www-data/TEsmall/hg38.tar.gz",
               "dm6": "http://labshare.cshl.edu/shares/mhammelllab/www-data/TEsmall/dm6.tar.gz",
               "mm10": "http://labshare.cshl.edu/shares/mhammelllab/www-data/TEsmall/mm10.tar.gz",
               }
        path = join(TESMALL, "genomes")
        filepath = join(path, basename(url[genome]))
        cmd = "wget --no-check-certificate {0} -O {1}".format(url[genome], filepath)
        subprocess.call(cmd, shell=True)

        tar = tarfile.open(filepath, "r:gz")
        tar.extractall(path)
        tar.close()
        os.remove(filepath)
