"""
Check the requirements
"""

# from contextlib import closing
import logging
import os
from os.path import basename, expanduser, isdir, isfile, join
import sys
# import shutil
import subprocess
import tarfile
import re
#from distutils.spawn import find_executable
# import urllib2

def check_requirements(folder,genome):
    WHOLE_GENOME = join(folder, "genomes/{0}/sequence")
    ANNOTATION = join(folder, "genomes/{0}/annotation")
    BOWTIE_INDEX = join(folder, "genomes/{0}/sequence/bowtie_index")

    logging.info("Checking if reference genome and annotation files exist...")
    if (isdir(folder) and
        isdir(WHOLE_GENOME.format(genome)) and
        isfile(join(WHOLE_GENOME.format(genome), "genome.fa")) and
        isfile(join(WHOLE_GENOME.format(genome), "genome.fa.fai")) and
        isfile(join(WHOLE_GENOME.format(genome), "rDNA.fa")) and
        isfile(join(WHOLE_GENOME.format(genome), "rDNA.fa.fai")) and
        isfile(join(WHOLE_GENOME.format(genome), "tDNA.fa")) and
        isfile(join(WHOLE_GENOME.format(genome), "tDNA.fa.fai")) and
        isdir(ANNOTATION.format(genome)) and
        isfile(join(ANNOTATION.format(genome), "structural_RNA.bed")) and
        isfile(join(ANNOTATION.format(genome), "miRNA.bed")) and
        isfile(join(ANNOTATION.format(genome), "hairpin.bed")) and
        isfile(join(ANNOTATION.format(genome), "exon.bed")) and
        isfile(join(ANNOTATION.format(genome), "TE.bed")) and
        isfile(join(ANNOTATION.format(genome), "intron.bed")) and
        isfile(join(ANNOTATION.format(genome), "piRNA_cluster.bed")) and
        isdir(BOWTIE_INDEX.format(genome)) and
        (isfile(join(BOWTIE_INDEX.format(genome), "genome.1.ebwt")) or (isfile(join(BOWTIE_INDEX.format(genome), "genome.1.ebwtl")))) and
        (isfile(join(BOWTIE_INDEX.format(genome), "genome.2.ebwt")) or (isfile(join(BOWTIE_INDEX.format(genome), "genome.2.ebwtl")))) and
        (isfile(join(BOWTIE_INDEX.format(genome), "genome.3.ebwt")) or (isfile(join(BOWTIE_INDEX.format(genome), "genome.3.ebwtl")))) and
        (isfile(join(BOWTIE_INDEX.format(genome), "genome.4.ebwt")) or (isfile(join(BOWTIE_INDEX.format(genome), "genome.4.ebwtl")))) and
        (isfile(join(BOWTIE_INDEX.format(genome), "genome.rev.1.ebwt")) or (isfile(join(BOWTIE_INDEX.format(genome), "genome.rev.1.ebwtl")))) and
        (isfile(join(BOWTIE_INDEX.format(genome), "genome.rev.2.ebwt")) or (isfile(join(BOWTIE_INDEX.format(genome), "genome.rev.2.ebwtl")))) and
        isfile(join(BOWTIE_INDEX.format(genome), "rDNA.1.ebwt")) and
        isfile(join(BOWTIE_INDEX.format(genome), "rDNA.2.ebwt")) and
        isfile(join(BOWTIE_INDEX.format(genome), "rDNA.3.ebwt")) and
        isfile(join(BOWTIE_INDEX.format(genome), "rDNA.4.ebwt")) and
        isfile(join(BOWTIE_INDEX.format(genome), "rDNA.rev.1.ebwt")) and
        isfile(join(BOWTIE_INDEX.format(genome), "rDNA.rev.2.ebwt")) and
        isfile(join(BOWTIE_INDEX.format(genome), "tDNA.1.ebwt")) and
        isfile(join(BOWTIE_INDEX.format(genome), "tDNA.2.ebwt")) and
        isfile(join(BOWTIE_INDEX.format(genome), "tDNA.3.ebwt")) and
        isfile(join(BOWTIE_INDEX.format(genome), "tDNA.4.ebwt")) and
        isfile(join(BOWTIE_INDEX.format(genome), "tDNA.rev.1.ebwt")) and
        isfile(join(BOWTIE_INDEX.format(genome), "tDNA.rev.2.ebwt"))) :
        logging.info("Genome and annotation files present")
        return True

    return False

def get_requirements(folder,genome):
    if (folder != "NULL"):
        TESMALL = folder
        TESMALL = TESMALL.replace('genomes','')
        TESMALL = re.sub('\/\/','/',TESMALL)
    else:
        HOME = expanduser("~")
        TESMALL = join(HOME, "TEsmall_db")

    if not check_requirements(TESMALL,genome):
        if not isdir(join(TESMALL, "genomes")):
            os.makedirs(join(TESMALL, "genomes"))
        logging.info("Downloading reference genome and annotation files...")
        url = {"hg19": "https://www.dropbox.com/scl/fo/nllcuktpoew8diu9kyf66/AIi9PSIc4cU0KTynNcDyYHg/hg19.tar.gz?rlkey=ucnup95d36kzvokcwiy0adb1l&dl=1",
#               "dm3": "http://labshare.cshl.edu/shares/mhammelllab/www-data/TEsmall/dm3.tar.gz",
#               "mm9": "http://labshare.cshl.edu/shares/mhammelllab/www-data/TEsmall/mm9.tar.gz",
               "hg38": "https://www.dropbox.com/scl/fo/nllcuktpoew8diu9kyf66/AMoleu7HsbKILZup2F6s1qo/hg38.tar.gz?rlkey=ucnup95d36kzvokcwiy0adb1l&dl=1",
               "dm6": "https://www.dropbox.com/scl/fo/nllcuktpoew8diu9kyf66/AACcdDFYG43-JhTg_01nwBo/dm6.tar.gz?rlkey=ucnup95d36kzvokcwiy0adb1l&dl=1",
               "mm10": "https://www.dropbox.com/scl/fo/nllcuktpoew8diu9kyf66/AEPIpJpR3CRb0rOa-TVHwMk/mm10.tar.gz?rlkey=ucnup95d36kzvokcwiy0adb1l&dl=1",
               "mm39": "https://www.dropbox.com/scl/fo/nllcuktpoew8diu9kyf66/AHY6vU2jl2T2XftO6bUcgDo/mm39.tar.gz?rlkey=ucnup95d36kzvokcwiy0adb1l&dl=1",
               "GRCz11": "https://www.dropbox.com/scl/fo/nllcuktpoew8diu9kyf66/AJerINAfNI7ZnmMuxeqfRtc/GRCz11.tar.gz?rlkey=ucnup95d36kzvokcwiy0adb1l&dl=1",
               "T2Tv2": "https://www.dropbox.com/scl/fo/nllcuktpoew8diu9kyf66/AKOaCKpfiYSxXjPgACXlQ8Y/T2Tv2.tar.gz?rlkey=ucnup95d36kzvokcwiy0adb1l&dl=1"
               }
        if genome not in url:
              logging.info("Reference genome not available for download")
              logging.info("Please contact mghcompbio@gmail.com for info")
              sys.exit(1)

        outname = (genome,"tar","gz")
        outfile=".".join(outname)
        path = join(TESMALL, "genomes")
        filepath = join(path, outfile)
        cmd = "wget --no-check-certificate \"{0}\" -O {1}".format(url[genome], filepath)
#        print(cmd)
        subprocess.call(cmd, shell=True)

        tar = tarfile.open(filepath, "r:gz")
        tar.extractall(path)
        tar.close()
        os.remove(filepath)

    WHOLE_GENOME = join(TESMALL, "genomes/{0}/sequence").format(genome)
    ANNOTATION = join(TESMALL, "genomes/{0}/annotation").format(genome)
    BOWTIE_INDEX = join(TESMALL, "genomes/{0}/sequence/bowtie_index").format(genome)

    return BOWTIE_INDEX, ANNOTATION
