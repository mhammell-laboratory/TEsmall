"""Alignment module."""
import logging
from os.path import splitext, basename
import subprocess

def trim_3p_adapters(fastq, prefix, adapter, minlen, maxlen):
    """Find and remove adapter sequences from the input FASTQ file."""
    outfile = "{0}.trimmed1.fastq".format(prefix)
    logfile = "{0}.cutadapt1.log".format(prefix)
    command = 'cutadapt -a "{0}" -m {1} -M {2} -o "{3}" "{4}" > "{5}"'
    command = command.format(adapter, minlen, maxlen, outfile, fastq, logfile)
    logging.info("Trimming 3' adapters...")
    subprocess.call(command, shell=True)
    return outfile

def trim_5p_adapters(fastq, prefix, adapter, minlen, maxlen):
    """Find and remove adapter sequences from the input FASTQ file."""
    outfile = "{0}.trimmed2.fastq".format(prefix)
    logfile = "{0}.cutadapt2.log".format(prefix)
    command = 'cutadapt -g "{0}" -m {1} -M {2} -o "{3}" "{4}" > "{5}"'
    command = command.format(adapter, minlen, maxlen, outfile, fastq, logfile)
    logging.info("Trimming 5' adapters...")
    subprocess.call(command, shell=True)
    return outfile

def map_multi_reads(fastq, ebwt, maxaln=100, mm=0):
    """Align short reads to the reference genome in multi mode."""
    prefix = splitext(basename(fastq))[0]
    prefix = splitext(prefix)[0]
    prefix = splitext(prefix)[0]
    logfile = "{0}.log".format(prefix)
    unfile = "{0}.unaligned.fastq".format(prefix)
    exfile = "{0}.exceeded.fastq".format(prefix)
    command = 'bowtie --chunkmbs 1024 -a -m {0} --best --strata -v {1} -S '\
              '--un "{2}" --max "{3}" "{4}" "{5}" 2> "{6}" | samtools view -F 4 -buS - '\
              '2>> "{6}" | samtools sort -n -o "{7}.multi.bam" - 2>> "{6}"'
    command = command.format(maxaln, mm, unfile, exfile, ebwt, fastq, logfile, prefix)
    logging.info("Aligning reads to reference sequences...")
    subprocess.call(command, shell=True)
    return "{0}.multi.bam".format(prefix)

def map_bestone_reads(fastq, ebwt, mm=0):
    """Align short reads to the reference genome in bestone mode."""
    prefix = splitext(basename(fastq))[0]
    prefix = splitext(prefix)[0]
    logfile = "{0}.rRNA.log".format(prefix)
    unfile = "{0}.rm_rRNA.fastq".format(prefix)
    command = 'bowtie --chunkmbs 1024 -k 1 --best -v {0} -S --un "{1}" "{2}" "{3}" 2> '\
              '"{4}" | samtools view -F 4 -buS - 2>> "{4}" | samtools sort -n '\
              '-o "{5}.bam" - 2>> "{4}"'
    command = command.format(mm, unfile, ebwt, fastq, logfile, prefix)
    logging.info("Removing rRNA-derived reads...")
    subprocess.call(command, shell=True)
    return "{0}.bam".format(prefix), unfile

