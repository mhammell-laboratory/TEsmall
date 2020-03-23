## TEsmall

A pipeline for profiling TE-derived small RNAs

### Install Miniconda (Linux)

```
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
```

### Setup channels

```
$ conda config --add channels conda-forge
$ conda config --add channels bioconda
```

### Install TEsmall

```
$ git clone git@github.com:wwliao/tesmall.git
$ cd tesmall
$ CONDA_RESTORE_FREE_CHANNEL=1 conda env create -f environment.txt -n tesmall
$ conda activate tesmall
$ python setup.py install
```

### How to run TEsmall

1. Before executing TEsmall, make sure you have activated the environment

	```
	$ conda activate tesmall
	```

2. For example, you would like to apply TEsmall on 2 FASTQ files: `Parental_1.fastq.gz` and `DroKO_1.fastq.gz`

	```
	$ tesmall -f Parental_1.fastq.gz DroKO_1.fastq.gz -l Parental DroKO
	```

3. When it's done, deactivate the environment

	```
	$ conda deactivate
	```
4. If you would like to specify the directory to which the genomes TEsmall uses for annotation are downloaded and read from please use the export command as follows
	
	```
	$ source activate tesmall
	$ export TESMALLROOT=/your/desired/directory
	$ tesmall -f Parental_1.fastq.gz DroKO_1.fastq.gz -l Parental DroKO
	$ source deactivate
	```
### For more information

```
$ tesmall -h
usage: TEsmall [-h] [-a STR] [-m INT] [-M INT] [-g STR] [--maxaln INT]
               [--mismatch INT] [-o STR [STR ...]] [-p INT] [-f STR [STR ...]]
               [-l STR [STR ...]] [--verbose INT] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -a STR, --adapter STR
                        Sequence of an adapter that was ligated to the 3' end.
                        The adapter itself and anything that follows is
                        trimmed. (default: TGGAATTCTCGGGTGCCAAGG)
  -m INT, --minlen INT  Discard trimmed reads that are shorter than INT. Reads
                        that are too short even before adapter removal are
                        also discarded. (default: 16)
  -M INT, --maxlen INT  Discard trimmed reads that are longer than INT. Reads
                        that are too long even before adapter removal are also
                        discarded. (default: 36)
  -g STR, --genome STR  Version of reference genome (default: hg19)
  --maxaln INT          Suppress all alignments for a particular read if more
                        than INT reportable alignments exist for it. (default:
                        100)
  --mismatch INT        Report alignments with at most INT mismatches.
                        (default: 0)
  -o STR [STR ...], --order STR [STR ...]
                        Annotation priority. (default: structural_RNA miRNA
                        hairpin exon TE intron piRNA_cluster)
  -p INT, --parallel INT
                        Parallel execute by INT CPUs. (default: 1)
  -f STR [STR ...], --fastq STR [STR ...]
                        Input in FASTQ format. Compressed input is supported
                        and auto-detected from the filename extension (.gz).
  -l STR [STR ...], --label STR [STR ...]
                        Unique label for each sample.
  --verbose INT         Set verbose level. 0: only show critical message, 1:
                        show additional warning message, 2: show process
                        information, 3: show debug messages. DEFAULT:2
  -v, --version         show program's version number and exit
```


