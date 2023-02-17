## TEsmall

Version 2.0.2

A pipeline for profiling TE-derived small RNAs.

Created by Wen-Wei Liao, Kat O'Neill & Molly Gale Hammell, March 2017

Contact: Oliver Tam (tam@cshl.edu) or Molly Gale Hammell (mhammell@cshl.edu)

### Install Miniconda 3 (Linux)

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
$ git clone https://github.com/mhammell-laboratory/TEsmall.git
$ cd TEsmall
$ conda env create -f environment.yaml -n TEsmall
$ conda activate TEsmall
$ python setup.py install
```

### How to run TEsmall

1. Before executing TEsmall, make sure you have activated the environment

	```
	$ conda activate TEsmall
	```

2. For example, you would like to apply TEsmall on 2 FASTQ files: `Parental_1.fastq.gz` and `DroKO_1.fastq.gz`

	```
	$ TEsmall -f Parental_1.fastq.gz DroKO_1.fastq.gz -l Parental DroKO
	```

3. When it's done, deactivate the environment

	```
	$ conda deactivate
	```
4. If you would like to specify the directory to which the genomes
   TEsmall uses for annotation are downloaded and read from, you can
   specify it at runtime using the `--dbfolder` parameter
	
	```
	$ TEsmall -f Parental_1.fastq.gz DroKO_1.fastq.gz -g hg19 -l
	Parental DroKO --dbfolder /path/to/another/folder/
	```
	The files used by TEsmall will be downloaded to/access from the
	`genomes` folder inside `/path/to/another/folder/`.
	
	The default location is `$HOME/TEsmall_db/`

### For more information

```
$ TEsmall -h
usage: TEsmall [-h] [-a STR] [-m INT] [-M INT] [-g STR] [--maxaln INT]
               [--mismatch INT] [-o STR [STR ...]] [-p INT] [-f STR [STR ...]]
               [-l STR [STR ...]] [--dbfolder STR] [--verbose INT] [-v]

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
  -g STR, --genome STR  Version of reference genome (hg38, hg19, mm10, mm39 dm6 or GRCz11; default: hg38)
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
  --dbfolder STR        Custom location of TEsmall database folder (containing the "genomes" folder). 
						DEFAULT: $HOME/TEsmall_db/

  --verbose INT         Set verbose level. 
                        0: only show critical message
						1: show additional warning message
						2: show process information
						3: show debug messages.
						DEFAULT: 2
  -v, --version         show program's version number and exit
```

### Copying & distribution

TEsmall is part of [TEToolkit suite](http://hammelllab.labsites.cshl.edu/software/).

TEsmall is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but *WITHOUT ANY WARRANTY*; without even the implied warranty of
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE*.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with TEsmall.  If not, see [this website](http://www.gnu.org/licenses/).

### Citation

If using the software in a publication, please cite the [following](https://pubmed.ncbi.nlm.nih.gov/30349559/):

O'Neill K, Liao WW, Patel A, Hammell MG. (2018) TEsmall Identifies Small RNAs Associated With Targeted Inhibitor Resistance in Melanoma. Front Genet. Oct 5;9:461.
