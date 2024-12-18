# Analyse reads at defined SNP positions and determine genotype abundances

## Installation

To get the script, just clone the git repository via `git clone https://github.com/m4rku5-5/IGR-caller.git`

It is written in python and has only two dependencies, namely `biopython` and `muscle` (http://www.drive5.com/muscle). 
It is recommended to use a `conda` ([Anaconda](https://anaconda.org/) or [miniconda](https://conda.io/miniconda.html)) environment for this. To install those dependencies issue the following commands:
```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -n IGR_caller python=3.6 biopython muscle
```

You can then activate the environment via `conda activate IGR_caller` to use the tool. When finished deactivate the environment by `conda deactivate` to get back to your normal shell.

## Usage

The tool has pretty basic usage: `python analyse_IGR_multithread.py merged_reads.fastq IGR_strains.fasta`
It takes in a fastq file with merged Illumina reads, while quality values are ignored (it is suggested to do read-merging with [PEAR](https://cme.h-its.org/exelixis/web/software/pear/)) as well as a fasta file with all reference genotype sequences. 
SNP positions are hard coded in this version for now, but this could easily be changed. They are expressed as relative nucleotide positions of the reference genotype sequences.
The tool is designed to leverage all multithreading capabilities of the system.

## How it works 

The tool iterates over each read in the fastq file and executes the following procedure for each one:
* align the read to all reference genotype sequences with muscle
* get the offset of the reference sequences, as they align somewhere in the middle of the read
* determine if there are any insertions or deletions in the alignment
	* if there is an insertion in the read, then shift the SNP positions to match up with the reference again
* compare all SNP positions of the read with all the reference sequences
* if there is an exact match to one of the reference sequences, then get its strain name and increment its abundance by 1

## Info

The tool is part of the following publication: Fuhren, J., Schwalbe, M., Boekhorst, J., Rösch, C., Schols, H. A., & Kleerebezem, M. (2024). Prebiotic utilisation provides Lactiplantibacillus plantarum a competitive advantage in vitro, but is not reflected by an increased intestinal fitness. Gut Microbes, 16(1), 2338946. https://doi.org/10.1080/19490976.2024.2338946

