# SNP genome placement
[![Build Status](https://travis-ci.org/CNuge/snp_genome_placement.svg?branch=master)](https://travis-ci.org/CNuge/snp_genome_placement)

## Overview
The goal of this program is to take a series of single nucleotide polymprphisms (snps) associated with short sequence reads, and determine their exact base pair position in the context of a larger genome. This can be useful for comparative study of snps identified in different studies/techniques as it orients the snps relative to a common reference. This can also be useful if you need more surrounding sequence information then provided by the origninal short sequence reads. By placing the SNP in the correct location on the genome, you are able to pull an exact amount of flanking sequence from either side of the snp for further analysis. This can be useful in compartive genome analyses (obtaining inputs for BLAST) or for things like designing genotyping arrays, were an exact amount of flanking sequence is needed on either side of the snp.

## Inputs
The program requires two inputs:

1. A list of short sequence reads with the information formatted as follows:
	* The name associated with the snp 
	* The alleles involved in the polymorphism (separated by a slash /) 
	* The base pair position of the snp. The first bp in the string is position 1.
	* The DNA sequence


```		
		SNP	Polymorphism	bp	Sequence
		CMN1211	G/T	24	TGCATATGGCTCATCACAAATACGCAGAAAAAATGTTGCAGGTGGAGCATCACATGCA
		CMN8988	A/C	51	TGCATATGGCTCTCCTATTCTTTGCCCAGTCATATTCAAGGTTAGAACTCATTTTCTAGGGTTC
```	
	
Where for the sequence 'CMN1211' there is a G at the 24th base pair of the sequence, and T is the alternate allele at this location.

2. A .sam file as input, aligning the short sequence reads from 1. to a genome (or simply larger contigs). The .sam file(s) can be generated using  [samtools](https://github.com/samtools/samtools) or any equivalent workflow of the user's choice that generates alignment data in the .sam format.

See the folder [example_data](https://github.com/CNuge/snp_genome_placement/tree/master/example_data) for example input files.

## Requirements
The program is written in [python3.6](https://www.python.org/downloads/) and utilizes [pandas](https://pandas.pydata.org/), which can be obtained by typing the following into the command line: `pip install pandas`

## Workflow

To use the program, you first filter the input sam file using `filter_sam_file.py` and then pass the output(s) from the filter step to `place_snps.py` along with the short sequence read information files.

### `filter_sam_file.py` - Filters a .sam file retaining only the required data
You can conduct an initial filtering of the .sam file, splitting the file into sequences that align to one location (default output name `one_location_alignments.sam`) and sequences that align to two+ locations (default output name `multiple_location_alignments.sam`. This is done because it is logical to treat these cases differently. Short sequence reads aligning to 2+ inherently introduce uncertainty into the exact bp location of snps, so these should be utilized with caution. 

`filter_sam.py` can be run directly from the command line using the follwing syntax

```		
python filter_sam.py input_sam_file.sam
```

The `input_sam_file.sam` data will be filtered and data is sent to the one location(`one_location_alignments.sam`) and two+ locations (`multiple_location_alignments.sam`) output files. Sequences aligning to no locations are discarded.

The output file names can be changed using the `-p` flag to indicate the one location output filename and the `-s` flag to indicate the  two+ locations output filename. For example:
```
python filter_sam.py input_sam_file.sam -p example_name_primary.sam -s example_name_secondary.sam
```

### `place_snps.py` - Place the snps from the short sequences into the correct position on the genome

Take a list of short sequence reads and filtered .sam file and accurately place the snps on the short reads in the context of the larger genome (the .sam file's subject sequences, the short reads are the .sam file's query sequences). Required inputs are one or more .sam files, one or more snp information files in .txt format.
```
python place_snps.py -s example_name_primary.sam -p short_snp_seqs.txt
```
The output is in .vcf format. Default output name is `placed_snps.vcf`, change this with the `-o` flag.
```
python place_snps.py -s example_name_primary.sam -p short_snp_seqs.txt -o example_output_name.vcf
```
Multiple files can be passed after the two input flags
```
python place_snps.py -s example_name_primary_a.sam example_name_primary_b.sam -p short_snp_seqs_a.txt short_snp_seqs_b.txt 
```
