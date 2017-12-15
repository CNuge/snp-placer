# SNP genome placement
[![Build Status](https://travis-ci.org/CNuge/snp_genome_placement.svg?branch=master)](https://travis-ci.org/CNuge/snp_genome_placement)

The goal of this program is to take a series of single nucleotide polymprphisms (snps) associated with short sequence reads, and determine their exact base pair position in the context of a larger genome. This can be useful for comparative study of snps identified in different studies/techniques as it orients the snps relative to a common reference. This can also be useful if you need more surrounding sequence information then provided by the origninal short sequence reads. By placing the SNP in the correct location on the genome, you are able to pull an exact amount of flanking sequence from either side of the snp for further analysis. This can be useful in compartive genome analyses (obtaining inputs for BLAST) or for things like designing genotyping arrays, were an exact amount of flanking sequence is needed on either side of the snp.

## Inputs:
The program requires two inputs:

1. A list of short sequence reads with the following information:
		- the name associated with the snp
		- the alleles of the polymorphism (separated by a slash /)
		- the base pair

		SNP	Polymorphism	bp	Sequence
		CMN1211	G/T	24	TGCATATGGCTCATCACAAATACGCAGAAAAAATGTTGCAGGTGGAGCATCACATGCA
		CMN8988	A/C	51	TGCATATGGCTCTCCTATTCTTTGCCCAGTCATATTCAAGGTTAGAACTCATTTTCTAGGGTTC
	
	Where for the sequence 'CMN1211' there is a G at the 24th base pair of the sequence, and T is the alternate allele at this location.

2. A .sam file as input, aligning the short sequence reads from 1. to a genome (or simply larger contigs). The .sam file(s) can be generated using  [samtools](https://github.com/samtools/samtools) or any equivalent workflow of the user's choice that generates alignment data in the .sam format.

See the folder [example_data](https://github.com/CNuge/snp_genome_placement/tree/master/example_data) for example input files.

## Requirements
The program is written in [python3.6](https://www.python.org/downloads/) and utilizes [pandas](https://pandas.pydata.org/), which can be obtained using the following command:

		pip install pandas

## Workflow

### `filter_sam_file.py` - Filters a .sam file retaining only the required data
You can conduct an initial filtering of the .sam file, splitting alignments to one location from alignments to two+ locations. This is because you may wish to treat these cases in different manners.

### `place_snps.py` - place snps from the short sequences into the genome

Take a list of short sequence reads and filtered .sam file and accurately place the snps on the short reads in the context of the larger genome (the .sam file's subject sequences, the short reads are the .sam file's query sequences).

The information is turned into a .vcf file for future use.

