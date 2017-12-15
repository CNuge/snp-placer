# SNP genome placement
[![Build Status](https://travis-ci.org/CNuge/snp_genome_placement.svg?branch=master)](https://travis-ci.org/CNuge/snp_genome_placement)

The goal of this program is to take a series of single nucleotide polymprphisms (snps) associated with short sequence reads, and determine their exact base pair position in the context of a larger genome. This can be useful for comparative study of snps identified in different studies/techniques as it orients the snps relative to a common reference. This can also be useful if you need more surrounding sequence information then your provided by the origninal short sequence reads. By placing the SNP in the correct location on the genome you are able to pull an exact amount of flanking sequence from either side of the snp. This can be useful in compartive genome analyses (inputs for BLAST) or for things like designing genotyping arrays, were an exact amount of flanking sequence on either side of the snp is required.

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

2. The program requires a .sam file as input, aligning the short sequence reads from 1. to a genome or 
(of either chromosomes, scaffolds or contigs). The .sam file(s) can be generated using  [samtools](https://github.com/samtools/samtools) or any equivalent workflow of the user's choice that generates alignment data in the .sam format.

See the folder [example_data](https://github.com/CNuge/snp_genome_placement/tree/master/example_data) for example input files.

## Requirements
The program is written in [python3.6](https://www.python.org/downloads/) and utilizes [pandas](https://pandas.pydata.org/), which can be obtained using the following command:

		pip install pandas

## Workflow


1. This program allows you to accurately place snps on a genome after you have conducted an alignment of short sequence reads against a genome. Assuming the position of the SNP on the short read is known, then you can use this to find the exact base pair on the reference chromosome/contig. For my use of this function the .sam alignment info was generated though a burrows-wheeler alignment. 

2. You can conduct an initial filtering of the .sam file, splitting alignments to one location from alignments to two+ locations. This is because you may wish to treat these cases in different manners.

3. The program will determine the bp position for the snp on the chromosome/contig of the reference genome you have aligned it to (using the .sam information)

4. The information is turned into a .vcf file for future use.

