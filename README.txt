
###################################################################
# this is a work in progress, appropriate at your own risk as it 
# has yet to be tested fully!
###################################################################


Goals
1. This set of programs will allow you to accurately place snps on a genome.

2. You can conduct an initial filtering of the .sam file, splitting alignments
to one location from alignments to two+ locations and discards alignments to 0 locations.
This is because you may wish to treat these cases in different manners.


3. Program will determine the bp position for the snp on the chromosome/contig of
the reference genome you have aligned it to (using the .sam information)

4. The output can be turned into a .vcf (I may make the program output this directly).



Inputs:
1. The program requires a .sam file, aligning short  query sequence reads that
contain snps (something like gbs or radseq data) to a subject genome
(of either chromosomes, scaffolds or contigs).

2. The program requires a file of all the query sequences, that contains tab
delimited information on their name, major/minor alleles, bp location and sequence
looks like the following:

SNP	Polymorphism_query_hit	bp_SNP_location	Sequence
TP10000	G/T	24	TGCATATGGCTCATCACAAATAGGCAGAAAAAATGTTGCAGGTGGAGCATCACATGCA
TP10002	A/C	51	TGCATATGGCTCTCCTATTCTTTGCCCAGTCATATTCAAGGTTAGAACTAAATTTCTAGGGTTC




