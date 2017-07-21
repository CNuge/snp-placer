

1. This set of programs allows you to accurately place snps on a genome.

2. You can conduct an initial filtering of the .sam file, splitting alignments
to one location from alignments to two+ locations and discards alignments to 0 locations.
This is because you may wish to treat these cases in different manners.

3. Program will determine the bp position for the snp on the chromosome/contig of
the reference genome you have aligned it to (using the .sam information)

4. The output can be turned into a .vcf (I may make the program output this directly).



Inputs:
1. The program requires a .sam file, aligning short query sequence reads that
contain snps (something like gbs or radseq data) to a subject genome
(of either chromosomes, scaffolds or contigs).

2. The program requires a file of all the query sequences, that contains tab
delimited information on their name, major/minor alleles, bp location and sequence
looks like the following:

SNP	Polymorphism_query_hit	bp_SNP_location	Sequence
TP10000	G/T	24	TGCATATGGCTCATCACAAATAGGCAGAAAAAATGTTGCAGGTGGAGCATCACATGCA
TP10002	A/C	51	TGCATATGGCTCTCCTATTCTTTGCCCAGTCATATTCAAGGTTAGAACTAAATTTCTAGGGTTC
Where for 'TP10000' there is a G at the 24th base pair of the sequence, and T
is the alternate allele at this location

What it does:

If one wants to place their snps onto a genome and make sure that they are in the
correct spot down to the base pair, then all the necessary information is in a .sam
file.

The .sam fields that act as 'keys' to this puzzle are the following:
	1. the Qname and Rname columns.
		- these are pretty clear, tell you which short sequence aligns to which reference location
	2. the POS column:
		- this give the bp number of the leftmost part of the alignment, letting us anchor the
		sequence in the right spot
	3. the flag column:
		- these are a little more cryptic 
		- a 0 or 256 indicated a forward alignment
		- a 16 or 272 indicates a reverse alignment 
		- if we have a reverse, our snp's base pair positon needs 
		to 'flip' along with the rest of the read!
	4. the CIGAR string:
		- this string of letters and numbers indicates the alignment details, insertions, 
		deltions, matches etc.
		examples:
			84M1S == 84 matches, 1 soft trim
			52M1D33M == 52 matches, 1 deletion, 33 matches
			85M == 85 matches
		- code in cigarParse applies the information in these strings to the alignment 
		in the correct manner

The program draws on the above data, to calculate the exact base pair on the contig where the 
snp falls, and adds this information to the input dataframe. 


Potential sources of error:
There are a few ways this script can fail that I've found, here I'll point them out and tell you the fix.

1. If your .sam file has extra columns on the right, you need to add these in to the sam_header list on line 75

2. If your snp names are numeric and not strings (i.e. 76 not CAM_SNP_76) then move the # from line 60 to line 61
to make the necessary type change.



