

"""
here we will read in the data from the .sam file
and also read the data from the SNP sequence files.

build a single dataframe with all of the SNP files in it
add new columns to the SNP dataframe
specifically : bp_on_contig,contig,alignment_type


- where bp_on_contig is determined using:
	- using the functions in cigar_string.py 
	+ the alignment start location
	+ the reverse flag column (0 or 256 ==  simple addition, 16 or 272 == do math to reverse the bp position with the read)

- where the alignment type is specified as either:
	'complete' or 'partial' based on: 
	- the .sam file's cigar string M+len(sequence) == complete



"""



