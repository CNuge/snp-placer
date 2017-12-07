import argparse
import cigarParse
import samParse
import pandas as pd
from pandas import Series, DataFrame
from itertools import groupby
import gc

def sam_subset(snp_names, sam_file_df):
	""" build a dictonary of relevant_sam_information
		for each snp, grab the following from the sam file:
		('Qname','Flag','Rname','Pos','MapQ','Cigar')
		if there are two+ hit locations, add an additional
		tuple to the dict value """
	sam_rows = sam_file_df.loc[sam_file_df['Qname'].isin(snp_names)]	
	return sam_rows


def sam_polymorphism_column_merger(sam_dataframe, snp_dataframe):
	"""grab the bp, and polymorphisms from the snp dataframe
	note that if two snps on one contig, there will be multiple rows for that snp"""
	return pd.merge(sam_dataframe, snp_dataframe, 
					how='left',left_on='Qname', right_on='SNP_name')


def calculate_new_bp_data(sam_dataframe):
	"""pass in sam df, this will call the relevant cigar functions
		determining the length of the alignment, and the adjusted
		position of the bp based on cigar alignment 
		apply cigarParse bp adjustment to each row"""
	sam_dataframe['adjusted_bp_SNP_location'] = sam_dataframe.apply(
		lambda x: cigarParse.cigar_string_change(x['Sequence'], 
													x['bp_SNP_location'], 
													x['Cigar']) , axis=1)
	""" count the alignment length for each row using cigar """
	sam_dataframe['alignment_length'] = sam_dataframe.apply(
		lambda x: cigarParse.alignment_length(x['Cigar']) , axis=1)
	return sam_dataframe


def snp_placement_dataframe(sam_dataframe):
	"""apply snp_contig_location across a dataframe """
	sam_dataframe['contig_location'] = sam_dataframe.apply(
		lambda x: samParse.snp_contig_location(x['Flag'], 
												x['Pos'], 
												x['adjusted_bp_SNP_location'], 
												x['alignment_length']), axis=1)
	return sam_dataframe


def output_to_vcf(output_df):
	""" need the following: #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""

#test this
	if type(output_df['SNP_name'][0]) == str:
		output_df['adj_name'] = output_df['SNP_name'] +'_' + output_df['Polymorphism'] + '_' + output_df['bp_SNP_location'].astype(str)		
	else:
		output_df['adj_name'] = output_df['SNP_name'].astype(str) +'_' + output_df['Polymorphism'] + '_' + output_df['bp_SNP_location'].astype(str)

	output_df['full_adj_name'] = output_df.apply(lambda x: samParse.compliment_name(x['adj_name'], x['Flag']), axis=1)

	vcf_out = output_df[['Rname','contig_location','full_adj_name', 'Polymorphism','MapQ','Flag']]
	vcf_out['FILTER'] = 'PASS'
	vcf_out['INFO'] = '.'
	vcf_out['REF_check'] = vcf_out['Polymorphism'].apply(lambda x: x.split('/')[0])
	vcf_out['ALT_a'] = vcf_out['Polymorphism'].apply(lambda x: x.split('/')[1:])
	vcf_out['ALT_check'] = vcf_out['ALT_a'].apply(lambda x: ','.join(x))

	vcf_out['REF'] = vcf_out.apply(lambda x: samParse.allele_comp_check(x['REF_check'] , x['Flag']), axis=1)
	vcf_out['ALT'] = vcf_out.apply(lambda x: samParse.allele_comp_check(x['ALT_check'] , x['Flag']), axis=1)

	vcf_out = vcf_out[['Rname','contig_location','full_adj_name','REF','ALT','MapQ','FILTER','INFO']]
	vcf_out.columns =['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
	return vcf_out


if __name__ == '__main__':


	parser = argparse.ArgumentParser()
	parser.add_argument('-s','--samfile', nargs = '+', required = True,
							help = 'The sam file(s) you wish to process. Pass in multiple files behind one flag. \
							i.e.  -s ex1.sam ex2.sam')
	parser.add_argument('-p', '--snpfile', nargs = '+', required = True,
							help = 'The snp file(s) you wish to process. Pass in multiple files behind one flag. \
							i.e.  -p ex1.txt ex2.txt')
	parser.add_argument('-o', '--output', default = 'placed_snps.vcf'
							help = 'The name of the output .vcf file. Default is placed_snps.vcf')
	args = parser.parse_args()


#clean this and come up with a way to do it dynamically
	#read in sam file
	sam_header = ['Qname','Flag','Rname','Pos','MapQ','Cigar','Rnext','Pnext', 'TLEN', 'SEQ', 'QUAL','tag','type','value']
	# if you have more columns, change this!
	#sam_header = ['Qname','Flag','Rname','Pos','MapQ','Cigar','Rnext','Pnext', 'TLEN', 'SEQ', 'QUAL','tag','type','value','bonus']


	if len(args.samfile) == 1:
		sam_dat = pd.read_table(args.samfile[0], sep='\t', names = sam_header, index_col=None)
	else:
		samfile_inputs = []
		for i in args.samfile:
			samfile_inputs.append(pd.read_table(i, sep='\t', index_col=None))
		sam_dat = pd.concat(samfile_inputs)
	
	#take the brackets out of the query section
	sam_dat['Qname'] = [x.split('(')[0] for x in sam_dat['Qname']]


	if len(args.snpfile) == 1:
		snp_input_dat= pd.read_table(args.snpfile, sep='\t', index_col=None)
	else:
		snpfile_inputs = []
		for i in args.snpfile:
			snpfile_inputs.append(pd.read_table(i, sep='\t', index_col=None))
		snp_input_dat = pd.concat(snpfile_inputs)

	samfile_inputs = [] #for garbage collection
	snpfile_inputs = [] #for garbage collection
	gc.collect()

"""
	#current tests:
	#args = parser.parse_args('-s ./example_data/numeric_ex.sam ./example_data/string_name_ex.sam -p ./example_data/numeric_ex.txt ./example_data/string_name_ex.txt'.split())
	#args = parser.parse_args('-s ./example_data/numeric_ex.sam  -p ./example_data/numeric_ex.txt ./example_data/string_name_ex.txt'.split())

	#example data for testing
	#read in SNP files
	sam_input_file1 = 'string_name_ex.sam'
	sam_input_file2 = 'numeric_ex.sam'

	snp_input_file1 = './example_data/numeric_ex.txt'
	snp_input_file2 = './example_data/string_name_ex.txt'

	snp_input_dat= pd.read_table(snp_input_file1, sep='\t', index_col=None)

	polymorphism_vcf.to_csv('example_data.vcf', sep='\t',index=False)

	Potential sources of error:
	There are a few ways this script can fail that I've found, here I point them out and tell you the fix.

	1. If your .sam file has extra columns on the right, you need to add these in to the sam_header list on line 75 (or delete the colums)

	2. If your snp names are numeric and not strings (i.e. 76 not CAM_SNP_76) then in the script 'place_snps.py' 
	move the # from line 60 to line 61 to make the necessary type change.

	3. add an argument parser to this 
"""

	#subset the sam alignments for rows matching the snp input
	snp_data_on_contigs = sam_subset(snp_input_dat['SNP_name'], sam_dat)

	all_polymorphism_data = sam_polymorphism_column_merger(snp_data_on_contigs, snp_input_dat)

	#calculate new bp data
	all_polymorphism_data = calculate_new_bp_data(all_polymorphism_data)

	all_polymorphism_data  = snp_placement_dataframe(all_polymorphism_data )


	polymorphism_vcf = output_to_vcf(all_polymorphism_data)

	polymorphism_vcf.to_csv(args.output, sep='\t',index=False)

