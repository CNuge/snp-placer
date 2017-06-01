import cigarParse
import pandas as pd
from pandas import Series, DataFrame
from itertools import groupby


def sam_subset(snp_names, sam_file_df):
	""" build a dictonary of relevant_sam_information
		for each snp, grab the following from the sam file:
		('Qname','Flag','Rname','Pos','MapQ','Cigar')
		if there are two+ hit locations, add an additional
		tuple to the dict value """
	sam_rows = sam_file_df.loc[sam_file_df['Qname'].isin(snp_names)]	
	return sam_rows


def sam_polymorphism_column_merger(sam_dataframe, snp_dataframe):
	"""grab the bp, and polymorphisms from the snp dataframe"""
	""" note that if two snps on one contig, there will be multiple rows for that snp"""
	return pd.merge(sam_dataframe, snp_dataframe, how='left',left_on='Qname', right_on='SNP_name')


def calculate_new_bp_data(sam_dataframe):
	"""pass in sam df, this will call the relevant sigar functions
		determining the length of the alignment, and the adjusted
		position of the bp based on cigar alignment """
	""" apply cigarParse bp adjustment to each row"""
	sam_dataframe['adjusted_bp_SNP_location'] = sam_dataframe.apply(
		lambda x: cigarParse.cigar_string_change(x['Sequence'], x['bp_SNP_location'], x['Cigar']),axis=1)
	""" count the alignment length for each row using cigar """
	sam_dataframe['alignment_length'] = sam_dataframe.apply(
		lambda x: cigarParse.alignment_length(x['Cigar']), axis=1)
	return sam_dataframe


def snp_contig_location(flag, pos, adjusted_bp_location, alignment_length):
	""" determine new bp position of the snp on the larger contig"""
	try:
		adjusted_bp_location / 1
	except:
		return '-' 
	if flag == 0 or flag == 256:
		""" forward aligment, add adj_bp to pos"""
		return (pos + adjusted_bp_location - 1 )
	elif flag == 16 or flag == 272:
		return (pos + alignment_length - adjusted_bp_location)
	else:
		return '-'


def snp_placement_dataframe(sam_dataframe):
	"""apply snp_contig_location across a dataframe """
	sam_dataframe['contig_location'] = sam_dataframe.apply(
		lambda x: snp_contig_location(x['Flag'], x['Pos'], x['adjusted_bp_SNP_location'], x['alignment_length']), axis=1)
	return sam_dataframe


if __name__ == '__main__':



	#read in sam file
	sam_header = ['Qname','Flag','Rname','Pos','MapQ','Cigar','Rnext','Pnext', 'TLEN', 'SEQ', 'QUAL','tag','type','value']

	sam_input_file = './sam_files/all_snps_samfile_one_location_alignments.sam'
	sam_dat = pd.read_table(sam_input_file, sep='\t', names = sam_header, index_col=None)

	#take the brackets out of the query section
	sam_dat['Qname'] = [x.split('(')[0] for x in sam_dat['Qname']]



	#read in SNP files

	snp_input_file1 = './snp_files/Stacks_SNP_info.txt'
	snp_dat1 = pd.read_table(snp_input_file1, sep='\t', index_col=None)

	snp_input_file2 = './snp_files/Tassel_PSV_info.txt'
	snp_dat2 = pd.read_table(snp_input_file2, sep='\t', index_col=None)

	snp_input_file3 = './snp_files/Tassel_SNP_info.txt'
	snp_dat3 = pd.read_table(snp_input_file3, sep='\t', index_col=None)

	frames=[snp_dat1,snp_dat2,snp_dat3]
	snp_input_dat = pd.concat(frames)

	#snp data now avaliable
	snp_input_dat.head()

	#subset the sam alignments for rows matching the snp input
	snp_data_on_contigs = sam_subset(snp_input_dat['SNP_name'], sam_dat)

	all_polymorphism_data = sam_polymorphism_column_merger(snp_data_on_contigs, snp_input_dat)

	#calculate new bp data
	all_polymorphism_data = calculate_new_bp_data(all_polymorphism_data)

	all_polymorphism_data  = snp_placement_dataframe(all_polymorphism_data )

	all_polymorphism_data.to_csv('snps_on_contigs_one_location.csv')










