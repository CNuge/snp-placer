#!/usr/bin/env python3
import gc
import argparse
import pandas as pd
from itertools import groupby
from pandas import Series, DataFrame

from parse import cigarParse
from parse import samParse


def read_input_files(list_of_inputs, names = None):
	""" read in a list of files (from argparse), turn them into dataframes
		and concatenate them if the length of the list exceeds 1."""
	if len(list_of_inputs) == 1:
		data = pd.read_table(list_of_inputs[0], sep='\t', names = names, index_col=None)

	else:
		data_files = []
		for i in list_of_inputs:
			data_files.append(pd.read_table(i, sep='\t',names = names, index_col=None))
		data = pd.concat(data_files)
	return data


def read_sam_files(list_of_inputs):
	""" this calls the read_input_files and stores the sam_header """
	sam_header = ['Qname','Flag','Rname','Pos','MapQ','Cigar','Rnext',
					'Pnext', 'TLEN', 'SEQ', 'QUAL','tag','type','value','value2']
	return read_input_files(list_of_inputs, sam_header)


def sam_subset(snp_names, sam_file_df):
	""" take the subset of rows from the sam file that contain snps found 
		in the snp list. Equivalent to an inner join"""
	sam_rows = sam_file_df.loc[sam_file_df['Qname'].isin(snp_names)]	
	return sam_rows


def sam_polymorphism_column_merger(sam_dataframe, snp_dataframe):
	""" grab the bp, and polymorphisms from the snp dataframe
		note that if two snps on one contig, there will be multiple 
		rows for that snp"""
	return pd.merge(sam_dataframe, snp_dataframe, 
					how='left', left_on='Qname', right_on='SNP')


def calculate_new_bp_data(sam_dataframe):
	""" pass in sam df, this will call the relevant cigar functions
		determining the length of the alignment, and the adjusted
		position of the bp based on cigar alignment 
		apply cigarParse bp adjustment to each row"""
	sam_dataframe['adjusted_bp_SNP_location'] = sam_dataframe.apply(
		lambda x: cigarParse.cigar_string_change(x['bp'], 
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
	""" take the informaiton in the dataframe and turn it
		into .vcf file format with the following header: 
		#CHROM POS ID REF ALT QUAL FILTER INFO"""
	if type(output_df['SNP'][0]) == str:
		output_df['adj_name'] = output_df['SNP'] + \
								'_' + \
								output_df['Polymorphism'] + \
								'_' + \
								output_df['bp'].astype(str)		
	else:
		output_df['adj_name'] = output_df['SNP'].astype(str) + \
								'_' + \
								output_df['Polymorphism'] + \
								'_' + \
								output_df['bp'].astype(str)

	output_df['full_adj_name'] = output_df.apply(
		lambda x: samParse.compliment_name(x['adj_name'], 
											x['Flag']), axis=1)

	vcf_out = output_df[['Rname',
						'contig_location',
						'full_adj_name', 
						'Polymorphism',
						'MapQ',
						'Flag']].copy()
	
	vcf_out['FILTER'] = '.'
	vcf_out['INFO'] = '.'
	vcf_out['QUAL'] = '.'
	
	vcf_out['REF_check'] = vcf_out['Polymorphism'].apply(
		lambda x: x.split('/')[0])

	vcf_out['ALT_a'] = vcf_out['Polymorphism'].apply(
		lambda x: x.split('/')[1:])
	
	vcf_out['ALT_check'] = vcf_out['ALT_a'].apply(
		lambda x: ','.join(x))

	vcf_out['REF'] = vcf_out.apply(
		lambda x: samParse.allele_comp_check(x['REF_check'],
											x['Flag']), axis=1)
	
	vcf_out['ALT'] = vcf_out.apply(
		lambda x: samParse.allele_comp_check(x['ALT_check'],
											x['Flag']), axis=1)

	vcf_out = vcf_out[['Rname','contig_location','full_adj_name','REF','ALT','MapQ','FILTER','INFO']]
	vcf_out.columns =['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
	return vcf_out


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-s','--samfile', nargs = '+', required = True,
							help = 'The sam file(s) you wish to process. Pass in multiple files behind one flag. \
							i.e.  -s ex1.sam ex2.sam')
	parser.add_argument('-p', '--snpfile', nargs = '+', required = True,
							help = 'The snp file(s) you wish to process. Pass in multiple files behind one flag. \
							i.e.  -p ex1.txt ex2.txt')
	parser.add_argument('-o', '--output', default = 'placed_snps.vcf',
							help = 'The name of the output .vcf file. Default is placed_snps.vcf')
	args = parser.parse_args()

	#note these must receive lists as first argument
	sam_dat =  read_sam_files(args.samfile)
	
	snp_input_dat = read_input_files(args.snpfile)

	#subset the sam alignments for rows matching the snp input
	sam_data_on_contigs = sam_subset(snp_input_dat['SNP'], sam_dat)

	#merge the snp and sam dataframes
	all_polymorphism_data = sam_polymorphism_column_merger(sam_data_on_contigs, snp_input_dat)

	#calculate new bp data
	all_polymorphism_data = calculate_new_bp_data(all_polymorphism_data)

	all_polymorphism_data  = snp_placement_dataframe(all_polymorphism_data )

	#default .vcf has name 'placed_snps.vcf'
	polymorphism_vcf = output_to_vcf(all_polymorphism_data)

	polymorphism_vcf.to_csv(args.output, sep='\t',index=False)

