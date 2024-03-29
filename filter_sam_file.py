#!/usr/bin/env python3
import argparse

""" this program will take a .sam input file (SAMTOOLs format)
	and separate the data into 3 categories:
	1. no alignments (discarded)
	2. alignments to one location 
	3. secondary aligments (2+ locations)

	For 2. and 3. new .sam files are produced with all 
	information retained from the alignment section

	default output names are:
	samfile_one_location_alignments.sam
	secondary_alignments.sam

	These names can be changed using the input flags
	-p for primary (one location) outputs and 
	-s for secondary (multiple location) outputs

	run from command line using:
	python filter_sam_file.py input_sam_name.sam
	"""

def filter_sam(input_sam_file, alignments_one_location, secondary_alignments_out, unaccounted_alignments_out):
	""" take an input .sam file, remove the header lines and split the results
		into alignments to one location and alignments to two+ locations """
	with open(input_sam_file) as file:
		for i, line in enumerate(file):
			if line[0] == "@":
				""" skip the header section"""
				continue
			line_split = line.split('\t')
			if line_split[1] == '4':
				""" a 4 in column 2 means the read was not aligned"""
				continue
			elif line_split[1] == '256' or line_split[1] == '272':
				""" a 256 or 272 in column 2 means the read aligned to 2+ locations """
				add_to = open(secondary_alignments_out,'a')
				add_to.write(line)
				add_to.close()
			elif line_split[1] == '0' or line_split[1] == '16':
				"""a 0 or 16 in column 2 means the read aligned to one location"""
				add_to = open(alignments_one_location,'a')
				add_to.write(line)
				add_to.close()
			else:
				#account for these!				
				""" there are other flags in this column, that do no appear in the dataset I'm working on"""
				print(f'line {i+1} not in the categories this program accounts for.\n')
				add_to = open(unaccounted_alignments_out,'a')
				add_to.write(line)
				add_to.close()

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('input', type = str, 
		help = 'The name of the .sam file you wish to filter.')
	parser.add_argument('-p', '--primary', type = str, default = 'one_location_alignments.sam',
		help = 'Optional: a name for the output .sam file containing the alignments to one location.\n\
		Default is: one_location_alignments.sam')
	parser.add_argument('-s', '--secondary', type = str, default = 'multiple_location_alignments.sam',
		help = 'Optional: a name for the output .sam file containing the alignments to two or more locations.\n\
		Default is: multiple_location_alignments.sam')
	parser.add_argument('-u', '--unaccounted', type = str, default = 'unaccounted_alignments.sam',
		help = 'Optional: a name for the output .sam file containing the alignments falling in neither of the previous categories.\n\
		Default is: unaccounted_alignments.sam')
	
	args = parser.parse_args()

	filter_sam(args.input, args.primary, args.secondary, args.unaccounted)

