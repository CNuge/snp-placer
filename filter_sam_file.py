
""" this program will take a .sam input file (SAMTOOLs format)
	and separate the data into 3 categories:
	1. no alignments (discarded)
	2. alignments to one location 
	3. secondary aligments (2+ locations)

	For 2. and 3. new .sam files are produced with all 
	information retained from the alignment section"""


input_sam_file = './example_data/unfiltered_sam_data.sam'

secondary_alignments_out = './example_data/example_secondary_alignments.sam'
alignments_one_location = './example_data/samfile_one_location_alignments.sam'


with open(input_sam_file) as file:
	for line in file:
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
			
			""" there are other flags in this column, that do no appear in the dataset I'm working on"""
			print('this line not in the categories this program accounts for:\n')
			print(line)