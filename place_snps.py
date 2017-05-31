



""" example cigar data cases """

SAC24863(Q30)	0	Contig1825	947267	30	85M	*	0	0	TGCATTCCGATTCTTTCATAGATTACAATGGACACATATTATGTGTGTAAAATATTTTTTTGAAGGTTGTACTGATTATGATGAG	*	RG:Z:group1	MD:Z:58G26	NM:i:1
SAC24863	T/G	59	TGCATTCCGATTCTTTCATAGATTACAATGGACACATATTATGTGTGTAAAATATTTTTTTGAAGGTTGTACTGATTATGATGAG


SAC32530(Q42)	16	Contig785	222749	42	85M	*	0	0	TCCAAATATTCTAAGTCAGAACCGTCCAGAGTAGTGATGCTGGGCGGGCGGGCAGGTGTTGGTAGGATCGGTTAAAAAGCATGCA	*	RG:Z:group1	MD:Z:85	NM:i:0
SAC32530	C/T	42	TGCATGCTTTTTAACCGATCCTACCAACACCTGCCCGCCCGCCCAGCATCACTACTCTGGACGGTTCTGACTTAGAATATTTGGA


SAC55327(Q30)	0	Contig1199	1693992	30	84M1S	*	0	0	TGCATGGAGCGAGACAGGAAACAGGCAGAACAAACTATGGCCTCTTCTCCCCAAAGCCACTTCATCCATCCAGTCTGTTATGTTG	*	RG:Z:group1	MD:Z:84	NM:i:0
SAC55327	G/T	85	TGCATGGAGCGAGACAGGAAACAGGCAGAACAAACTATGGCCTCTTCTCCCCAAAGCCACTTCATCCATCCAGTCTGTTATGTTG

SAC200463(Q30)	16	Contig1290	5387428	30	84M1S	*	0	0	TCACACTCTCCGTAACTGGTTGCATCACTGCCTGGAATGGCAATTGCTCGGCCTCCGACCGCAAGGCACTACAGAGGGTAATGCA	*	RG:Z:group1	MD:Z:4G79	NM:i:1
SAC200463	T/C	81	TGCATTACCCTCTGTAGTGCCTTGCGGTCGGAGGCCGAGCAATTGCCATTCCAGGCAGTGATGCAACCAGTTACGGAGAGTGTGA


SAC162239(Q30)	0	Contig967	159258	30	74M11S	*	0	0	TGCATTACATAAAGACTTACAGTATACAGAAATCCAAAAGACCGCAAAACCCTGCTACTTGCAGATAGTATGCAAGATCGGAAGA	*	RG:Z:group1	MD:Z:74	NM:i:0
SAC162239	A/C	63	TGCATTACATAAAGACTTACAGTATACAGAAATCCAAAAGACCGCAAAACCCTGCTACTTGCAGATAGTATGCAAGATCGGAAGA


SAC45959(Q0)	0	Contig0	422443	0	52M1D33M	*	0	0	TGCATATTGTACGATCTAGAATAGAGTACGAGGCAGTTTAATTTGGGCACGAGTTTTCCAAAGTGGAAACAGCGCCTCCATATTG	*	RG:Z:group1	MD:Z:52^T0T23C8	NM:i:3
SAC45959	T/C	77	TGCATATTGTACGATCTAGAATAGAGTACGAGGCAGTTTAATTTGGGCACGAGTTTTCCAAAGTGGAAACAGCGCCTCCATATTG



sequence_string='TGCATATTGTACGATCTAGAATAGAGTACGAGGCAGTTTAATTTGGGCACGAGTTTTCCAAAGTGGAAACAGCGCCTCCATATTG'
bp_of_snp=77
cigar_string='52M1D33M'
cigar='52M1D33M'


###########################################
# code starts here
###########################################

from itertools import groupby


def cigar_cutter(cigar):
	"""split a cigar string into tuples of bp# and cigar identifier """
	list_of_data = [''.join(g) for _, g in groupby(cigar, str.isalpha)]
	cigar_tuples = list(zip([int(x) for x in list_of_data[::2]],list_of_data[1::2]))
	return cigar_tuples

def adjust_bp(bp_of_snp, cigar_dat):
	""" scan the cigar dat, making front trims, insertions, and deletions"""
	change_to_bp = 0
	bp_scan = 0
	for cigar_bit in cigar_dat:
		if cigar_bit[1] == 'S':
			"""find the soft clipping strings, subtract from bp location"""
			change_to_bp -= cigar_bit[0]
			bp_scan += cigar_bit[0]
		elif cigar_bit[1] == 'M':
			""" count the matches towards the scan, no change to location"""
			bp_scan += cigar_bit[0]
		
		elif cigar_bit[1] == 'D':
			"""minus one from location, move bp scan count up"""
			change_to_bp -= cigar_bit[0]
			bp_scan += cigar_bit[0]
		
		elif cigar_bit[1] == 'I':
			"""add one to location, no change to scan count"""
			change_to_bp += cigar_bit[0]

		""" if the scan has passed the snp, we can return the result,"""
		"""	as no more changes will happen"""
		if bp_scan >= cigar_bit[0]:
				return (bp_of_snp + change_to_bp)

def fringe_snp_check(bp_of_snp, cigar_dat, sequence_string):
	""" look at the first and last tuples, if snp falls outside aligned region, return True"""
	if cigar_dat[0][1] == 'S':
		if cigar_dat[0][0] > bp_of_snp:
			return True
	elif cigar_dat[-1][1] == "S":
		if (len(sequence_string) - cigar_dat[-1][0]) < bp_of_snp:
			return True
	return False


def cigar_string_change(sequence_string, bp_of_snp, cigar_string):
	""" take in the original string, and the snp location, adjust location based on
		cigar data, returns a new bp integer that can be used relative to the start
		of the sequence's alignment to place the bp of the snp
		NOTE: both the input and output string are NOT zero indexed"""
	cigar_dat = cigar_cutter(cigar_string)
	#first, identify the snps with no indels or font trimming
	if cigar_dat[0][1] == 'M' and cigar_dat[0][0] > bp_of_snp:
		return bp_of_snp
	else:
		new_bp = adjust_bp(bp_of_snp, cigar_dat)
		if fringe_snp_check(bp_of_snp, cigar_dat, sequence_string) == True:
			new_bp = 'snp_outside_aligned_region'
		return new_bp














