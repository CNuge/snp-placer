from itertools import groupby
import unittest

def cigar_cutter(cigar):
	"""split a cigar string into tuples of bp# and cigar identifier """
	list_of_data = [''.join(g) for _, g in groupby(cigar, str.isalpha)]
	cigar_tuples = list(zip([int(x) for x in list_of_data[::2]],list_of_data[1::2]))
	return cigar_tuples

def adjust_bp(bp_of_snp, cigar_dat):
	""" scan the cigar data, making front trims, insertions, and deletions"""
	change_to_bp = 0
	bp_scan = 0
	for x, cigar_bit in enumerate(cigar_dat):
		if cigar_bit[1] == 'M':
			""" count the matches towards the scan, no change to location"""
			bp_scan += cigar_bit[0]		
		elif cigar_bit[1] == 'D':
			"""minus one from location, move bp scan count up"""
			change_to_bp -= cigar_bit[0]
			bp_scan += cigar_bit[0]
			#exception: if the deleted bp housed the SNP
			if bp_scan == bp_of_snp:
				return 'snp_outside_aligned_region'		
		elif cigar_bit[1] == 'I':
			"""add one to location, no change to scan count"""
			change_to_bp += cigar_bit[0]
		elif cigar_bit[1] == 'S' and (x != (len(cigar_dat) - 1)):
			"""find the soft clipping strings, subtract from bp location"""
			change_to_bp -= cigar_bit[0]
			bp_scan += cigar_bit[0]
		elif cigar_bit[1] == 'S':
			bp_scan += cigar_bit[0]

		""" if the scan has passed the snp, we can return the result,"""
		"""	as no more changes will happen"""
		if bp_scan >= bp_of_snp:
				return (bp_of_snp + change_to_bp)


def fringe_snp_check(bp_of_snp, cigar_dat):
	""" look at the first and last tuples, if snp falls outside aligned region, return True"""
	sequence_len = 0
	for x in cigar_dat:
		sequence_len += x[0]
	if cigar_dat[0][1] == 'S':
		if cigar_dat[0][0] > bp_of_snp:
			return True
	if cigar_dat[-1][1] == 'S':
		if (sequence_len - cigar_dat[-1][0]) < bp_of_snp:
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


def alignment_length(cigar_string):
	""" Take a list of cigar data tuples count total length of alignment
		Note this is the length of sequence traversed on the genome, 
		not the length on the short read!"""
	cigar_cutter_output = cigar_cutter(cigar_string)
	align_length = 0
	for pair in cigar_cutter_output:
		if pair[1] == 'M' or pair[1] == 'D':
			align_length += pair[0]
		elif pair[1] == 'S' or pair[1] == 'I':
			align_length -= pair[0]
	return align_length


class CigarTests(unittest.TestCase):
	def test_cigar_cutter(self):
		self.assertEqual(
			cigar_cutter('52M1D33M'), 
			[(52, 'M'), (1, 'D'), (33, 'M')])
		self.assertEqual(
			cigar_cutter('84M1S'), 
			[(84, 'M'), (1, 'S')])
		self.assertEqual(
			cigar_cutter('74M11S'), 
			[(74, 'M'), (11,'S')])

	def test_adjust_bp(self):
		self.assertEqual(
			adjust_bp(70 ,[(52, 'M'), (1, 'D'), (33, 'M')]),
			69)
		self.assertEqual(
			adjust_bp(70 ,[(52, 'M'), (5, 'I'), (33, 'M')]),
			75)
		self.assertEqual(
			adjust_bp(33, [(45,'M'),(23,'S')]),
			33)
		self.assertEqual(
			adjust_bp(53, [(52, 'M'), (1, 'D'), (33, 'M')]),
			'snp_outside_aligned_region')

	def test_fringe_snp_check(self):
		self.assertTrue(
			fringe_snp_check(75, [(74, 'M'), (11,'S')]))
		self.assertFalse(
			fringe_snp_check(74 , [(74, 'M'), (11,'S')]))
		self.assertTrue(
			fringe_snp_check(5, [(12, 'S'),(84, 'M')]))
		self.assertFalse(
			fringe_snp_check(5, [(1, 'S'),(84, 'M')]))

"""
TODO

'52M1D33M'
[(52, 'M'), (1, 'D'), (33, 'M')]

'84M1S'


'74M11S'
[(74, 'M'), (11,'S')]


70 [(52, 'M'), (1, 'D'), (33, 'M')]
69

[(74, 'M'), (11,'S')]

31M8D46M8S


74 , [(74, 'M'), (11,'S')]
False

75, [(74, 'M'), (11,'S')]
True

5, [(1, 'S'),(84, 'M')]
False

5, [(12, 'S'),(84, 'M')]
True

"""
if __name__ == "__main__":

	unittest.main()





