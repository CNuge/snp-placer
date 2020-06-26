#!/usr/bin/env python
import unittest
from itertools import groupby

def cigar_cutter(cigar):
	""" split a cigar string into tuples of bp# and cigar identifier """
	list_of_data = [''.join(g) for _, g in groupby(cigar, str.isalpha)]
	cigar_tuples = list(zip([int(x) for x in list_of_data[::2]],list_of_data[1::2]))
	return cigar_tuples


def adjust_bp(bp_of_snp, cigar_dat, flag = 0):
	""" scan the cigar data, making front trims, insertions, and deletions
		This makes the changes relative to the REFERENCE GENOME's base pairs
		Thereby orienting the SNPS correctly
		i.e. a sequence with a SNP at 20, who's cigar is 10M5D10M would have the 
		SNP 5 positions furhter right then where it is indicated on the short read at 25
		Similarly, a 10M5I10M would have the SNP 5bp left of the indiction on the short read
		as five base pairs are skipped over and not used in the reference genome"""
	change_to_bp = 0
	bp_scan = 0

	#bug caught by Eric Rondeau
	#But for reads aligned in reverse complement, 
	#the bp_of_snp is provided with respect to the original read orientation. 
	#Therefore the cigar string needs to be read from right to left
	if flag == 16 or flag == 272:
		cigar_dat = cigar_dat[::-1]		

	for x, cigar_bit in enumerate(cigar_dat):
		if cigar_bit[1] == 'M':
			""" count the matches towards the scan, no change to location"""
			bp_scan += cigar_bit[0]		
		elif cigar_bit[1] == 'D':
			"""add one to location, no change to scan count"""
			change_to_bp += cigar_bit[0]
		elif cigar_bit[1] == 'I':
			"""minus one from location, move bp scan count up"""
			change_to_bp -= cigar_bit[0]
			bp_scan += cigar_bit[0]	
			#exception: if the inserted bp housed the SNP
			if bp_scan == bp_of_snp:
				return 'snp_outside_aligned_region'
		elif cigar_bit[1] == 'S' and (x != (len(cigar_dat) - 1)):
			"""find the soft clipping strings, subtract from bp location"""
			change_to_bp -= cigar_bit[0]
			bp_scan += cigar_bit[0]
		elif cigar_bit[1] == 'S':
			bp_scan += cigar_bit[0]

		""" if the scan has passed the snp, return the result
			as no more changes are needed"""
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


def cigar_string_change(bp_of_snp, cigar_string, flag = 0):
	""" take in the original string, and the snp location, adjust location based on
		cigar data, returns a new bp integer that can be used relative to the start
		of the sequence's alignment to place the bp of the snp
		NOTE: both the input and output string are NOT zero indexed """
	cigar_dat = cigar_cutter(cigar_string)
	#first, identify the snps with no indels or font trimming
	if flag != 16 and flag != 272 and cigar_dat[0][1] == 'M' and cigar_dat[0][0] > bp_of_snp:
		return bp_of_snp
	else:
		new_bp = adjust_bp(bp_of_snp, cigar_dat, flag)
		if fringe_snp_check(bp_of_snp, cigar_dat) == True:
			new_bp = 'snp_outside_aligned_region'
		return new_bp


def alignment_length(cigar_string):
	""" take a list of cigar data tuples count total length of alignment: 
		'M' is match to the reference so these are counted
		'D' is deletion from  the reference, so these are counted in the length
		
		'I' is more sequence on the short read not on the referece, thefore
		don't add to length of sequence covered on the reference
		
		'S' is trimmed so don't add to length of sequence covered on the reference """
	align_length = 0
	for pair in cigar_cutter(cigar_string):
		if pair[1] == 'M' or pair[1] == 'D':
			align_length += pair[0]

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
			71)
		self.assertEqual(
			adjust_bp(70 ,[(52, 'M'), (5, 'I'), (33, 'M')]),
			65)
		self.assertEqual(
			adjust_bp(33, [(45,'M'),(23,'S')]),
			33)
		self.assertEqual(
			adjust_bp(53, [(52, 'M'), (1, 'I'), (33, 'M')]),
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

	def test_alignment_length(self):
		self.assertEqual(
			alignment_length('74M11S'),
			74)
		self.assertEqual(
			alignment_length('17M9D68M'),
			94)
		self.assertEqual(
			alignment_length('21M6I58M'),
			79)
		self.assertEqual(
			alignment_length('85M'),
			85)
		self.assertEqual(
			alignment_length('31M8D2I46M8S'),
			85)

	def test_cigar_string_change(self):
		self.assertEqual(
			cigar_string_change(70 ,'52M1D33M'),
			71)
		self.assertEqual(
			cigar_string_change(15 ,'52M1D33M'),
			15)
		self.assertEqual(
			cigar_string_change(70 ,'52M5I33M'),
			65)
		self.assertEqual(
			cigar_string_change(33, '45M23S'),
			33)
		self.assertEqual(
			cigar_string_change(53, '52M1I33M'),
			'snp_outside_aligned_region')
		self.assertEqual(
			cigar_string_change(46, '45M23S'),
			'snp_outside_aligned_region')		


if __name__ == "__main__":

	unittest.main()





