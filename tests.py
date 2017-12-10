#!/usr/bin/env python3

import os
import unittest
import filter_sam_file
import place_snps

class SamFilterTests(unittest.TestCase):
	""" Test the functionality of the filter sam file module """
	@classmethod
	def setUpClass(self):
		""" read in the data nexessary for the tests """
		primary_example_file = open('example_data/samfile_one_location_alignments.sam','r')
		self._primary_example = primary_example_file.read()
		primary_example_file.close()

		primary_test_file = open('example_data/unittest_out1.sam' , 'r')
		self._primary_test = primary_test_file.read()
		primary_test_file.close()

		secondary_example_file = open('example_data/secondary_alignments.sam','r')
		self._secondary_example = secondary_example_file.read()
		secondary_example_file.close()

		secondary_test_file = open('example_data/unittest_out2.sam' , 'r')
		self._secondary_test = secondary_test_file.read()
		secondary_test_file.close()

	@classmethod
	def tearDown(self):
		""" once the unittest is run, remove the temporary test outputs"""
		try:
			os.remove('./example_data/unittest_out1.sam')
			os.remove('./example_data/unittest_out2.sam')
		except OSError:
			pass
	def test_primary_output(self):
		self.assertEqual(self._primary_example, self._primary_test)
	def test_secondary_output(self):
		self.assertEqual(self._secondary_example, self._secondary_test)

"""
class PlaceSnpsTests(unittest.TestCase):
	#test all of the functions for placing snps in the genome 
	@classmethod
	def setUpClass(self):


	@classmethod
	def tearDown(self):

TODO
	def test_read_input_files
	def test_sam_subset
	def test_sam_polymorphism_column_merger
	def test_calculate_new_bp_data
	def test_snp_placement_dataframe
	def test_output_to_vcf

also test the argument parsers?

		#example data for testing
		#read in SNP files
		sam_input_file1 = './example_data/string_name_ex.sam'
		sam_input_file2 = './example_data/numeric_ex.sam'

		snp_input_file1 = './example_data/numeric_ex.txt'
		snp_input_file2 = './example_data/string_name_ex.txt'

		polymorphism_vcf.to_csv('example_data.vcf', sep='\t',index=False)

"""

if __name__ == '__main__':

	filter_sam_file.filter_sam('example_data/unfiltered_sam_data.sam', 
								'example_data/unittest_out1.sam', 
								'example_data/unittest_out2.sam')

	unittest.main()






