#!/usr/bin/env python3

import os
import unittest
import filter_sam_file
import place_snps
import pandas as pd

class SamFilterTests(unittest.TestCase):
	""" Test the functionality of the filter sam file module
		after the function is run, the outputs are compared 
		to the expected files in the example_data folder """
	@classmethod
	def setUpClass(self):
		""" read in the data necessary for the tests """
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


class PlaceSnpsTests(unittest.TestCase):
	"""test all of the functions for placing snps in the genome"""
	@classmethod
	def setUpClass(self):
		"""load in the data using the read input functions 
			list of inputs is to simulate the argument parser output.
			This test will fail if the files do not load, later tests
			will catch errors in the formatting of the files """
		self._sam_data = place_snps.read_sam_files(['example_data/numeric_ex.sam',
													'example_data/string_name_ex.sam'])
		self._snp_data = place_snps.read_input_files(['example_data/numeric_ex_snps.txt',
													'example_data/string_name_ex_snps.txt'])
			
	@classmethod
	def tearDown(self):
		"""once the unittest is run, remove the temporary test outputs"""
		try:
			os.remove('./example_data/temp.vcf')
		except OSError:
			pass

	def test_pipline(self):
		"""run the place snps pipeline and record the result in a temp .vcf
			the output is then evaluated against the expected data in the
			example_output.vcf"""
		self._filtered_sam = place_snps.sam_subset(self._snp_data['SNP'], 
													self._sam_data)

		self._all_data = place_snps.sam_polymorphism_column_merger(self._filtered_sam, 
																		self._snp_data)
		
		self._all_data = place_snps.calculate_new_bp_data(self._all_data)

		self._all_data = place_snps.snp_placement_dataframe(self._all_data)

		self._polymorphism_vcf = place_snps.output_to_vcf(self._all_data)

		self._polymorphism_vcf.to_csv('./example_data/temp.vcf', sep='\t', index=False)

		#then read in the example and test.vcfs as strings, compare the strings.
		template_vcf = open('example_data/example_output.vcf','r')
		self._template_vcf = template_vcf.read()
		template_vcf.close()

		test_gen_vcf = open('example_data/temp.vcf' , 'r')
		self._primary_test = test_gen_vcf.read()
		test_gen_vcf.close()

		self.assertEqual(self._template_vcf, self._primary_test)

if __name__ == '__main__':

	filter_sam_file.filter_sam('example_data/unfiltered_sam_data.sam', 
								'example_data/unittest_out1.sam', 
								'example_data/unittest_out2.sam')

	unittest.main()






