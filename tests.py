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


class PlaceSnpsTests(unittest.TestCase):
	"""test all of the functions for placing snps in the genome """
	def test_read_input_files(self):
		""" load in the data using the read input functions 
			list of inputs is to simulate the argument parser output.
			This test will fail if the files do not load, later tests
			will catch errors in the formatting of the files """
		self._sam_data = place_snps.read_sam_files(['example_data/string_name_ex.sam',
													'example_data/numeric_ex.sam'])
		self._snp_data = place_snps.read_input_files(['example_data/numeric_ex_snps.txt',
													'example_data/string_name_ex_snps.txt'])

	def test_sam_subset(self):
		""" Dataframe here is small example. Returns rows containing ids from list"""
		_ex_df = pd.DataFrame([['Name1',1,2],
								['Name9',3,4],
								['Name2',5,6]], 
								columns = ['Qname', 'Contig','Pos'])
		_test_ids = ['Name1', 'Name2', 'Name3']
		#what it will be compared against
		_comp_df = pd.DataFrame([['Name1',1,2],
								['Name2',5,6]], 
								columns = ['Qname', 'Contig','Pos'])
		_test_out = place_snps.sam_subset(_test_ids, _ex_df)
		self.assertEqual(list(_test_out['Qname']),
						list(_comp_df['Qname']))

"""

TODO
	def test_sam_polymorphism_column_merger
	def test_calculate_new_bp_data
	def test_snp_placement_dataframe
	def test_output_to_vcf

ALT one big test of pipeline and compare output to hand made .vcf:

class PlaceSnpsTests(unittest.TestCase):
	#test all of the functions for placing snps in the genome

	@classmethod
	def setUpClass(self):
		# load in the data using the read input functions 
		#	list of inputs is to simulate the argument parser output.
		#	This test will fail if the files do not load, later tests
		#	will catch errors in the formatting of the files 
		self._sam_data = place_snps.read_sam_files(['example_data/string_name_ex.sam',
													'example_data/numeric_ex.sam'])
		self._snp_data = place_snps.read_input_files(['example_data/numeric_ex_snps.txt',
			
	@classmethod
	def tearDown(self):
		#once the unittest is run, remove the temporary test outputs
		try:
			os.remove('./example_data/temp.vcf')
		except OSError:
			pass

	def test_pipline(self):
		#this could be moved out to the ifmain? or bring the other one into its unittest
		#either way consistency would be good.
		self._filtered_sam = place_snps.sam_subset(self._snp_data['SNPs'], 
													self._sam_data)

		self._all_data = place_snps.sam_polymorphism_column_merger(self._filtered_sam, 
																		self._sam_data)
		
		self._all_data = place_snps.calculate_new_bp_data(self._all_data)

		self._all_data = place_snps.snp_placement_dataframe(self._all_data)

		self._polymorphism_vcf = place_snps.output_to_vcf(self._all_data)

		self._polymorphism_vcf.to_csv('./example_data/temp.vcf', sep='\t', index=False)

		#then read in the example and test.vcfs as strings, compare the strings.

"""

if __name__ == '__main__':

	filter_sam_file.filter_sam('example_data/unfiltered_sam_data.sam', 
								'example_data/unittest_out1.sam', 
								'example_data/unittest_out2.sam')

	unittest.main()






