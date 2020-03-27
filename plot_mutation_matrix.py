"""Plotting"""
__author__ = "Vivekananda Sarangi"
__email__ = "sarangi.vivekananda@mayo.edu"
__status__ = "Development"

import os
import argparse
import sys
import Util
import pandas as pd
import pickle
import seaborn as sns
import matplotlib.pyplot as plt

def argument_parse():
	'''Parses the command line arguments'''
	parser=argparse.ArgumentParser(description='Preprocessing of vcf file')
	parser.add_argument("-p","--ALL2_output_path",help="Path to output folder of ALL2",required="True",
										type=Util.DirectoryValidator)
	parser.add_argument("-o","--output_dir",help="Path to directory where results will be written",required="True")
	parser.add_argument("-m","--mutation",help="Mutation format should be chr_pos_alt_ref(eg. chr1_12553232_A_C)",
											required="True",action='append')
	return parser

def main():
	# Extracting passed arguments
	parser = argument_parse()
	arg = parser.parse_args()
	
	# Assigning values to variable
	ALL2_output = arg.ALL2_output_path
	output_dir = arg.output_dir
	mutation_list = arg.mutation

	Util.ensure_dir(output_dir)
		
	# loading pickle file
	print("Loading pickle file")
	mutation_matrix_file = os.path.join( ALL2_output,"mutation_matrix.pkl" )	
	mutation_matrix_file_fh=open( mutation_matrix_file, "rb" )
	mutation_matrix_dict = pickle.load( mutation_matrix_file_fh )
	
	# ploting mutation
	print("Plotting mutation matrix")
	for mutation in mutation_list:
		print(mutation)
		mutation_df=mutation_matrix_dict[mutation]
		plt.figure()
		sns.heatmap(mutation_df,cmap="Blues",cbar=False)
		plt.tight_layout()
		plt.title(mutation)
		plt.show()
		plt.savefig(os.path.join(output_dir,mutation+".png"))
	
	
		
if __name__ == "__main__":
	Util.OsEnviron(os.environ)
	main()
