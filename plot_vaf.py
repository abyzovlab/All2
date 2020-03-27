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
	return parser

def main():
	# Extracting passed arguments
	parser = argument_parse()
	arg = parser.parse_args()
	
	# Assigning values to variable
	ALL2_output = arg.ALL2_output_path
	output_dir = arg.output_dir

	Util.ensure_dir(output_dir)
	
	es_file=os.path.join(ALL2_output,"explaination_score.txt")
	vaf_file=os.path.join(ALL2_output,"variant_vaf_per_pair.txt")
	vaf_ef_file=os.path.join(ALL2_output,"variant_vaf_per_pair_es.txt")
	vaf_ef_file_fh=open(vaf_ef_file,'w')
	es_dict={}
	vaf_dict={}
	
	for i in open(es_file):
		line=i.strip().split("\t")
		if i.startswith("#"):
			continue
		mutation="_".join(line[:4])
		germline_score=line[5]
		mosaic_score=line[4]
		germline_samples=line[7]
		mosaic_samples=line[6]
		es_dict[mutation]=[(germline_score,mosaic_score),(germline_samples,mosaic_samples)]
	for i in open(vaf_file):
		line=i.strip().split("\t")
		if i.startswith("Case"):
			vaf_ef_file_fh.write(i.strip()+"\tMosaic_score\tGermline_score\tMosaic_samples\tGermline_Samples\n")
			continue	
		case=line[0]
		control=line[1]
		mutation="_".join(line[2:6])
		mosaic_score=es_dict[mutation][0][1]
		germline_score=	es_dict[mutation][0][0]
		mosaic_samples=es_dict[mutation][1][1]
		germline_samples=es_dict[mutation][1][0]
		outline="\t".join([i.strip(),mosaic_score,germline_score,mosaic_samples, germline_samples])
		vaf_ef_file_fh.write(outline+"\n")
	vaf_ef_file_fh.close()
		
		
	
	
		
if __name__ == "__main__":
	Util.OsEnviron(os.environ)
	main()
