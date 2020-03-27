"""Takes a manifest file, and does all to all comparision to detect mosaic mutations"""
__author__ = "Vivekananda Sarangi"
__email__ = "sarangi.vivekananda@mayo.edu"
__status__ = "Development"

import os
import argparse
import sys
import Util
import gzip
import pandas as pd
import numpy as np 
from math import sqrt
import pickle
import seaborn as sns
import matplotlib.pyplot as plt

def argument_parse():
	'''Parses the command line arguments'''
	parser=argparse.ArgumentParser(description='Preprocessing of vcf file')
	parser.add_argument("-m","--manifest_file",help="Path to manifest file",required="True",type=Util.FileValidator)
	parser.add_argument("-o","--output_dir",help="Path to directory where results will be written",required="True")
	parser.add_argument("-a","--af_cutoff",help="Allele frequency cut-off for varaints",default="0.35")
	return parser

def extract_mutation_information(manifest_file, output_dir):
	output_file=os.path.join(output_dir,"variant_vaf_per_pair.txt")
	output_file_fh=open(output_file,'w')
	variant_dict={}
	pairs_vaf_dict={}
	head={}
	for i in open(manifest_file):
		line=i.strip().split("\t")
		if i.startswith("#"):
			for n,j in enumerate(line):
				head[j.replace("#","")] = n
			continue
		case=line[head["Case"]]
		control=line[head["Control"]]
		filename=line[head["Filename"]]
		pair=(case,control)
		if filename.endswith("vcf.gz"):
			filename_fh=gzip.open(filename)
		elif filename.endswith("vcf"):
			filename_fh=open(filename)
		else:
			print("Please make sure you provide a valid vcf file")
			exit()
		variant_head={}
	        # Reading variants from the vcf file
		for variant in filename_fh:
			line=variant.strip().split("\t")
			if variant.startswith("#"):
				if variant.startswith("#CHROM"):
        	                	for n,j in enumerate(line):
	                                	variant_head[j]=n
				continue
			chrm=line[variant_head["#CHROM"]]
			pos=line[variant_head["POS"]]
			ref=line[variant_head["REF"]]
			all_alt=line[variant_head["ALT"]]
			# Getting AD and DP field for case
			case_format=line[variant_head[case]-1].split(":")
			case_genotype=line[variant_head[case]].split(":")
			ad_index=case_format.index("AD")
			depth_index=case_format.index("DP")
			if depth_index>-1:
				 case_genotype_depth=case_genotype[depth_index]
			else:
				case_genotype_depth="Absent"
			if ad_index>-1:
				case_genotype_ad=case_genotype[ad_index]
			else:
				case_genotype_ad="Absent"
                        # making sure to take care of multiallelic locations
			for n,alt in enumerate(all_alt.split(",")):
				mutation="_".join([chrm,pos,ref,alt])
                                # storing variants
				if mutation in variant_dict:
					variant_dict[mutation].append(pair)
				else:
					variant_dict[mutation]=[pair]
				# storing VAFs
				if case_genotype_ad is not "Absent" and case_genotype_depth is not "Absent":
					alt_supporting_read=case_genotype_ad.split(",")[n+1]
					vaf=str(float(float(alt_supporting_read)/float(case_genotype_depth)))
				else:
					 vaf="Absent"
                                # creating a dictionary of dictionary ({pairs:{mutation:[vaf]}})
				if pair in pairs_vaf_dict.keys():
					if mutation in pairs_vaf_dict[pair].keys():
						pairs_vaf_dict[pair][mutation].append(vaf)
					else:
						pairs_vaf_dict[pair][mutation]=[vaf]
				else:
					pairs_vaf_dict[pair]={mutation:[vaf]}
	filename_fh.close()
	list_of_samples=[]
	output_file_fh.write("Case\tControl\tChrom\tPos\tRef\tAlt\tVAF\n")
	for pairs in pairs_vaf_dict:
		if pairs[0] not in list_of_samples:
			list_of_samples.append(pairs[0])
		case,control=pairs
		for mutation in pairs_vaf_dict[pairs]:
			for vaf in pairs_vaf_dict[pairs][mutation]:
				out_line="\t".join([case,control,"\t".join(mutation.split("_")),vaf])
				output_file_fh.write(out_line+"\n")
	output_file_fh.close()
	return variant_dict,pairs_vaf_dict,list_of_samples	


def explaination_score(variant_dict,list_of_samples,output_dir):
		
	output_file=os.path.join(output_dir,"explaination_score.txt")
	output_file_fh=open(output_file,'w')
	vaf_file=os.path.join(output_dir,"variant_vaf_per_pair.txt")
	vaf_ef_file=os.path.join(output_dir,"variant_vaf_per_pair_es.txt")
	vaf_ef_file_fh=open(vaf_ef_file,'w')
	es_dict={}
	#output_file_fh.write("#Chrom\tPos\tRef\tAlt\tMosaic_score\tGermline_score\tSample_with_Mosaic\tSamples_without_Germline\n")
	output_file_fh.write("#Chrom\tPos\tRef\tAlt\tMosaic_score\tGermline_score\tSamples_with_mutation\n")
	# number_of_cells_N is the total number of cells in the experiment
	total_number_of_cells_N = len(list_of_samples)

	# Preparing to store the data matrix for each mutation
	mutation_matrix_dict={}	
	mutation_matrix_file=os.path.join(output_dir,"mutation_matrix.pkl")
	mutation_matrix_file_fh=open(mutation_matrix_file,'wb')

	for mutation in variant_dict:
		# pairs_list is a list ofpairs the mutation was called in
		pairs_list = variant_dict[mutation]
		# pairs_list_n is the number of pairs the mutation was called in
		pairs_list_n = len(pairs_list)	
		# cell_fraction_f is the fraction of cells carrying the mutatin
		cell_fraction_f = 1/2 - sqrt(1/4 - pairs_list_n/total_number_of_cells_N**2)
		# is the number of cells carrying the mutation
		cells_carrying_mutation_Nv = round(cell_fraction_f*total_number_of_cells_N)
		#Creating an 'zero' data frame/matrix and updating mutation specific dataframe/matrix
		mutation_df = pd.DataFrame(np.zeros((total_number_of_cells_N,total_number_of_cells_N )),
						index=list_of_samples,columns=list_of_samples)
		for case,control in pairs_list:
			mutation_df.loc[case,control] = 1
		
		mutation_matrix_dict[mutation]=mutation_df	
	
		# calculating explaination score
		ordered_col_sum = mutation_df.sum(axis=1).sort_values(ascending=False)
		ordered_row_sum = mutation_df.sum(axis=0).sort_values(ascending=False)
		explained_call_n_mosaic = ordered_col_sum[:cells_carrying_mutation_Nv].sum()
		explained_call_n_germ = ordered_row_sum[:cells_carrying_mutation_Nv].sum()
		explanation_score_mosaic = explained_call_n_mosaic / pairs_list_n
		explanation_score_germ=explained_call_n_germ / pairs_list_n
		cells_with_mosaic=""
		cells_without_germline=""
		cells_with_mutation="-"
		if explanation_score_mosaic==1:
			cells_with_mosaic=",".join(ordered_col_sum[:cells_carrying_mutation_Nv].index)
		else:
			cells_with_mosaic="-"
		if explanation_score_germ==1:
			cells_without_germline=",".join(ordered_row_sum[:cells_carrying_mutation_Nv].index)
		else:
			cells_without_germline="-"
		print()
		cells_with_mutation=",".join(ordered_col_sum[:cells_carrying_mutation_Nv].index)
		output_line="\t".join(["\t".join(mutation.split("_")),str(explanation_score_mosaic),
						str(explanation_score_germ),cells_with_mutation])
		output_file_fh.write(output_line+"\n")
		es_dict[mutation]=[(explanation_score_germ,explanation_score_mosaic),(cells_without_germline,cells_with_mosaic)]
		
	pickle.dump(mutation_matrix_dict,mutation_matrix_file_fh)
	mutation_matrix_file_fh.close()
	output_file_fh.close()	
	for i in open(vaf_file):
		line=i.strip().split("\t")
		if i.startswith("Case"):
			vaf_ef_file_fh.write(i.strip()+"\tMosaic_score\tGermline_score\tMosaic_samples\tGermline_Samples\n")
			continue
		mutation="_".join(line[2:6])
		mosaic_score=es_dict[mutation][0][1]
		germline_score= es_dict[mutation][0][0]
		mosaic_samples=es_dict[mutation][1][1]
		germline_samples=es_dict[mutation][1][0]
		outline="\t".join([i.strip(),str(mosaic_score),str(germline_score),mosaic_samples, germline_samples])
		vaf_ef_file_fh.write(outline+"\n")
	vaf_ef_file_fh.close()
	
	
#def plotting(output_dir,list_of_samples,af_cutoff):
def plotting(output_dir,af_cutoff):
	explaination_score=os.path.join(output_dir,"explaination_score.txt")
	vaf_es_file=os.path.join(output_dir,"variant_vaf_per_pair_es.txt")
	pass
	
	# Plotting germline versus mosaic score heatmap
	df_es=pd.read_table(explaination_score)
	x=df_es["Mosaic_score"]
	y=df_es["Germline_score"]
	size_dict={}
	for n,i in enumerate(x):
		score_pair=str(i)+"_"+str(y[n])
		if score_pair in size_dict:
			size_dict[score_pair]+=1
		else:
			size_dict[score_pair]=1
	size=[]
	for n,i in enumerate(x):
		score_pair=str(i)+"_"+str(y[n])
		size.append(float(size_dict[score_pair]))	
	df_es.plot.scatter("Mosaic_score","Germline_score",c=size,s=[float(s / 10.0) for s in size],cmap='tab10')
	plt.tight_layout()
	plt.savefig(os.path.join(output_dir,"Explaination_score_scatter.png"))
	
	# Plotting per sample mosaic VAF plot
	
	# Plotting per sample germline VAF plot

	# 

	

def main():
	# Extracting passed arguments
	parser=argument_parse()
	arg=parser.parse_args()

	# Assigning values to variable
	manifest_file=arg.manifest_file
	output_dir=arg.output_dir
	af_cutoff=float(arg.af_cutoff)
	Util.ensure_dir(output_dir)

	# Extracting variant information from the manifest file.
	print("Extracting variant information")
	variant_dict,pairs_vaf_dict,list_of_samples=extract_mutation_information(manifest_file,output_dir)
	# variant_dict={mutation:[(case,control)]}
        # pairs_vaf_dict={pairs:{mutation:[vaf]}}
	# list_of_samples = list of all samples in the analysis

	# Generating explaination score
	print("Generating explaination scores")
	explaination_score(variant_dict,list_of_samples,output_dir)	

	# Plotting
	print("Plotting")
	#plotting(output_dir,list_of_samples,af_cutoff)
	plotting(output_dir,af_cutoff)
	

if __name__ == "__main__":
        Util.OsEnviron(os.environ)
        main()
