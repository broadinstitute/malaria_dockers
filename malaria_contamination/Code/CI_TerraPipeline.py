#!/usr/bin/env python

import os
import sys
import argparse
import json
import subprocess
import csv
import gzip
import pandas as pd

from Bio import SeqIO

import ci_detection as ci

def main():
	"""
	Implementation of the amplicon contamination detection pipeline for use and 
	distribution in TERRA v=1.0

	Usage: python Code/CI_TerraPipeline.py --config config.json ...
	
	Returns:
	None. However, see the functions documentation.
	"""

	### LOAD ARGUMENTS

	#Parse command line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('--config', help="Path to config.json file.", required =True)
	parser.add_argument('--terra', action="store_true", help="Specify whether the pipeline is being run in the Terra platform.")
	parser.add_argument('--meta', action="store_true", help="Specify if metadata file must be created. This flag runs the pipeline from the beginning.")
	parser.add_argument('--repo', action="store_true", help="Specify if the reports must be created.")
	parser.add_argument('--contamination', action="store_true", help="Specify if contamination detection is performed.")
	parser.add_argument('--adaptor_removal', action="store_true", help="Specify if adaptor removal needed.")

	args = parser.parse_args()

	#Check minimum arguments and contradicting flags
	if args.terra:
		print("Pipeline is running in Terra. Adjusted paths will be used.")
	else:
		print("Pipeline not running in Terra. Default paths will be used.")
				
	#Configuration aguments will be parsed from config.json
	with open(args.config, 'r') as config_file:
		config_inputs = json.load(config_file)
		path_to_fq = config_inputs['path_to_fq']
		path_to_flist = config_inputs['path_to_flist']
		if 'pattern_fw' in config_inputs.keys(): pattern_fw = config_inputs['pattern_fw']
		if 'pattern_rv' in config_inputs.keys(): pattern_rv = config_inputs['pattern_rv']

	### PREPARE OUTPUT DIRECTORIES
	#Generate the path to the Results directory.
	global res_dir
	global rep_dir
	res_dir = os.path.abspath(os.path.join("Results"))
	rep_dir = os.path.abspath(os.path.join("Report"))

	#Restart Results and Report directory if the workflow runs form the start
	if args.meta:
		os.system("rm -rf " + res_dir)
		os.mkdir(res_dir)

	if args.repo:
		os.system("rm -rf " + rep_dir)
		os.mkdir(rep_dir)
	
	#Create metadata files
	if args.meta:
		print("Creating meta files")
		ci.flush_dir(res_dir, "Fq_metadata")
		ci.create_meta(path_to_fq, res_dir, "Fq_metadata", "rawfilelist.tsv", pattern_fw, pattern_rv)

		#List missing file from the complete list of files
		with open(path_to_flist, 'r') as input_file:
			next(input_file)
			samples = [line.split(',')[0] for line in input_file]

		with open('Results/Fq_metadata/rawfilelist.tsv', 'r') as raw_files:
			raw_file_samples = [line.split('\t')[0] for line in raw_files]

		missing_samples = [sample.strip() for sample in samples if sample not in raw_file_samples]

		with open('Results/missing_files.tsv', 'w') as output_file:
			output_file.write('\n'.join(missing_samples))

	### EXECUTE

	#Remove adaptors
	#Most sequences must be adaptor free; just in case, run this step to eliminate any lingering adaptors.
	if args.adaptor_removal:
		print("Removing adaptors")
		ci.flush_dir(res_dir, "AdaptorRem")
		meta = open(os.path.join(res_dir, "Fq_metadata", "rawfilelist.tsv"), 'r')
		samples = meta.readlines()
		
		for sample in samples: #[TODO: This can be parallelized - how do we keep the order of information]
			slist = sample.split()
			ci.adaptor_rem(slist[0], slist[1], slist[2], res_dir, "AdaptorRem")
	
		ci.create_meta(os.path.join(res_dir, "AdaptorRem"), res_dir, "AdaptorRem", "adaptorrem_meta.tsv",
			pattern_fw="*_val_1.fq.gz", pattern_rv="*_val_2.fq.gz")

	#Merge forward and reverse reads with bbmerge. Only for reads that overlap and
	#Process merge report and generate RStudio plots
	if args.contamination:
		print("Concatenating reads for contamination detection")
		ci.flush_dir(res_dir, "Merge")
		meta = open(os.path.join(res_dir, "AdaptorRem", "adaptorrem_meta.tsv") , 'r')
		samples = meta.readlines()
		for sample in samples:
			slist = sample.split()
			ci.mergereads(slist[0], slist[1], slist[2], res_dir, "Merge")

		ci.flush_dir(rep_dir, "Merge")
		meta = open(os.path.join(res_dir, "Merge", "merge_meta.tsv"), 'r')
		samples = meta.readlines()
		
		print("Extracting fields from bbmerge reports")
		for sample in samples:
			slist = sample.split()
			ci.extract_bbmergefields(slist[0], slist[1], slist[3], path_to_flist, res_dir, rep_dir, "Merge", args.terra)
		
if __name__ == "__main__":
	main()
