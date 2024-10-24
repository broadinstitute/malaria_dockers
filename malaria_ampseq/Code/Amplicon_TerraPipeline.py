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

import amplseq_functionalities as af
import asv_to_cigar as ac

def main():
	"""
	Implementation of the amplicon decontamination pipeline for use and 
	distribution in TERRA v=1.0

	Usage: python Code/Amplicon_TerraPipeline.py --config config.json ...
	
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
	parser.add_argument('--adaptor_removal', action="store_true", help="Specify if adaptor removal needed.")
	parser.add_argument('--demultiplex', action="store_true", help="Specify if reads must be separated by size.")
	parser.add_argument('--demultiplexed', action="store_true", help="Specify if reads have been demultiplexed.")
	parser.add_argument('--primer_removal', action="store_true", help="Specify if primer removal needed.")
	parser.add_argument('--dada2', action="store_true", help="Specifiy if standard preprocess merge with DADA2 is performed.")
	parser.add_argument('--postproc_dada2', action="store_true", help="Specify if postProcess of DADA2 results is perfomed.")
	parser.add_argument('--asv_to_cigar', action="store_true", help="Specify if the ASV to CIGAR transformation is perfomed.")

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
		forward_primers_file = config_inputs['forward_primers_file']
		reverse_primers_file = config_inputs['reverse_primers_file']
		reference_amplicons = config_inputs['reference_amplicons']
		if 'reference2'	in config_inputs.keys(): reference2 = config_inputs['reference2']
		if 'path_to_snv' in config_inputs.keys(): path_to_snv = config_inputs['path_to_snv']
		if 'pattern_fw' in config_inputs.keys(): pattern_fw = config_inputs['pattern_fw']
		if 'pattern_rv' in config_inputs.keys(): pattern_rv = config_inputs['pattern_rv']
		if 'Class' in config_inputs.keys(): Class = config_inputs['Class']
		if 'maxEE' in config_inputs.keys(): maxEE = config_inputs['maxEE']
		if 'trimRight' in config_inputs.keys(): trimRight = config_inputs['trimRight']
		if 'minLen' in config_inputs.keys(): minLen = config_inputs['minLen']
		if 'truncQ' in config_inputs.keys(): truncQ = config_inputs['truncQ']
		if 'matchIDs' in config_inputs.keys(): matchIDs = config_inputs['matchIDs']
		if 'max_consist' in config_inputs.keys(): max_consist = config_inputs['max_consist']
		if 'omegaA' in config_inputs.keys(): omegaA = config_inputs['omegaA']
		if 'saveRdata' in config_inputs.keys(): saveRdata = config_inputs['saveRdata']
		if 'justConcatenate' in config_inputs.keys(): justConcatenate = config_inputs['justConcatenate']
		if 'maxMismatch' in config_inputs.keys(): maxMismatch = config_inputs['maxMismatch']
		if 'adjust_mode' in config_inputs.keys(): adjust_mode = config_inputs['adjust_mode']
		if 'no_ref' in config_inputs.keys(): no_ref = config_inputs['no_ref']
		if 'strain' in config_inputs.keys(): strain = config_inputs['strain']
		if 'strain2' in config_inputs.keys(): strain2 = config_inputs['strain2']
		if 'polyN' in config_inputs.keys(): polyN = int(config_inputs['polyN'])
		if 'min_reads' in config_inputs.keys(): min_reads = int(config_inputs['min_reads'])
		if 'min_samples' in config_inputs.keys(): min_samples = int(config_inputs['min_samples'])
		if 'max_snv_dist' in config_inputs.keys(): max_snv_dist = int(config_inputs['max_snv_dist'])
		if 'max_indel_dist' in config_inputs.keys(): max_indel_dist = int(config_inputs['max_indel_dist'])
		if 'include_failed' in config_inputs.keys(): include_failed = eval(config_inputs['include_failed'])
		if 'exclude_bimeras' in config_inputs.keys(): exclude_bimeras = eval(config_inputs['exclude_bimeras'])
		if 'verbose' in config_inputs.keys(): verbose = eval(config_inputs['verbose'])

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
		af.flush_dir(res_dir, "Fq_metadata")
		af.create_meta(path_to_fq, res_dir, "Fq_metadata", "rawfilelist.tsv", pattern_fw, pattern_rv)

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
		af.flush_dir(res_dir, "AdaptorRem")
		meta = open(os.path.join(res_dir, "Fq_metadata", "rawfilelist.tsv"), 'r')
		samples = meta.readlines()
		
		for sample in samples: #[TODO: This can be parallelized - how do we keep the order of information]
			slist = sample.split()
			af.adaptor_rem(slist[0], slist[1], slist[2], res_dir, "AdaptorRem")
	
		af.create_meta(os.path.join(res_dir, "AdaptorRem"), res_dir, "AdaptorRem", "adaptorrem_meta.tsv",
			pattern_fw="*_val_1.fq.gz", pattern_rv="*_val_2.fq.gz")

	if args.demultiplex:
		print("Entering demultiplexing algorithm")
		meta = open(os.path.join(res_dir, "AdaptorRem", "adaptorrem_meta.tsv"), 'r')
		samples = meta.readlines()
		af.flush_dir(res_dir, "Demultiplex_by_Size")

		#Get the reads size
		print("Getting the read size")
		lengths_fw = []
		lengths_rv = []
		for sample in samples:
			longest_sequence_length_fw = af.find_longest_sequence_length(sample.split()[1])
			longest_sequence_length_rv= af.find_longest_sequence_length(sample.split()[2])
			lengths_fw.append(longest_sequence_length_fw)
			lengths_rv.append(longest_sequence_length_rv)

		percentile = lambda lengths: sorted(lengths)[int(len(lengths) * 0.95):]
		p_fw = percentile(lengths_fw)
		p_rv = percentile(lengths_rv)

		all_equal_fw = len(set(p_fw)) == 1
		all_equal_rv = len(set(p_rv)) == 1
		if all_equal_fw and all_equal_rv:
			print("Are all values in the top 95% percentile of read size equal?", all_equal_fw and all_equal_rv)
			print("The following values will be used as the standard read sizes in this run:")
		else:
			print("Are all values in the top 95% percentile of read size equal?", all_equal_fw and all_equal_rv)
			print("Largest values found for the forward and reverse read will be used. These values are:")
		read_size_fw = sorted(lengths_fw)[-1]
		read_size_rv = sorted(lengths_rv)[-1]
		print("Forward read:", read_size_fw)
		print("Reverse read:", read_size_rv)
		#print("If these sizes do not match the expected size of your technology, consider rerurring the pipeline after manually providing the size of your reads")

		#Get and remove the adapter sequence
		print("Removing the adapter sequence from the primers")
		adapter_fw = af.find_common_subsequence(forward_primers_file)
		adapter_rv = af.find_common_subsequence(reverse_primers_file)
		print("Adapter forward:", adapter_fw)
		print("Adapter reverse:", adapter_rv)
		af.remove_adapter(forward_primers_file, adapter_fw, 'primer_fw_no_adapter.fasta')
		af.remove_adapter(reverse_primers_file, adapter_rv, 'primer_rv_no_adapter.fasta')

		#Get the size of the reference ASVs
		print("Getting the size of the reference ASV")
		if os.path.exists(reference_amplicons):
			ref_files = [reference_amplicons]
			if os.path.exists(reference2):
				ref_files.append(reference2)
		else:
			sys.exit('At the least one reference is necessary to separate reads according to ASV target size.')

		#Get the ASV sizes
		asv_lengths = {}

		for ref_file in ref_files:
			for record in SeqIO.parse(ref_file, "fasta"):
				if record.id not in asv_lengths:
					asv_lengths[record.id] = len(record.seq)
				elif asv_lengths[record.id] < len(record.seq):
					#When two reference genomes are provided, if an ASV is already present, always reference the longest version.
					asv_lengths[record.id] = len(record.seq)

		with open(os.path.join(res_dir, "asv_lengths.tsv"), 'w', newline='') as csvfile:
			writer = csv.writer(csvfile)
			writer.writerow(['ASV', 'Length'])
			for seq_id, length in asv_lengths.items():
				writer.writerow([seq_id, length])

		print("Demultiplexing reads by size of reads according to their target amplicon")
		#Make dictionary of wells
		sample_number = [f"S{i}" for i in range(1, 193)]
		illumina_well = [f"{char}{num}" for char in list("ABCDEFGH") for num in range(1, 13)]*2
		sample_dict = dict(zip(sample_number, illumina_well))

		print(sample_dict)

		for sample in samples:
			slist = sample.split()
			print(slist)
			af.demultiplex_per_size(slist[0], slist[1], slist[2], 'primer_fw_no_adapter.fasta', 'primer_rv_no_adapter.fasta', res_dir, "Demultiplex_by_Size", read_size_fw, read_size_rv, asv_lengths, args.contamination, sample_dict)
		
		#Create Metafile for reads with no overlap
		af.create_meta(os.path.join(res_dir, "Demultiplex_by_Size"), res_dir, "Demultiplex_by_Size", "demux_nop_meta.tsv",
			pattern_fw="*_nop_L001_R1_001.fastq.gz", pattern_rv="*_nop_L001_R2_001.fastq.gz")
		af.create_meta(os.path.join(res_dir, "Demultiplex_by_Size"), res_dir, "Demultiplex_by_Size", "demux_op_meta.tsv",
			pattern_fw="*_op_L001_R1_001.fastq.gz", pattern_rv="*_op_L001_R2_001.fastq.gz")
	else:
		print("Skipping demultiplexing algorithm")
		meta = open(os.path.join(res_dir, "AdaptorRem", "adaptorrem_meta.tsv"), 'r')
		samples = meta.readlines()
		af.flush_dir(res_dir, "Demultiplex_by_Size")

		#Get and remove the adapter sequence
		print("Removing the adapter sequence from the primers")
		adapter_fw = af.find_common_subsequence(forward_primers_file)
		adapter_rv = af.find_common_subsequence(reverse_primers_file)
		print("Adapter forward:", adapter_fw)
		print("Adapter reverse:", adapter_rv)
		af.remove_adapter(forward_primers_file, adapter_fw, 'primer_fw_no_adapter.fasta')
		af.remove_adapter(reverse_primers_file, adapter_rv, 'primer_rv_no_adapter.fasta')

		for sample in samples:
			slist = sample.split()
			output_fastq_fw_op = os.path.join(res_dir, "Demultiplex_by_Size", f"{slist[0]}_op_L001_R1_001.fastq.gz")
			output_fastq_rv_op = os.path.join(res_dir, "Demultiplex_by_Size", f"{slist[0]}_op_L001_R2_001.fastq.gz")
			os.system(f"cp {slist[1]} {output_fastq_fw_op}")
			os.system(f"cp {slist[2]} {output_fastq_rv_op}")
			af.create_meta(os.path.join(res_dir, "Demultiplex_by_Size"), res_dir, "Demultiplex_by_Size", "demux_op_meta.tsv",
				pattern_fw="*_op_L001_R1_001.fastq.gz", pattern_rv="*_op_L001_R2_001.fastq.gz")
		
	#Remove primers
	#For a set where all reads have overlap
	if args.primer_removal: # [TODO: Need an indicator whether primers contain combinatorial_indices or not]
		# if args.combinatorial_indices:...
		# else:

		print("Removing primers")
		#Extract primer for the target without amplicons
		fw = 'primer_fw_no_adapter.fasta'
		rv = 'primer_rv_no_adapter.fasta'
			
		records_fw = SeqIO.parse(fw, 'fasta')
		common_subsequences = af.find_common_subsequences(records_fw)
		af.write_common_subsequences_to_fasta(common_subsequences, 'amp_primer_fw.fasta')
		records_rv = SeqIO.parse(rv, 'fasta')
		common_subsequences = af.find_common_subsequences(records_rv)
		af.write_common_subsequences_to_fasta(common_subsequences, 'amp_primer_rv.fasta')

		af.flush_dir(res_dir, "PrimerRem")

		#Trim primers off non-overlapping targets
		if args.demultiplex:
			meta = open(os.path.join(res_dir, "Demultiplex_by_Size", "demux_nop_meta.tsv"), 'r')
			samples = meta.readlines()
			for sample in samples:
				slist = sample.split()
				# af.trim_primer(slist[0], slist[1], slist[2], res_dir, "PrimerRem", "amp_primer_fw.fasta", "amp_primer_rv.fasta", "mixed_nop")
				# af.trim_primer(slist[0], slist[1], slist[2], res_dir, "PrimerRem", fw, rv, "mixed_nop")
				af.trim_primer(slist[0], slist[1], slist[2], res_dir, "PrimerRem", forward_primers_file, reverse_primers_file, "mixed_op")

			#Metafile for trimmed non-op target reads
			af.create_meta(os.path.join(res_dir, "PrimerRem"), res_dir, "PrimerRem", "mixed_nop_prim_meta.tsv", 
				pattern_fw="*_mixed_nop_1.fq.gz", pattern_rv="*_mixed_nop_2.fq.gz")

		#Trim primers off overlapping targets
		meta = open(os.path.join(res_dir, "Demultiplex_by_Size", "demux_op_meta.tsv"), 'r')
		samples = meta.readlines()
		for sample in samples:
			slist = sample.split()
			# af.trim_primer(slist[0], slist[1], slist[2], res_dir, "PrimerRem", "amp_primer_fw.fasta", "amp_primer_rv.fasta", "mixed_op")
			# af.trim_primer(slist[0], slist[1], slist[2], res_dir, "PrimerRem", fw, rv, "mixed_op")
			af.trim_primer(slist[0], slist[1], slist[2], res_dir, "PrimerRem", forward_primers_file, reverse_primers_file, "mixed_op")

		#Metafile for trimmed overlapping target reads
		af.create_meta(os.path.join(res_dir, "PrimerRem"), res_dir, "PrimerRem", "mixed_op_prim_meta.tsv",
			pattern_fw="*_mixed_op_1.fq.gz", pattern_rv="*_mixed_op_2.fq.gz")

#	#For a set that mixes reads with and without overlap
	if args.dada2:
		print("Running DADA2")
		#Run DADA2 on op targets
		af.flush_dir(res_dir, "DADA2_OP", "QProfile")
		path_to_meta = os.path.join(res_dir, "PrimerRem", "mixed_op_prim_meta.tsv")
		justConcatenate=0	
		af.run_dada2(path_to_meta, path_to_fq, path_to_flist, Class, maxEE, trimRight, minLen, truncQ, matchIDs, max_consist, omegaA, justConcatenate, maxMismatch, saveRdata, res_dir, "DADA2_OP", args.terra)
		seqtab_op = os.path.join(res_dir, 'DADA2_OP', 'seqtab.tsv')
		bimera_op = os.path.join(res_dir, 'DADA2_OP', 'ASVBimeras.txt')

		if args.demultiplexed:
			#Run DADA2 on non-op targets
			af.flush_dir(res_dir, "DADA2_NOP", "QProfile")
			path_to_meta = os.path.join(res_dir, "PrimerRem", "mixed_nop_prim_meta.tsv")
			justConcatenate=1	
			af.run_dada2(path_to_meta, path_to_fq, path_to_flist, Class, maxEE, trimRight, minLen, truncQ, matchIDs, max_consist, omegaA, justConcatenate, maxMismatch,saveRdata, res_dir, "DADA2_NOP", args.terra)
			seqtab_nop = os.path.join(res_dir, 'DADA2_NOP', 'seqtab.tsv')
			bimera_nop = os.path.join(res_dir, 'DADA2_NOP', 'ASVBimeras.txt')

		#ASV modification block for non-op targets and merge two ASV tables
		if os.path.exists(reference_amplicons) and args.demultiplexed:
			if args.terra:
				path_to_program = os.path.join("/", "Code/adjustASV.R")
			else:
				path_to_program = os.path.join("Code/adjustASV.R")
			adjASV = ['Rscript', path_to_program, '-s', seqtab_nop, '-ref', reference_amplicons,
			'-dist', adjust_mode,
			'-o', os.path.join(res_dir, 'DADA2_NOP', 'correctedASV.txt')]
			print(adjASV)
			procASV = subprocess.Popen(adjASV)
			procASV.wait()
			seqtab_corrected = os.path.join(res_dir, 'DADA2_NOP', 'seqtab_corrected.tsv')
			seqtab = af.merge_seqtab(seqtab_op, seqtab_corrected)
			bimera = af.merge_bimeras(bimera_op, bimera_nop)
			seqtab.fillna(0).to_csv(os.path.join(res_dir, 'seqtab.tsv'), sep = "\t", index=False)
			bimera.to_csv(os.path.join(res_dir, 'ASVBimeras.txt'), sep = "\t", index=False)
		elif args.demultiplexed:
			print('--reference_amplicons file not found. skipping ASV correction..')
			seqtab = af.merge_seqtab(seqtab_op, seqtab_nop)
			bimera = af.merge_bimeras(bimera_op, bimera_nop)
			seqtab.fillna(0).to_csv(os.path.join(res_dir, 'seqtab.tsv'), sep = "\t", index=False)
			bimera.to_csv(os.path.join(res_dir, 'ASVBimeras.txt'), sep = "\t", index=False)
		else:
			print("No non-op targets found. Skipping ASV merging")
			print("Creating placeholder files for non-op targets")
			seqtab_final = os.path.join(res_dir, 'seqtab.tsv')
			bimera_final = os.path.join(res_dir, 'ASVBimeras.txt')
			os.system(f"cp {seqtab_op} {seqtab_final}")
			os.system(f"cp {bimera_op} {bimera_final}")

	#Perform postprocessing (postproc_dada2)
	if args.postproc_dada2:
		print("Performing post processing")
		af.flush_dir(res_dir, "PostProc_DADA2")
		
		path_to_seqtab = os.path.join(res_dir, 'seqtab.tsv')

		if args.terra:
			path_to_program = os.path.join("/", "Code/postProc_dada2.R")
		else:
			path_to_program = os.path.join("Code/postProc_dada2.R")

		if "path_to_snv" not in locals():
			path_to_snv = "No_File"
		if "reference2" not in locals():
			reference2 = "No_File"
		
		postProc = ['Rscript', path_to_program,
				'-s', path_to_seqtab,
				'-b', os.path.join(res_dir, 'ASVBimeras.txt'),
				'-snv', os.path.join(path_to_snv),
				'--indel_filter', '0.01',
				'-o', os.path.join(res_dir, 'PostProc_DADA2', 'ASVTable.txt'),
				'--fasta']
		
		if no_ref == 'True':
			postProc.extend(['-no_ref'])
		else:
			postProc.extend(['--reference_amplicons', reference_amplicons, '--strain', strain])
			if os.path.exists(reference2):
				postProc.extend(['--reference2', reference2, '--strain2', strain2])
		print(postProc)
		procASV = subprocess.Popen(postProc)
		procASV.wait()

	#ASV to CIGAR
	#Convert ASVs from DADA2 pipeline to pseudo-CIGAR strings.
	if args.asv_to_cigar:
		print("Converting ASVs to CIGARs")
		af.flush_dir(res_dir, "ASV_to_CIGAR", "alignments")

		path_to_seqtab = os.path.join(res_dir, 'seqtab.tsv')
		path_to_fasta = os.path.join(res_dir, "PostProc_DADA2", "ASVSeqs.fasta") #Fasta file of ASV sequences from DADA2 pipeline"
		path_to_table = os.path.join(res_dir, "PostProc_DADA2", "ASVTable.txt") #ASV table from DADA2 pipeline
		path_to_out = os.path.join(res_dir, "CIGARVariants_Bfilter.out.tsv") #Output seqtab tsv file with amplicon/variant counts
		path_asv_to_cigar = os.path.join(res_dir, "ASV_to_CIGAR", "ASV_to_CIGAR.out.txt") #Output file for ASV -> CIGAR string table 
		path_to_zero_read_samples = os.path.join(res_dir, "ASV_to_CIGAR", "ZeroReadsSampleList.txt") #Output file for 
		path_to_amp_db = reference_amplicons #Amplicon sequence fasta file
		path_to_alignments = os.path.join(res_dir, "ASV_to_CIGAR", "alignments") #Directory to store ASV alignment files

		print(f"INFO: Loading {path_to_amp_db}")
		amplicons = ac.parse_amp_db(path_to_amp_db)
		if not amplicons:
			print(f"ERROR: No amplicons in {path_to_amp_db}")
			sys.exit(1)

		if os.path.exists("amp_mask.txt"):
			print(f"INFO: Loading amp_mask.txt")
			mask = ac.parse_dustmasker("amp_mask.txt")
		else:
			print(f"INFO: No mask data specified.")
			mask = {}

		print(f"INFO: Loading {path_to_fasta}")
		asvs = ac.get_asv_seqs(path_to_fasta)
		if not asvs:
			print(f"ERROR: No ASV sequences in {path_to_fasta}")
			sys.exit(1)

		print(f"INFO: Parsing {path_to_table} with total reads >= {min_reads}, samples >= {min_samples}, snv_dist <= {max_snv_dist}, indel_dist <= {max_indel_dist}")

		if include_failed:
			print("WARNING: Including ASVs that failed post-DADA2 filters! This is not recommended.")
		else:
			print("INFO: Excluding ASVs that failed post-DADA2 filters.")

		if exclude_bimeras:
			print("INFO: Excluding ASVs that DADA2 marked as bimeras.")

		bins = ac.parse_asv_table(path_to_table, min_reads=min_reads, min_samples=min_samples, max_snv_dist=max_snv_dist, max_indel_dist=max_indel_dist, include_failed=include_failed, exclude_bimeras=exclude_bimeras) #This function only matches to the first strain.
		if not bins:
			print(f"ERROR: No useable data in {path_to_table}")
			sys.exit(1)

		print(f"INFO: Writing amplicon fasta files to {path_to_alignments}")
		ac.write_amplicon_fastas(asvs, bins, amplicons, outdir=path_to_alignments)

		print("INFO: Running MUSCLE aligner on amplicon fasta files. Please wait...")
		ac.run_muscle(bins, outdir=path_to_alignments)

		print("INFO: Parsing alignments to CIGAR strings")
		cigars = ac.parse_alignments(bins, mask=mask, min_homopolymer_length=polyN, outdir=path_to_alignments, verbose=False)
		if not cigars:
			print("ERROR: could not determine CIGAR strings")
			sys.exit(1)

		if path_asv_to_cigar:
			ac.write_cigar_strings(cigars, path_asv_to_cigar)
			print(f"INFO: Wrote ASV->CIGAR table to {path_asv_to_cigar}")

		print(f"INFO: Converting DADA2 seqtab file {path_to_seqtab} to {path_to_out}")
		if ac.convert_seqtab(path_to_seqtab, cigars, path_to_out):
			print("INFO: Completed conversion of seqtab to CIGAR variants successfully!")

			if ac.get_zero_reads_samples(path_to_out, path_to_zero_read_samples):
				print("INFO: Obtained samples with zero reads successfully!")
		
if __name__ == "__main__":
	main()
