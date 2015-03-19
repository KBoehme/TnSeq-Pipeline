#Author: Kevin Boehme
#Date: 12-20-2014
#This program is part of a TnSeq analysis pipeline, designed to take raw fastq reads, and produce tabulated data on hop count occurance.

import sys
sys.dont_write_bytecode = True

import ConfigParser

import logging
import gzip
import os
import collections
import glob
import subprocess
import re
from objects import *
from pprint import pprint

from time import time
from datetime import datetime

class hops_pipeline(object):
	"""docstring for hops_pipeline"""
	def __init__(self):

		#Inputs
		self.input_files = []
		self.ref   = ""
		self.ptt   = []
		self.out   = ""

		#Parameters
		self.transposon = ""
		self.mismatches = 0
		self.gene_trim = 0
		self.read_length = 0
		#bools
		self.debug = False
		self.normalize = True

		#Intermediate files for output
		self.int_prefix = []
		self.int_trimmed = [] # list of files for the intermediate trimmed read files.
		self.int_sam = [] # list of sam output files.
		self.tabulated_filename = ""

		self.num_conditions = 0

		#information variables
		self.removed_notn = 0
		self.removed_tooshort = 0
		self.kept_reads = 0
		self.original_read_count = 0
		self.starttime = time()

		# Variables for tabulating hop hits.
		self.gene_info = []
		self.normalization_coefficients = []
		

		

############ Read Config File and Run Pipeline ###############################
	def read_config(self, config_path):
		cp = ConfigParser.RawConfigParser()
		try:
			cp.read(config_path)
		except:
			sys.exit("Error reading config file.")
		# input paths
		self.input_files = glob.glob(cp.get('input', 'Reads'))
		self.ref = cp.get('input', 'BowtieReference')
		self.ptt = glob.glob(cp.get('input', 'Ptt'))
		self.out  = cp.get('input', 'Out')

		re.match
		if len(self.input_files) == 0:
			sys.exit('Error with input Reads parameter.')
		elif self.ref == None:
			sys.exit('Error with BowtieReference parameter.')
		elif len(self.ptt) == 0:
			sys.exit('Error with Ptt parameter.')
		elif self.out == None:
			sys.exit("Error with Out parameter.")

		#Parameters
		try:
			self.transposon = cp.get('parameters', 'Transposon')
		except:
			sys.exit('Error with Transposon parameter.')


		if self.transposon == None:
			sys.exit('Error with Transposon parameter.')
		elif not (re.match('^[ACGTacgt]+$',self.transposon)):
			sys.exit('Error with Transposon parameter (Make sure it only contains [ATCG]).')
	
		try:
			self.mismatches = int ( cp.get('parameters', 'Mismatches') )
		except:
			sys.exit('Error with Mismatches parameter (Not an integer).')

		if self.mismatches < 0:
			sys.exit('Mismatches parameter is negative.')
		elif self.mismatches >= len(self.transposon):
			logging.info("Mismatches parameter is same length or greater than transposon sequence. Be aware that all reads will be mapped because of this.")

		try:
			self.gene_trim = int ( cp.get('parameters', 'GeneTrim') )
		except:
			sys.exit('Error with GeneTrim parameter (Not an integer).')


		if self.gene_trim < 0:
			sys.exit('GeneTrim parameter is negative.')
		elif self.gene_trim > 99:
			sys.exit('Error with GeneTrim parameter (Must be 99 or smaller).')


		try:
			self.read_length = int ( cp.get('parameters', 'ReadLength') )
		except:
			sys.exit('Error with ReadLength parameter (Not an integer).')



		if self.read_length < 0:
			sys.exit('ReadLength parameter is negative.')

		#Options
		try:
			self.debug = cp.getboolean('options', 'Debug')
		except:
			sys.exit('Error with Debug parameter (Not True/False).')

		try:
			self.normalize = cp.getboolean('options','Normalize')
		except:
			sys.exit('Error with Normalize parameter (Not True/False)')


		# Generate other variables
		self.num_conditions = len(self.input_files)
		self.output_directory = os.path.dirname(config_path)
		if self.output_directory == "":
			self.output_directory = "./"
		else:
			self.output_directory += "/"
		self.tabulated_filename = self.output_directory + "output_files/" + self.out + "-HOPS.txt"
		self.gene_tabulated_filename = self.output_directory + "output_files/" + self.out + "-GENE.txt"
	

		if not os.path.exists( self.output_directory + "output_files/"):
			subprocess.check_output(["mkdir", self.output_directory + "output_files"])
		
		if not os.path.exists( self.output_directory + "output_files/"):
			subprocess.check_output(["mkdir", self.output_directory + "output_files"])
		
		self.set_up_logger(self.output_directory + "output_files/" + self.out + ".log")

		with open(config_path,'r') as f:
			logging.info("\n\n=============== Config File Settings ===============\n\n" + ''.join(f.readlines()) + "\n--------------------------------------\n\n")

		if not os.path.exists( self.output_directory + "intermediate_files"):
			subprocess.check_output(["mkdir", self.output_directory + "intermediate_files"])

		for f in self.input_files:	
			file_prefix = os.path.basename(f).split('.')[0]
			self.int_prefix.append(file_prefix)

			#create both trimmed and sam intermediate names for this particular input file.
			trimmed = ""
			fi = None
			trimmed = self.output_directory + "intermediate_files/" + file_prefix + "-trimmed.fasta"
			fi = open(trimmed,'w')

			self.int_trimmed.append(fi)
			self.int_sam.append(self.output_directory + "intermediate_files/" + file_prefix + ".sam")

	def run_pipeline(self):
		self.process_reads() #pass in number corresponding to location of input file in list.
		self.call_bowtie2()
		self.process_sam()
		self.print_time_output("Total run time,", self.starttime)


############ Useful Functions ###############################
	def set_up_logger(self, log_file_name):
		logging.basicConfig( filename=log_file_name , filemode='w', level=logging.INFO,format='%(message)s' )
		logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
		
		logging.info("=========================================")
		logging.info("             Tn-Seq Pipeline")
		logging.info("The MIT License (MIT) Copyright (c) 2014 Kevin")
		logging.info("     Date: "+ str(datetime.now()))
		logging.info("           Author: Kevin Boehme")
		logging.info("      Email: kevinlboehme@gmail.com")
		logging.info("=========================================")

	def debugger(self, *text):
		if self.debug:
			sys.stdout.write("[DEBUG] ")
			for i in text:
				sys.stdout.write(str(i) + " ")
			print

	def print_time_output(self, command, start):
		time_to_run = time()-start
		if time_to_run > 60:
			time_to_run = time_to_run/60
			logging.info(command + " ran in " + "%.2f" % (time_to_run) + " minutes." + "\n")
		else:
			logging.info(command + " ran in " + "%.2f" % (time_to_run) +" seconds." + "\n")
	
	def fuzzy_match_beginning(self, pattern, genome, mismatches):
		for i in range(3):
			chunk = genome[i : i + len(pattern)]
			# now compare chunk with pattern to see if they match with at least mismatches.
			if(self.compareChunks(pattern, chunk, mismatches)):
				return i+len(pattern)
		return -1

	def compareChunks(self, pattern, chunk, mismatches):
		misses = 0
		for i in range(len(pattern)):
			if(pattern[i] != chunk[i]):
				misses = misses + 1
				if(misses > mismatches):
					return False
		return True


############ Process Reads ###############################
	def print_summary_stats(self, start_process_reads_time, out_file_num):
		prob = float(self.removed_notn)/float(self.original_read_count)
	
		#Done processing file, print out stats.
		self.print_time_output("Processing " + str(self.original_read_count) + " reads ", start_process_reads_time)

		logging.info("Removed " + str(self.removed_notn) 
			+ " reads with no detected transposon sequence (" 
			+ str(self.removed_notn) + "/"+ str(self.original_read_count)
			+ ") = " + "%.2f%%" % (prob * 100) + ".")
		logging.info("Removed " + str(self.removed_tooshort) 
			+ " reads that were too short (Less than " + str(self.read_length) + " bp after transposon trimming).")
		logging.info("")
		logging.info("Kept " + str(self.kept_reads) + " reads written to " + self.int_trimmed[out_file_num].name+"\n")
		logging.info("         ---------------\n")

	def process_reads(self):

		logging.info("\n\n=============== Read pre-processing ===============\n")
		logging.info("         ---------------\n")		
		for out_file_num in range(self.num_conditions):
			start_process_reads_time = time()
			filepath = self.input_files[out_file_num]

			self.removed_notn = 0
			self.removed_tooshort = 0
			self.kept_reads = 0
			self.original_read_count = 0

			isfastq = None

			if filepath.find("fasta") != -1:
				isfastq = False
			elif filepath.find("fastq") != -1:
				isfastq = True
			else:
				logging.error("Didn't find the string fasta or fastq in input files. Please make sure your data has a fasta or fastq file extension (Could also be gzipped).")
				sys.exit('')

			#try:
			logging.info( "Input fastq = " + filepath + "." )
			f = None
			if (filepath[-3:] == ".gz"):
				f = gzip.open(filepath)
			else:
				f = open(filepath)

			if isfastq:
				self.read_fastq(f, out_file_num)
			else:
				self.read_fasta(f, out_file_num)

			self.print_summary_stats(start_process_reads_time, out_file_num)

			#except:
			#	logging.error("Failed to process the reads.")
		for f in self.int_trimmed:
			f.close()
		logging.info("--------------------------------------\n\n")

	def read_fastq(self, f, out_file_num):
		t0 = time()
		while True:
			if self.original_read_count != 0 and self.original_read_count % 1000000 == 0:
				logging.info('Processed ' + str(self.original_read_count) + ' reads.')
			name = f.readline().strip()[1:] #Strip the @ sign infront of fastq names.
			seq = f.readline().strip().upper()
			plus = f.readline().strip()
			score = f.readline().strip()
			if not name or not seq or not plus or not score:
				break #We are done, lets break out of the loop.

			#we have all the contents of one read, now lets look for the transposon.
			self.process_read(name, seq, out_file_num)
			self.original_read_count += 1

		else: #Its a fasta file.
			while True:

				if self.original_read_count != 0 and self.original_read_count % 1000000 == 0:
					logging.info('Processed ' + str(self.original_read_count) + ' reads.')
				name = f.readline().strip()[1:] #Strip the @ sign infront of fastq names.
				seq = f.readline().strip().upper()
				if not name or not seq: break #We are done, lets break out of the loop.

				#we have all the contents of one read, now lets look for the transposon.
				print "processing read",seq
				self.process_read(name, seq, out_file_num)
				self.original_read_count += 1

	def read_fasta(self, f, out_file_num):
		t0 = time()
		while True:
			print "looping through fasta."

			if self.original_read_count != 0 and self.original_read_count % 1000000 == 0:
				logging.info('Processed ' + str(self.original_read_count) + ' reads.')
			name = f.readline().strip()[1:] #Strip the @ sign infront of fastq names.
			seq = f.readline().strip().upper()
			if not name or not seq: break #We are done, lets break out of the loop.

			#we have all the contents of one read, now lets look for the transposon.
			print "processing read",seq
			self.process_read(name, seq, out_file_num)
			self.original_read_count += 1

	def process_read(self, name, seq, out_file_num):
			tn_trimmed_seq  = self.remove_transposon(name, seq)
			if tn_trimmed_seq: # True if sequence had a transposon. False otherwise.
				prepped_seq = self.qc_genomic_region(tn_trimmed_seq)
				if prepped_seq:
					#write to fasta file.
					self.kept_reads += 1
					self.write_to_output(name, prepped_seq, self.int_trimmed[out_file_num])
				else:
					self.removed_tooshort += 1
					pass
			else: #no match.
				self.removed_notn += 1
				pass

	def write_to_output(self, name, seq, f):
		f.write(">"+name+"\n")
		f.write(seq+"\n")

	def remove_transposon(self, name, seq):

		num = self.fuzzy_match_beginning( self.transposon, seq, self.mismatches)
		#capture only the bacterial genome region.
		tn_trimmed_seq = seq[num:]

		if num != -1:
			return tn_trimmed_seq
		else:
			return None

	def qc_genomic_region(self, seq):
		#noNread = seq.replace("N","") #remove all N.
		if len(seq) < self.read_length:
			return None
		else:
			return seq[0:self.read_length]



############ Bowtie2 ###############################
	def call_bowtie2(self):

		logging.info("\n\n=============== Bowtie2 Mapping ===============\n")
		logging.info("         ---------------\n")		
		for out_file_num in range(self.num_conditions):
			start_time = time()
			
			#bowtie_command = ["bowtie2", "-x", self.ref,"--phred33",
			#                    "-f",self.int_trimmed[out_file_num].name, "-D" , "25" , "-R" ,"3", "-L" , "10" , "-i" , "S,1,0.50" ,
			#                    "-S",self.int_sam[out_file_num], "--no-hd"]

			bowtie_command = ["bowtie2", "-x", self.ref,"--phred33",
								"-f",self.int_trimmed[out_file_num].name,
								"-S",self.int_sam[out_file_num],"--no-hd"]
								

			logging.info("Bowtie Command Used: " + ' '.join(bowtie_command)+"\n\n")


			logging.info("Writing output to = " + self.int_sam[out_file_num]+"\n")

			try:
				logging.info(subprocess.check_output(bowtie_command,stderr=subprocess.STDOUT))#,shell=True))
			except:
				logging.error("Bowtie2 doesn't seem to be installed. Make sure it can be run with the command: bowtie2")
				sys.exit('Exiting')

			self.print_time_output("Bowtie2",start_time)
			logging.info("Zipping up trimmed file.\n")
			subprocess.check_output(["gzip","-f", self.int_trimmed[out_file_num].name])
			logging.info("         ---------------\n")

		logging.info("--------------------------------------\n\n")
		return True



############ Tabulate Sams ###############################
	def process_sam(self):
		start_time = time()
		self.prepare_gene_info()
		logging.info("\n\n=============== Process SAM File ===============\n")
		logging.info("         ---------------\n")
		self.read_sam_file()
		logging.info("         ---------------\n")		
		self.write_output()
		logging.info("         ---------------\n")

		self.print_time_output("Done processing SAM files,", start_time)
		logging.info("--------------------------------------\n\n")

	def prepare_gene_info(self):
		start_time = time()
		for i,file_name in enumerate(self.ptt):
			name = os.path.basename(file_name).split('.')[0]
			with open(file_name, 'r') as f:
				new_chrom = Chromsome(name)
				# Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
				title = f.readline() # Title
				num_prot = f.readline() # number of proteins
				f.readline() # header


				new_chrom.fill_info(title, num_prot)
				self.gene_info.append(new_chrom)
				our_genes = []
				for line in f:
					ptt_entry = line.split()

					new_gene = Gene()
					new_gene.create_from_ptt_entry(ptt_entry,self.gene_trim, self.num_conditions)

					our_genes.append(new_gene)
				self.gene_info[i].set_gene_list(sorted(our_genes))

		#This function finishes up the ptt info by adding another column "Order" as well as generating intergenic regions.
		self.add_intergenic_regions_and_order_column()

	def add_intergenic_regions_and_order_column(self):
		#Lets add in the order column and intergenic regions.
		for j,chrom in enumerate(self.gene_info):
			self.debugger("On replicon = " + str(chrom))
			example_gene = chrom.gene_list[0]
			self.debugger("Example gene from this replicon = " + str(example_gene))
			gene_name = example_gene.synonym
			prefix = "".join(re.findall("[a-zA-Z]+", gene_name))
			self.debugger("Prefix = " + prefix)
			num_genes = len(chrom.gene_list)
			self.debugger("Number of genes on replicon = " + str(num_genes))
			zfill_digits = len(str(num_genes))
			self.debugger("Digits to use in zfill = " + str(zfill_digits))


			intergenic_genes = []
			for i,cur_gene in enumerate(chrom.gene_list,start=1):
				ordered_name = prefix + str(i).zfill(zfill_digits)
				cur_gene.order_code = ordered_name

				# Now lets add intergenic regions.
				new_gene = Gene()
				if i == 1: # On the first gene.
					if cur_gene.start > 1: # We have some room to capture.
						new_gene.create_intergenic_region(1, cur_gene.start - 1, "int_BEG-"+cur_gene.synonym, self.num_conditions)
						intergenic_genes.append(new_gene)
					if cur_gene.end + 1 < chrom.gene_list[i].start: # We have some room to capture.
						new_gene.create_intergenic_region(cur_gene.end + 1, 
							chrom.gene_list[i].start - 1, 
							"int_" + cur_gene.synonym + "-" + chrom.gene_list[i].synonym, self.num_conditions)
						intergenic_genes.append(new_gene)
				elif i == chrom.num_proteins: # On the last gene.
					if cur_gene.end < chrom.end: # We have some room to capture.
						new_gene.create_intergenic_region(cur_gene.end + 1, chrom.end, "int_"+cur_gene.synonym+"-END", self.num_conditions)
						intergenic_genes.append(new_gene)
				else: # Make sure the end of the current gene and the beginning of the next gene have a space.
					if cur_gene.end + 1 < chrom.gene_list[i].start: # We have some room to capture.
						new_gene.create_intergenic_region(cur_gene.end + 1, 
							chrom.gene_list[i].start - 1, 
							"int_" + cur_gene.synonym + "-" + chrom.gene_list[i].synonym, self.num_conditions )
						intergenic_genes.append(new_gene)
			#Now lets mash those intergenic regions with the current genes.
			self.gene_info[j].gene_list = sorted(self.gene_info[j].gene_list + intergenic_genes)

	def update_progress(self,progress):
			sys.stdout.write('\r[{0}] {1}%'.format('#'*(progress/10), progress))

	def read_sam_file(self):
		sam_file_contents = {}
		for i in self.gene_info:
			sam_file_contents[str(i)] = []

		for i,sam_file in enumerate(self.int_sam):
			start_time = time()
			treatment = self.int_prefix[i]
			logging.info("Reading file = " + sam_file +".")
			with open(sam_file) as f:
				for line in f:
					if line[0] == "@":
						pass
					else:
						sam_entry = line.split()
						code = sam_entry[1]
						pos = int ( sam_entry[3] )

						if code == "4": # unmapped
							pass
						elif code == "0" or code == "16":
							ref_name = sam_entry[2].split('|')[3][:-2]
							strand = '+'
							if code == "16":
								strand = '-'

							#Lets see if this position already exists.
							try:
								index = sam_file_contents[ref_name].index(pos)
								sam_file_contents[ref_name][index].increment_hop_count(strand,i)
								print "It already exists"
							except:
								new_hop = HopSite(pos, strand, self.num_conditions)
								new_hop.increment_hop_count(i)
								sam_file_contents[ref_name].append(new_hop)
								#print "Creating new hop,",new_hop
			subprocess.check_output(["gzip", "-f", sam_file])
		for ref in sam_file_contents:
			sam_file_contents[ref] = sorted(sam_file_contents[ref])
		self.tabulate_gene_hits(sam_file_contents)

	def get_hops_in_gene(self, it, sam_contents_from_ref, gene):
		#print gene
		beg = gene.start_trunc
		end = gene.end_trunc
		gather_hops = []
		while sam_contents_from_ref[it].position < beg:
			it += 1
			if it >= len(sam_contents_from_ref):
				return -1
		while beg <= sam_contents_from_ref[it].position <= end:
			gather_hops.append(sam_contents_from_ref[it])
			it += 1
			if it >= len(sam_contents_from_ref):
				return -1

		#print "cur it = ",it
		#if it != -1:
		#	print "cur pos = ",sam_contents_from_ref[it].position
		gene.hop_list = gather_hops
		#if len(gene.hop_list) != 0:
		#	print "Not empty:",gene.hop_list
		#else:
		#	print "No hops in this gene."

		it = it - 50
		if it < 0:
			it = 0
		return it

	def tabulate_gene_hits(self, sam_file_contents):
		logging.info("         ---------------\n")
		logging.info("Begin tabulating gene hits...\n")
		start_time = time()
		for chrom in self.gene_info:
			ref = str(chrom)
			it = 0
			logging.info("Working on reference = " + ref)
			for gene in chrom.gene_list:
				# We are on the first gene of the first reference.
				# Lets get all the hops within the boundaries of this gene.
				it = self.get_hops_in_gene(it, sam_file_contents[ref], gene)
				if it == -1: # We are done with these hops.
					break

		self.get_normalized_coefficients()
		self.print_time_output("Done tabulating gene hits,",start_time)	
	
	def get_normalized_coefficients(self):
		logging.info("\nBegin Normalization Steps.\n")
		# Lets collect the total hits in intergenic regions to see what is going on.
		intergenic_totals = [0] * self.num_conditions		
		for chrom in self.gene_info:
			for gene in chrom.gene_list:
				if gene.is_intergenic:
					intergenic_totals = [x + y for x, y in zip(intergenic_totals, gene.hop_totals())]


		logging.info("Total intergenic hops observed is: " + str(sum(intergenic_totals)))
		for i,total in enumerate(intergenic_totals):
			logging.info(self.int_prefix[i] + " has " + str(total) + " intergenic hops observed.")

		minimum = min(intergenic_totals)

		self.debugger("min = ",minimum)
		if minimum <= 0:
			logging.error("Normalization couldn't be completed. It appears a condition has no hop hits.")
			sys.exit('Exiting')
		for i,totals in enumerate(intergenic_totals):
			self.normalization_coefficients.append(float(minimum)/float(totals))
		logging.info('Normalization coefficients used:')
		for i,condition in enumerate(self.int_prefix):
			logging.info(condition + " multiplied by " + str(self.normalization_coefficients[i]))

	def write_output(self):
		logging.info("Begin calculating gene totals and writing to output...")
		start_time = time()
		
		with open(self.tabulated_filename, 'w') as hf, open(self.gene_tabulated_filename, 'w') as gf, open(os.path.splitext(self.tabulated_filename)[0]+"-intergenic.txt", 'w') as intf:
			hops_header = ["Num","GeneID"]
			hops_header.extend(self.int_prefix)
			hops_header.extend([s + "(Normalized)" for s in self.int_prefix])
			hops_header.extend(["Start","Stop","Order","Strand","Length","PID","Gene","Function"])
			hf.write("\t".join(hops_header)+"\n")
			gf.write("\t".join(hops_header)+"\n")
			intf.write("\t".join(hops_header)+"\n")
			count = 1

			for chrom in self.gene_info:
				for gene in chrom.gene_list:
					if not gene.is_intergenic:
						hf.write(gene.write_hops(count,self.normalization_coefficients)+"\n")
						gf.write(gene.write_gene()+"\n")
						count += 1
					else: # We can write the intergenic stuff to a file as well for fun.
						intf.write(gene.write_hops(count,self.normalization_coefficients)+"\n")

		self.print_time_output("Done calculating totals and writing to output,", start_time)



################################
############# Main #############
################################

def main():

	#our main object
	hp = hops_pipeline()
	config = ""
	try:
		config = sys.argv[1]
	except:
		sys.exit("\nUSAGE: python TnSeq-Pipeline.py example.config\n")

	hp.read_config(config)

	hp.run_pipeline()

if __name__ == "__main__":
	main()
