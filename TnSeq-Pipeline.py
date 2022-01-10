#Author: Kevin Boehme
#Date: 12-20-2014
#This program is part of a TnSeq analysis pipeline, designed to take raw fastq or fasta reads, and produce tabulated data on hop count occurrence.

import sys
sys.dont_write_bytecode = True

from configparser import RawConfigParser

import logging
import gzip
import os
import glob
import subprocess
import re
from objects import *

from time import time
from datetime import datetime

class hops_pipeline(object):
	"""docstring for hops_pipeline"""
	def __init__(self):

		########### From config ###########
		#Inputs
		self.input_files = []
		self.ref   = ""
		self.ptt   = []
		self.out   = ""

		#Parameters
		self.transposon = ""
		self.mismatches = 0
		self.minbaseoffset = 0
		self.maxbaseoffset = 5
		self.gene_trim = 0
		self.read_length = 0
		self.minimum_hop_count = 0

		#bools
		self.debug = False
		self.normalize = True
		self.delete_intermediate_files = True
		self.check_transposon = True
		self.reverse_complement_reads = False
		self.igv_normalize = False
		self.negateIGV = False
		########### End From Config ###########


		#Intermediate files for output
		self.int_prefix = []
		self.int_trimmed = [] # list of files for the intermediate trimmed read files.
		self.int_sam = [] # list of sam output files.
		self.tabulated_filename = ""
		self.gene_tabulated_filename = ""
		self.intergenic_filename = "" 
		self.igv_filenames = {} # key: ref_name | value: file

		self.num_conditions = 0

		#information variables
		self.removed_notn = 0
		self.removed_tooshort = 0
		self.kept_reads = 0
		self.original_read_count = 0
		self.starttime = time()

		# Variables for tabulating hop hits.
		self.chromosomes = {}
		self.sam_file_contents = {}
		self.normalization_coefficients = []	

	


############ Read Config File and Run Pipeline ###############################
	def read_config(self, config_path):
		cp = RawConfigParser()
		try:
			cp.read(config_path)
		except:
			sys.exit("Error reading config file.")
		# input paths
		self.input_files = glob.glob(cp.get('input', 'Reads'))
		self.ref = cp.get('input', 'BowtieReference')
		self.ptt = glob.glob(cp.get('input', 'Ptt'))
		self.out  = cp.get('input', 'Out')

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
			self.minbaseoffset = int(cp.get('parameters', 'MinBaseOffset'))
		except:
			sys.exit('Error with MinBaseOffset parameter')
		try:
			self.maxbaseoffset = int(cp.get('parameters', 'MaxBaseOffset'))
		except:
			sys.exit('Error with MaxBaseOffset parameter')


		try:
			self.transposon = cp.get('parameters', 'Transposon')
		except:
			sys.exit('Error with Transposon parameter.')

		if self.transposon == "":
			print("Transposon parameter was empty. This means no check will be made for a transposon sequence and all reads will move to the mapping stage.")
			self.check_transposon = False
		elif not (re.match('^[ACGTacgt]+$',self.transposon)):
			sys.exit('Error with Transposon parameter (Make sure it only contains [ATCG]).')
	
		try:
			self.mismatches = int ( cp.get('parameters', 'Mismatches') )
		except:
			sys.exit('Error with Mismatches parameter (Not an integer).')

		if self.mismatches < 0:
			sys.exit('Mismatches parameter is negative.')
		elif self.check_transposon and self.mismatches >= len(self.transposon):
			sys.exit("Mismatches parameter is same length or greater than transposon sequence (Note if you dont want to check for a transposon sequence, leave it blank in the config file.")

		try:
			self.gene_trim = int ( cp.get('parameters', 'GeneTrim') )
		except:
			sys.exit('Error with GeneTrim parameter (Not an integer).')


		if self.gene_trim < 0:
			sys.exit('GeneTrim parameter is negative.')
		elif self.gene_trim > 49:
			sys.exit('Error with GeneTrim parameter (Must be 49 or smaller).')


		try:
			self.read_length = int ( cp.get('parameters', 'ReadLength') )
		except:
			sys.exit('Error with ReadLength parameter (Not an integer).')

		if self.read_length < 0:
			sys.exit('ReadLength parameter is negative.')


		try:
			self.minimum_hop_count = int ( cp.get('parameters', 'MinimumHopCount'))
		except:
			sys.exit('Error with MinimumHopCount parameter (Not an integer).')

		if self.minimum_hop_count < 0:
			sys.exit('MinimumHopCount parameter is negative.')


		#Options
		try:
			self.debug = cp.getboolean('options', 'Debug')
		except:
			sys.exit('Error with Debug parameter (Not True/False).')

		try:
			self.normalize = cp.get('options','Normalize')
			if self.normalize != "Intergenic" and self.normalize != "Total":
				sys.exit('Normalize parameter not one of the following [Intergenic, Total].')

		except:
			sys.exit('Error with Normalize parameter.')

		try:
			self.delete_intermediate_files = cp.getboolean('options','DeleteIntermediateFiles')
		except:
			sys.exit('Error with DeleteIntermediateFiles parameter (Not True/False)')


		try:
			self.reverse_complement_reads = cp.getboolean('options','ReverseComplementReads')
		except:
			sys.exit('Error with ReverseComplementReads parameter')

		try:
			self.igv_normalize = cp.getboolean('options','IGVNormalize')
		except:
			sys.exit('Error with IGVNormalize parameter')

		try:
			self.negateIGV = cp.getboolean('options','IGVNegateNegStrand')
		except:
			sys.exit('Error with IGVNormalize parameter')


		# Generate other variables
		self.num_conditions = len(self.input_files)
		self.output_directory = os.path.dirname(config_path)
		if self.output_directory == "":
			self.output_directory = "./"
		else:
			self.output_directory += "/"
		self.tabulated_filename = self.output_directory + "output_files/" + self.out + "-HOPS.txt"
		self.gene_tabulated_filename = self.output_directory + "output_files/" + self.out + "-GENE.txt"
		self.intergenic_filename = self.output_directory + "output_files/" + self.out + "-INTERGENIC.txt"

		if not os.path.exists( self.output_directory + "output_files/"):
			subprocess.check_output(["mkdir", self.output_directory + "output_files"])

		self.set_up_logger(self.output_directory + "output_files/" + self.out + ".log")

		with open(config_path,'r') as f:
			logging.info( bcolors.HEADER + "\n\n=============== Config File Settings ===============\n\n"  
				+ bcolors.ENDC + ''.join(f.readlines()) +  bcolors.HEADER 
				+ "\n--------------------------------------\n\n" + bcolors.ENDC)

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
		self.process_reads()
		self.call_bowtie2()
		self.process_sam()
		self.print_time_output("Total run time,", self.starttime)



############ Useful Functions ###############################
	def set_up_logger(self, log_file_name):
		logging.basicConfig( filename=log_file_name , filemode='w', level=logging.INFO,format='%(message)s' )
		logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

		logging.info(bcolors.OKGREEN + "=========================================")
		logging.info("             Tn-Seq Pipeline")
		logging.info("The MIT License (MIT) Copyright (c) 2014 Kevin")
		logging.info("     Date: "+ str(datetime.now()))
		logging.info("           Author: Kevin Boehme")
		logging.info("      Email: kevinlboehme@gmail.com")
		logging.info("=========================================" + bcolors.ENDC)

	def debugger(self, *text):
		if self.debug:
			sys.stdout.write("[DEBUG] ")
			for i in text:
				sys.stdout.write(bcolors.OKBLUE + str(i) + " " + bcolors.ENDC)
			print

	def print_time_output(self, command, start):
		time_to_run = time()-start
		if time_to_run > 60:
			time_to_run = time_to_run/60
			logging.info(command + " ran in " + "%.2f" % (time_to_run) + " minutes." + "\n")
		else:
			logging.info(command + " ran in " + "%.2f" % (time_to_run) +" seconds." + "\n")
	
	def fuzzy_match_beginning(self, pattern, genome, mismatches):
		for i in range(self.minbaseoffset, self.maxbaseoffset+1):
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

	def reverse_complement(self, sequence):
		reverse_complement = []
		for index in range(len(sequence) - 1, -1 , -1):
				if(sequence[index] == "A"):
						reverse_complement.append("T")
				if(sequence[index] == "T"):
						reverse_complement.append("A")
				if(sequence[index] == "C"):
						reverse_complement.append("G")
				if(sequence[index] == "G"):
						reverse_complement.append("C")

		return ''.join(reverse_complement)



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
		logging.info(bcolors.WARNING + "         ---------------\n" + bcolors.ENDC)

	def process_reads(self):
		logging.info( bcolors.HEADER + "\n\n=============== Read pre-processing ===============\n" + bcolors.ENDC)
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
				logging.error("Didn't find the string fasta or fastq in input files. Please make sure your data has a fasta or fastq file extension (Could also be 2ped).")
				sys.exit('')

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

		for f in self.int_trimmed:
			f.close()
		logging.info( bcolors.HEADER + "--------------------------------------\n\n" + bcolors.ENDC)

	def read_fastq(self, f, out_file_num):
		t0 = time()
		while True:
			if self.original_read_count != 0 and self.original_read_count % 1000000 == 0:
				sys.stdout.write('\rProcessed ' + str(self.original_read_count) + ' reads.')
				sys.stdout.flush()
			#	logging.info('Processed ' + str(self.original_read_count) + ' reads.')
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
					sys.stdout.write('\rProcessed ' + str(self.original_read_count) + ' reads.')
				#	logging.info('Processed ' + str(self.original_read_count) + ' reads.')
				name = f.readline().strip()[1:] #Strip the @ sign infront of fastq names.
				seq = f.readline().strip().upper()
				if not name or not seq: break #We are done, lets break out of the loop.

				#we have all the contents of one read, now lets look for the transposon.
				self.process_read(name, seq, out_file_num)
				self.original_read_count += 1

	def read_fasta(self, f, out_file_num):
		t0 = time()
		while True:

			if self.original_read_count != 0 and self.original_read_count % 1000000 == 0:
				sys.stdout.write('\rProcessed ' + str(self.original_read_count) + ' reads.')
			#	logging.info('Processed ' + str(self.original_read_count) + ' reads.')
			name = f.readline().strip()[1:] #Strip the @ sign infront of fastq names.
			seq = f.readline().strip().upper()
			if not name or not seq: break #We are done, lets break out of the loop.

			self.process_read(name, seq, out_file_num)
			self.original_read_count += 1

	def process_read(self, name, seq, out_file_num):
		if self.reverse_complement_reads:
			seq = self.reverse_complement(seq)
		tn_trimmed_seq = ""
		if self.check_transposon:
			tn_trimmed_seq  = self.remove_transposon(name, seq)
		else:
			tn_trimmed_seq = seq

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
		if len(seq) < self.read_length:
			return None
		else:
			return seq[0:self.read_length]



############ Bowtie2 ###############################
	def call_bowtie2(self):

		logging.info( bcolors.HEADER + "\n\n=============== Bowtie2 Mapping ===============\n" + bcolors.ENDC)
		for out_file_num in range(self.num_conditions):
			start_time = time()
			
			#bowtie_command = ["bowtie2", "-x", self.ref,"--phred33",
			#                    "-f",self.int_trimmed[out_file_num].name, "-D" , "25" , "-R" ,"3", "-L" , "10" , "-i" , "S,1,0.50" ,
			#                    "-S",self.int_sam[out_file_num], "--no-hd"]

			bowtie_command = ["bowtie2", "-x", self.ref,"--phred33",
								"-f",self.int_trimmed[out_file_num].name,
								"-S",self.int_sam[out_file_num]]
								

			logging.info("Bowtie Command Used: " + ' '.join(bowtie_command)+"\n\n")


			logging.info("Writing output to = " + self.int_sam[out_file_num]+"\n")

			try:
				logging.info(subprocess.check_output(bowtie_command,stderr=subprocess.STDOUT))#,shell=True))
			except:
				logging.error("Bowtie2 doesn't seem to be installed. Make sure it can be run with the command: bowtie2")
				sys.exit('Exiting')

			self.print_time_output("Bowtie2",start_time)
			if self.delete_intermediate_files:		
				logging.info("Deleted trimmed Fasta file.")
				subprocess.check_output(["rm",self.int_trimmed[out_file_num].name])
			else:
				logging.info("Zipping up trimmed file.\n")
				subprocess.check_output(["gzip","-f", self.int_trimmed[out_file_num].name])
			logging.info(bcolors.WARNING + "         ---------------\n" + bcolors.ENDC)

		logging.info( bcolors.HEADER + "--------------------------------------\n\n" + bcolors.ENDC)
		return True



############ Tabulate Sams ###############################
	def process_sam(self):
		start_time = time()
		logging.info( bcolors.HEADER + "\n\n=============== Process SAM File ===============\n" + bcolors.ENDC)
		self.prepare_gene_info()
		self.add_intergenic_regions_and_order_column()
		self.read_sam_file()
		self.tabulate_gene_hits()
		self.get_normalized_coefficients()
		self.prepare_igv_files()
		self.write_output()
		logging.info(bcolors.WARNING + "         ---------------\n" + bcolors.ENDC)
		self.print_time_output("Done processing SAM files,", start_time)
		logging.info( bcolors.HEADER + "--------------------------------------\n\n" + bcolors.ENDC)

	def prepare_gene_info(self):
		self.debugger("On function: prepare_gene_info")
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
				self.chromosomes[name] = new_chrom
				our_genes = []
				for line in f:
					ptt_entry = line.split()

					new_gene = Gene()
					new_gene.create_from_ptt_entry(ptt_entry,self.gene_trim, self.num_conditions)

					our_genes.append(new_gene)
				self.chromosomes[name].set_gene_list(sorted(our_genes))

	def add_intergenic_regions_and_order_column(self):
		self.debugger("On function: add_intergenic_regions_and_order_column")
		for j, (ref_name, chrom) in enumerate(self.chromosomes.iteritems()):
			self.debugger("On replicon = " + ref_name)
			example_gene = chrom.gene_list[0]
			self.debugger("Example gene from this replicon = " + str(example_gene))
			gene_name = example_gene.synonym
			prefix = "".join(re.findall("[a-zA-Z]+", gene_name))
			self.debugger("Prefix = " + prefix)
			num_genes = len(chrom.gene_list)
			self.debugger("Number of genes on replicon = " + str(num_genes))
			zfill_digits = len(str(num_genes))
			self.debugger("Digits to use in zfill = " + str(zfill_digits))

			intergenic_genes = set()
			for i,cur_gene in enumerate(chrom.gene_list,start=1):
				ordered_name = prefix + str(i).zfill(zfill_digits)
				cur_gene.order_code = ordered_name

				# Now lets add intergenic regions.
				new_gene = Gene()
				if i == 1: # On the first gene.
					if cur_gene.start > 1: # We have some room to capture.
						new_gene.create_intergenic_region(1, cur_gene.start - 1, "int_BEG-"+cur_gene.synonym, self.num_conditions)
						intergenic_genes.add(new_gene)
					if cur_gene.end + 1 < chrom.gene_list[i].start: # We have some room to capture.
						new_gene.create_intergenic_region(cur_gene.end + 1, 
							chrom.gene_list[i].start - 1, 
							"int_" + cur_gene.synonym + "-" + chrom.gene_list[i].synonym, self.num_conditions)
						intergenic_genes.add(new_gene)
				elif i == chrom.num_proteins: # On the last gene.
					if cur_gene.end < chrom.end: # We have some room to capture.
						new_gene.create_intergenic_region(cur_gene.end + 1, chrom.end, "int_"+cur_gene.synonym+"-END", self.num_conditions)
						intergenic_genes.add(new_gene)
				else: # Make sure the end of the current gene and the beginning of the next gene have a space.
					if cur_gene.end + 1 < chrom.gene_list[i].start: # We have some room to capture.
						new_gene.create_intergenic_region(cur_gene.end + 1, 
							chrom.gene_list[i].start - 1, 
							"int_" + cur_gene.synonym + "-" + chrom.gene_list[i].synonym, self.num_conditions )
						intergenic_genes.add(new_gene)
			#Now lets smash those intergenic regions with the current genes.
			self.chromosomes[ref_name].gene_list = sorted(self.chromosomes[ref_name].gene_list + list(intergenic_genes))
		logging.info(bcolors.WARNING + "         ---------------\n" + bcolors.ENDC)

	def read_sam_file(self):
		self.debugger("On function: read_sam_file")
		for i in self.chromosomes:
			self.sam_file_contents[str(i)] = {}

		for i,sam_file in enumerate(self.int_sam):
			start_time = time()
			treatment = self.int_prefix[i]
			logging.info("Reading file = " + sam_file +".")
			with open(sam_file) as f:
				#num_lines = float(sum(1 for line in f))
				#f.seek(0)
				for j,line in enumerate(f):
					#self.update_progress(float(j)/num_lines)
					if line[0] == "@": #Pass the headers
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
							hop_exists = False
							if pos in self.sam_file_contents[ref_name]:
								self.sam_file_contents[ref_name][pos].increment_hop_count(i)
							else:
								new_hop = HopSite(pos, self.negateIGV, strand, self.num_conditions)
								new_hop.increment_hop_count(i)
								self.sam_file_contents[ref_name][pos] = new_hop
				logging.info("")
				self.print_time_output("Reading file", start_time)
			
			if self.delete_intermediate_files:
				logging.info("Deleting SAM file.")
				subprocess.check_output(["rm",sam_file])
			else:
				logging.info("Zipping up SAM file.")
				subprocess.check_output(["gzip", "-f", sam_file])

		self.filter_on_min_hops()

		for ref, value in self.sam_file_contents.iteritems():
			temp_dict = {}
			for i,item in enumerate(sorted(value.iteritems())):
				temp_dict[i] = item[1]
			self.sam_file_contents[ref] = temp_dict	
		logging.info(bcolors.WARNING + "         ---------------\n" + bcolors.ENDC)

	def filter_on_min_hops(self):
		for ref, positions in self.sam_file_contents.iteritems():
			for index,hop in positions.items():
				if hop.total_hops() < self.minimum_hop_count:
					del positions[index]
		
	def tabulate_gene_hits(self):
		self.debugger("On function: tabulate_gene_hits")
		logging.info("Begin tabulating gene hits...\n")
		for ref, pos in self.sam_file_contents.iteritems():
			if len(pos) == 0:
				logging.info("Reference " + ref + " has no reads which mapped to it.")
			else:
				start_time = time()
				it = 0
				logging.info("Working on reference = " + ref)
				chrom = self.chromosomes[ref]
				num_genes = float(len(chrom.gene_list))
				for i,gene in enumerate(chrom.gene_list,start=1):
					self.update_progress(float(i)/num_genes)
					beg = gene.start_trunc
					end = gene.end_trunc

					curr_hop = self.update_current_hop(it, self.sam_file_contents[ref])
					while curr_hop.position < beg:
						it += 1
						curr_hop = self.update_current_hop(it, self.sam_file_contents[ref])
						if not curr_hop:
							break
					if not curr_hop:
						break
					while beg <= curr_hop.position <= end:
						gene.hop_list.append(curr_hop)
						it += 1
						curr_hop = self.update_current_hop(it, self.sam_file_contents[ref])
						if not curr_hop: # We ran out of hop hits, so we are done.
							break

					if not curr_hop:
						break
					else:
						it -= 50
						if it < 0:
							it = 0
				self.update_progress(1)
				self.print_time_output(" Done tabulating gene hits,",start_time)
		logging.info(bcolors.WARNING + "         ---------------\n" + bcolors.ENDC)

	def get_normalized_coefficients(self):
		self.debugger("On function: get_normalized_coefficients")
		logging.info("\nBegin Normalization Steps.\n")

		total_counted_hops = [0] * self.num_conditions	
		if self.normalize == "Total":
			for ref_name,chrom in self.chromosomes.iteritems():
				for gene in chrom.gene_list:
					total_counted_hops = [x + y for x, y in zip(total_counted_hops, gene.hop_totals())]
			logging.info("Total hops observed is: " + str(sum(total_counted_hops)))
			for i,total in enumerate(total_counted_hops):
				logging.info(self.int_prefix[i] + " has " + str(total) + " [total] hops observed.")

		elif self.normalize == "Intergenic":
			for ref_name,chrom in self.chromosomes.iteritems():
				for gene in chrom.gene_list:
					if gene.is_intergenic:
						total_counted_hops = [x + y for x, y in zip(total_counted_hops, gene.hop_totals())]
			logging.info("Total [intergenic] hops observed is: " + str(sum(total_counted_hops)))
			for i,total in enumerate(total_counted_hops):
				logging.info(self.int_prefix[i] + " has " + str(total) + " [intergenic] hops observed.")
		else:
			logging.error("Error with normalization.")
			sys.exit()



		minimum = min(total_counted_hops)

		self.debugger("min = ",minimum)
		if minimum <= 0:
			logging.error("Normalization couldn't be completed. It appears a condition has no hop hits.")
			self.normalization_coefficients = [1] * self.num_conditions
			return
			#sys.exit('Exiting')
		for i,totals in enumerate(total_counted_hops):
			self.normalization_coefficients.append(float(minimum)/float(totals))
		logging.info('Normalization coefficients used:')
		for i,condition in enumerate(self.int_prefix):
			logging.info(condition + " multiplied by " + str(self.normalization_coefficients[i]))

	def write_output(self):
		self.debugger("On function: write_output")
		logging.info("Begin calculating gene totals and writing to output...")
		start_time = time()

		if self.delete_intermediate_files:
			logging.info("Deleting Intermediate Folder.")
			subprocess.check_output(["rmdir",self.output_directory + "intermediate_files/"])

		with open(self.tabulated_filename, 'w') as hf, open(self.gene_tabulated_filename, 'w') as gf, open(self.intergenic_filename, 'w') as intf:
			hops_header = ["Num","GeneID"]
			hops_header.extend(self.int_prefix)
			hops_header.extend([s + "(Normalized)" for s in self.int_prefix])
			hops_header.extend(["Start","Stop","Order","Strand","Length","PID","Gene","Function"])
			hf.write("\t".join(hops_header)+"\n")
			gf.write("\t".join(hops_header)+"\n")
			intf.write("\t".join(hops_header)+"\n")
			count = 1

			for ref_name, chrom in self.chromosomes.iteritems():
				for gene in chrom.gene_list:
					if len(gene.hop_list) > 0: # If the gene even has hops in it to write.
						self.igv_filenames[ref_name].write(gene.write_igv(ref_name, self.igv_normalize, self.normalization_coefficients) + "\n")
					else:
						pass

					if not gene.is_intergenic:
						hf.write(gene.write_hops(count,self.normalization_coefficients)+"\n")
						gf.write(gene.write_gene()+"\n")
						count += 1
					else: # We can write the intergenic stuff to a file as well for fun.
						intf.write(gene.write_hops(count,self.normalization_coefficients)+"\n")

		self.print_time_output("Done calculating totals and writing to output,", start_time)

	def update_progress(self,progress):
		progress = int(round(progress * 100.0))
		sys.stdout.write('\r[{0}] {1}%'.format('#'*(progress/5), progress))
		sys.stdout.flush()
		if progress == 100:
			sys.stdout.write(' ')

	def update_current_hop(self, it, sam_file_contents):
		if it > len(sam_file_contents)-1:
			return None
		hop = sam_file_contents[it]
		return hop

	def prepare_igv_files(self):
		igv_path = self.output_directory + "output_files/IGV/"
		if not os.path.exists( self.output_directory + "output_files/IGV/"):
			subprocess.check_output(["mkdir", igv_path])
		for ref_name in self.chromosomes.keys():
			filename = igv_path + ref_name + ".igv"
			self.igv_filenames[ref_name] = open(filename, 'w+')

			#Now lets give them a header.
			self.igv_filenames[ref_name].write("#Transpon"+ "\n")
			header = ["Chromosome",	"Start", "End", "Feature"]
			for i in self.int_prefix:
				header.append(i)
			self.igv_filenames[ref_name].write('\t'.join(header) + "\n")


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
		sys.exit("\nUSAGE: python TnSeq-Pipeline.py pathtoconfig.config\n")

	hp.read_config(config)

	hp.run_pipeline()

if __name__ == "__main__":
	main()
