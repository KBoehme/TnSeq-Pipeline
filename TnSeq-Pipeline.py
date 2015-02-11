#Author: Kevin Boehme
#Date: 12-20-2014
#This program is part of a TnSeq analysis pipeline, designed to take raw fastq reads, and produce tabulated data on hop count occurance.

import sys
import ConfigParser

import logging
import gzip
import os
import collections
import glob
import subprocess

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
		self.ptt_info = {}
		self.gene_info = {}
		self.gene_totals = {}

		self.reference_names = []
		self.gene_keys = {}
		self.current_gene = {}
		self.current_gene_tuple = None

		self.total_counts = []
		self.normalization_coefficients = []

############ Read Config File and Run Pipeline ###############################
	def read_config(self, config_path):
		cp = ConfigParser.RawConfigParser()

 		cp.read(config_path)

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
 		self.transposon = cp.get('parameters', 'Transposon')
 		if self.transposon == None:
 			sys.exit('Error with Transposon parameter.')

 		try:
  			self.mismatches = int ( cp.get('parameters', 'Mismatches') )
  		except:
  			sys.exit('Error with Mismatches parameter (Not an integer).')

  		try:
 			self.gene_trim = int ( cp.get('parameters', 'GeneTrim') )
 		except:
  			sys.exit('Error with GeneTrim parameter (Not an integer).')

  		try:
 			self.read_length = int ( cp.get('parameters', 'ReadLength') )
 		except:
  			sys.exit('Error with ReadLength parameter (Not an integer).')




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
		self.total_counts = [0] * self.num_conditions
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
			logging.info("\n\n+============== Config File Settings ===============\n|\n|" + '| '.join(f.readlines()) + "|\n+-------------------------------------\n\n")

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
			
			bowtie_command = ["bowtie2", "-x", self.ref,"--phred33",
                                "-f",self.int_trimmed[out_file_num].name, "-D" , "25" , "-R" ,"3", "-L" , "10" , "-i" , "S,1,0.50" ,
                                "-S",self.int_sam[out_file_num], "--no-hd"]


			logging.info("Bowtie Command Used: " + ' '.join(bowtie_command)+"\n\n")


			logging.info("Writing output to = " + self.int_sam[out_file_num]+"\n")

			try:
				logging.info(subprocess.check_output(bowtie_command,stderr=subprocess.STDOUT))
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
		self.post_process_gene_info()
		logging.info("         ---------------\n")

		self.print_time_output("Done processing SAM files,", start_time)
		logging.info("--------------------------------------\n\n")

	def prepare_gene_info(self):
		start_time = time()
		for file_name in self.ptt:
			name = os.path.basename(file_name).split('.')[0]
			self.reference_names.append(name)
			with open(file_name, 'r') as f:
				
				gene_dict = {}
				gene_tuple_list = []
				

				# Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
				f.readline() # information title.
				f.readline() # number of proteins
				f.readline() # header
				for line in f:

					ptt_entry = line.split()
					location = ptt_entry[0].split('..')

					# Lets remove around the gene.
					beg = int ( location[0])
					end = int ( location[1])
					length = end - beg
					if self.gene_trim != 0:
						user_percentage = int ( length/self.gene_trim )
						if user_percentage > 0:
							beg += user_percentage
							end -= user_percentage
					
					sm = ptt_entry[5]
					#ptt_entry[0] = str(location[0]) + "(" + str(beg) + ")" +".."+str(location[1])+"("+str(end)+")"
					#ptt_entry[2] = str(ptt_entry[2]) + "(" + str(end-beg) + ")"
					self.ptt_info[sm] = ptt_entry # get a dictionary relating pid of the gene to its info.
					self.gene_totals[sm] = [0] * len(self.input_files)

					tup_key = (sm, beg , end )
					gene_key_item = (beg, end , sm )

					gene_tuple_list.append(gene_key_item)

					gene_dict[tup_key] = {}
					gene_dict[tup_key]['+'] = {}
					gene_dict[tup_key]['-'] = {}

			gene_key_dict = {}

			self.gene_keys[name] = gene_tuple_list
			self.gene_info[name] = collections.OrderedDict(sorted(gene_dict.items())) #key=lambda x: x[2]))

	def read_sam_file(self):
		logging.info("Begin reading SAM files into memory...\n")
		sam_file_contents = {}
		for i in self.reference_names:
			sam_file_contents[i] = {}
			
		for i,sam_file in enumerate(self.int_sam):
			start_time = time()
			treatment = self.int_prefix[i]
			logging.info("Reading file = " + sam_file +".")
			with open(sam_file) as f:
				num_lines = sum(1 for line in f)
				print num_lines,"lines."
				f.seek(0)
				ten_percent = num_lines/10
				current_percent = 1
				current_line = 0
				for line in f:
					current_line += 1
					if current_line == ten_percent*current_percent:
						print "Finished " + str(current_percent*10) +"%."
						current_percent += 1
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
							shortcut = sam_file_contents[ref_name]
							if (pos,strand) in shortcut:
								sam_file_contents[ref_name][(pos,strand)][i] += 1
							else:
								sam_file_contents[ref_name][(pos,strand)] = [0] * self.num_conditions
								sam_file_contents[ref_name][(pos,strand)][i] += 1
						else:
							pass
			self.print_time_output("Done reading SAM file into memory,",start_time)
			logging.info("Zipping up SAM file.\n")
			subprocess.check_output(["gzip", "-f", sam_file])
		for ref, value in sam_file_contents.iteritems():
			temp_dict = {}
			for i,item in enumerate(sorted(value.iteritems())):
				temp_dict[i] = item
			sam_file_contents[ref] = temp_dict	
		self.tabulate_gene_hits(sam_file_contents)

	def update_pos_info(self, it, sam_file_contents):
		if it > len(sam_file_contents)-1:
			return None, None
		pos_info = sam_file_contents[it]
		position = pos_info[0][0]
		return pos_info,position

	def tabulate_gene_hits(self, sam_file_contents):
		logging.info("         ---------------\n")
		logging.info("Begin tabulating gene hits...\n")
		start_time = time()
		for ref, pos in sam_file_contents.iteritems():
			it = 0
			logging.info("Working on reference = " + ref)
			for gene_tuple in self.gene_keys[ref]:
				beg = gene_tuple[0]
				end = gene_tuple[1]
				sm = gene_tuple[2]
				gene_tuple_key = (sm,beg,end)
				pos_info, position = self.update_pos_info(it, sam_file_contents[ref])
				while position < beg:
					it += 1
					pos_info , position = self.update_pos_info(it, sam_file_contents[ref])
					if pos_info == None:
						break
				if pos_info == None:
					break
				while beg<=position<=end:
					for i,item in enumerate(pos_info[1]):
						self.total_counts[i] += item
					it += 1
					self.gene_info[ref][gene_tuple_key][pos_info[0][1]][position] = pos_info[1]
					pos_info, position = self.update_pos_info(it,sam_file_contents[ref])
					if pos_info == None: #We ran out of hop hits, so we are done.
						break
				if pos_info == None:
					break
				else:
					it -= 50
					if it < 0:
						it = 0

		self.get_normalized_coefficients()
		self.print_time_output("Done tabulating gene hits,",start_time)	
	
	def get_normalized_coefficients(self):
		minimum = min(self.total_counts)
		print minimum
		for i,totals in enumerate(self.total_counts):
			self.normalization_coefficients.append(float(minimum)/float(totals))
		print self.normalization_coefficients

	def post_process_gene_info(self):
		
		logging.info("Begin calculating gene totals and writing to output...")
		start_time = time()
		
		with open(self.tabulated_filename, 'w') as hf, open(self.gene_tabulated_filename, 'w') as gf:
			gene_header = ["Num","GeneID"]
			hops_header = gene_header
			if self.normalize:
				for i in self.int_prefix:
					gene_header.append(i+"(NORMALIZED)")
			else:
				gene_header.extend(self.int_prefix)

			hops_header.extend(self.int_prefix)
			gene_header.extend(["Start(Truncated)","Stop(Truncated)","Strand","Length","PID","Gene","Function"])
			hops_header.extend(["Start(Truncated)","Stop(Truncated)","Strand","Length","PID","Gene","Function"])

			hf.write("\t".join(hops_header)+"\n")
			gf.write("\t".join(gene_header)+"\n")
			count = 1

			for ref_name, gene in self.gene_info.iteritems():
				for i, (gene_info, strand) in enumerate(gene.iteritems()):
					collect_to_print = []
					
					sm_key = gene_info[0] #sm
					gene_count = [0] * self.num_conditions
					for pos_neg, pos in strand.iteritems():
						for position, hops in pos.iteritems():
							hop_line = ["",sm_key]
				#			print "position",position
				#			print "hops",hops
							for i,hop_count in enumerate(hops):
								if hop_count == 1:
									#hops[i] = 0
									gene_count[i] += hops[i]
								else:
									gene_count[i] += hops[i]
								hop_line.append(hops[i])
							hop_line.extend([position,"",pos_neg])
							collect_to_print.append(hop_line)
					self.gene_totals[sm_key] = gene_count
					total_line = []
					ptt_entry = self.ptt_info[sm_key]
					loc = ptt_entry[0].split('..')
					total_line.append(count)
					total_line.append(sm_key)
					total_line.extend(self.gene_totals[sm_key])
					copy_total_line = total_line
					if self.normalize:
						for num, l in enumerate( self.gene_totals[sm_key] ):
							total_line.append(l*self.normalization_coefficients[num])
					else:
						total_line.extend(self.gene_totals[sm_key])
					copy_total_line.extend(self.gene_totals[sm_key])

					total_line.append(str(loc[0])+" ("+str(gene_info[1])+")")
					total_line.append(str(loc[1])+" ("+str(gene_info[2])+")")
					total_line.append(ptt_entry[1])
					total_line.append(gene_info[2] - gene_info[1])
					total_line.append(ptt_entry[3])
					total_line.append(ptt_entry[4])
					total_line.append(' '.join(ptt_entry[8:]))
					collect_to_print = sorted(collect_to_print, key=lambda x: x[2+self.num_conditions])
					copy_total_line.extend(total_line[2+self.num_conditions:])
					collect_to_print.insert(0,copy_total_line)
					gf.write('\t'.join(map(str,total_line))+"\n") #write gene only file
					for item in collect_to_print:
						hf.write('\t'.join(map(str,item))+ "\n")
					count += 1
		self.print_time_output("Done calculating totals and writing to output,", start_time)


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
