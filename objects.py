

#This bcolors stuff will make the output nice looking. For now I comment out the fun stuff because not all terminals are compatible with it. 
#If your feeling adventurus you might try to comment out the bottom one and uncomment out the top one...
'''
class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'
'''
class bcolors:
	HEADER = ''
	OKBLUE = ''
	OKGREEN = ''
	WARNING = ''
	FAIL = ''
	ENDC = ''
	BOLD = ''
	UNDERLINE = ''

# This chromosome object represents each ptt file read in.
class Chromsome(object):
	"""docstring for Chromsome"""
	def __init__(self, name):
		self.name = name
		self.gene_list = []
		self.start = 1
		self.end = -1
		self.num_proteins = -1

	def __str__(self):
		return self.name

	def __repr__(self):
		return self.name

	def __cmp__(self, other):
		if hasattr(other, 'file_name'):
			return self.name.__cmp__(other.name)

	def __eq__(self, other):
		return self.name == other

	def set_gene_list(self, gene_list):
		self.gene_list = gene_list

	def fill_info(self, title, num_prot):
		self.end = int ( title.split('..')[1] )
		self.num_proteins = int ( num_prot.split()[0] )


# This Gene class will represent an entry from the ptt file.
class Gene(object):
	"""docstring for Gene"""
	def __init__(self):
		#Ptt info
		self.start = -1
		self.end = -1
		self.strand = None
		self.length = -1
		self.pid = ""
		self.gene = ""
		self.synonym = ""
		self.code = ""
		self.cog = ""
		self.product = ""

		#Unique to script info.
		self.start_trunc = -1
		self.end_trunc = -1
		self.order_code = ""
		self.is_intergenic = False

		self.num_conditions = -1
		self.hop_list = [] #This will be a list of HopSite objects.
		self.gene_total_line = []

	def __repr__(self):
		return '{}: {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(
			self.__class__.__name__,
			self.start,
			self.end,
			self.strand ,
			self.length,
			self.pid ,
			self.gene ,
			self.synonym ,
			self.code ,
			self.cog ,
			self.product,
			self.start_trunc,
			self.end_trunc,
			self.order_code,
			self.hop_list
			)

	def __cmp__(self, other):
		if hasattr(other, 'start'):
			return self.start.__cmp__(other.start)
			
	def __eq__(self, other):
		return self.synonym == other

	def create_from_ptt_entry(self, ptt_entry, trim, num_conditions):
		location = ptt_entry[0].split('..')
		self.start = int( location[0] )
		self.end = int ( location[1] )
		self.strand = ptt_entry[1]
		self.length = self.end - self.start + 1
		#self.length = int ( ptt_entry[2] ) # The ptt file lengths were a little wonky.
		self.pid = ptt_entry[3]
		self.gene = ptt_entry[4]
		self.synonym = ptt_entry[5]
		self.code = ptt_entry[6]
		self.cog = ptt_entry[7]
		self.function = ' '.join(ptt_entry[8:])
		self.num_conditions = num_conditions
		self.start_trunc = self.start
		self.end_trunc = self.end

		try:
			if trim != 0:
				trim_percent = float(trim)/100.0
				gene_length = self.end - self.start
				user_percentage = int ( gene_length * trim_percent )
				if user_percentage > 0:
					self.start_trunc = self.start + user_percentage
					self.end_trunc = self.end - user_percentage
		except:
			print "Seems like the trim length is giving me some problems. Setting it to 0 and continuing."
			
	def create_intergenic_region(self, start, end, synonym, num_conditions):
		#print "Creating intergenic region"
		#print "start:",start
		#print "end:",end
		#print "synonym:",synonym

		self.start = start
		self.end = end
		self.strand = "na"
		self.length = self.end - self.start + 1
		self.pid = "na"
		self.gene = "na"
		self.synonym = synonym
		self.code = "na"
		self.cog = "na"
		self.function = "na"
		self.is_intergenic = True
		self.num_conditions = num_conditions

		self.start_trunc = start
		self.end_trunc = end
		self.order_code = "na"

	def hop_totals(self):
		total = [0] * self.num_conditions
		for hop in self.hop_list:
			total = [x + y for x, y in zip(total, hop.hops)]
		return total

	def write_hops(self, i, norm_coef):
		hop_entry = []
		raw_gene_totals = [0] * self.num_conditions
		normalized_gene_totals = [0] * self.num_conditions
		for hop in self.hop_list:
			new_line, normalized_hops = hop.print_hop_line( self.synonym, norm_coef)

			raw_gene_totals = [x + y for x, y in zip(raw_gene_totals, hop.hops)]
			normalized_gene_totals = [x + y for x, y in zip(normalized_gene_totals, normalized_hops)]
			hop_entry.append(new_line)

		self.generate_gene_total_line(i, raw_gene_totals, normalized_gene_totals)

		hop_entry.insert(0, self.gene_total_line)
		for i,entry in enumerate(hop_entry):
			hop_entry[i] = '\t'.join(map(str,entry))
		#print "writing:",'\n'.join(hop_entry)
		return '\n'.join(hop_entry)

	def generate_gene_total_line(self, i, raw_gene_totals, normalized_gene_totals):
		self.gene_total_line.append(i)
		self.gene_total_line.append(self.synonym)
		self.gene_total_line.extend(raw_gene_totals)
		self.gene_total_line.extend(normalized_gene_totals)
		self.gene_total_line.append(self.start)
		self.gene_total_line.append(self.end)
		self.gene_total_line.append(self.order_code)
		self.gene_total_line.append(self.strand)
		self.gene_total_line.append(self.length)
		self.gene_total_line.append(self.pid)
		self.gene_total_line.append(self.gene)
		self.gene_total_line.append(self.function)

	def write_gene(self):
		return '\t'.join(map(str,self.gene_total_line))

	def write_igv(self, ref_name, igv_normalize, norm_coef):
		#Chromosome      Start   End     Feature condition1      condition2      condition3
		igv_lines = []
		for hop in self.hop_list:
			igv_lines.append(hop.get_hops_for_igv(self.synonym, ref_name, igv_normalize, norm_coef))
		return '\n'.join(igv_lines)


# This HopSite object represents a specific hop site
class HopSite(object):
	"""docstring for HopSite"""
	def __init__(self, position, strand, num_conditions):
		self.position = position
		#self.hops = {}
		self.strand = strand
		self.hops = [0] * num_conditions

	def __str__(self):
		return '{}: {} {}\n'.format(
			self.__class__.__name__,
			self.position,
			self.strand,
			self.hops
			)

	def __repr__(self):
		return '{}: {} {}\n'.format(
			self.__class__.__name__,
			self.position,
			self.strand,
			self.hops
			)

	def __cmp__(self, other):
		if hasattr(other, 'position'):
			return self.position.__cmp__(other.position)

	def __eq__(self, other):
		return self.position == position

	def total_hops(self):
		return sum(self.hops)

	def increment_hop_count(self, condition):
		self.hops[condition] += 1

	def print_hop_line(self, syn, norm_coef):
		new_line = [""]
		new_line.append(syn)
		new_line.extend(self.hops) # raw hops
		normalized_hops = [int(round(a*b)) for a,b in zip(self.hops,norm_coef)]
		new_line.extend(normalized_hops) #normalized hops
		new_line.append(self.position)
		new_line.append("")
		new_line.append("")
		new_line.append(self.strand)
		return new_line, normalized_hops

	def get_hops_for_igv(self, feature, ref_name, igv_normalize, norm_coef):
		igv_line = []
		igv_line.append(ref_name)
		igv_line.append(str(self.position))
		igv_line.append(str( self.position + 1))
		igv_line.append(feature)
		if igv_normalize:
			igv_line.extend(str(int(round(a*b))) for a,b in zip(self.hops,norm_coef))
		else:
			igv_line.extend(list(str(i) for i in self.hops))
		return '\t'.join(igv_line)

