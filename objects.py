

# This chromosome object represents each ptt file read in.
class Chromsome(object):
	"""docstring for Chromsome"""
	def __init__(self, file_name):
		self.file_name = file_name
		self.gene_list = []
		self.start = 1
		self.end = -1
		self.num_proteins = -1

	def __str__(self):
		return self.file_name

	def __repr__(self):
		return self.file_name

	def __cmp__(self, other):
		if hasattr(other, 'file_name'):
			return self.file_name.__cmp__(other.file_name)

	def __eq__(self, other):
		return self.file_name == other

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
		self.is_intergenic = False

		#Unique to script info.
		self.start_trunc = -1
		self.end_trunc = -1
		self.order_code = ""

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
		self.length = int ( ptt_entry[2] )
		self.pid = ptt_entry[3]
		self.gene = ptt_entry[4]
		self.synonym = ptt_entry[5]
		self.code = ptt_entry[6]
		self.cog = ptt_entry[7]
		self.function = ''.join(ptt_entry[8:])
		self.num_conditions = num_conditions

		if trim != 0:
			user_percentage = int ( ( self.end - self.start )/trim )
			if user_percentage > 0:
				self.start_trunc = self.start + user_percentage
				self.end_trunc = self.end - user_percentage

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
		#	SM_b20002	92	79	101	72	78	57	5585		-	
		new_line = [""]
		raw_gene_totals = [0] * self.num_conditions
		normalized_gene_totals = [0] * self.num_conditions
		for hop in self.hop_list:
			new_line.append(self.synonym)
			if hop.hops: 
				new_line.extend(hop.hops) # raw hops
				normalized_hops = [int(round(a*b)) for a,b in zip(hop.hops,norm_coef)]
				new_line.extend(normalized_hops) #normalized hops
				new_line.append(hop.position)
				new_line.append("")
				new_line.append("")
				new_line.append(hop.strand)
				raw_gene_totals = [x + y for x, y in zip(raw_gene_totals, hop.hops)]
				normalized_gene_totals = [x + y for x, y in zip(normalized_gene_totals, normalized_hops)]
			else:
				logging.error("Seems a hop was created but has no hops inside it.")
				pass
			hop_entry.append(new_line)
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
		hop_entry.insert(0, self.gene_total_line)
		for i,entry in enumerate(hop_entry):
			hop_entry[i] = '\t'.join(map(str,entry))
		#print "writing:",'\n'.join(hop_entry)
		return '\n'.join(hop_entry)

	def write_gene(self):
		return '\t'.join(map(str,self.gene_total_line))

# This HopSite object represents a specific hop site
class HopSite(object):
	"""docstring for HopSite"""
	def __init__(self, position, strand, num_conditions):
		self.position = position
		self.hops = {}
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

	def increment_hop_count(self, condition):
		self.hops[condition] += 1