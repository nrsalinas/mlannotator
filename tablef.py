from numpy import polyfit

def parse_supertable(table_file):
	"""
	Retrieves domains and SCOP (Super ID) from superfamily table. Columns should
	be in the following order: Domains, Proteins, Architectures, Families,
	SCOP-Superfamily_ID, Superfamily_name.
	"""
	lineCount = 1
	table = []

	with open(table_file, 'r') as fhandle:
		for line in fhandle:
			if lineCount > 1:
				bits = line.split(',')
				table.append((int(bits[4]),int(bits[0])))
			lineCount += 1

	table.sort(key = lambda x : x[1], reverse = True)

	return table

def sf_ratio(supertable, superfamilies = 40, max_reads = 20, min_reads_prop = 0.5):
	"""
	Returns vector of short read counts per superfamily to simulate.
	Input table should be sorted by number of domains.
	"""
	ratios = []
	total_domains = float(sum([x[1] for x in supertable]))
	regr_func = polyfit([supertable[-1][1], supertable[0][1]], \
				[(max_reads * min_reads_prop), max_reads], 1)
	intercept = regr_func[1]
	coefficient = regr_func[0]

	for x in supertable[:superfamilies]:
		#reads = int(round((x[1] / total_domains) * const))
		reads = int(round(coefficient * x[1] + intercept)) 
		ratios.append((x[0], reads))

	return ratios
