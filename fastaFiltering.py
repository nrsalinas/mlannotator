import os

infile = 'GCA_000364345.1_Macaca_fascicularis_5.0_genomic.fna'
outfile = 'nice_chromosomes.fasta'

try:
	os.remove(outfile)
except OSError:
	pass

bffr = ''

good_seq = 0

with open(infile,'r') as file_handle:
	for line in file_handle:
		if len(bffr) > 3000:
			with open(outfile,'a') as out_handle:
				out_handle.write(bffr)
			bffr = ''
		if line.startswith('>CM'):
			bffr += line
			good_seq = 1
			continue
		elif line.startswith('>'):
			good_seq = 0
			continue
		elif good_seq:
			bffr += line
			continue

with open(outfile,'a') as out_handle:
	out_handle.write(bffr)
bffr = ''
exit()
