import re
sam_file = "030.sam"
fasta_file = "030_from_sam.fasta"
bffr = ""

compl = {'A':'T', 'G':'C', 'C':'G', 'T':'A', 'Y':'R', 'R':'Y', 'W':'W', 'S':'S', \
'K':'M', 'M':'K', 'D':'H', 'V':'B', 'H':'D', 'B':'V', 'X':'X', 'N':'N', '-':'-'}

with open(sam_file, "r") as fhandle:
	for line in fhandle:
		if not line.startswith("@"):
			bits = re.split("\s+", line)
			if bits[1] == '0':
				bffr += ">{0}\n{1}\n".format(bits[0], bits[9])
			elif bits[1] == '16':
				reverse = bits[9][::-1]
				co = "".join([compl[x] for x in reverse.upper()])
				bffr += ">{0}\n{1}\n".format(bits[0], co)
		if len(bffr) > 10000:
			with open(fasta_file, "a") as ihandle:
				ihandle.write(bffr)
			bffr = ""
with open(fasta_file, "a") as ihandle:
	ihandle.write(bffr)

exit()
