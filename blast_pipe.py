from subprocess import Popen, PIPE
from os import remove

outcsv = "030_annotations.csv"
try:
	remove(outcsv)
except OSError:
	pass

infasta = "030_from_sam.fasta"
#infasta = "/data/GroupOne/superprotdb/query.fa"
bach_size = 1000
batches = -1 # if set to -1 will loop through the whole file

myout = ""
myerr = ""
myin = ""

seqcount = 0
batchcount = 0
with open(infasta,"r") as fhandle:
	for line in fhandle:
		if line.startswith(">"):
			seqcount += 1
		if seqcount > bach_size:
			subp = Popen(["blastx","-db","/data/GroupOne/superprotdb/super_family_prots.fasta", \
"-evalue","1e-6","-outfmt","10 qseqid sacc stitle bitscore evalue","-max_target_seqs", \
"1","-query","-"], stdout = PIPE, stderr = PIPE, stdin = PIPE)
			myout, myerr = subp.communicate(myin)
			with open(outcsv, "a") as ohandle:
				ohandle.write(myout)
			myin = ""
			seqcount = 0
			batchcount += 1
			if batches > 0 and batchcount > batches:
				break

		myin += line

if len(myin) > 0:
	subp = Popen(["blastx","-db","/data/GroupOne/superprotdb/super_family_prots.fasta", \
"-evalue","1e-6","-outfmt","10 qseqid sacc stitle bitscore evalue","-max_target_seqs", \
"1","-query","-"], stdout = PIPE, stderr = PIPE, stdin = PIPE)
	myout, myerr = subp.communicate(myin)
	with open(outcsv, "a") as ohandle:
		ohandle.write(myout)
	

exit()

