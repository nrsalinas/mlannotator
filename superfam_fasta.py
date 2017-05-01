import sqlite3
import re

outfile = "super_family_prots.fasta"
bfrout = ""
db_file = 'superprot.sqlite'

con = sqlite3.connect(db_file)
cursor = con.cursor()
cursor.execute("SELECT SuperID, SuperName FROM Superfamilies")
supfam = cursor.fetchall()

for sup in supfam:
	title = sup[1]
	#title = re.sub("\s+","_",sup[1])
	#title = re.sub(",","_",title)
	#print sup[0], sup[1]

	cursor.execute("SELECT Protein, Start, End FROM Assignments LEFT JOIN Superfamilies ON SuperID = Super WHERE SuperID = ? " , (sup[0],))
	myassignments = cursor.fetchall()
	myassignments.sort(key = lambda x: x[0])
	#print myassignments[:20]

	curr_prot = 0
	dat = []
	curr_seq = []
	for (prot, start, end) in myassignments:
		if prot != curr_prot:
			curr_prot = prot
			cursor.execute("SELECT Sequence FROM Proteins WHERE ProteinID = ?", (prot,) )
			curr_seq = cursor.fetchall()
			curr_seq = curr_seq[0][0]

		protseq = curr_seq[start:end]
		bfrout += ">Domain:{0}_{1}-{2} SuperfamilyID:{3} SuperfamilyDescr:{4}\n".format(prot, start, end, sup[0], title)
		bfrout += "{0}\n".format(protseq)

		if len(bfrout) > 10000:
			with open(outfile, "a") as fhandle:
				fhandle.write(bfrout)
			bfrout = ""

with open(outfile, "a") as fhandle:
	fhandle.write(bfrout)

con.close()
