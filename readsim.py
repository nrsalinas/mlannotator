import sqlite3
import tablef
import vectorf
import random


bffr = ''

db_file = 'superprot.sqlite'
csv_file = 'supertable.csv'
peptide_length = 33

table = tablef.parse_supertable(csv_file)
reads_sf = tablef.sf_ratio(table, superfamilies = 5, max_reads = 50)

con = sqlite3.connect(db_file)
cursor = con.cursor()

for sfid,num in reads_sf:
	cursor.execute("SELECT Protein FROM Assignments LEFT JOIN Superfamilies ON SuperID = Super WHERE SuperID = ?" , (sfid,))
	protein_ids = cursor.fetchall()
	random.shuffle(protein_ids)
	selected_proteins = [x[0] for x in protein_ids[:num]]
	#print selected_proteins
	cursor.execute("SELECT Sequence, Start, End FROM Proteins LEFT JOIN Assignments ON ProteinID = Protein WHERE ProteinID IN ({0}) AND End - Start > {1}".format(', '.join('?' for _ in selected_proteins), peptide_length), selected_proteins)
	dat = cursor.fetchall()

	for line in dat:
		i = line[1] + 1
		f = line[2] + 1
		protseq = line[0][i:f]
		if "X" in protseq:
			continue
		else:
			print protseq
			kmer_freqs = vectorf.kmerme(protseq)

		# Store kmer freqs in Scipy sparse matrix
		# Parse matrix to scikit-learn

con.close()
