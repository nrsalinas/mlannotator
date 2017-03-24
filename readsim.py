import sqlite3
import tablef
import vectorf
import random


bffr = ''

db_file = 'superprot.sqlite'
csv_file = 'supertable.csv'

table = tablef.parse_supertable(csv_file)
reads_sf = tablef.sf_ratios(table, superfamilies = 5, max_reads = 20)

con = sqlite3.connect(db_file)
cursor = con.cursor()

for sfid,num in reads_sf:
	cursor.execute('''SELECT Protein FROM Assignments LEFT JOIN Superfamilies ON SuperID = Super WHERE SuperID = ?''' , (sfid,))
	protein_ids = cursor.fetchall()
	random.shuffle(protein_ids)
	selected_proteins = protein_ids[:num]

	cursor.executemany('''SELECT Sequence, Start, End FROM Proteins WHERE ProteinID = ?''' , selected_proteins)
	dat = cursor.fetchall()
	i = dat[1] + 1
	f = dat[2] + 1
	protseq = dat[0][i,f]

	vectorf.kmerme(protseq)
