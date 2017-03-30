import sqlite3
import tablef
import vectorf
import random
from scipy.sparse import dok_matrix
from sklearn.naive_bayes import MultinomialNB

bffr = ''
db_file = 'superprot.sqlite'
csv_file = 'supertable.csv'
sim_seq_length = 33
kmer_size = 2
classes = 5 # Number of most common superfamilies of proteins to select for training
most_repres = 50 # Maximum number of simulated reads to generate from the most common superfamily

table = tablef.parse_supertable(csv_file)
reads_sf = tablef.sf_ratio(table, superfamilies = classes, max_reads = most_repres)
no_train_samples = classes * most_repres

con = sqlite3.connect(db_file)
cursor = con.cursor()

samples_indx = 0
spa_mat = dok_matrix((no_train_samples, (20 ** kmer_size + 1)))

for sfid,num in reads_sf:
	print "sfid:",sfid," - num:",num
	cursor.execute("SELECT Protein, Start, End FROM Assignments LEFT JOIN Superfamilies ON SuperID = Super WHERE SuperID = ? and End - Start > ?" , (sfid, sim_seq_length))
	myassignments = cursor.fetchall()
	random.shuffle(myassignments)
	myassignments = myassignments[:num]
	myassignments.sort(key = lambda x: x[0])

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
		if "X" in protseq:
			continue
		else:
			#print protseq
			kmer_freqs = vectorf.kmerme(protseq, kmer_length = kmer_size, peptide_length = sim_seq_length)

		# Store kmer freqs in Scipy sparse matrix

		for kmer_indx in kmer_freqs:
			#print "samples_indx:", samples_indx, "kmer_indx:", kmer_indx
			spa_mat[samples_indx, kmer_indx] = kmer_freqs[kmer_indx]
		spa_mat[samples_indx, (20 ** kmer_size)] = sfid
		samples_indx += 1

con.close()

# Reshape the matrix
max_row = (max(spa_mat.keys())[0] + 1)
_,cols = spa_mat.shape
#print "len(spa_mat.keys())", len(spa_mat.keys())
spa_mat.resize((max_row, cols))
#print "len(spa_mat.keys())", len(spa_mat.keys())
X = spa_mat.tocsc()[:,:-1]
Y = spa_mat.tocsc()[:,-1]

# Parse matrix to scikit-learn

model = MultinomialNB()
model.fit(X,Y)
