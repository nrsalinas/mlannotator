import sqlite3
import tablef
import vectorf
from scipy.sparse import dok_matrix
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.model_selection import train_test_split

db_file = 'superprot.sqlite'
csv_file = 'supertable.csv'

def get_arrays(kmer_size = 2, classes = 30, sim_seq_length = 33, reads_per_domain = 2):
	"""
	classes (int): Number of the most common superfamilies of proteins to be
	selected for training.

	reads_per_domain (int): Reads to be simulated per domain assignment.

	sim_seq_length(int): Length of simulated aminoacid sequences.
	"""

	bffr = "Kmer size: {0}\nSuperfamilies included: {1}\nShort reads simulated per superfamily: {2}\n\n".format(kmer_size, classes, reads_per_domain)

	table = tablef.parse_supertable(csv_file)
	table = table[:classes]

	no_train_samples = classes * table[0][1] * reads_per_domain

	con = sqlite3.connect(db_file)
	cursor = con.cursor()

	samples_indx = 0
	spa_mat = dok_matrix((no_train_samples , (20 ** kmer_size + 1)))

	bffr += "Superprotein\tSeqs attempted\tSeqs retrieved\tMean seq length\n"
	for sfid,num in table:
		bffr += "{0}\t\t{1}\t\t".format(sfid, num)
		seqs_retrieved = num
		mean_length = 0
		cursor.execute("SELECT Protein, Start, End FROM Assignments LEFT JOIN Superfamilies ON SuperID = Super WHERE SuperID = ? and End - Start > ?" , (sfid, sim_seq_length))
		myassignments = cursor.fetchall()
		myassignments.sort(key = lambda x: x[0])

		curr_prot = 0
		dat = []
		curr_seq = []
		for (prot, start, end) in myassignments:
			mean_length += (end - start)
			if prot != curr_prot:
				curr_prot = prot
				cursor.execute("SELECT Sequence FROM Proteins WHERE ProteinID = ?", (prot,) )
				curr_seq = cursor.fetchall()
				curr_seq = curr_seq[0][0]

			protseq = curr_seq[start:end]
			if "X" in protseq:
				seqs_retrieved -= (1 * reads_per_domain)
				continue
			else:
				#print protseq
				for x in xrange(reads_per_domain):
					kmer_freqs = vectorf.kmerme(protseq, kmer_length = kmer_size, peptide_length = sim_seq_length)

					# Store kmer freqs in Scipy sparse matrix
					for kmer_indx in kmer_freqs:
						#print "samples_indx:", samples_indx, "kmer_indx:", kmer_indx
						spa_mat[samples_indx, kmer_indx] = kmer_freqs[kmer_indx]
					spa_mat[samples_indx, (20 ** kmer_size)] = sfid
					samples_indx += 1
		mean_length /= seqs_retrieved
		bffr += "{0}\t\t{1}\n".format(seqs_retrieved, mean_length)
	con.close()
	print bffr

	# Reshape the matrix
	max_row = (max(spa_mat.keys())[0] + 1)
	_,cols = spa_mat.shape
	#print "len(spa_mat.keys())", len(spa_mat.keys())
	spa_mat.resize((max_row, cols))
	#print "len(spa_mat.keys())", len(spa_mat.keys())
	X = spa_mat.tocsc()[:,:-1]
	Y = spa_mat.tocsc()[:,-1]
	Y = Y.toarray().flatten()

	tfidrer = TfidfTransformer()
	Xtrans = tfidrer.fit_transform(X)
	#X_train, X_test, Y_train, Y_test = train_test_split(Xtrans, Y, test_size = 0.1, random_state = 0)
	return (Xtrans, Y)
