from random import random
from numpy import array


aaindx = {"A":0, "C":1, "D":2, "E":3, "F":4, "G":5, "H":6, "I":7, "K":8, "L":9,
"M":10, "N":11, "P":12, "Q":13, "R":14, "S":15, "T":16, "V":17, "W":18, "Y":19}

def kmerme(aa_seq, kmer_length = 2, peptide_length = 33):
	"""
	Returns a sparse vector of kmer frequencies from a peptide sequence. Requires
	a kmer length (default = 2) and a target subsequence length (default = 33).
	"""
	vector = {}
	norm_const = float(peptide_length - kmer_length + 1)

	if len(aa_seq) < peptide_length:
		print len(aa_seq), ":" , aa_seq
		raise ValueError('Input protein sequence is shorter than target subsequence')

	else:

		overhang = len(aa_seq) - peptide_length
		st = int(overhang * random())
		subseq = aa_seq[st : (st + peptide_length)]
		#print "subseq",subseq

		for ia in xrange(peptide_length - kmer_length + 1):
			iz = ia + kmer_length
			kmer = subseq[ia:iz]
			#print "kmer",kmer
			index = 0
			for ik,k in enumerate(kmer):
				#print "k:",k
				index += aaindx[k] * 20 ** ik

			if index in vector:
				vector[index] += 1
			else:
				vector[index] = 1

		# Normalize vector
		vector = {ind : vector[ind] / norm_const for ind in vector}

	return vector
