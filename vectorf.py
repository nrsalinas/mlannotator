from random import random
from numpy import array


aaindx = {"A":0, "C":1, "D":2, "E":3, "F":4, "G":5, "H":6, "I":7, "K":8, "L":9,
"M":10, "N":11, "P":12, "Q":13, "R":14, "S":15, "T":16, "V":17, "W":18, "Y":19}

codons = {
"ATT":"I", "ATC":"I", "ATA":"I",
"CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L", "TTA":"L", "TTG":"L",
"GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
"TTT":"F", "TTC":"F",
"ATG": "M",
"TGT":"C", "TGC":"C",
"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
"GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
"CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
"ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
"TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
"TAT":"Y", "TAC":"Y",
"TGG":"W",
"CAA":"Q", "CAG":"Q",
"AAT":"N", "AAC":"N",
"CAT":"H", "CAC":"H",
"GAA":"E", "GAG":"E",
"GAT":"D", "GAC":"D",
"AAA":"K", "AAG":"K",
"CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
"TAA":"0", "TAG":"0", "TGA":"0"}


def kmerme(aa_seq, kmer_length = 2, peptide_length = 33, subsequencing = True):
	"""
	Returns a sparse vector of kmer frequencies from a peptide sequence. Requires
	a kmer length (default = 2) and a target subsequence length (default = 33).
	"""
	if subsequencing == False:
		peptide_length = len(aa_seq)

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
			for ik,k in enumerate(kmer[::-1]):
				#print "k:",k
				index += aaindx[k] * 20 ** ik
			#print index
			if index in vector:
				vector[index] += 1
			else:
				vector[index] = 1

		# Normalize vector
		vector = {ind : vector[ind] / norm_const for ind in vector}

	return vector

def column_cal(kmer_length = 2):
	"""Estimates the size of a kmer vector"""
	col = 20 ** kmer_length

def dna2aa(inseq):
	inseq = inseq.upper()
	seq0 = ""
	seq1 = ""
	seq2 = ""
	seqlen = len(inseq)
	mymod = seqlen % 3
	if mymod == 0:
		len0 = seqlen
		len1 = seqlen - 2
		len2 = seqlen - 1
	if mymod == 1:
		len0 = seqlen - 1
		len1 = seqlen
		len2 = seqlen - 2
	if mymod == 2:
		len0 = seqlen - 2
		len1 = seqlen - 1
		len2 = seqlen

	#First orf
	#print "First"
	for i in xrange(0,len0,3):
		#print inseq[i:(i+3)]
		aa = codons[inseq[i:(i+3)]]
		if aa == '0':
			break
		else:
			seq0 += aa

	#Second orf
	#print "Second"
	for i in xrange(1,len1,3):
		#print inseq[i:(i+3)]
		aa = codons[inseq[i:(i+3)]]
		if aa == '0':
			break
		else:
			seq1 += aa

	#Third orf
	#print "Third"
	for i in xrange(2,len2,3):
		#print inseq[i:(i+3)]
		aa = codons[inseq[i:(i+3)]]
		if aa == '0':
			break
		else:
			seq2 += aa

	return (seq0, seq1, seq2)
