import readsim
import tablef
import time
import vectorf
import re
import pickle
from scipy.sparse import dok_matrix
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.externals import joblib
from sklearn.feature_extraction.text import TfidfTransformer
from numpy import array, concatenate
from os import remove

train = False
predict = True
model_file = "gbc_model_no_tfidf.pkl"
accuracy_file = "accuracy_nt.txt"

annotation_table = "../030_annotations.csv"
sam_file = "../030.sam"
super_ref = 'supertable.csv'
max_predict_batches = 200

# Get table superfamilies abundances and their id numbers
table = tablef.parse_supertable(super_ref)
table = table[:30]
mySPIDs = [x[0] for x in table]

t0 = time.time()

model = None

if train:
	X, Y = readsim.get_arrays(kmer_size = 2, classes = 30, sim_seq_length = 33, \
		reads_per_domain = 10)

	model = GradientBoostingClassifier(min_samples_split = 7,min_samples_leaf = 20,\
	max_depth = 8, learning_rate = 0.1, n_estimators = 100, subsample = 1.0, \
	max_features = 'sqrt')

	model.fit(X, Y)
	joblib.dump(model, model_file)

else:
	model = joblib.load(model_file)

if predict:
#	target_reads_in_csv = 0
	reads_processed_in_sam = 0
	reads_parsed_pred = 0

	try:
		remove("predictions.pkl")
	except OSError:
		pass
	try:
		remove(accuracy_file)
	except OSError:
		pass
	spa_mat = dok_matrix((1100, (20 ** 2)))
	predCount = 0
	preds = array([])
	reads_annotated = []
	reads_predicted = []
	final_ref = {}
	target_id = []
	with open(annotation_table, "r") as csvhandle:
		for ann in csvhandle:
			ann = ann.rstrip()
			abits = ann.split(",")
			thisSPID = int(re.search("SuperfamilyID\:(\d+)\s",ann).group(1))
			if thisSPID in mySPIDs:
				#target_reads_in_csv += 1
				reads_annotated.append(abits[0])
				final_ref[abits[0]] = [thisSPID]

	with open(sam_file, "r") as samhandle:
		classif_attemp = 0
		for alg in samhandle:
			dnaseq = ""
			orientation = ""
			bbits = re.split("\s+", alg)
			to_remove = None
			for indr, read in enumerate(reads_annotated):
				if read == bbits[0]:
					orientation = bbits[1]
					dnaseq = bbits[9]
					reads_processed_in_sam += 1
					to_remove = indr
					break
			if to_remove is not None:
				reads_annotated.pop(to_remove)
			if orientation == "0" and not re.search('[^acgtACGT]',dnaseq):
				aaseqs = vectorf.dna2aa(dnaseq)
				topredict = filter(lambda x: len(x) > 10, aaseqs)
				if topredict:
					classif_attemp += 1
					reads_parsed_pred += 1
					reads_predicted += [bbits[0] for x in xrange(len(topredict))]
					veclist = [vectorf.kmerme(x, subsequencing = False) for x in topredict]
					for ind,dic in enumerate(veclist):
						for key in dic:
							spa_mat[(predCount + ind),key] = dic[key]
					predCount += len(topredict)
					if predCount > 1000:
						spa_mat.resize((predCount, (20 ** 2)))
						tfidrer = TfidfTransformer()
						matTrans = tfidrer.fit_transform(spa_mat.toarray())
						thisPred = model.predict(matTrans.toarray())
						#preds = concatenate((preds, thisPred))
						max_predict_batches -= 1
						hits = 0
						for myre,pr in zip(reads_predicted, thisPred):
							if int(final_ref[myre][0]) == int(pr):
								hits += 1
						accuracy = float(hits) / float(classif_attemp)
						with open(accuracy_file,"a") as fhandle:
							fhandle.write("{0}\n".format(accuracy))
						reads_predicted = []
						classif_attemp = 0
						spa_mat.clear()
						spa_mat.resize((1100, (20 ** 2)))
						predCount = 0

			if max_predict_batches < 1:
				break

	if predCount > 0:
		spa_mat.resize((predCount, (20 ** 2)))
		tfidrer = TfidfTransformer()
		matTrans = tfidrer.fit_transform(spa_mat.toarray())
		thisPred = model.predict(matTrans.toarray())
		#preds = concatenate((preds, thisPred))
		max_predict_batches -= 1
		hits = 0
		for myre,pr in zip(reads_predicted, thisPred):
			if int(final_ref[myre][0]) == int(pr):
				hits += 1
		accuracy = float(hits) / float(classif_attemp)
		with open(accuracy_file,"a") as fhandle:
			fhandle.write("{0}\n".format(accuracy))

	print "target reads in csv:", len(final_ref)
	print "reads parsed for prediction:" , reads_parsed_pred

	#pickle.dump(final_ref, open("predictions.pkl","w"))


print ((time.time() - t0) / 60.0), "mins"

# process 28552
