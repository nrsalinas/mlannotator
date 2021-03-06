import readsim
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import GridSearchCV

bffr = ''
X, Y = readsim.get_arrays(kmer_size = 2, classes = 30, sim_seq_length = 33, reads_per_domain = 10)

#model = MultinomialNB()
#model = KNeighborsClassifier()
#model = GradientBoostingClassifier(min_samples_split = 5, min_samples_leaf = 5, max_depth = 5, learning_rate = 0.1,  n_estimators = 100, subsample = 1.0)
#model.fit(X_train, Y_train)
#myscore = model.score(X_test.toarray(), Y_test)
#print "Score :",myscore

##############################################
# Grid search GradientBoostingClassifier()
##############################################

param2test = {"min_samples_split" : [3, 5, 7],
	"min_samples_leaf" : [8, 10, 15, 20],
	"max_depth" : [6, 7, 8],
	"max_features" : ["sqrt"]
	}

model = GradientBoostingClassifier()
grid_search = GridSearchCV(estimator = model, param_grid = param2test, cv = 10)
bffr += "Training matrix shape: {0}\n".format(X.shape)
bffr += "Shape of training label vector: {0}\n".format(Y.shape)
grid_search.fit(X.toarray(), Y)

bffr += "Best score: {0}\n".format(grid_search.best_score_)
bffr += "Best parameter set: {0}\n".format(grid_search.best_params_)
bffr += "\n\nResults DataFrame:\n{0}\n".format(grid_search.cv_results_)
print bffr
