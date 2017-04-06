import readsim
from sklearn.svm import LinearSVC
from sklearn.model_selection import GridSearchCV

bffr = ''
X, Y = readsim.get_arrays(kmer_size = 2, classes = 30, sim_seq_length = 33, reads_per_domain = 10)

##############################################
# Grid search GradientBoostingClassifier()
##############################################

param2test = {"C" : [0.5, 1.0, 2.0, 5.0],
	"loss" : ["hinge", "squared_hinge"]
#	"penalty" : ["l1" , "l2"]
	}

model = LinearSVC()
grid_search = GridSearchCV(estimator = model, param_grid = param2test, cv = 10)
bffr += "Training matrix shape: {0}\n".format(X.shape)
bffr += "Shape of training label vector: {0}\n".format(Y.shape)
grid_search.fit(X, Y)

bffr += "Best score: {0}\n".format(grid_search.best_score_)
bffr += "Best parameter set: {0}\n".format(grid_search.best_params_)
bffr += "\n\nResults DataFrame:\n{0}\n".format(grid_search.cv_results_)
print bffr
