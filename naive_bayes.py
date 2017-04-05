import readsim
from sklearn.naive_bayes import MultinomialNB
from sklearn.model_selection import GridSearchCV, train_test_split

bffr = ''

param2test = {'alpha': [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]}
model = MultinomialNB()
grid_search = GridSearchCV(estimator = model, param_grid = param2test, cv = 10)
bffr += "Training matrix shape: {0}\n".format(readsim.Xtrans.shape)
bffr += "Shape of training label vector: {0}\n".format(readsim.Y.shape)
#grid_search.fit(readsim.X_train.toarray(), readsim.Y_train)
grid_search.fit(readsim.Xtrans, readsim.Y)

bffr += "Best score: {0}\n".format(grid_search.best_score_)
bffr += "Best parameter set: {0}\n".format(grid_search.best_params_)
bffr += "\n\nResults DataFrame:\n{0}\n".format(grid_search.cv_results_)
print bffr
