#! /usr/bin/env python
'''
python train_TPOT.py --training_data path/to/training_data.csv \
--outfile tpot_output.py \
--classifier_subset sklearn.ensemble.ExtraTreesClassifier sklearn.ensemble.RandomForestClassifier \
--generations 50 --population_size 50 --n_jobs 10 --temp_dir ./tpot_temp
'''
import gc
import argparse
import numpy as np
import pandas as pd
from tpot import TPOTClassifier
from sklearn.model_selection import train_test_split
import matplotlib
matplotlib.use('Agg')     # Issues warning on spyder - don't worry abt it
import matplotlib.pyplot as plt
import re
font = {'family' : 'normal',
        'size'   : 20}

matplotlib.rc('font', **font)
#matplotlib.rcParams.update({'font.size': 20})

from deap import creator
from sklearn.model_selection import cross_val_score

from tpot.export_utils import generate_pipeline_code
from sklearn.pipeline import make_pipeline,make_union
import sklearn.metrics
'''
                from sklearn.preprocessing import *
                from copy import copy
                from sklearn.ensemble import *
                from sklearn.svm import *
                from sklearn.tree import *
                from tpot.builtins import *
                
                from sklearn.feature_selection import *
                from sklearn.linear_model import *
                from sklearn.naive_bayes import *
'''
# Config dicts from https://github.com/EpistasisLab/tpot/blob/master/tpot/config/classifier.py

Classifiers = {

    'sklearn.naive_bayes.GaussianNB': {
    },

    'sklearn.naive_bayes.BernoulliNB': {
        'alpha': [1e-3, 1e-2, 1e-1, 1., 10., 100.],
        'fit_prior': [True, False]
    },

    'sklearn.naive_bayes.MultinomialNB': {
        'alpha': [1e-3, 1e-2, 1e-1, 1., 10., 100.],
        'fit_prior': [True, False]
    },

    'sklearn.tree.DecisionTreeClassifier': {
        'criterion': ["gini", "entropy"],
        'max_depth': range(1, 11),
        'min_samples_split': range(2, 21),
        'min_samples_leaf': range(1, 21)
    },
    'sklearn.ensemble.ExtraTreesClassifier': {
        'n_estimators': [100],
        'criterion': ["gini", "entropy"],
        'max_features': np.arange(0.05, 1.01, 0.05),
        'min_samples_split': range(2, 21),
        'min_samples_leaf': range(1, 21),
        'bootstrap': [True, False]
    },

    'sklearn.ensemble.RandomForestClassifier': {
        'n_estimators': [100],
        'criterion': ["gini", "entropy"],
        'max_features': np.arange(0.05, 1.01, 0.05),
        'min_samples_split': range(2, 21),
        'min_samples_leaf':  range(1, 21),
        'bootstrap': [True, False]
    },

    'sklearn.ensemble.GradientBoostingClassifier': {
        'n_estimators': [100],
        'learning_rate': [1e-3, 1e-2, 1e-1, 0.5, 1.],
        'max_depth': range(1, 11),
        'min_samples_split': range(2, 21),
        'min_samples_leaf': range(1, 21),
        'subsample': np.arange(0.05, 1.01, 0.05),
        'max_features': np.arange(0.05, 1.01, 0.05)
    },

    'sklearn.neighbors.KNeighborsClassifier': {
        'n_neighbors': range(1, 101),
        'weights': ["uniform", "distance"],
        'p': [1, 2]
    },

    'sklearn.svm.LinearSVC': {
        'penalty': ["l1", "l2"],
        'loss': ["hinge", "squared_hinge"],
        'dual': [True, False],
        'tol': [1e-5, 1e-4, 1e-3, 1e-2, 1e-1],
        'C': [1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1., 5., 10., 15., 20., 25.]
    },

    'sklearn.linear_model.LogisticRegression': {
        'penalty': ["l1", "l2"],
        'C': [1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1., 5., 10., 15., 20., 25.],
        'dual': [True, False]
    },

    'xgboost.XGBClassifier': {
        'n_estimators': [100],
        'max_depth': range(1, 11),
        'learning_rate': [1e-3, 1e-2, 1e-1, 0.5, 1.],
        'subsample': np.arange(0.05, 1.01, 0.05),
        'min_child_weight': range(1, 21),
        'nthread': [1]
    },
}

Preprocessors = {
    'sklearn.preprocessing.Binarizer': {
        'threshold': np.arange(0.0, 1.01, 0.05)
    },

    'sklearn.decomposition.FastICA': {
        'tol': np.arange(0.0, 1.01, 0.05)
    },

    'sklearn.cluster.FeatureAgglomeration': {
        'linkage': ['ward', 'complete', 'average'],
        'affinity': ['euclidean', 'l1', 'l2', 'manhattan', 'cosine']
    },

    'sklearn.preprocessing.MaxAbsScaler': {
    },

    'sklearn.preprocessing.MinMaxScaler': {
    },

    'sklearn.preprocessing.Normalizer': {
        'norm': ['l1', 'l2', 'max']
    },

    'sklearn.kernel_approximation.Nystroem': {
        'kernel': ['rbf', 'cosine', 'chi2', 'laplacian', 'polynomial', 'poly', 'linear', 'additive_chi2', 'sigmoid'],
        'gamma': np.arange(0.0, 1.01, 0.05),
        'n_components': range(1, 11)
    },

    'sklearn.decomposition.PCA': {
        'svd_solver': ['randomized'],
        'iterated_power': range(1, 11)
    },

    'sklearn.preprocessing.PolynomialFeatures': {
        'degree': [2],
        'include_bias': [False],
        'interaction_only': [False]
    },

    'sklearn.kernel_approximation.RBFSampler': {
        'gamma': np.arange(0.0, 1.01, 0.05)
    },

    'sklearn.preprocessing.RobustScaler': {
    },

    'sklearn.preprocessing.StandardScaler': {
    },

    'tpot.builtins.ZeroCount': {
    },

    'tpot.builtins.OneHotEncoder': {
        'minimum_fraction': [0.05, 0.1, 0.15, 0.2, 0.25],
        'sparse': [False]
    },

    # Selectors
    'sklearn.feature_selection.SelectFwe': {
        'alpha': np.arange(0, 0.05, 0.001),
        'score_func': {
            'sklearn.feature_selection.f_classif': None
        }
    },

    'sklearn.feature_selection.SelectPercentile': {
        'percentile': range(1, 100),
        'score_func': {
            'sklearn.feature_selection.f_classif': None
        }
    },

    'sklearn.feature_selection.VarianceThreshold': {
        'threshold': [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2]
    },

    'sklearn.feature_selection.RFE': {
        'step': np.arange(0.05, 1.01, 0.05),
        'estimator': {
            'sklearn.ensemble.ExtraTreesClassifier': {
                'n_estimators': [100],
                'criterion': ['gini', 'entropy'],
                'max_features': np.arange(0.05, 1.01, 0.05)
            }
        }
    },

    'sklearn.feature_selection.SelectFromModel': {
        'threshold': np.arange(0, 1.01, 0.05),
        'estimator': {
            'sklearn.ensemble.ExtraTreesClassifier': {
                'n_estimators': [100],
                'criterion': ['gini', 'entropy'],
                'max_features': np.arange(0.05, 1.01, 0.05)
            }
        }
    }

}

parser = argparse.ArgumentParser("Run TPOT to find a good machine learning pipeline on training data")
parser.add_argument("--training_data", required=True, help="Features with training labels. Columns: feat1,feat2...featN,label")
parser.add_argument("--testing_data", required=True, help="Features with testing labels. Columns: feat1,feat2...featN,label")
parser.add_argument("--outfile", required=True, help="File name to write the output pipeline to")
parser.add_argument("--outfile_classifiers", required=True, help="File name to write the output pipeline to")
parser.add_argument("--outfile_fig", required=True, help="File name to write the output pipeline to")
parser.add_argument("--classifier_subset", default=None, nargs="+", choices=Classifiers.keys(), help="Use a subset of sklearn's classifiers in search")
parser.add_argument("--score", default="average_precision", help="Which scoring function to use, default=average_precision")
parser.add_argument("--generations", type=int, default=100, help="How many generations to run, default=100")
parser.add_argument("--population_size", type=int, default=100, help="Size of the recombining population, default=100")
parser.add_argument("--n_jobs", type=int, default=1, help="How many jobs to run in parallel. Warning, n_jobs>1 is crashy")
parser.add_argument("--labels", type=int, nargs="+", default=[0,1], help="Which labels to retain, default=[0,1]")
parser.add_argument("--delimiter", default=",", help="Delimiter of training data, default=','")
parser.add_argument("--temp_dir", default="tpot_tmp", help="Temporary directory to stash intermediate results")
parser.add_argument("--warm_start", action='store_true', help="Flag: Whether to re-start TPOT from the results in the temp_dir")
args = parser.parse_args()


tpot_config = Preprocessors.copy()
if args.classifier_subset != None:
    classifiers = {i:Classifiers[i] for i in args.classifier_subset}
    tpot_config.update(classifiers)
else:
    tpot_config.update(Classifiers) # use all
    
    
print("Loading data")
df = pd.read_csv(args.training_data, sep=args.delimiter)
label_name = df.columns[-1]
print("Using '{}' as label column".format(label_name))

print("Dropping unlabeled rows")
df = df[df[label_name].isin(args.labels)]

labels = df.pop(label_name)
data = df.values

print("Loading test data")
df = pd.read_csv(args.testing_data, sep=args.delimiter)
label_name = df.columns[-1]
print("Using '{}' as label column".format(label_name))

print("Dropping unlabeled rows")
df = df[df[label_name].isin(args.labels)]

labels_test = df.pop(label_name)
data_test = df.values
        
print("Running TPOT")
#cv=5 , i.e by default it does 5-fold crossvalidation on the training set
tpot = TPOTClassifier(verbosity=2, scoring=args.score, config_dict=tpot_config,
                        generations=args.generations, population_size=args.population_size,
                        memory=args.temp_dir, n_jobs=args.n_jobs, warm_start=args.warm_start)
tpot.fit(data, labels)

tpot.export(args.outfile)

my_dict = dict(list(tpot.evaluated_individuals_.items()))
'''
# print a pipeline and its values
pipeline_str = list(tpot.evaluated_individuals_.keys())[0]
print(pipeline_str)
#print(tpot.evaluated_individuals_[pipeline_str])
# convert pipeline string to scikit-learn pipeline object
optimized_pipeline = creator.Individual.from_string(pipeline_str, tpot._pset) # deap object
fitted_pipeline = tpot._toolbox.compile(expr=optimized_pipeline ) # scikit-learn pipeline object
# print scikit-learn pipeline object
print(fitted_pipeline)
# Fix random state when the operator allows  (optional) just for get consistent CV score 
#tpot._set_param_recursive(fitted_pipeline.steps, 'random_state', 42)
# CV scores from scikit-learn
scores = cross_val_score(fitted_pipeline, data, labels, cv=5, scoring='accuracy', verbose=0)
print(np.mean(scores))
#print(tpot._evaluated_individuals[pipeline_str][1])
'''
final_dict = {}
for key,val in my_dict.items():
    cur_classifier = str(key)
    main_classifier = cur_classifier.split("(")[0]
    cur_score = val['internal_cv_score']
    if main_classifier in final_dict.keys():
        final_dict[main_classifier].append([cur_score,cur_classifier])
    else:
        final_dict[main_classifier] = []

with open(args.outfile_classifiers,"w") as f:
    f.write("Classifier\tScore\tPipeline\n")
    for key,val in final_dict.items():
        main_classifier = key
        val_scores = [item[0] for item in val]
        pipelines = [item[1] for item in val]
        if val_scores:
            cur_score = max(val_scores)
            cur_pipeline = pipelines[np.argmax(val_scores)]
            f.write(main_classifier + "\t" + str(cur_score) + "\t" + str(cur_pipeline) + "\n")

with open(args.outfile_classifiers) as f:
    raw_lines = f.readlines()
    words = [line.rstrip("\n").split("\t") for line in raw_lines[1:]]
    classifiers = [word[0] for word in words]
    pipelines = [word[2] for word in words]

fig = plt.figure(figsize=(8,6),dpi=300)   
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])    
plt.title('Precision-Recall curve')          
for i,classi in enumerate(classifiers):
    # if model is linearsvc then convert to svc
    pipeline_string = pipelines[i]
    #print(pipeline_string)
    # convert pipeline string to scikit-learn pipeline object
    deap_pipeline = creator.Individual.from_string(pipeline_string, tpot._pset)
    clf = tpot._toolbox.compile(expr=deap_pipeline)		
    if classi == "LinearSVC":
        print(clf)
        n = len(clf.steps)
        linsvc = str(clf.steps.pop(n-1))
        print(linsvc)
        print(clf)
        match = re.search(r"C=(\d*.\d*)",linsvc)
        C_val = float(match.group(1))
        print(C_val)
        from sklearn.svm import SVC
        clf.steps.append(('svc',SVC(kernel='linear',probability=True,C=C_val,tol=1e-05)))
        print(clf)
    clf.fit(data,labels)    
    test_fit_probs = clf.predict_proba(data_test)[:,1]
    test_aps = sklearn.metrics.average_precision_score(labels_test,test_fit_probs)
    test_p, test_r, _ = sklearn.metrics.precision_recall_curve(labels_test, test_fit_probs)
    plt.plot(test_r,test_p,label = classi+" - AP: "+str(round(test_aps,3)),linewidth=2.5)

plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
plt.savefig(args.outfile_fig,bbox_inches='tight')
#plt.show()            # Does not work on pod 
plt.close(fig)
