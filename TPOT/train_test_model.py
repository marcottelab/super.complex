import os
import numpy as np
import pandas as pd
import sklearn.metrics
from sklearn.externals import joblib

### BJL: paste necessary imports from TPOT here
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from tpot.builtins import ZeroCount

### BJL: training and test data are csv files with first two columns as IDs and label column called "label"
### all other columns should be features.
training_infile = "/project/bjl786/Synaptosome/Mouse/features/Mouse_features_07142018.train_labeled.csv"
test_infile = "/project/bjl786/Synaptosome/Mouse/features/Mouse_features_07142018.test_labeled.csv"

### BJL: outfiles
serialized_trained_model_outfile = "tpot_07142018_fitted_model.p"
pr_curve_outfile = "tpot_07142018_test_PRC.csv"
results_df_outfile = "tpot_07142018_test_resultsDF.csv"

##### BJL: paste TPOT pipeline here

exported_pipeline = make_pipeline(
    ZeroCount(),
    VarianceThreshold(threshold=0.0005),
    ExtraTreesClassifier(bootstrap=False, criterion="entropy", max_features=0.1, min_samples_leaf=3, min_samples_split=5, n_estimators=100)
)

##### Don't change anything below here

assert os.path.exists(training_infile), "{} not found".format(training_infile)
assert os.path.exists(test_infile), "{} not found".format(test_infile)

print "Reading training data"
train = pd.read_csv(training_infile, index_col=[0,1])
train_label = train.pop("label").values
train_data = train.values
    
print "Training model"
exported_pipeline.fit(train_data, train_label)
joblib.dump(exported_pipeline, serialized_trained_model_outfile)

train_fit_probs = exported_pipeline.predict_proba(train_data)[:,1]
train_aps = sklearn.metrics.average_precision_score(train_label,train_fit_probs)
print "Training set average precision score: {}".format(train_aps)

del train
del train_data

print "Reading test data"
test = pd.read_csv(test_infile, index_col=[0,1])
test_label = test.pop("label")
test_data = test.values

test_probs = exported_pipeline.predict_proba(test_data)[:,1]

test_aps = sklearn.metrics.average_precision_score(test_label, test_probs)
print "Test set average precision score: {}".format(test_aps)

test_p, test_r, _ = sklearn.metrics.precision_recall_curve(test_label, test_probs)

test_PRC = pd.DataFrame({"precision": test_p, "recall": test_r}).sort_values("recall")
test_PRC.to_csv(pr_curve_outfile,index=False)

test_DF = pd.DataFrame({"label":test_label,"P_1":test_probs}, index=test.index).sort_values("P_1", ascending=False)
test_DF["FDR"] = 1 - (test_DF.label.cumsum() / (np.arange(test_DF.shape[0]) + 1))
test_DF.to_csv(results_df_outfile)