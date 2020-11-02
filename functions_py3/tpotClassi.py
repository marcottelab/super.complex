# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 19:06:04 2019

@author: Meghana
"""
from tpot import TPOTClassifier
from deap import creator
from re import search as re_search


def train_model(pipeline_string, classi):
    tpot = TPOTClassifier()
    # if model is linearsvc then convert to svc
    # convert pipeline string to scikit-learn pipeline object
    deap_pipeline = creator.Individual.from_string(pipeline_string, tpot._pset)
    clf = tpot._toolbox.compile(expr=deap_pipeline)
    if classi == "LinearSVC":
        n = len(clf.steps)
        linsvc = str(clf.steps.pop(n - 1))
        match = re_search(r"C=(\d*.\d*)", linsvc)
        C_val = float(match.group(1))
        from sklearn.svm import SVC
        clf.steps.append(('svc', SVC(kernel='linear', probability=True, C=C_val, tol=1e-05)))
    return clf


def tpot_classi(inputs):
    classifier_file = inputs['classifier_file']
    clf = None
    with open(classifier_file) as f:
        raw_lines = f.readlines()

    words = [line.rstrip("\n").split("\t") for line in raw_lines[1:]]
    if 'best_classifier' in inputs:
        select_classifier = inputs['best_classifier']
        classifiers = [word[0] for word in words]
        pipelines = [word[2] for word in words]
        try:
            ind = classifiers.index(select_classifier)
            pipeline_string = pipelines[ind]
            clf = train_model(pipeline_string, select_classifier)
        except ValueError:
            print("Classifier not in list")
    else:
        best_pipeline = max(words, key=lambda x: float(x[1]))
        classi = best_pipeline[0]
        out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
        with open(out_comp_nm + '_metrics.out', "a") as fid:
            print("Model trained = Best classifier = ", classi, file=fid)
        pipeline_string = best_pipeline[2]
        clf = train_model(pipeline_string, classi)
    return clf
