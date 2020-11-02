# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 15:20:02 2020

@author: Meghana
"""

from tpot import TPOTClassifier
from deap import creator
import re
    
def tpotClassi2(inputs):
    classifier_file = inputs['classifier_file']

    best_classifier = inputs['best_classifier']       
    tpot = TPOTClassifier()
    with open(classifier_file) as f:
        raw_lines = f.readlines()
        words = [line.rstrip("\n").split("\t") for line in raw_lines[1:]]
        classifiers = [word[0] for word in words]
        pipelines = [word[2] for word in words]
        
    for i,classi in enumerate(classifiers):
        if classi == best_classifier:
            # if model is linearsvc then convert to svc
            pipeline_string = pipelines[i]
            #print(pipeline_string)
            # convert pipeline string to scikit-learn pipeline object
            deap_pipeline = creator.Individual.from_string(pipeline_string, tpot._pset)
            clf = tpot._toolbox.compile(expr=deap_pipeline)		
            if classi == "LinearSVC":
                #print(clf)
                n = len(clf.steps)
                linsvc = str(clf.steps.pop(n-1))
                #print(linsvc)
                #print(clf)
                match = re.search(r"C=(\d*.\d*)",linsvc)
                C_val = float(match.group(1))
                #print(C_val)
                from sklearn.svm import SVC
                clf.steps.append(('svc',SVC(kernel='linear',probability=True,C=C_val,tol=1e-05)))
    return clf
