# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 15:15:30 2020

@author: Meghana
"""
from __future__ import print_function
import numpy as np
import sklearn.metrics
import matplotlib
matplotlib.use('Agg')     # Issues warning on spyder - don't worry abt it
import matplotlib.pyplot as plt
import logging

#@profile
def testClassi2(model,scaler,inputs,X_pos_test,X_neg_test,test_complex_graphs,X_test,y_test):
    logging.info("Evaluating test complexes...")  
    model_type = inputs['model_type']
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']      
    
    if(model_type == "tpot"):
        res_pos = model.predict(X_pos_test)
        res = model.predict(X_neg_test)            
        
        if hasattr(model, 'decision_function'):            
            score = model.decision_function(X_pos_test)
            np.savetxt(out_comp_nm+'_test_pos_score.out',score)
            #print("Scores for positive complexes are",score)
            score = model.decision_function(X_neg_test)
            np.savetxt(out_comp_nm+'_test_neg_score.out',score)                
            
        # Write the else case 
        
    elif(model_type == "NN"):

        X_pos_test = scaler.transform(X_pos_test)        
        
        preds = model.predict(X_pos_test)
        res_pos = [np.argmax(pred) for pred in preds]
        score = np.array([pred[1] for pred in preds])
        np.savetxt(out_comp_nm+'_test_pos_score.out',score)
        
        X_neg_test = scaler.transform(X_neg_test)             
        preds = model.predict(X_neg_test)
        res = [np.argmax(pred) for pred in preds]            
        # Score of being negative !!
        score = np.array([pred[0] for pred in preds])
        np.savetxt(out_comp_nm+'_test_pos_score.out',score)  
    #print("Scores for negative complexes are",score)
        
    n_pos = len(test_complex_graphs)
    n_neg = len(X_neg_test)
               
    TP = sum(res_pos) # assuming negatives are 0s
    FN = n_pos - TP
    acc = TP/float(n_pos) 
    
    TN = sum([res[ind] == 0 for ind in range(len(res))])
    FP = n_neg - TN
    acc_neg = TN /float(n_neg) # assuming negatives are 0s
    
    Recall = float(TP)/(TP+FN) # Just accuracy of test positives
    Precision = float(TP)/(TP + FP)
    F1_score = 2*Precision*Recall/(Precision + Recall)
    
    if(model_type == "tpot"):
        test_fit_probs = model.predict_proba(X_test)[:,1]
        
    elif(model_type == "NN"):
        X_test = scaler.transform(X_test)             
        preds = model.predict(X_test)
        test_fit_probs = np.array([pred[1] for pred in preds])           
        
    test_aps = sklearn.metrics.average_precision_score(y_test,test_fit_probs)
    with open(out_comp_nm+'_metrics.out',"a") as fid:          
        print("Training set average precision score = ",test_aps,file = fid)    
    
    test_p, test_r, _ = sklearn.metrics.precision_recall_curve(y_test, test_fit_probs)
        

    fig = plt.figure()        
    plt.plot(test_r,test_p)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title('Precision-Recall curve: AP={0:0.2f}'.format(test_aps))
    plt.savefig(out_comp_nm+'_pr_curve.png')
    #plt.show()            # Does not work on pod 
    plt.close(fig)        

    with open(out_comp_nm+'_metrics.out',"a") as fid:        
        print("Accuracy for test positive complexes = ",acc,file = fid)                   
        print("Accuracy for test negative complexes = ",acc_neg,file = fid) #Really just tells you complex or not for random graphs    
        print("Test Precision = ",Precision,file = fid)
        print("Test Recall = ",Recall,file = fid)
        print("Test F1 score = ",F1_score,file = fid)
    

    logging.info("Finished evaluating test complexes.")        
    
