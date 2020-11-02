# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 14:57:46 2020

@author: Meghana
"""
from __future__ import print_function
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import logging

def eval_complex2(rf=0,rf_nm = 0,inputs={},known_complex_nodes_list=[],prot_list=[],G=None,fin_list_graphs=[]):
    logging.info("Evaluating complexes...")     
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    
    if(rf == 1):
        if(rf_nm == 0):
            rf_nm = out_comp_nm+'_pred.out'
        with open(rf_nm) as fn:
            complex_list = [line.rstrip('\n').split() for line in fn ]# Space separated text only
            fin_list_graphs = [nx.Graph(G.subgraph(complex)) for complex in complex_list]                    
            # Just list of list of nodes 
            #fin_list_graphs = [nx.Graph(complex) for complex in complex_list]                    
    
    sizes_orig = [len(comp) for comp in fin_list_graphs]  
    
    p = inputs["eval_p"]
    
    N_pred_comp = len(fin_list_graphs)
    with open(out_comp_nm+'_metrics.out',"a") as fid:                
        print("No. of predicted complexes = ",N_pred_comp,file = fid)
    # Remove all proteins in predicted complexes that are not present in known complex protein list
    
    comp_remove=[]
    for comp in fin_list_graphs:
        to_remove=[]
        for node in comp.nodes():
            if(node not in prot_list):
                to_remove.append(node)
        comp.remove_nodes_from(to_remove)
        
        if(len(comp) <= 1): # Removing complexes with only one node or none 
            comp_remove.append(comp)
    
    fin_list_graphs = [graph for graph in fin_list_graphs if graph not in comp_remove]

    sizes_known = [len(comp) for comp in known_complex_nodes_list]
    # Size distributions 
    sizes_new = [len(comp) for comp in fin_list_graphs]        
    fig = plt.figure()
    sns.distplot(sizes_known,hist=False,label="known")
    sns.distplot(sizes_orig,hist=False,label="predicted")     
    sns.distplot(sizes_new,hist=False,label="predicted_known_prots")  
    
    plt.savefig(out_comp_nm+'_size_dists_known_pred.png')        
    #plt.show()
    plt.close(fig)            
    

    N_test_comp = len(known_complex_nodes_list)
    N_pred_comp = len(fin_list_graphs)
    with open(out_comp_nm+'_metrics.out',"a") as fid:        
        print("No. of predicted complexes after removing non-gold std proteins = ",N_pred_comp,file = fid)        

    N_matches_test = 0
    
    Metric = np.zeros((N_test_comp,N_pred_comp))
    
    for i,test_complex in enumerate(known_complex_nodes_list):
        N_match_pred=0
        for j,pred_complex in enumerate(fin_list_graphs):
            T = set(test_complex)
            P = set(pred_complex.nodes())
            C = len(T.intersection(P))
            A = len(P.difference(T))
            B = len(T.difference(P))
            
            if(float(C)/(A+C) > p and float(C)/(B+C) > p):
                Metric[i,j] = 1; 
                N_match_pred=N_match_pred+1
                
        if(N_match_pred > 0):
            N_matches_test=N_matches_test+1
            
   
    Recall = float(N_matches_test)/N_test_comp
    
    N_matches_pred = np.count_nonzero(np.sum(Metric,axis=0))
    Precision = float(N_matches_pred)/N_pred_comp
    
    with open(out_comp_nm+'_metrics.out',"a") as fid:                
        print("No. of known complexes = ",N_test_comp,file = fid)         
        print("Prediction Precision = ",Precision,file = fid)
        print("Prediction Recall = ",Recall,file = fid)        
        try:
            F1_score = 2*Precision*Recall/(Precision + Recall)
            print("Prediction F1 score = ",F1_score,file = fid)            
        except:
            print("Error in calculating F1 score - likely divide by 0")
            
                    
    logging.info("Finished Evaluating complexes.")        
   
