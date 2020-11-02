# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 12:42:40 2020

@author: Meg_94
"""
from __future__ import print_function
import numpy as np
import pandas as pd
from create_feat_mat2 import create_feat_mat2
from construct_neg_comps2 import construct_neg_comps2
import networkx as nx
import logging

def feature_extract2(inputs,complex_graphs,test_complex_graphs,G):
    n_feats = inputs['feats']
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    mode = inputs['mode']
    #mode = "non_gen" # Change to gen if you want to generate matrices
 
    #n_pos = len(complex_graphs)
    sizes = [len(complex) for complex in complex_graphs]
    min_size = min(sizes)
    max_size = max(sizes)
    max_size_train = max_size

    #n_pos_test = len(test_complex_graphs)
    sizes = [len(complex) for complex in test_complex_graphs]
    min_size_test = min(sizes)
    max_size_test = max(sizes)
    max_size_test = max_size_test        

    if inputs['model_type'] == "tpot" and mode == "non_gen": # CHANGE X_POS, Y_POS later !!!!
        # Read X,y from file
        df_full = pd.read_csv(inputs['train_feat_mat'])
        y = df_full.pop('complex')
        X = df_full
        
        neg_start_ind = y[y == 0].index[0]
        X_pos = X.iloc[0:neg_start_ind]
        y_pos = y[0:neg_start_ind]
        X_neg = X.iloc[neg_start_ind:]
        y_neg = y[neg_start_ind:]
        
        df_full = pd.read_csv(inputs['test_feat_mat'])
        y_test = df_full.pop('complex')
        X_test = df_full
        
        neg_start_ind = y_test[y_test == 0].index[0]
        X_pos_test = X_test.iloc[0:neg_start_ind]
        y_pos_test = y_test[0:neg_start_ind]
        X_neg_test = X_test.iloc[neg_start_ind:]
        y_neg_test = y_test[neg_start_ind:]
    else:               

        logging.info("Feature extraction...")        
        
        X_pos = create_feat_mat2(complex_graphs,n_feats)
       
        n_pos = len(X_pos)
        
        dims = X_pos.shape
        n_feats = dims[1]
        with open(out_comp_nm+'_metrics.out',"a") as fid:        
            print("No. of features = ",n_feats,file = fid)          
            print("No. of train positive complexes = ",n_pos,file = fid) 
            
        logging.info("Constructing negative complexes...") 
        

        neg_comp_list = construct_neg_comps2(min_size,max_size,n_pos,inputs,G) 
        logging.info("Finished constructing negative complexes")                  

        X_neg = create_feat_mat2(neg_comp_list,n_feats)        
         
        #Removing negative feature rows that exactly match any row in positives
        for ind in range(n_pos):
            matching_inds = np.where((X_neg == X_pos[ind]).all(axis=1))
            X_neg = np.delete(X_neg,matching_inds,axis=0)
            for index in sorted(list(matching_inds[0]), reverse=True):
                del neg_comp_list[index]             

        n_neg = len(X_neg)
        #print(n_neg)
        # HHANDLE CASE WHEN n_neg = 0 !!!!!
        with open(out_comp_nm+'_metrics.out',"a") as fid:                
            print("No. of train negative complexes = ",n_neg,file = fid)
        
        with open(out_comp_nm+'_neg_train.out',"w") as fn:
            with open(out_comp_nm+'_neg_train_edges.out',"wb") as f_edges:
                for index in range(len(neg_comp_list)):
                    for node in neg_comp_list[index].nodes():
                        fn.write("%s " % node)            
                    nx.write_weighted_edgelist(neg_comp_list[index],f_edges)
                    f_edges.write("\n".encode())
                    fn.write("\n")
             
         
        X=np.vstack((X_pos,X_neg))
        
        y_pos = [1 for i in range(n_pos)]
        y_neg = [0 for i in range(n_neg)]            
        y = y_pos + y_neg 
        y = np.array(y) 
        y_pos = np.array(y_pos) 
        y_neg = np.array(y_neg)   
        
        # Writing raw training data to csv in tpot format 
        feat_list = ["dens","nodes","degree_max","degree_mean","degree_median","degree_var","CC_max","CC_mean","CC_var","edge_wt_mean","edge_wt_max","edge_wt_var","DC_mean","DC_var","DC_max","sv1","sv2","sv3","complex"]

        dat = np.hstack((X,y[:,None]))
        df = pd.DataFrame(dat)
        df.to_csv(path_or_buf=out_comp_nm+"_train_dat.csv",index=False,header=feat_list)
        
        X_pos_test = create_feat_mat2(test_complex_graphs,n_feats)
                
        n_pos = len(X_pos_test)
        # Negative tests
        
        #N_neg = len(test_complex_graphs)
        #neg_test_comp = [random_graph(n) for n in range(N_neg)]
        
            
        logging.info("Constructing test negative complexes...")        
        neg_test_comp = construct_neg_comps2(min_size_test,max_size_test,n_pos,inputs,G) 
        logging.info("Finished constructing test negative complexes.")     
    
        
        X_neg_test = create_feat_mat2(neg_test_comp,n_feats)
        
        X_allpos=np.vstack((X_pos,X_pos_test))
        n_allpos = len(X_allpos)
        
        #Removing negative feature rows that exactly match any row in positives
        for ind in range(n_allpos):
            matching_inds = np.where((X_neg_test == X_allpos[ind]).all(axis=1))
            X_neg_test = np.delete(X_neg_test,matching_inds,axis=0) 
            for index in sorted(list(matching_inds[0]), reverse=True):
                del neg_test_comp[index] 
                
        with open(out_comp_nm+'_neg_test.out',"w") as fn:
            with open(out_comp_nm+'_neg_test_edges.out',"wb") as f_edges:
                for index in range(len(neg_test_comp)):
                    for node in neg_test_comp[index].nodes():
                        fn.write("%s " % node)            
                    nx.write_weighted_edgelist(neg_test_comp[index],f_edges)
                    fn.write("\n")
                    f_edges.write("\n".encode())                    


        n_neg = len(X_neg_test)
        
        y_pos = [1 for i in range(n_pos)]
        y_neg = [0 for i in range(n_neg)]            
        y = y_pos + y_neg 
        y_test = np.array(y)         
        
        X_test=np.vstack((X_pos_test,X_neg_test))    
        X_pos_test = X_pos_test
        X_neg_test = X_neg_test
        
        # Writing raw test data to csv in tpot format 
        dat = np.hstack((X_test,y_test[:,None]))
        df = pd.DataFrame(dat)
        df.to_csv(path_or_buf=out_comp_nm+"_test_dat.csv",index=False,header=feat_list)

        with open(out_comp_nm+'_metrics.out',"a") as fid:        
            print("No. of test positive complexes = ",n_pos,file = fid)     
            print("No. of test negative complexes = ",n_neg,file = fid)  
            
    logging.info("Finished Feature extraction") 
    return max_size_train,max_size_test,X_pos_test,X_neg_test,X_test,y_test,X_pos,y_pos,X,y,X_neg,y_neg
