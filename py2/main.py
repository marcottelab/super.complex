# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 21:38:39 2018

@author: Meg_94
"""
from __future__ import print_function
from sup_graph import read_graphs2
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, './functions_py2/')
import yaml
import argparse
import os
#from test import test_class
import time 
start_time = time.time()
from testClassi2 import testClassi2
from trainClassi2 import trainClassi2
from feat_extract2 import feature_extract2
import logging

import joblib

#import joblib


# Write an input file parser for all the parameters 

def main():
    parser = argparse.ArgumentParser("Input parameters")
    parser.add_argument("--input_file_name", default="input_toy.yaml", help="Input parameters file name")
    parser.add_argument("--out_dir_name", default="/results", help="Output directory name")	
    args = parser.parse_args()
    
    with open(args.input_file_name, 'r') as f:
        inputs = yaml.load(f, yaml.Loader)

	# Override output directory name if same as gen 
    if(inputs['out_comp_nm'] == "/results/res"):
        if not os.path.exists(inputs['dir_nm']+args.out_dir_name):
            os.mkdir(inputs['dir_nm']+args.out_dir_name)
        inputs['out_comp_nm'] = args.out_dir_name + "/res"

    with open(inputs['dir_nm']+inputs['out_comp_nm']+"_input.yaml", 'w') as outfile:
        yaml.dump(inputs, outfile, default_flow_style=False)	
    
    logging.basicConfig(filename=inputs['dir_nm']+inputs['out_comp_nm']+"_logs.yaml",level=logging.INFO)
    start_time_read = time.time()                
    known_complex_nodes_list,prot_list,myGraph,complex_graphs,test_complex_graphs = read_graphs2(inputs)
    read_time = time.time() - start_time_read
    
    protlistfname = inputs['dir_nm']+"/res_protlist"
    known_complex_nodes_listfname = inputs['dir_nm']+"/res_known_complex_nodes_list"
    myGraphName = inputs['dir_nm']+"/res_myGraph"
    joblib.dump(prot_list, protlistfname)                
    joblib.dump(known_complex_nodes_list, known_complex_nodes_listfname) 
    joblib.dump(myGraph,myGraphName)
    
    if inputs['split_flag'] == 0:
        start_time_feat = time.time()            
        max_size_train,max_size_test,X_pos_test,X_neg_test,X_test,y_test,X_pos,y_pos,X,y,X_neg,y_neg=feature_extract2(inputs,complex_graphs,test_complex_graphs,myGraph)
        feat_time = time.time() - start_time_feat  
        
        max_size_trainF = inputs['dir_nm']+"/res_max_size_train"
        max_size_testF = inputs['dir_nm']+"/res_max_size_test"
        
        joblib.dump(max_size_train,max_size_trainF)
        joblib.dump(max_size_test,max_size_testF)
        
        
        if inputs['mode'] == 'non_gen':
            start_time_train = time.time()        
            model,scaler=trainClassi2(inputs['model_name'],inputs,X_pos,y_pos,X,y,X_neg,y_neg)
            train_time = time.time() - start_time_train 
            
            modelfname = inputs['dir_nm']+"/res_model"
            scalerfname = inputs['dir_nm']+"/res_scaler"
            joblib.dump(model, modelfname)                
            joblib.dump(scaler, scalerfname)
            
            start_time_test = time.time()      
            testClassi2(model,scaler,inputs,X_pos_test,X_neg_test,test_complex_graphs,X_test,y_test)
            test_time = time.time() - start_time_test                 
            
            tot_time = time.time() - start_time
            
            out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
            # Write to yaml file instead 
            with open(out_comp_nm+'_metrics.out',"a") as fid:
                print("--- Runtime performance ---",file = fid)
                print("Read network, complexes time (s) = ",read_time,"[",round(100*float(read_time)/tot_time,2),"%]",file = fid)                 
                print("Feature extraction time (s) = ",feat_time,"[",round(100*float(feat_time)/tot_time,2),"%]",file = fid)         
                print("Train time (s) = ",train_time,"[",round(100*float(train_time)/tot_time,2),"%]",file = fid) 
                print("Test time (s) = ",test_time,"[",round(100*float(test_time)/tot_time,2),"%]",file = fid)                 
                print("Total time (s) = ",tot_time,file = fid)

if __name__ == '__main__':
    main()
