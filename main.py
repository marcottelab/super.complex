# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 21:38:39 2018

@author: Meg_94
"""
from __future__ import print_function
from sup_graph import Sup_graph
import yaml
import argparse
import os
#from test import test_class
import time 
start_time = time.time()
# Write an input file parser for all the parameters 

def main():
    parser = argparse.ArgumentParser("Input parameters")
    parser.add_argument("--input_file_name", default="input_toy.yaml", help="Input parameters file name")
    parser.add_argument("--out_dir_name", default="/results", help="Output directory name")	
    args = parser.parse_args()
    
    with open(args.input_file_name, 'r') as f:
        try:
            inputs = yaml.load(f)
        except yaml.YAMLError as exc:
            print(exc)    

	# Override output directory name if same as gen 
    if(inputs['out_comp_nm'] == "/results/res"):
        if not os.path.exists(inputs['dir_nm']+args.out_dir_name):
            os.mkdir(inputs['dir_nm']+args.out_dir_name)
        inputs['out_comp_nm'] = args.out_dir_name + "/res"

    with open(inputs['dir_nm']+inputs['out_comp_nm']+"_input.yaml", 'w') as outfile:
        yaml.dump(inputs, outfile, default_flow_style=False)	
    # Initializing object 
    my_sup_graph = Sup_graph(inputs)
    #my_sup_graph_child = test_class()
    #my_sup_graph_child.blah_blah()
    
    start_time_read = time.time()                
    my_sup_graph.read_graphs()
    read_time = time.time() - start_time_read
    
    rf = 0 # 1 when you want to read from file fro evaluation
    rf_nm = "yeast/dip_res/results_clustereps/result_clustereps_tap_train_default.txt"
    
    if rf == 1:
        my_sup_graph.eval_complex(rf,rf_nm)
    else:
        if inputs['split_flag'] == 0:
            start_time_feat = time.time()            
            my_sup_graph.feature_extract(inputs['feats'])
            feat_time = time.time() - start_time_feat  
            
            if inputs['mode'] == 'non_gen':
                start_time_train = time.time()        
                my_sup_graph.train_model(inputs['model_name'])
                train_time = time.time() - start_time_train    
                
                start_time_test = time.time()      
                my_sup_graph.test()
                test_time = time.time() - start_time_test        
             
                start_time_sample = time.time()    
                my_sup_graph.sample(inputs['num_comp'],inputs['run_mode'],inputs['seed_mode'])
                sample_time = time.time() - start_time_sample
            
                start_time_eval = time.time()        
                my_sup_graph.eval_complex(rf,rf_nm)
                eval_time = time.time() - start_time_eval
                
                tot_time = time.time() - start_time
                
                out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
                # Write to yaml file instead 
                with open(out_comp_nm+'_metrics.out',"a") as fid:
                    print("--- Runtime performance ---",file = fid)
                    print("Read network, complexes time (s) = ",read_time,"[",round(100*float(read_time)/tot_time,2),"%]",file = fid)                 
                    print("Feature extraction time (s) = ",feat_time,"[",round(100*float(feat_time)/tot_time,2),"%]",file = fid)         
                    print("Train time (s) = ",train_time,"[",round(100*float(train_time)/tot_time,2),"%]",file = fid) 
                    print("Test time (s) = ",test_time,"[",round(100*float(test_time)/tot_time,2),"%]",file = fid)                 
                    print("Sample time (s) = ",sample_time,"[",round(100*float(sample_time)/tot_time,2),"%]",file = fid)       
                    print("Evaluate complex time (s) = ",eval_time,"[",round(100*float(eval_time)/tot_time,2),"%]",file = fid)           
                    print("Total time (s) = ",tot_time,file = fid)
                 

main()
