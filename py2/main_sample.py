# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 20:45:17 2020

@author: Meg_94
"""

from __future__ import print_function
from time import time as time_time
start_time = time_time()

from matplotlib import use as mpl_use
mpl_use('Agg')     # Issues warning on spyder - don't worry abt it
from sys import path as sys_path
# insert at 1, 0 is the script path (or '' in REPL)
sys_path.insert(1, './functions_py2/')
from yaml import load as yaml_load,dump as yaml_dump, Loader as yaml_Loader
from argparse import ArgumentParser as argparse_ArgumentParser
from os import path as os_path, mkdir as os_mkdir
from eval_complex2 import eval_complex2
from sample2 import sample2
from postprocess2 import postprocess2
#from random_walk_control import control

from logging import basicConfig as logging_basicConfig, INFO as logging_INFO, DEBUG as logging_DEBUG
from cProfile import Profile as cProfile_Profile
from pstats import Stats as pstats_Stats
from joblib import load as joblib_load

def main():
    parser = argparse_ArgumentParser("Input parameters")
    parser.add_argument("--input_file_name", default="input_toy.yaml", help="Input parameters file name")
    parser.add_argument("--out_dir_name", default="/results", help="Output directory name")	
    parser.add_argument("--seed_mode", help= "Seed mode - specify 'cliques' for the cliques algo")
    parser.add_argument("--search_method", help = "Sampling algorithm")
    args = parser.parse_args()
    
    with open(args.input_file_name, 'r') as f:
        inputs = yaml_load(f, yaml_Loader)

    if(args.seed_mode):
        inputs['seed_mode'] = args.seed_mode
    if(args.search_method):
        inputs['search_method'] = args.search_method
        
	# Override output directory name if same as gen 
    if(inputs['out_comp_nm'] == "/results/res"):
        if not os_path.exists(inputs['dir_nm']+args.out_dir_name):
            os_mkdir(inputs['dir_nm']+args.out_dir_name)
        inputs['out_comp_nm'] = args.out_dir_name + "/res" 
        
    with open(inputs['dir_nm']+inputs['out_comp_nm']+"_input_sample.yaml", 'w') as outfile:
        yaml_dump(inputs, outfile, default_flow_style=False)	
        
    logging_basicConfig(filename=inputs['dir_nm']+inputs['out_comp_nm']+"_logs.yaml",level=logging_INFO)
    #fin_list_graphs = control(myGraph,inputs,n=50)

    rf = 0 # 1 when you want to read from file for evaluation
    rf_nm = "humap/humap_2stage_clustering_res.txt"

    #eval_complex(rf,rf_nm,inputs,known_complex_nodes_list,prot_list,myGraph,fin_list_graphs)
    modelfname = inputs['dir_nm']+"/res_model"
    scalerfname = inputs['dir_nm']+"/res_scaler"
    protlistfname = inputs['dir_nm']+"/res_protlist"
    known_complex_nodes_listfname = inputs['dir_nm']+"/res_known_complex_nodes_list"
    myGraphName = inputs['dir_nm']+"/res_myGraph"
    max_size_trainF = inputs['dir_nm']+"/res_max_size_train"
    max_size_testF = inputs['dir_nm']+"/res_max_size_test"
    
    model = joblib_load(modelfname)                
    scaler = joblib_load(scalerfname)
    prot_list = joblib_load(protlistfname)
    known_complex_nodes_list = joblib_load(known_complex_nodes_listfname)
    myGraph = joblib_load(myGraphName)
    
    max_size_train = joblib_load(max_size_trainF)
    max_size_test = joblib_load(max_size_testF)
    
    if rf == 1:
        eval_complex2(rf,rf_nm,inputs,known_complex_nodes_list,prot_list,myGraph)
    else:
        start_time_sample = time_time() 
        prof = cProfile_Profile()
        pred_comp_list,num_comp = prof.runcall(sample2,inputs,myGraph,model,scaler,max_size_train,max_size_test,prot_list)
        
       # fin_list_graphs,num_comp = sample(inputs['num_comp'],inputs['run_mode'],inputs['seed_mode'],inputs,myGraph,model,scaler,max_size_train,max_size_test,prot_list,complex_graphs,test_complex_graphs,known_complex_nodes_list)
        sample_time = time_time() - start_time_sample
        sample_time_avg = sample_time/num_comp
        
        start_time_pp = time_time() 
        fin_list_graphs = postprocess2(pred_comp_list,model,scaler,inputs,myGraph)
        pp_time = time_time() - start_time_pp
        
        start_time_eval = time_time()        
        eval_complex2(rf,rf_nm,inputs,known_complex_nodes_list,prot_list,myGraph,fin_list_graphs)
        eval_time = time_time() - start_time_eval
        
        tot_time = time_time() - start_time
        
        out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
        # Write to yaml file instead 
        with open(out_comp_nm+'_metrics.out',"a") as fid:
            print("--- Runtime performance ---",file = fid)
            print("Sample time (s) = ",sample_time,"[",round(100*float(sample_time)/tot_time,2),"%]",file = fid)       
            print("Average sample time (s) = ",sample_time_avg,"[",round(100*float(sample_time)/tot_time,2),"%]",file = fid)               
            print("Post processing complex time (s) = ",pp_time,"[",round(100*float(eval_time)/tot_time,2),"%]",file = fid)                       
            print("Evaluate complex time (s) = ",eval_time,"[",round(100*float(eval_time)/tot_time,2),"%]",file = fid)           
            print("Total time (s) = ",tot_time,file = fid)
            print("Sampling breakup:",file = fid)
            ps = pstats_Stats(prof, stream=fid).sort_stats('cumulative')
            ps.print_stats()
if __name__ == '__main__':
    main()
        
