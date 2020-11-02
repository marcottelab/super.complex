# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 17:36:07 2020

@author: Meg_94
"""
from time import time as time_time
start_time = time_time()

from matplotlib import use as mpl_use
mpl_use('Agg')     # Issues warning on spyder - don't worry abt it
from sys import path as sys_path
# insert at 1, 0 is the script path (or '' in REPL)
sys_path.insert(1, './functions_py3/')
from yaml import load as yaml_load,dump as yaml_dump, Loader as yaml_Loader
from argparse import ArgumentParser as argparse_ArgumentParser
from os import path as os_path, mkdir as os_mkdir
from eval_complex import eval_complex
from postprocess import postprocess
from pickle import load as pickle_load
from logging import basicConfig as logging_basicConfig, INFO as logging_INFO, DEBUG as logging_DEBUG


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


    #eval_complex(rf,rf_nm,inputs,known_complex_nodes_list,prot_list,myGraph,fin_list_graphs)
    modelfname = inputs['dir_nm']+"/res_model"
    scalerfname = inputs['dir_nm']+"/res_scaler"
    protlistfname = inputs['dir_nm']+"/res_protlist"
    known_complex_nodes_listfname = inputs['dir_nm']+"/res_known_complex_nodes_list"
    myGraphName = inputs['dir_nm']+"/res_myGraph"
    
    #with open(modelfname,'rb') as f:                
        #model = pickle_load(f)   
    with open(scalerfname,'rb') as f:                                         
        scaler = pickle_load(f)
    with open(protlistfname,'rb') as f:    
        prot_list = pickle_load(f)
    with open(known_complex_nodes_listfname,'rb') as f:
        known_complex_nodes_list = pickle_load(f)
    with open(myGraphName,'rb') as f:
        myGraph = pickle_load(f)
            
    with open("./humap/results_eval_fixed_isa_o0.9/res_pred_complexes_covid_ids.out") as f:
        dat = f.readlines()
    
    pred_comp_list = []    
    for line in dat:
        words = line.strip().split(" ")
        prots = words[:-1]
        score = words[-1]
        pred_comp_list.append((prots,score))        
        
    start_time_pp = time_time() 
    fin_list_graphs = postprocess(pred_comp_list,modelfname,scaler,inputs,myGraph)
    pp_time = time_time() - start_time_pp
    
    start_time_eval = time_time()        
    eval_complex(0,0,inputs,known_complex_nodes_list,prot_list,myGraph,fin_list_graphs)
    eval_time = time_time() - start_time_eval
    
    tot_time = time_time() - start_time
    
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    # Write to yaml file instead 
    with open(out_comp_nm+'_runtime_performance.out',"a") as fid:
        print("--- Runtime performance ---",file = fid)
        print("Post processing complex time (s) = ",pp_time,"[",round(100*float(pp_time)/tot_time,2),"%]",file = fid)                       
        print("Evaluate complex time (s) = ",eval_time,"[",round(100*float(eval_time)/tot_time,2),"%]",file = fid)           
        print("Total time (s) = ",tot_time,file = fid)
        
if __name__ == '__main__':
    main() 
