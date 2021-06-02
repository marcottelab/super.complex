# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 17:36:07 2020

@author: Meg_94
"""
from time import time as time_time

start_time = time_time()

from matplotlib import use as mpl_use

mpl_use('Agg')  # Issues warning on spyder - don't worry abt it
from os import path as os_path, mkdir as os_mkdir, chdir as os_chdir

os_chdir(os_path.dirname(os_path.abspath(__file__)))
from sys import path as sys_path

# insert at 1, 0 is the script path (or '' in REPL)
sys_path.insert(1, './functions_py3/')
from yaml import load as yaml_load, dump as yaml_dump, Loader as yaml_Loader
from argparse import ArgumentParser as argparse_ArgumentParser
from postprocess import postprocess
from pickle import load as pickle_load, dump as pickle_dump
from logging import basicConfig as logging_basicConfig, INFO as logging_INFO, DEBUG as logging_DEBUG
from logging import info as logging_info, debug as logging_debug
from glob import glob
from calc_max_overlap_train import get_overlap_threshold, get_overlap_threshold_qi


def get_prot_list(test_complex_path):
    with open(test_complex_path, 'rb') as f:
        test_complex_list = pickle_load(f)

    test_complex_nodes = [item for sublist in test_complex_list for item in sublist]
    test_prot_list = set(test_complex_nodes)

    with open(test_complex_path + "_prot_list",'wb') as f:
        pickle_dump(test_prot_list, f)
    return test_prot_list


def main():
    parser = argparse_ArgumentParser("Input parameters")
    parser.add_argument("--input_file_name", default="input_toy.yaml", help="Input parameters file name")
    parser.add_argument("--out_dir_name", default="/results", help="Output directory name")
    parser.add_argument("--train_test_files_dir", default="", help="Train test file path")
    parser.add_argument("--n_pts", default=1, help="number of partitions (computers)")
    parser.add_argument("--over_t", help="Overlap threshold")
    parser.add_argument("--model_dir", help="Directory containing model")
    parser.add_argument("--sample_dir", help="Sample files dir + /res")
    parser.add_argument("--sample_folders_prefix", help="Input parameters file name /results..")    
    parser.add_argument("--sample_folders_prefix_final", help="Input file name to use final merged results Use as /results..")    
    parser.add_argument("--sample_folders_list",  nargs='+',help="Input parameters file name /results.. separated by commas")    
    parser.add_argument("--graph_files_dir", default="", help="Graph files' folder path")
    
    parser.add_argument("--overlap_method",default=1,help="Overlap method option: qi, default: jaccard")    
    parser.add_argument("--infer_overlap_threshold",default='n',help="y or n")    

    args = parser.parse_args()
    with open(args.input_file_name, 'r') as f:
        inputs = yaml_load(f, yaml_Loader)
        
    if args.overlap_method:
        inputs['overlap_method'] = args.overlap_method
    if args.over_t:
        inputs['over_t'] = float(args.over_t)
    if args.sample_dir:
        inputs['sample_dir'] = args.sample_dir
    if args.model_dir:
        inputs['model_dir'] = args.model_dir
    if args.infer_overlap_threshold:
        inputs['infer_overlap_threshold'] = args.infer_overlap_threshold        
 
        
    # Override output directory name if same as gen
    if args.out_dir_name or inputs['out_comp_nm'] == "/results/res":
        if not os_path.exists(inputs['dir_nm'] + args.out_dir_name):
            os_mkdir(inputs['dir_nm'] + args.out_dir_name)
        inputs['out_comp_nm'] = args.out_dir_name + "/res"
                
    inputs['train_test_files_dir'] = ''
    if args.train_test_files_dir:
        if not os_path.exists(inputs['dir_nm'] + args.train_test_files_dir):
            os_mkdir(inputs['dir_nm'] + args.train_test_files_dir)
        inputs['train_test_files_dir'] = args.train_test_files_dir     
        
    inputs['graph_files_dir'] = ''
    if args.graph_files_dir:
        if not os_path.exists(inputs['dir_nm'] + args.graph_files_dir):
            os_mkdir(inputs['dir_nm'] + args.graph_files_dir)
        inputs['graph_files_dir'] = args.graph_files_dir               

    logging_basicConfig(filename=inputs['dir_nm'] + inputs['out_comp_nm'] + "_logs.yaml", level=logging_INFO)
    # fin_list_graphs = control(myGraph,inputs,n=50)

    if "sample_dir" not in inputs:
        inputs['sample_dir']=inputs['out_comp_nm']
        
    myGraphName = inputs['dir_nm']+ inputs['graph_files_dir'] + "/res_myGraph"
    with open(myGraphName, 'rb') as f:
        myGraph = pickle_load(f)
        
    if 'infer_overlap_threshold' in inputs:
        if inputs['infer_overlap_threshold'] == 'y':
            pp_flag = 0
            if inputs['dir_nm'] == 'yeast':
                pp_flag = 1
            if 'overlap_method' in inputs:
                if inputs['overlap_method'] == 'qi':
                    inputs['over_t'] = get_overlap_threshold_qi(inputs,pp_flag,myGraph) 
                else:
                    inputs['over_t'] = get_overlap_threshold(inputs,pp_flag,myGraph) 
                    
    with open(inputs['dir_nm'] + inputs['out_comp_nm'] + "_input_pp.yaml", 'w') as outfile:
        yaml_dump(inputs, outfile, default_flow_style=False)                    
        
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    out_comp_nm_sample = inputs['dir_nm'] + inputs['sample_dir']
    out_comp_nm_model = inputs['dir_nm'] + inputs['model_dir']

    modelfname = out_comp_nm_model + "_model"
    scalerfname = out_comp_nm_model + "_scaler"
        
    with open(scalerfname, 'rb') as f:
        scaler = pickle_load(f)

    pred_comp_list = []
    sdndap = pred_comp_list.append

    if args.sample_folders_list:
        for folder in args.sample_folders_list:
            allfiles = './' + inputs['dir_nm'] + folder + '/res_pred_comp_list*'
            for fname in glob(allfiles, recursive=True):                
                with open(fname, 'rb') as f:
                    pred_comp_tmp = pickle_load(f)
                for snode in pred_comp_tmp:
                    sdndap(snode)      
    elif args.sample_folders_prefix_final:
        allsubd = './' + inputs['dir_nm'] + args.sample_folders_prefix_final + '*/res_pred.out'
        for fname in glob(allsubd, recursive=True):
            with open(fname) as f:
                complexes_score = [line.rstrip().split() for line in f.readlines()]
                pred_comp_tmp = [(frozenset(comp[:-1]),float(comp[-1])) for comp in complexes_score]
                
            for snode in pred_comp_tmp:
                sdndap(snode) 
    elif args.sample_folders_prefix:
        allsubd = './' + inputs['dir_nm'] + args.sample_folders_prefix + '*/res_pred_comp_list*'
        for fname in glob(allsubd, recursive=True):
            with open(fname, 'rb') as f:
                pred_comp_tmp = pickle_load(f)
            for snode in pred_comp_tmp:
                sdndap(snode)            
    else:
        for i in range(int(args.n_pts)):
            with open(out_comp_nm_sample + "_pred_comp_list" + str(i), 'rb') as f:
                pred_comp_tmp = pickle_load(f)
            for snode in pred_comp_tmp:
                sdndap(snode)
    len_pred_comp_list = 'No. of complexes before pp = ' + str(len(pred_comp_list))
    logging_info(len_pred_comp_list)
    test_complex_path = inputs['dir_nm'] + inputs['train_test_files_dir']+ "/res_test_known_complex_nodes_list"
    test_prot_list = get_prot_list(test_complex_path)

    train_complex_path = inputs['dir_nm'] + inputs['train_test_files_dir']+ "/res_train_known_complex_nodes_list"
    train_prot_list = get_prot_list(train_complex_path)

    protlistfname = inputs['dir_nm'] + inputs['train_test_files_dir']+ "/res_protlist"
    with open(protlistfname, 'rb') as f:
        prot_list = pickle_load(f)

    start_time_pp = time_time()
    fin_list_graphs = postprocess(pred_comp_list, modelfname, scaler, inputs, myGraph, prot_list, train_prot_list, test_prot_list)
    pp_time = time_time() - start_time_pp

    tot_time = time_time() - start_time

    # Write to yaml file instead 
    with open(out_comp_nm + '_runtime_performance.out', "a") as fid:
        print("--- Runtime performance ---", file=fid)
        print("Post processing complex time (s) = ", pp_time, "[", round(100 * float(pp_time) / tot_time, 2), "%]",file=fid)
        print("Total time (s) = ", tot_time, file=fid)


if __name__ == '__main__':
    main()
