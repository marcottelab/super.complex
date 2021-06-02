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
from pickle import load as pickle_load, dump as pickle_dump
from logging import basicConfig as logging_basicConfig, INFO as logging_INFO, DEBUG as logging_DEBUG
from numpy.random import permutation as rand_perm
from networkx import find_cliques as nx_find_cliques
from os import listdir as os_listdir, path as os_path
from numpy import percentile as np_percentile
from math import ceil as math_ceil

def main():
    parser = argparse_ArgumentParser("Input parameters")
    parser.add_argument("--input_file_name", default="input_toy.yaml", help="Input parameters file name")
    parser.add_argument("--out_dir_name", default="/results", help="Output directory name")
    parser.add_argument("--train_test_files_dir", default="", help="Train test file path")    
    parser.add_argument("--graph_files_dir", default="", help="Graph files' folder path") 
    parser.add_argument("--seed_mode", help="Seed mode - specify 'cliques' for the cliques algo")
    parser.add_argument("--max_size_thres", help="Max size threshold")    
    parser.add_argument("--n_pts", default=1, help="number of partitions (computers)")
    args = parser.parse_args()

    with open(args.input_file_name, 'r') as f:
        inputs = yaml_load(f, yaml_Loader)

    if args.seed_mode:
        inputs['seed_mode'] = args.seed_mode
    if args.max_size_thres:
        inputs['max_size_thres'] = int(args.max_size_thres)        

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

    with open(inputs['dir_nm'] + inputs['out_comp_nm'] + "_input_sample_partition.yaml", 'w') as outfile:
        yaml_dump(inputs, outfile, default_flow_style=False)

    logging_basicConfig(filename=inputs['dir_nm'] + inputs['out_comp_nm'] + "_logs.yaml", level=logging_INFO)
        
    neig_dicts_folder = inputs['dir_nm'] +inputs['graph_files_dir']+ "/neig_dicts"

    num_comp = inputs['num_comp']
    max_size_thres = inputs['max_size_thres']
    max_size_trainF = inputs['dir_nm'] + inputs['train_test_files_dir']+ "/res_max_size_train"
    with open(max_size_trainF, 'rb') as f:
        max_size_train = pickle_load(f)

    max_size = max_size_train
    
    max_sizeF_feat = inputs['dir_nm'] + inputs['train_test_files_dir']+ "/res_max_size_search"  
    if os_path.exists(max_sizeF_feat):
        with open(max_sizeF_feat, 'rb') as f:
            max_size = pickle_load(f)
    else:            
        with open(inputs['dir_nm'] + inputs['comf_nm']) as f:
            sizes = [len(line.rstrip().split()) for line in f.readlines()]    
        max_size = max(sizes)
        q1 = np_percentile(sizes, 25)
        q3 = np_percentile(sizes, 75)
        max_wo_outliers = math_ceil(q3 + 4.5*(q3-q1))  # Maximum after removing outliers    
        max_size = min(max_size,max_wo_outliers)
        
        
    if max_size >= max_size_thres:
        max_size = max_size_thres
        
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']

    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("Max number of steps for complex growth = ", max_size, file=fid)  # NOT actual max size since you merge later
    
    max_sizeF = inputs['dir_nm'] + inputs['train_test_files_dir']+ "/res_max_size_search_par"
    
    with open(max_sizeF, 'wb') as f:
        pickle_dump(max_size, f)

    seed_mode = inputs['seed_mode']

    if seed_mode == "all_nodes":
        #graph_nodes = list(myGraph.nodes())
        seed_nodes = rand_perm(os_listdir(neig_dicts_folder))
    elif seed_mode == "n_nodes":
        seed_nodes = rand_perm(os_listdir(neig_dicts_folder))[:num_comp]
    elif seed_mode == "all_nodes_known_comp":
        protlistfname = inputs['dir_nm']+ inputs['train_test_files_dir'] + "/res_protlist"
        with open(protlistfname, 'rb') as f:
            prot_list = pickle_load(f)        
        seed_nodes = list(prot_list)
    elif seed_mode == "cliques":
        myGraphName = inputs['dir_nm'] + inputs['graph_files_dir']+ "/res_myGraph"
        with open(myGraphName, 'rb') as f:
            myGraph = pickle_load(f)        
        clique_list = list(nx_find_cliques(myGraph))
        to_rem = []
        # Removing 2 node and big complexes
        for comp in clique_list:
            if len(comp) <= 2 or len(comp) >= max_size:
                to_rem.append(comp)

        for comp in to_rem:
            clique_list.remove(comp)

        seed_nodes = clique_list  # Remove duplicates later.

    # partition
    ptns = int(args.n_pts)

    nc = len(seed_nodes)
    if seed_mode == 'n_nodes':
        seed_nodes_F = out_comp_nm + "_seed_nodes"
        each_ptn = nc // ptns
        for i in range(ptns - 1):
            with open(seed_nodes_F + str(i), 'wb') as f:
                pickle_dump(seed_nodes[i * each_ptn:(i + 1) * each_ptn], f)
        with open(seed_nodes_F + str(ptns - 1), 'wb') as f:
            pickle_dump(seed_nodes[(ptns - 1) * each_ptn:], f)
    else:
        seed_nodes_dir =  inputs['dir_nm'] + inputs['graph_files_dir']+ "/" + seed_mode + "_n_pts_" + str(ptns)

        if not os_path.exists(seed_nodes_dir):
            os_mkdir(seed_nodes_dir)
            seed_nodes_F = seed_nodes_dir + "/res_seed_nodes"
            each_ptn = nc // ptns
            for i in range(ptns - 1):
                with open(seed_nodes_F + str(i), 'wb') as f:
                    pickle_dump(seed_nodes[i * each_ptn:(i + 1) * each_ptn], f)

            with open(seed_nodes_F + str(ptns - 1), 'wb') as f:
                pickle_dump(seed_nodes[(ptns - 1) * each_ptn:], f)


if __name__ == '__main__':
    main()
