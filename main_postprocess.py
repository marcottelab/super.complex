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
    parser.add_argument("--n_pts", default=1, help="number of partitions (computers)")
    parser.add_argument("--over_t", help="Overlap threshold")
    parser.add_argument("--model_dir", help="Directory containing model")
    parser.add_argument("--sample_dir", help="Sample files dir + /res")

    args = parser.parse_args()

    with open(args.input_file_name, 'r') as f:
        inputs = yaml_load(f, yaml_Loader)

    if args.over_t:
        inputs['over_t'] = float(args.over_t)
    if args.sample_dir:
        inputs['sample_dir'] = args.sample_dir
    if args.model_dir:
        inputs['model_dir'] = args.model_dir
    # Override output directory name if same as gen
    if inputs['out_comp_nm'] == "/results/res":
        if not os_path.exists(inputs['dir_nm'] + args.out_dir_name):
            os_mkdir(inputs['dir_nm'] + args.out_dir_name)
        inputs['out_comp_nm'] = args.out_dir_name + "/res"

    with open(inputs['dir_nm'] + inputs['out_comp_nm'] + "_input_pp.yaml", 'w') as outfile:
        yaml_dump(inputs, outfile, default_flow_style=False)

    logging_basicConfig(filename=inputs['dir_nm'] + inputs['out_comp_nm'] + "_logs.yaml", level=logging_INFO)
    # fin_list_graphs = control(myGraph,inputs,n=50)

    if "sample_dir" not in inputs:
        inputs['sample_dir']=inputs['out_comp_nm']
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    out_comp_nm_sample = inputs['dir_nm'] + inputs['sample_dir']
    out_comp_nm_model = inputs['dir_nm'] + inputs['model_dir']

    modelfname = out_comp_nm_model + "_model"
    scalerfname = out_comp_nm_model + "_scaler"
    myGraphName = inputs['dir_nm'] + "/res_myGraph"

    with open(scalerfname, 'rb') as f:
        scaler = pickle_load(f)
    with open(myGraphName, 'rb') as f:
        myGraph = pickle_load(f)

    pred_comp_list = []
    sdndap = pred_comp_list.append

    for i in range(int(args.n_pts)):
        with open(out_comp_nm_sample + "_pred_comp_list" + str(i), 'rb') as f:
            pred_comp_tmp = pickle_load(f)
        for snode in pred_comp_tmp:
            sdndap(snode)


    test_complex_path = inputs['dir_nm'] + "/res_test_known_complex_nodes_list"
    test_prot_list = get_prot_list(test_complex_path)

    train_complex_path = inputs['dir_nm'] + "/res_train_known_complex_nodes_list"
    train_prot_list = get_prot_list(train_complex_path)

    protlistfname = inputs['dir_nm'] + "/res_protlist"
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
