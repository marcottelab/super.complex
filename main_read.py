# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 21:38:39 2018

@author: Meg_94
"""
from time import time as time_time

start_time = time_time()
from sup_graph import read_graphs
from os import path as os_path, mkdir as os_mkdir, chdir as os_chdir

os_chdir(os_path.dirname(os_path.abspath(__file__)))
from sys import path as sys_path

# insert at 1, 0 is the script path (or '' in REPL)
sys_path.insert(1, './functions_py3/')
from yaml import load as yaml_load, dump as yaml_dump, Loader as yaml_Loader
from argparse import ArgumentParser as argparse_ArgumentParser
from logging import basicConfig as logging_basicConfig, INFO as logging_INFO, DEBUG as logging_DEBUG
from pickle import dump as pickle_dump


def main():
    parser = argparse_ArgumentParser("Input parameters")
    parser.add_argument("--input_file_name", default="input_toy.yaml", help="Input parameters file name")
    parser.add_argument("--out_dir_name", default="/results", help="Output directory name")
    args = parser.parse_args()
    with open(args.input_file_name, 'r') as f:
        inputs = yaml_load(f, yaml_Loader)

    # Override output directory name if same as gen
    if inputs['out_comp_nm'] == "/results/res":
        if not os_path.exists(inputs['dir_nm'] + args.out_dir_name):
            os_mkdir(inputs['dir_nm'] + args.out_dir_name)
        inputs['out_comp_nm'] = args.out_dir_name + "/res"

    with open(inputs['dir_nm'] + inputs['out_comp_nm'] + "_input.yaml", 'w') as outfile:
        yaml_dump(inputs, outfile, default_flow_style=False)

    logging_basicConfig(filename=inputs['dir_nm'] + inputs['out_comp_nm'] + "_logs.yaml", level=logging_INFO)
    start_time_read = time_time()
    myGraph = read_graphs(inputs)
    read_time = time_time() - start_time_read

    myGraphName = inputs['dir_nm'] + "/res_myGraph"
    with open(myGraphName, 'wb') as f:
        pickle_dump(myGraph, f)

    tot_time = time_time() - start_time

    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    # Write to yaml file instead 
    with open(out_comp_nm + '_runtime_performance.out', "a") as fid:
        print("Read network time (s) = ", read_time, "[", round(100 * float(read_time) / tot_time, 2), "%]", file=fid)
        print("Total time (s) = ", tot_time, file=fid)


if __name__ == '__main__':
    main()
