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
from sample import sample
# from random_walk_control import control

from logging import basicConfig as logging_basicConfig, INFO as logging_INFO, DEBUG as logging_DEBUG
from pickle import load as pickle_load, dump as pickle_dump


def main():
    parser = argparse_ArgumentParser("Input parameters")
    parser.add_argument("--input_file_name", default="input_toy.yaml", help="Input parameters file name")
    parser.add_argument("--out_dir_name", default="/results", help="Output directory name")
    parser.add_argument("--seed_mode", help="Seed mode - specify 'cliques' for the cliques algo")
    parser.add_argument("--search_method", help="Sampling algorithm")
    parser.add_argument("--model_dir", help="Directory containing model")
    parser.add_argument("--ptnum", default='0', help="partition number")
    parser.add_argument("--explore_prob", default=0.01, help="probability of exploring")
    parser.add_argument("--prob_metropolis", default=0.1, help="metropolis probability")
    parser.add_argument("--T0", default=0.88, help="isa T0")
    parser.add_argument("--alpha", default=1.8, help="isa alpha")
    parser.add_argument("--transfer2tmp",default=True,help="Transfer to tmp folder")

    args = parser.parse_args()

    with open(args.input_file_name, 'r') as f:
        inputs = yaml_load(f, yaml_Loader)

    if args.seed_mode:
        inputs['seed_mode'] = args.seed_mode
    if args.search_method:
        inputs['search_method'] = args.search_method
    if args.model_dir:
        inputs['model_dir'] = args.model_dir
    if args.explore_prob:
        inputs['explore_prob'] = float(args.explore_prob)
    if args.prob_metropolis:
        inputs['prob_metropolis'] = float(args.prob_metropolis)
    if args.T0:
        inputs['T0'] = float(args.T0)
    if args.alpha:
        inputs['alpha'] = float(args.alpha)

    # Override output directory name if same as gen
    if inputs['out_comp_nm'] == "/results/res":
        if not os_path.exists(inputs['dir_nm'] + args.out_dir_name):
            os_mkdir(inputs['dir_nm'] + args.out_dir_name)
        inputs['out_comp_nm'] = args.out_dir_name + "/res"

    with open(inputs['dir_nm'] + inputs['out_comp_nm'] + "_input_sample.yaml", 'w') as outfile:
        yaml_dump(inputs, outfile, default_flow_style=False)

    logging_basicConfig(filename=inputs['dir_nm'] + inputs['out_comp_nm'] + "_logs.yaml", level=logging_INFO)
    # fin_list_graphs = control(myGraph,inputs,n=50)
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    out_comp_nm_model = inputs['dir_nm'] + inputs['model_dir']

    modelfname = out_comp_nm_model + "_model"
    scalerfname = out_comp_nm_model + "_scaler"

    max_sizeF = inputs['dir_nm'] + "/res_max_size"
    with open(max_sizeF, 'rb') as f:
        max_size = pickle_load(f)

    with open(scalerfname, 'rb') as f:
        scaler = pickle_load(f)

    myGraph = None
    if inputs['seed_mode'] == "cliques":
        myGraphName = inputs['dir_nm'] + "/res_myGraph"
        with open(myGraphName, 'rb') as f:
            myGraph = pickle_load(f)

    seed_nodes_F = out_comp_nm + "_seed_nodes" + args.ptnum
    with open(seed_nodes_F, 'rb') as f:
        seed_nodes = pickle_load(f)

    start_time_sample = time_time()
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']

    num_comp = sample(inputs, myGraph, modelfname, scaler, seed_nodes, max_size,args.transfer2tmp)

    sample_time = time_time() - start_time_sample
    sample_time_avg = sample_time / num_comp
    folNm_out = "/tmp/" + out_comp_nm + "_orig_comps"

    pred_comp_list = [pickle_load(open(folNm_out + "/" + seed_node, 'rb')) for seed_node in seed_nodes if
                      os_path.exists(folNm_out + "/" + seed_node)]

    with open(out_comp_nm + "_pred_comp_list" + args.ptnum, "wb") as f:
        pickle_dump(pred_comp_list, f)
    tot_time = time_time() - start_time

    with open(out_comp_nm + '_runtime_performance.out', "a") as fid:
        print("--- Runtime performance ---", file=fid)
        print("Sample time (s) = ", sample_time, "[", round(100 * float(sample_time) / tot_time, 2), "%]", file=fid)
        print("Average sample time (s) = ", sample_time_avg, file=fid)
        print("Total time (s) = ", tot_time, file=fid)


if __name__ == '__main__':
    main()
