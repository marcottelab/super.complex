# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 21:38:39 2018

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
from testClassi import test_classi
from trainClassi import train_classi
from feat_extract import feature_extract
from read_complexes import read_complexes

from logging import basicConfig as logging_basicConfig, INFO as logging_INFO, DEBUG as logging_DEBUG
from pickle import dump as pickle_dump, load as pickle_load


def main():
    parser = argparse_ArgumentParser("Input parameters")
    parser.add_argument("--input_file_name", default="input_toy.yaml", help="Input parameters file name")
    parser.add_argument("--out_dir_name", default="/results", help="Output directory name")
    parser.add_argument("--split_flag", help="Train test split to do")
    parser.add_argument("--classifier_file", help="classifier file")
    parser.add_argument("--train_feat_mat", help="Train feat mat")
    parser.add_argument("--test_feat_mat", help="Test feat mat")
    parser.add_argument("--scale_factor", help="No. of times negatives are greater than positives")
    parser.add_argument("--neg_sample_meth", help="Method for sampling negatives - uniform / same")
    parser.add_argument("--mode", help="Generate feature matrices or not")
    args = parser.parse_args()

    with open(args.input_file_name, 'r') as f:
        inputs = yaml_load(f, yaml_Loader)

    # Override output directory name if same as gen
    if inputs['out_comp_nm'] == "/results/res":
        if not os_path.exists(inputs['dir_nm'] + args.out_dir_name):
            os_mkdir(inputs['dir_nm'] + args.out_dir_name)
        inputs['out_comp_nm'] = args.out_dir_name + "/res"

    # Override split flag and mode if present
    if args.split_flag:
        inputs['split_flag'] = int(args.split_flag)
    if args.mode:
        inputs['mode'] = args.mode
    if args.classifier_file:
        inputs['classifier_file'] = args.classifier_file
    if args.train_feat_mat:
        inputs['train_feat_mat'] = args.train_feat_mat
    if args.test_feat_mat:
        inputs['test_feat_mat'] = args.test_feat_mat
    if args.scale_factor:
        inputs['scale_factor'] = int(args.scale_factor)
    if args.neg_sample_meth:
        inputs['neg_sample_method'] = args.neg_sample_meth

    with open(inputs['dir_nm'] + inputs['out_comp_nm'] + "_input.yaml", 'w') as outfile:
        yaml_dump(inputs, outfile, default_flow_style=False)

    logging_basicConfig(filename=inputs['dir_nm'] + inputs['out_comp_nm'] + "_logs.yaml", level=logging_INFO)

    myGraphName = inputs['dir_nm'] + "/res_myGraph"
    with open(myGraphName, 'rb') as f:
        myGraph = pickle_load(f)

    start_time_read_c = time_time()
    known_complex_nodes_list, complex_graphs, test_complex_graphs, prot_list, test_known_complex_nodes_list, train_known_complex_nodes_list = read_complexes(inputs, myGraph)
    read_time_c = time_time() - start_time_read_c

    known_complex_nodes_listfname = inputs['dir_nm'] + "/res_known_complex_nodes_list"
    with open(known_complex_nodes_listfname, 'wb') as f:
        pickle_dump(known_complex_nodes_list, f)

    train_known_complex_nodes_listfname = inputs['dir_nm'] + "/res_train_known_complex_nodes_list"
    with open(train_known_complex_nodes_listfname, 'wb') as f:
        pickle_dump(train_known_complex_nodes_list, f)
    
    test_known_complex_nodes_listfname = inputs['dir_nm'] + "/res_test_known_complex_nodes_list"
    with open(test_known_complex_nodes_listfname, 'wb') as f:
        pickle_dump(test_known_complex_nodes_list, f)
   
    protlistfname = inputs['dir_nm'] + "/res_protlist"
    with open(protlistfname, 'wb') as f:
        pickle_dump(prot_list, f)

    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    if inputs['split_flag'] == 0:
        start_time_feat = time_time()
        max_size_train, max_size_test, X_pos_test, X_neg_test, X_test, y_test, X_pos, y_pos, X, y, X_neg, y_neg = feature_extract(
            inputs, complex_graphs, test_complex_graphs, myGraph)
        feat_time = time_time() - start_time_feat

        max_size_trainF = inputs['dir_nm'] + "/res_max_size_train"
        max_size_testF = inputs['dir_nm'] + "/res_max_size_test"

        with open(max_size_trainF, 'wb') as f:
            pickle_dump(max_size_train, f)
        with open(max_size_testF, 'wb') as f:
            pickle_dump(max_size_test, f)

        if inputs['mode'] == 'non_gen':
            start_time_train = time_time()
            model, scaler = train_classi(inputs['model_name'], inputs, X_pos, y_pos, X, y, X_neg, y_neg)
            train_time = time_time() - start_time_train

            modelfname = out_comp_nm + "_model"
            scalerfname = out_comp_nm + "_scaler"
            with open(modelfname, 'wb') as f:
                pickle_dump(model, f)
            with open(scalerfname, 'wb') as f:
                pickle_dump(scaler, f)

            start_time_test = time_time()
            test_classi(model, scaler, inputs, X_pos_test, X_neg_test, test_complex_graphs, X_test, y_test)
            test_time = time_time() - start_time_test

            tot_time = time_time() - start_time

            # Write to yaml file instead
            with open(out_comp_nm + '_runtime_performance.out', "a") as fid:
                print("--- Runtime performance ---", file=fid)
                print("Read complexes time (s) = ", read_time_c, "[", round(100 * float(read_time_c) / tot_time, 2),
                      "%]", file=fid)
                print("Feature extraction time (s) = ", feat_time, "[", round(100 * float(feat_time) / tot_time, 2),
                      "%]", file=fid)
                print("Train time (s) = ", train_time, "[", round(100 * float(train_time) / tot_time, 2), "%]",
                      file=fid)
                print("Test time (s) = ", test_time, "[", round(100 * float(test_time) / tot_time, 2), "%]", file=fid)
                print("Total time (s) = ", tot_time, file=fid)


if __name__ == '__main__':
    main()
