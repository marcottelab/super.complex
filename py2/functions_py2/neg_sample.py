# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 22:11:32 2020

@author: Meg_94
"""
import argparse
import numpy as np

def main():
    parser = argparse.ArgumentParser("Input parameters")
    parser.add_argument("--neg_edge_list", default="../humap/results_50step/results_50step_isa_best/res_neg_comp_edges_list_all.out", help="Known edge list")
    parser.add_argument("--neg_edge_fin", default="../humap/results_50step/results_50step_isa_best/res_neg_comp_edges_list_fin.out", help="Predicted complex list")	    
 
    args = parser.parse_args()

    with open(args.neg_edge_list) as f:
        edg_all = f.readlines()
        
    fin = np.random.permutation(edg_all)[:10000]
    
    with open(args.neg_edge_fin,'w') as f:
        f.writelines(fin)
    
main()