# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 21:42:41 2018
Copyright (C) - Unpublished work of Meghana Venkata Palukuri
All Rights Reserved
Do not redistribute without the express permission of the author 
@author: Meg_94
"""
from matplotlib import use as mpl_use
mpl_use('Agg')  # Issues warning on spyder - don't worry abt it
from re import findall as re_findall
from logging import info as logging_info
from os import mkdir as os_mkdir, path as os_path
from joblib import dump as joblib_dump
from numpy import mean as np_mean
import networkx as nx


def write_graph_stats_neig_lists(G,inputs):
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    n_nodes_net = nx.number_of_nodes(G)
    n_edges_net = nx.number_of_edges(G)

    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("No. of nodes in network = ", n_nodes_net, file=fid)
        print("No. of edges in network = ", n_edges_net, file=fid)

    myGraphdict = nx.to_dict_of_dicts(G)

    folNm = inputs['dir_nm'] + inputs['graph_files_dir']+ "/neig_dicts"
    if not os_path.exists(folNm):
        os_mkdir(folNm)
    neig_lens = []
    for node, val in myGraphdict.items():
        with open(folNm + "/" + node, 'wb') as f:
            joblib_dump(val, f)
        neig_lens.append(len(val))
    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("Max number of neighbors = ", max(neig_lens), file=fid)
        print("Avg number of neighbors = %.2f " % np_mean(neig_lens), file=fid)
    logging_info("Finished writing neighbor lists.")    
    
    
def read_graphs(inputs):
    network_path = inputs['dir_nm'] + inputs['netf_nm']
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    use_full = inputs['use_full']

    with open(out_comp_nm + '_metrics.out', "w") as fid:
        print("Network:", network_path, file=fid)

    logging_info("Reading network...")
    G = nx.read_weighted_edgelist(network_path, nodetype=str)

    logging_info("Finished Reading network")

    # Removing self loops 
    G.remove_edges_from(nx.selfloop_edges(G))
    
    # CHECK IF DUPLICATE NODES AND EDGES ARE REMOVED

    if use_full == 0:
        logging_info("Extracting sub network with proteins from complexes...")
        # Extracting subnetwork of the original PPIN containing proteins in the known complexes
        complex_path = inputs['dir_nm'] + inputs['comf_nm']
        test_complex_path = inputs['dir_nm'] + inputs['comf_test_nm']
        with open(complex_path) as f:
            data = f.read()
            id_list_train = re_findall(r"[\w']+", data)

        with open(test_complex_path) as f:
            data = f.read()
            id_list_test = re_findall(r"[\w']+", data)

        prot_list = set(id_list_train + id_list_test)
        G = G.subgraph(prot_list)
  
    if use_full == 1:
        write_graph_stats_neig_lists(G,inputs)        

    return G
