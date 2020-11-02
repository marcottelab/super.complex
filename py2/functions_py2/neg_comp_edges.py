# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 21:21:58 2020

@author: Meg_94
"""
import argparse
import networkx as nx
def main():
    parser = argparse.ArgumentParser("Input parameters")
    parser.add_argument("--prot_list", default="../known_prot_list.txt", help="Known prot list")
    parser.add_argument("--network", default="../humap/test_graph.txt", help="Network")	    
    parser.add_argument("--reduced_network", default="../humap/test_graph_reduced_known.txt", help="Predicted edge list")	
    parser.add_argument("--neg_edge_list_out", default="../humap/results_50step/results_50step_isa_best/res_neg_comp_edges_list.out", help="Negative complex edge list")	
    parser.add_argument("--known_edge_list", default="../humap/results_50step/results_50step_isa_best/res_known_edges_list.out", help="Known edge list")	
    
    args = parser.parse_args()
    
    G_full = nx.read_weighted_edgelist(args.network,nodetype=str);

    
    # Removing self loops 
    G_full.remove_edges_from(nx.selfloop_edges(G_full))
    
    with open(args.known_edge_list) as f:
        known_edges = [tuple(line.strip("\n").split("\t")) for line in f]    
    G_full.remove_edges_from(known_edges)
    
    '''
    with open(args.prot_list) as f:
        prot_list = [line.strip("\n") for line in f]
        
    G_full =   G_full.subgraph(prot_list)    
    '''
    nx.write_edgelist(G_full, args.neg_edge_list_out,delimiter='\t',data=False)
        
        
    
     

main()