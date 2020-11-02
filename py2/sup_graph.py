# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 21:42:41 2018
Copyright (C) - Unpublished work of Meghana Venkata Palukuri
All Rights Reserved
Do not redistribute without the express permission of the author 
@author: Meg_94
"""
from __future__ import print_function
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, './functions_py2/')
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')     # Issues warning on spyder - don't worry abt it
import matplotlib.pyplot as plt
import seaborn as sns
from jaccard_coeff2 import jaccard_coeff2
import re
import logging
from os import mkdir as os_mkdir, path as os_path
from joblib import dump as joblib_dump


#from profilehooks import profile

def split_train_test_complexes(inputs,G):
    sep = inputs['sep']
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    
    all_complexes_path = inputs['dir_nm'] + inputs["comf_nm_all"]    
    with open(all_complexes_path) as f:
        raw_lines = f.readlines()
        
    lines1_orig = [line.rstrip("\r\n").split(sep) for line in raw_lines if len(line.rstrip("\r\n").split(sep)) > 2]
        
    # Remove Nones 
    lines1 = []
    for line in lines1_orig:
        line_noNone = [item for item in line if item != "None"]
        if line_noNone:
            lines1.append(line_noNone)
        
    lines1=np.random.permutation(lines1)   
    
    # Removing redundancy from list: 
    # merge complexes with jaccard similarity greater than 0.6 
    # Loop since it is a forward merging algorithm 
    n_merge = 1
    while (n_merge != 0):
        n_merge = 0
        complex_list = []
        
        lines2 = list(lines1) # Assign as new variable to avoid call by ref
        for line1 in lines1:
            lines2_rem = list(lines2)
            lines2_rem.remove(line1)
            for line2 in lines2_rem:
                jc = jaccard_coeff2(line1,line2)
                
                if jc >= 0.6:
                    # merge and store in line1
                    n_merge +=1 
                    line1 = line1 + line2
            
            complex_list.append(line1)
        lines1 = list(complex_list)
        logging.debug("No. of merges=",n_merge)  
    
    complexes = [G.subgraph(complex) for complex in complex_list if len(G.subgraph(complex).nodes())>2 and len(G.subgraph(complex).edges())>=2 ]                    
    perm_lines=np.random.permutation(complexes)        
    fact = inputs['fact'] #0.99
    train_list = [line for line in perm_lines[0:int(round(len(perm_lines)*fact))]]
    test_list = [line for line in perm_lines[int(round(len(perm_lines)*fact)):]]            

    # Start with something that has a biased size distribution !! 
    
    sizes = [len(line) for line in train_list]
    train_mean = np.mean(sizes)
    
    # Transferring some of the smaller complexes to the test list 
    train_list_lower_mean = [line for line in train_list if len(line) < train_mean]
    perc_transfer = inputs['perc_transfer'] #0.3 # You can optimize these parameters !
    to_transfer = train_list_lower_mean[:int(round(len(train_list_lower_mean)*perc_transfer))]
    test_list = test_list + to_transfer
    
    # Now remove from train set 
    for line in to_transfer:
        train_list.remove(line)        
        
    # Finding complexes in train that share an edge with a complex in test 
    com_comp = 10
    while (com_comp !=0): # Do until train and test sets are completely separated
    
        # Removing super huge complexes also (nodes >30 ) from test set         
        test_list = [line for line in test_list if len(line)<30]    
        
        # REMOVE OVERLAP B/W TRAIN AND TEST DATA 
        # Remove complexes from train set sharing two proteins with test set 
        train_rem = []
        for train_line in train_list:
            pres = 0
            for test_line in test_list:
                common = len(list(set(train_line.edges()).intersection(set(test_line.edges))))
                if common >= 1:
                    pres += 1
            if pres > 0:
                train_rem.append(train_line)
                    
        com_comp = len(train_rem)
        logging.debug("No. of train complexes transferred = ",com_comp)
                
        for t_line in train_rem:        
            test_list.append(t_line)
    
        for line in train_rem:
            train_list.remove(line)      
            
    complex_graphs =  train_list  
    test_complex_graphs = test_list

    
    with open(out_comp_nm+"_train_complexes_wo2s.txt","w") as f:
        for line in train_list:
            f.write(sep.join(line)+"\n")
        
    with open(out_comp_nm+"_test_complexes_wo2s.txt","w") as f:
        for line in test_list:
            f.write(sep.join(line)+"\n")   
    with open(out_comp_nm+'_metrics.out',"a") as fid:
        print("No. of training complexes = ",len(train_list),file = fid)
        print("No. of test complexes = ",len(test_list),file = fid)            
        print("Initial train_test split = ",fact,file = fid)    
        print("Percentage of low sizes transferred from train to test = ",perc_transfer,file = fid)                

    return (complex_graphs,test_complex_graphs)
 #   @profile
def read_graphs2(inputs):
    
    logging.info("Initializing...")
    network_path = inputs['dir_nm'] + inputs['netf_nm']
    complex_path = inputs['dir_nm'] + inputs['comf_nm']
    test_complex_path = inputs['dir_nm'] + inputs['comf_test_nm']
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    split_flag = inputs['split_flag']
    use_full = inputs['use_full']
    sep = inputs['sep']
    
    with open(out_comp_nm+'_metrics.out',"w") as fid:
        print("Network:",network_path,file = fid)
    logging.info("Finished Initializing")        
    # Extracting subnetwork of the original PPIN containing proteins in the known complexes
    with open(complex_path) as f:
        data = f.read()
        id_list_train = re.findall(r"[\w']+", data)
        
    with open(test_complex_path) as f:
        data = f.read()
        id_list_test = re.findall(r"[\w']+", data)
    
    prot_list_temp = id_list_train + id_list_test
    
    prot_list = []
    for prot in prot_list_temp:
        if prot not in prot_list:
            prot_list.append(prot)
 
    logging.info("Reading network...")        
    G_full = nx.read_weighted_edgelist(network_path,nodetype=str);

    logging.info("Finished Reading network")  
    
    # Removing self loops 
    G_full.remove_edges_from(nx.selfloop_edges(G_full))

    logging.info("Extracting sub network with proteins from complexes...")          

    if use_full == 1:
        G = G_full
    else:
        G =   G_full.subgraph(prot_list)      

    n_nodes_net = nx.number_of_nodes(G)
    n_edges_net = nx.number_of_edges(G)
    
    with open(out_comp_nm+'_metrics.out',"a") as fid:
        print("No. of nodes in network = ",n_nodes_net,file = fid)
        print("No. of edges in network = ",n_edges_net,file = fid)
    
    logging.info("Reading complexes...")
    
    split_flag = split_flag
    if (split_flag == 1):
        split_train_test_complexes()
    else:
        with open(complex_path) as f:
            complex_list = [line.rstrip('\n').split(sep) for line in f ]# Space separated text only
            
            #complex_graphs = [G.subgraph(complex) for complex in complex_list if G.subgraph(complex).nodes() ]                    
            
            # Keeping only complexes with nodes greater than 2 FOR TOPOLOGICAL FEATURES ONLY CASE 

            complex_graphs = [G.subgraph(complex) for complex in complex_list if len(G.subgraph(complex).nodes())>2 and len(G.subgraph(complex).edges())>=2 ]                    
        
        with open(test_complex_path) as f:
            complex_list = [line.rstrip('\n').split(sep) for line in f ]# Space separated text only
            
            #test_complex_graphs = [G.subgraph(complex) for complex in complex_list if G.subgraph(complex).nodes() ]                    
            
            # Keeping only complexes with nodes greater than 2 FOR TOPOLOGICAL FEATURES ONLY CASE 
            test_complex_graphs = [G.subgraph(complex) for complex in complex_list if len(G.subgraph(complex).nodes())>2 and len(G.subgraph(complex).edges())>=2]                    
        # Change later to extract features of complexes not present in network from CORUM edge predictor 
    logging.info("Finished Reading complexes...")  
    
    known_complexes = complex_graphs + test_complex_graphs
    
    known_complex_nodes_list = [list(complex.nodes()) for complex in known_complexes]
    
    known_complex_nodes = [item for sublist in known_complex_nodes_list for item in sublist]
    prot_list=[]
    for protein in known_complex_nodes:
        if protein not in prot_list:
            prot_list.append(protein)         
    # Plotting size distributions of complexes in the network 
    train_sizes = [len(comp) for comp in complex_graphs]
    test_sizes = [len(comp) for comp in test_complex_graphs]        
    fig = plt.figure()
    sns.distplot(train_sizes,hist=False,label="train")
    sns.distplot(test_sizes,hist=False,label="test")
    plt.savefig(out_comp_nm+'_size_dists_train_test.png')        
    plt.close(fig)        
    
    # T-test for sizes 
    #tval, pval = ttest_ind(train_sizes, test_sizes)
    #print("pval",pval)
    #print("tval",tval)
    
    logging.info("Writing known complexes.")        
    known_complexes = complex_graphs + test_complex_graphs
    with open(out_comp_nm+'_known.out',"w") as fn:
        with open(out_comp_nm+'_known_edges.out',"wb") as f_edges:
            for index in range(len(known_complex_nodes_list)):
                for node in known_complex_nodes_list[index]:
                    fn.write("%s " % node)
                nx.write_weighted_edgelist(known_complexes[index],f_edges)
                fn.write("\n")
                f_edges.write("\n".encode())

    logging.info("Finished writing known complexes.")     
    logging.info("Writing neighbor lists...")     
    
    myGraphdict = nx.to_dict_of_dicts(G)
    
    folNm = inputs['dir_nm']+"/neig_dicts"
    if not os_path.exists(folNm):
        os_mkdir(folNm)

    for node,val in myGraphdict.items():        
        joblib_dump(val,folNm+"/"+node)    
    logging.info("Finished writing neighbor lists.")     
        
    return (known_complex_nodes_list,prot_list,G,complex_graphs,test_complex_graphs)
    
