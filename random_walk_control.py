# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 16:43:00 2020

@author: Meg_94
"""

import networkx as nx
from random_walk import random_walk
from jaccard_coeff import jaccard_coeff
import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed
import multiprocessing 


def merge_filter_overlapped(list_comp,inputs,G):
    print("Filtering complexes...")        

    fin_list = list(list_comp)
     
    n = len(fin_list)
    while True:
        n_changes = 0
        ind = 0
        while ind < n:
            comp = fin_list[ind]
            temp_list = list(fin_list)
            OS_comp = []            
            temp_list.remove(comp)
            
            if(len(temp_list)):
                for comp2 in temp_list:
                    
                    #print("1:",comp.nodes())
                    #print("2:",comp2.nodes())
                    Overlap = jaccard_coeff(comp.nodes(),comp2.nodes())
                    #print("Overlap:",Overlap)
                    OS_comp.append(Overlap)
                
                OS_max = max(OS_comp)
                OS_max_ind = np.argmax(OS_comp)
                max_over_comp = fin_list[OS_max_ind]
                
                #print("OSmax=",OS_max)

                if( OS_max >= inputs['over_t']):
                    n_changes += 1                    
                    # Merge 
                    merge_comp_nodes = list(set(comp.nodes()).union(set(max_over_comp.nodes())))
                    merge_comp = nx.Graph(G.subgraph(merge_comp_nodes),comp_score=0)
                                     
                    fin_list.append(merge_comp)
                    fin_list.remove(comp)
                    del fin_list[OS_max_ind]
                    if OS_max_ind <= ind:
                        ind -= 1                                
                else:
                    ind += 1
                    
            else:
                print("temp_list is empty")
                ind += 1
            n = len(fin_list)

        print(n_changes)
        if n_changes == 0:
            break
    
    print("Finished filtering complexes.")          

    return fin_list

def control(G,inputs={},n=50):
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    seed_nodes = list(G.nodes()) # all nodes
    pred_comp_list = []

    num_cores = multiprocessing.cpu_count()
    print("Sampling complexes by random walks...")
    pred_comp_list = Parallel(n_jobs=num_cores,prefer="threads")(delayed(random_walk)(node,n,G) for node in tqdm(seed_nodes))
    '''
    for node in seed_nodes:
        rand_walk = random_walk(node,n,G)
        pred_comp_list.append(rand_walk)
    '''    
    # Removing complexes with only two nodes 
    pred_comp_list = [comp for comp in pred_comp_list if len(comp)>2]
    pred_comp_node_list = [list(comp.nodes()) for comp in pred_comp_list]
    # Finding unique ccomplexes     

    
    fin_list=[]
    indices=[]
    for j,comp in enumerate(pred_comp_node_list):
        if comp not in fin_list:
            fin_list.append(comp)        
            indices.append(j)

    fin_list_graphs = [pred_comp_list[i] for i in indices]
    

    print("Finished sampling complexes.")  
 
    # Filtering complexes with high overlap with bigger complexes 
    fin_list_graphs = merge_filter_overlapped(fin_list_graphs,inputs,G)
    
    print("Writing predicted complexes.")        
    
    with open(out_comp_nm+'_random_walk_control.out',"w") as fn:
        with open(out_comp_nm+'_random_walk_control_edges.out',"wb") as f_edges:
            for index in range(len(fin_list_graphs)):
                for node in fin_list_graphs[index].nodes():
                    fn.write("%s " % node)     
                
                nx.write_weighted_edgelist(fin_list_graphs[index],f_edges)
                fn.write("\n")
                f_edges.write("\n".encode())

    print("Finished writing predicted complexes.")    
    
    complex_sizes = [len(gr.nodes()) for gr in fin_list_graphs]
    biggest_complex_size = max(complex_sizes)
    avg_complex_size = np.mean(complex_sizes)
    smallest_complex_size = min(complex_sizes)
    
    with open(out_comp_nm+'_metrics.out',"a") as fid:                
       print("No. of cores = ",num_cores,file = fid)  
       print("No. of steps taken = ",n,file = fid)  
       print("Predicted biggest_complex_size = ",biggest_complex_size,file = fid)
       print("Predicted avg_complex_size = ",avg_complex_size,file = fid) 
       print("Predicted smallest_complex_size = ",smallest_complex_size,file = fid) 

    return fin_list_graphs

   
    
        
        
