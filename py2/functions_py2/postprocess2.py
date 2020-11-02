# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 17:22:47 2020

@author: Meghana
"""
from logging import info as logging_info,debug as logging_debug
import networkx as nx
from jaccard_coeff2 import jaccard_coeff2
from numpy import argmax as np_argmax, argsort as np_argsort
from sample2 import get_score


def filter_overlapped(list_comp,inputs):
    logging_info("Filtering complexes...")        

    # Sort by size 
    sizes = [nx.number_of_nodes(comp) for comp in list_comp]
    
    sorted_ind = np_argsort(sizes) # ascending order.
    
    list_comp = [list_comp[i] for i in sorted_ind]

    fin_list = list(list_comp)
    list_comp2 = list(list_comp)
    
    #print(len(list_comp))
    # Ensure ascending order 
    for comp in list_comp:
        OS_comp = []            
        list_comp2.remove(comp)
        
        if(len(list_comp2)):
            for comp2 in list_comp2:
                
                Overlap = jaccard_coeff2(comp.nodes(),comp2.nodes())
                
                OS_comp.append(Overlap)
            
            OS_max = max(OS_comp)
            #print(OS_max)
            
            if( OS_max > inputs['over_t']):
                fin_list.remove(comp)
    
    logging_info("Finished filtering complexes.")        

    return fin_list

def merge_filter_overlapped_score(list_comp,model,scaler,inputs,G):
    logging_info("Filtering complexes...")        

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
                    Overlap = jaccard_coeff2(comp.nodes(),comp2.nodes())
                    #print("Overlap:",Overlap)
                    OS_comp.append(Overlap)
                    
                OS_max_ind = np_argmax(OS_comp)                
                OS_max = OS_comp[OS_max_ind]
                max_over_comp = temp_list[OS_max_ind]
                OS_max_ind_fin = fin_list.index(max_over_comp)
                
                
                #print("OSmax=",OS_max)

                if( OS_max >= inputs['over_t']):
                    n_changes += 1                    
                    # Merge and find score. If score is higher than individual complexes 
                    # Keep as new complex                                        
                    merge_comp_nodes = list(set(comp.nodes()).union(set(max_over_comp.nodes())))
                    merge_comp = nx.Graph(G.subgraph(merge_comp_nodes),comp_score=0)
                    
                    (score_merge,comp_bool) = get_score(merge_comp,model,scaler,inputs)
                    merge_comp.graph['comp_score'] = score_merge
                    sc1 = comp.graph['comp_score']
                    sc2 = max_over_comp.graph['comp_score']                    
                    if(score_merge > sc1 and score_merge > sc2):
                        fin_list.append(merge_comp)
                        fin_list.remove(comp)
                        fin_list.remove(max_over_comp)
                        if OS_max_ind_fin <= ind:
                            ind -= 1                                

                    # Otherwise: remove lower scoring complex
                    elif sc1 <= sc2:
                        fin_list.remove(comp)
                    else:
                        fin_list.remove(max_over_comp)
                        if OS_max_ind_fin > ind:
                            ind += 1
                else:
                    ind += 1
                    
            else:
                logging_debug("temp_list is empty")
                ind += 1
            n = len(fin_list)

        logging_info("No. of changes = %s", str(n_changes))
        if n_changes == 0:
            break
    
    logging_info("Finished filtering complexes.")          

    return fin_list

def postprocess2(pred_comp_list,model,scaler,inputs,G):
    # Removing complexes with only two nodes 
    pred_comp_list = [comp for comp in pred_comp_list if len(comp)>2]
    pred_comp_node_list = [list(comp.nodes()) for comp in pred_comp_list]
    # Finding unique complexes     

    fin_list=[]
    indices=[]
    for j,comp in enumerate(pred_comp_node_list):
        if comp not in fin_list:
            fin_list.append(comp)        
            indices.append(j)

    fin_list_graphs = [pred_comp_list[i] for i in indices]
    

    logging_info("Finished sampling complexes.")  
 
    # Filtering complexes with high overlap with bigger complexes 
    fin_list_graphs = merge_filter_overlapped_score(fin_list_graphs,model,scaler,inputs,G)
    
    logging_info("Writing predicted complexes.")        
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    
    with open(out_comp_nm+'_pred.out',"w") as fn:
        with open(out_comp_nm+'_pred_edges.out',"wb") as f_edges:
            for index in range(len(fin_list_graphs)):
                for node in fin_list_graphs[index].nodes():
                    fn.write("%s " % node)     
                
                fn.write("%.3f" % fin_list_graphs[index].graph['comp_score'])                        
                nx.write_weighted_edgelist(fin_list_graphs[index],f_edges)
                fn.write("\n")
                f_edges.write("\n".encode())

    logging_info("Finished writing predicted complexes.")    
    

    return fin_list_graphs    