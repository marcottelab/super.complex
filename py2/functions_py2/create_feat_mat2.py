# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 17:02:04 2020

@author: Meghana
"""
import networkx as nx
import numpy as np

def create_feat_mat2(graph_list,n_feats):
    dens_pos = [nx.density(graph) for graph in graph_list]
    nodes_pos = [nx.number_of_nodes(graph) for graph in graph_list]
    
    #CC statistics - mean and max  - faster to use a big loop mostly
    CC_mean=[]
    CC_max = []
    CC_var = []   
    # Degree correlation - avg degree of the neighborhood     
    DC_mean=[]
    DC_max = []
    DC_var = []       
    #Degree statistics 
    degree_mean=[]
    degree_max = []
    degree_median = []
    degree_var = []    
    # Edge weight statistics 
    edge_wt_mean=[]
    edge_wt_max = []
    edge_wt_var=[]    
    # First 3 singular values 
    sv1=[]
    sv2=[]
    sv3=[]    
    for graph in graph_list:
        
        CCs = list(nx.clustering(graph).values())
        CC_max.append(max(CCs))
        CC_mean.append(np.mean(CCs))
        CC_var.append(np.var(CCs))
        
        DCs = list(nx.average_neighbor_degree(graph).values())
        DC_max.append(max(DCs))
        DC_mean.append(np.mean(DCs))
        DC_var.append(np.var(DCs))        
        
        degrees = [tup[1] for tup in graph.degree()]
        degree_mean.append(np.mean(degrees))
        degree_median.append(np.median(degrees))
        degree_max.append(max(degrees))
        degree_var.append(np.var(degrees))
    
        edge_wts = [tup[2] for tup in graph.edges.data('weight')]
        edge_wt_mean.append(np.mean(edge_wts))
        edge_wt_var.append(np.var(edge_wts))            
        edge_wt_max.append(max(edge_wts))        
        
        A_mat = nx.to_numpy_matrix(graph)
        svs = np.linalg.svd(A_mat,full_matrices=False,compute_uv=False)
        
        if(len(svs) >=3):
            sv1.append(svs[0])
            sv2.append(svs[1])
            sv3.append(svs[2])
        elif(len(svs) >=2):
            sv1.append(svs[0])
            sv2.append(svs[1])
            sv3.append(0)                
        else:
            sv1.append(svs[0])
            sv2.append(0)
            sv3.append(0)          
            
    #feat_mat = np.vstack((dens_pos,nodes_pos)).T
    feat_mat = np.vstack((dens_pos,nodes_pos,degree_max,degree_mean,degree_median,degree_var,CC_max,CC_mean,CC_var,edge_wt_mean,edge_wt_max,edge_wt_var,DC_mean,DC_var,DC_max,sv1,sv2,sv3)).T
    
    #feat_mat = np.vstack((dens_pos,nodes_pos,degree_max,degree_mean,degree_var,CC_max,CC_mean,CC_var,edge_wt_mean,edge_wt_max,edge_wt_var,DC_mean,DC_var,sv2,sv3)).T
    if(n_feats == 1):
        feat_mat = np.array(dens_pos).reshape(-1,1)
        
    return feat_mat
 