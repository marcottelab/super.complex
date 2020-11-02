# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 12:58:38 2020

@author: Meg_94
"""
import numpy as np
from random_graph2 import random_graph2
def construct_neg_comps2(min_size,max_size,n_pos,inputs,G):
    n_each = int(np.ceil( n_pos/(max_size - min_size + 1))) # Uniform distribution 
    n_each=int(np.ceil(n_each*inputs['scale_factor']))    
        
    # Negative complexes - random walks from random seeds with incremental steps from min_size to max_size
    neg_comp_list_init = []
    for n in range(min_size,max_size+1): # Generate with same distribution as sizes later
        for k in range(n_each):
            # Pick a random seed node 
            neg_comp = random_graph2(n,G)
            neg_comp_list_init.append(neg_comp)  
        
    # Remove random walks with 1 or 2 nodes 
    neg_comp_list = []
    for walk in neg_comp_list_init:
        if len(walk) >= 3:
            neg_comp_list.append(walk)
    return neg_comp_list 
