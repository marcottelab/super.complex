# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 17:28:40 2020

@author: Meghana
"""

import networkx as nx
import numpy as np
import logging

def random_walk2(seed_node,n,G): # in graph G 
    neg_comp = nx.Graph()
    neg_comp.add_node(seed_node)
    node_num = 1
    pres_node = seed_node
    while node_num<n:
        neig_list = G[pres_node]
        if not neig_list:
            logging.debug("No neighbor")
            break
        if(len(neig_list) != 1):
        	new_node = np.random.choice(neig_list)
        else:
            temp = neig_list.keys()
            new_node = list(temp)[0]
        wt = neig_list[new_node]
        wt_edge = wt['weight']
        neg_comp.add_edge(pres_node,new_node,weight=wt_edge)
        pres_node = new_node
        node_num = node_num + 1
	    #print(pres_node) 
    return neg_comp