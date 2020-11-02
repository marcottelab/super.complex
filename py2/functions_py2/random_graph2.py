# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 12:55:12 2020

@author: Meg_94
"""
import numpy as np
from random_walk2 import random_walk2
def random_graph2(n,G): # random seed 
    seed_node_rand = np.random.choice(G.nodes())
    rand_graph = random_walk2(seed_node_rand,n,G)
    return rand_graph