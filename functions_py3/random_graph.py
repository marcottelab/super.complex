# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 10:38:13 2020

@author: Meg_94
"""
from numpy.random import choice as rand_choice
from networkx import Graph as nx_Graph
from logging import debug as logging_debug
from pickle import load as pickle_load


# More efficient to have a global set of random graphs, and sample again for missing number of steps. Start from a high number of steps and go lower, recycling the smaller random walks.
def random_walk(seed_node, n, folNm):  # in graph G
    neg_comp = nx_Graph()
    neg_comp.add_node(seed_node)
    node_num = 1
    pres_node = seed_node
    extra = 10
    while node_num < n + extra:
        with open(folNm + "/" + pres_node, 'rb') as f:
            neig_list = pickle_load(f)
        if not neig_list:
            logging_debug("No neighbours")
            break
        if len(neig_list) != 1:
            new_node = rand_choice(list(neig_list.keys()))
        else:
            new_node = list(neig_list.keys())[0]
        wt = neig_list[new_node]
        wt_edge = wt['weight']
        neg_comp.add_edge(pres_node, new_node, weight=wt_edge)
        if len(neg_comp.nodes()) == n:
            break
        pres_node = new_node
        node_num = node_num + 1
        # print(pres_node)
    return neg_comp


def random_graph(n, G_nodes, folNm):  # random seed
    attempt_thres = 10
    attempt_ct = 0
    while attempt_ct < attempt_thres:
        seed_node_rand = rand_choice(G_nodes)
        rand_graph = random_walk(seed_node_rand, n, folNm)
        if len(rand_graph.nodes()) == n:
            break
        attempt_ct += 1
    return rand_graph
