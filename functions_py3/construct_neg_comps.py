# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 10:35:30 2020

@author: Meg_94
"""
from numpy import ceil as np_ceil
from random_graph import random_graph


def uniform_dist(min_size, max_size, n_pos, scale_fac, G_nodes, folNm):
    neg_comp_list_init = []
    n_each = int(np_ceil(float(n_pos) / (max_size - min_size + 1)))  # Uniform distribution
    n_each = int(np_ceil(n_each * scale_fac))

    neg_comp_list_init_append = neg_comp_list_init.append
    # Negative complexes - random walks from random seeds with incremental steps from min_size to max_size
    for n in range(min_size, max_size + 1):
        for k in range(n_each):
            # Pick a random seed node
            neg_comp_list_init_append(random_graph(n, G_nodes, folNm))
    return neg_comp_list_init


def same_dist(scale_fac, G_nodes, sizes, folNm):
    neg_comp_list_init = []
    ct_sizes = {}
    for sz in sizes:
        if sz in ct_sizes:
            ct_sizes[sz] += 1
        else:
            ct_sizes[sz] = 1
    neg_comp_list_init_append = neg_comp_list_init.append
    wrong_size_ct = 0
    for sz in ct_sizes:
        nk = int(np_ceil(ct_sizes[sz] * scale_fac))
        for k in range(nk):
            neg_comp = random_graph(sz, G_nodes, folNm)
            neg_comp_sz = len(neg_comp.nodes())
            if neg_comp_sz in ct_sizes:
                if ct_sizes[neg_comp_sz] > 0:
                    ct_sizes[neg_comp_sz] -= 1
                else:
                   wrong_size_ct += 1
            else:
                wrong_size_ct += 1
            neg_comp_list_init_append(neg_comp)

    print("Wrong sizes count ", wrong_size_ct)
    return neg_comp_list_init


def same_dist_exact(scale_fac, G_nodes, sizes, folNm):
    thres_ct = round(len(sizes)*scale_fac)
    neg_comp_dict = {}
    for sz in sizes:
        if sz in neg_comp_dict:
            neg_comp_dict[sz]['count'] += 1
        else:
            neg_comp_dict[sz]= {}
            neg_comp_dict[sz]['count'] = 1
            neg_comp_dict[sz]['complexes'] = []

    for sz in neg_comp_dict:
        for k in range(neg_comp_dict[sz]['count']):
            graph_rand = random_graph(sz, G_nodes, folNm)
            sz_graph_rand = len(graph_rand.nodes())
            if sz_graph_rand in neg_comp_dict:
                if neg_comp_dict[sz]['count'] > 0:
                    neg_comp_dict[sz]['complexes'].append(graph_rand)
                    neg_comp_dict[sz]['count'] -= 1
    return neg_comp_list_init


def construct_neg_comps(max_size, n_pos, scale_fac, G_nodes, sizes, dist="uniform", folNm = ""):
    neg_comp_list_init = []
    print("No. of positive comps = ", len(sizes))
    if dist == "same":
        neg_comp_list_init = same_dist(scale_fac, G_nodes, sizes, folNm)
    else:
        min_size = min(sizes)
        neg_comp_list_init = uniform_dist(min_size, max_size, n_pos, scale_fac, G_nodes, folNm)

    print("No. of negative complexes = ", len(neg_comp_list_init))
    # Remove random walks with 1 or 2 nodes
    neg_comp_list = [walk for walk in neg_comp_list_init if len(walk) >= 3]
    print("No. of negative complexes with sizes greater than 2 = ", len(neg_comp_list))
    return neg_comp_list
