# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 16:47:00 2019

@author: Meghana
"""

from networkx import clustering as nx_clustering, average_neighbor_degree as nx_average_neighbor_degree, \
    to_numpy_matrix as nx_to_numpy_matrix, density as nx_density, number_of_nodes as nx_number_of_nodes
from numpy import mean as np_mean, var as np_var, vstack as np_vstack, median as np_median
from numpy.linalg import svd as np_linalg_svd


def create_feat_mat_1(graph):
    CCs = list(nx_clustering(graph).values())

    DCs = list(nx_average_neighbor_degree(graph).values())

    degrees = [tup[1] for tup in graph.degree()]

    edge_wts = [tup[2] for tup in graph.edges.data('weight')]

    A_mat = nx_to_numpy_matrix(graph)
    svs = np_linalg_svd(A_mat, full_matrices=False, compute_uv=False)

    if len(svs) >= 3:
        sv1 = svs[0]
        sv2 = svs[1]
        sv3 = svs[2]
    elif len(svs) >= 2:
        sv1 = svs[0]
        sv2 = svs[1]
        sv3 = 0
    else:
        sv1 = svs[0]
        sv2 = sv3 = 0

    feat_mat = np_vstack((nx_density(graph), nx_number_of_nodes(graph), max(degrees),
                          np_mean(degrees), np_median(degrees), np_var(degrees),
                          max(CCs), np_mean(CCs), np_var(CCs), np_mean(edge_wts),
                          max(edge_wts), np_var(edge_wts), np_mean(DCs), np_var(DCs), max(DCs), sv1, sv2, sv3)).T

    return feat_mat
