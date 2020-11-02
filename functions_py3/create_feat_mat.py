# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 16:47:00 2019

@author: Meghana
"""

from networkx import clustering as nx_clustering, average_neighbor_degree as nx_average_neighbor_degree, \
    to_numpy_matrix as nx_to_numpy_matrix, density as nx_density, number_of_nodes as nx_number_of_nodes
from numpy import mean as np_mean, var as np_var, array as np_array, vstack as np_vstack, median as np_median
from numpy.linalg import svd as np_linalg_svd


def create_feat_mat(graph_list, n_feats):
    dens_pos = [nx_density(graph) for graph in graph_list]
    nodes_pos = [nx_number_of_nodes(graph) for graph in graph_list]

    # CC statistics - mean and max  - faster to use a big loop mostly
    CC_mean = []
    CC_mean_append = CC_mean.append
    CC_max = []
    CC_max_append = CC_max.append
    CC_var = []
    CC_var_append = CC_var.append
    # Degree correlation - avg degree of the neighborhood     
    DC_mean = []
    DC_mean_append = DC_mean.append
    DC_max = []
    DC_max_append = DC_max.append
    DC_var = []
    DC_var_append = DC_var.append
    # Degree statistics
    degree_mean = []
    degree_mean_append = degree_mean.append
    degree_max = []
    degree_max_append = degree_max.append
    degree_median = []
    degree_median_append = degree_median.append
    degree_var = []
    degree_var_append = degree_var.append
    # Edge weight statistics 
    edge_wt_mean = []
    edge_wt_mean_append = edge_wt_mean.append
    edge_wt_max = []
    edge_wt_max_append = edge_wt_max.append
    edge_wt_var = []
    edge_wt_var_append = edge_wt_var.append
    # First 3 singular values 
    sv1 = []
    sv1_append = sv1.append
    sv2 = []
    sv2_append = sv2.append
    sv3 = []
    sv3_append = sv3.append
    for graph in graph_list:

        CCs = list(nx_clustering(graph).values())
        CC_max_append(max(CCs))
        CC_mean_append(np_mean(CCs))
        CC_var_append(np_var(CCs))

        DCs = list(nx_average_neighbor_degree(graph).values())
        DC_max_append(max(DCs))
        DC_mean_append(np_mean(DCs))
        DC_var_append(np_var(DCs))

        degrees = [tup[1] for tup in graph.degree()]
        degree_mean_append(np_mean(degrees))
        degree_median_append(np_median(degrees))
        degree_max_append(max(degrees))
        degree_var_append(np_var(degrees))

        edge_wts = [tup[2] for tup in graph.edges.data('weight')]
        edge_wt_mean_append(np_mean(edge_wts))
        edge_wt_var_append(np_var(edge_wts))
        edge_wt_max_append(max(edge_wts))

        A_mat = nx_to_numpy_matrix(graph)
        svs = np_linalg_svd(A_mat, full_matrices=False, compute_uv=False)

        if len(svs) >= 3:
            sv1_append(svs[0])
            sv2_append(svs[1])
            sv3_append(svs[2])
        elif len(svs) >= 2:
            sv1_append(svs[0])
            sv2_append(svs[1])
            sv3_append(0)
        else:
            sv1_append(svs[0])
            sv2_append(0)
            sv3_append(0)

    feat_mat = np_vstack((dens_pos, nodes_pos, degree_max, degree_mean, degree_median, degree_var, CC_max, CC_mean,
                          CC_var, edge_wt_mean, edge_wt_max, edge_wt_var, DC_mean, DC_var, DC_max, sv1, sv2, sv3)).T

    if n_feats == 1:
        feat_mat = np_array(dens_pos).reshape(-1, 1)

    return feat_mat
