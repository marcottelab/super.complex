# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 10:18:32 2020

@author: Meg_94
"""
from numpy import vstack as np_vstack, hstack as np_hstack, where as np_where, array as np_array, delete as np_delete, percentile as np_percentile
from pandas import DataFrame as pd_DataFrame, read_csv as pd_read_csv
from create_feat_mat import create_feat_mat
from construct_neg_comps import construct_neg_comps
from networkx import write_weighted_edgelist as nx_write_weighted_edgelist
from logging import info as logging_info
from seaborn import boxplot as sns_boxplot
from math import ceil as math_ceil
from pickle import dump as pickle_dump

import matplotlib.pyplot as plt

def extract_features(out_comp_nm, split_type, max_size, inputs, G_nodes, feat_list, X_pos, X_allpos, n_allpos,
                     sizes):
    n_pos = len(X_pos)

    folNm = inputs['dir_nm']+ inputs['graph_files_dir'] + "/neig_dicts"
    dims = X_pos.shape
    n_feats = dims[1]
    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("No. of " + split_type + " features = ", n_feats, file=fid)
        print("No. of " + split_type + " positive complexes = ", n_pos, file=fid)

    logging_info("Constructing " + split_type + " negative complexes...")
    if "neg_sample_method" not in inputs:
        inputs["neg_sample_method"] = "uniform"
    neg_comp_list = construct_neg_comps(max_size, n_pos, inputs['scale_factor'], G_nodes, sizes,inputs["neg_sample_method"],folNm)
    logging_info("Finished constructing " + split_type + " negative complexes")

    X_neg = create_feat_mat(neg_comp_list, n_feats)

    X_neg, neg_comp_list, n_neg = remove_same_rows(n_allpos, X_neg, X_allpos, neg_comp_list)

    # print(n_neg)
    # HHANDLE CASE WHEN n_neg = 0 !!!!!
    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("No. of " + split_type + " negative complexes = ", n_neg, file=fid)

    write_neg2out(out_comp_nm + '_neg_' + split_type + '.out', out_comp_nm + '_neg_' + split_type + '_edges.out',
                  neg_comp_list)

    X = np_vstack((X_pos, X_neg))

    y_pos = [1] * n_pos
    y_neg = [0] * n_neg
    y = y_pos + y_neg
    y = np_array(y)
    y_pos = np_array(y_pos)
    y_neg = np_array(y_neg)

    # Writing raw training data to csv in tpot format
    write2csv_tpot(X, y, out_comp_nm + "_" + split_type + "_dat.csv", feat_list)
    return y, X, X_pos, y_pos, X_neg, y_neg


def write_neg2out(neg, neg_edges, neg_test_comp):
    with open(neg, "w") as fn:
        with open(neg_edges, "wb") as f_edges:
            for index in range(len(neg_test_comp)):
                for node in neg_test_comp[index].nodes():
                    fn.write("%s " % node)
                nx_write_weighted_edgelist(neg_test_comp[index], f_edges)
                fn.write("\n")
                f_edges.write("\n".encode())


def write2csv_tpot(X, y, outfName, feat_list):
    dat = np_hstack((X, y[:, None]))
    df = pd_DataFrame(dat)
    df.to_csv(path_or_buf=outfName, index=False, header=feat_list)


def remove_same_rows(n_pos, X_neg, X_pos, neg_comp_list):
    # Removing negative feature rows that exactly match any row in positives
    cout = 0
    for ind in range(n_pos):
        matching_inds = np_where((X_neg == X_pos[ind]).all(axis=1))
        X_neg = np_delete(X_neg, matching_inds, axis=0)
        for index in sorted(list(matching_inds[0]), reverse=True):
            cout += 1
            del neg_comp_list[index]
    print("No. of negs removed due to same feature vector = ", cout)
    n_neg = len(X_neg)
    return X_neg, neg_comp_list, n_neg


def read_from_csv(fileName):
    df_full = pd_read_csv(fileName)
    y = df_full.pop('complex')
    X = df_full

    neg_start_ind = y[y == 0].index[0]
    X_pos = X.iloc[0:neg_start_ind]
    y_pos = y[0:neg_start_ind]
    X_neg = X.iloc[neg_start_ind:]
    y_neg = y[neg_start_ind:]

    return y, X, X_pos, y_pos, X_neg, y_neg


def feature_extract(inputs, complex_graphs, test_complex_graphs, G):
    G_nodes = G.nodes()
    n_feats = inputs['feats']
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    mode = inputs['mode']
    # mode = "non_gen" # Change to gen if you want to generate matrices

    # n_pos = len(complex_graphs)
    sizes = [len(comp) for comp in complex_graphs]

    # get quartiles
    q1 = np_percentile(sizes, 25)
    q3 = np_percentile(sizes, 75)
    max_wo_outliers = math_ceil(q3 + 4.5*(q3-q1)) # Maximum after removing outliers

    max_size_train = max(sizes)
    recommended_max_size = min(max_size_train,max_wo_outliers)
    
    max_sizeF = inputs['dir_nm'] + inputs['train_test_files_dir']+ "/res_max_size_search"
    with open(max_sizeF, 'wb') as f:
        pickle_dump(recommended_max_size, f)
    
    # n_pos_test = len(test_complex_graphs)
    sizes_test = [len(comp) for comp in test_complex_graphs]
    max_size_test = max(sizes_test)

    fig = plt.figure()
    # Plot box plot of sizes to know the outliers (for setting step size in sampling)
    sns_boxplot(sizes)
    plt.xlabel("Size")
    plt.title("Size distribution of training complexes")
    plt.savefig(out_comp_nm + "_known_train_size_dist_box_plot")
    plt.close(fig)

    fig = plt.figure()
    # Plot box plot of sizes to know the outliers (for setting step size in sampling)
    sns_boxplot(sizes + sizes_test)
    plt.xlabel("Size")
    plt.title("Size distribution of known complexes")
    plt.savefig(out_comp_nm + "_known_size_dist_box_plot")
    plt.close(fig)

    if inputs['model_type'] == "tpot" and mode == "non_gen":  # CHANGE X_POS, Y_POS later !!!!
        logging_info("Reading labeled feature matrix from file...")
        # Read X,y from csv file

        y, X, X_pos, y_pos, X_neg, y_neg = read_from_csv(inputs['train_feat_mat'])

        y_test, X_test, X_pos_test, y_pos_test, X_neg_test, y_neg_test = read_from_csv(inputs['test_feat_mat'])

        logging_info("Finished reading feature matrix")
    else:

        logging_info("Feature extraction...")

        feat_list = ["dens", "nodes", "degree_max", "degree_mean", "degree_median", "degree_var", "CC_max", "CC_mean",
                     "CC_var",
                     "edge_wt_mean", "edge_wt_max", "edge_wt_var", "DC_mean", "DC_var", "DC_max", "sv1", "sv2", "sv3",
                     "complex"]

        X_pos = create_feat_mat(complex_graphs, n_feats)
        X_pos_test = create_feat_mat(test_complex_graphs, n_feats)

        X_allpos = np_vstack((X_pos, X_pos_test))
        n_allpos = len(X_allpos)
        y, X, X_pos, y_pos, X_neg, y_neg = extract_features(out_comp_nm, 'train', max_size_train, inputs, G_nodes,
                                                            feat_list, X_pos, X_allpos, n_allpos, sizes)
        y_test, X_test, X_pos_test, y_pos_test, X_neg_test, y_neg_test = extract_features(out_comp_nm, 'test',
                                                                                          max_size_test, inputs, G_nodes,
                                                                                          feat_list, X_pos_test,
                                                                                          X_allpos,
                                                                                          n_allpos, sizes_test)

        logging_info("Finished Feature extraction")
    return max_size_train, max_size_test, X_pos_test, X_neg_test, X_test, y_test, X_pos, y_pos, X, y, X_neg, y_neg
