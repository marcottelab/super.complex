# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 16:25:40 2020

@author: Meg_94
"""
from math import pi
from matplotlib.pyplot import figure as plt_figure, savefig as plt_savefig, close as plt_close
from seaborn import distplot as sns_distplot
from jaccard_coeff import jaccard_coeff
from numpy import mean as np_mean, argmax as np_argmax, argsort as np_argsort, var as np_var, sqrt as sqrt, \
    exp as np_exp
from numpy.random import permutation as rand_perm, choice as rand_choice
from logging import info as logging_info
from networkx import write_weighted_edgelist as nx_write_weighted_edgelist, is_connected as nx_is_connected
from scipy.stats import norm as norm_dist
from convert_humap_ids2names import convert2names


def write_known_comps(known_complex_nodes_list,known_complexes,pref, out_comp_nm, dir_nm):
    total_unique_edges = set()
    if dir_nm == "humap":
        convert2names(known_complex_nodes_list, out_comp_nm + pref + '_known_names.out',
                      known_complexes, out_comp_nm + pref + '_known_edge_names.out')

    with open(out_comp_nm + pref + '_known.out', "w") as fn:
        with open(out_comp_nm + pref + '_known_edges.out', "wb") as f_edges:
            fn_write = fn.write
            f_edges_write = f_edges.write

            for index in range(len(known_complex_nodes_list)):
                for node in known_complex_nodes_list[index]:
                    fn_write("%s " % node)
                nx_write_weighted_edgelist(known_complexes[index], f_edges)
                complex_edges = known_complexes[index].edges()
                for edge in complex_edges:
                    total_unique_edges.add(frozenset([edge[0], edge[1]]))
                fn_write("\n")
                f_edges_write("\n".encode())
    with open(out_comp_nm + pref + '_tot_known_edges_unique.out', "w") as f:
        for edge in total_unique_edges:
            edge = list(edge)
            f.write(edge[0] + "\t" + edge[1] + "\n")

def plot_size_dists(complex_graphs, test_complex_graphs, out_comp_nm):
    train_sizes = [len(comp) for comp in complex_graphs]
    test_sizes = [len(comp) for comp in test_complex_graphs]
    try:
        fig = plt_figure()
        sns_distplot(train_sizes, hist=False, label="train")
        sns_distplot(test_sizes, hist=False, label="test")
        plt_savefig(out_comp_nm + '_size_dists_train_test.png')
        plt_close(fig)
    except:
        print("Cannot plot")


def check_independence(test_list, train_list):
    for line in test_list:
        pres = check_overlap(train_list, line)
        if pres == 1:
            return "Not independent"
    return "Independent"


def nonindependence_num(test_list, train_list):
    ct = 0
    for line in test_list:
        pres = check_overlap(train_list, line)
        if pres == 1:
            ct += 1
    return ct


def check_overlap(train_list, test_line):
    pres = 0
    for train_line in train_list:
        common = len(set(test_line.edges()).intersection(set(train_line.edges)))
        if common >= 1:
            pres = 1
            break
    return pres


def transfer_same_dist(test_list, train_list, com_comp, test_rem):
    if len(test_rem) == 0:
        return test_list, train_list, com_comp

    sizes = [len(line) for line in test_rem]
    mean_test_size = np_mean(sizes)
    sd = sqrt(np_var(sizes))
    if sd != 0:
        test_rem_dist = norm_dist(mean_test_size, sd)
        p_dist = [test_rem_dist.pdf(len(line)) for line in train_list]
        norm_ct = sum(p_dist)
        if norm_ct != 0:
            p_dist = [val / norm_ct for val in p_dist]
        train_rem = rand_choice(train_list, size=com_comp, replace=False, p=p_dist)
    else:
        train_rem = [line for line in train_list if len(line) == mean_test_size][:com_comp]
    test_list = test_list + train_rem
    for line in train_rem:
        train_list.remove(line)
    return test_list, train_list


def transfer_common(test_list, train_list):  # Test to train
    test_rem = []
    test_rem_append = test_rem.append
    com_comp = 0
    for test_line in test_list:
        pres = check_overlap(train_list, test_line)
        if pres == 1:
            com_comp += 1
            test_rem_append(test_line)

    logging_info("No. of transfers = %s", str(com_comp))
    train_list = train_list + test_rem
    for t_line in test_rem:
        test_list.remove(t_line)
    train_list = list(rand_perm(train_list))
    test_list = list(rand_perm(test_list))
    return test_list, train_list, com_comp, test_rem


def transfer_final(test_list, train_list, n_transfers):
    i = 0
    tr_transfers = 0
    while i < len(train_list):
        if tr_transfers == n_transfers:
            break
        line = train_list[i]
        mean_test_size = np_mean([len(line) for line in test_list])
        if len(line) > mean_test_size:
            temp_train_list = list(train_list)
            temp_train_list.remove(line)
            pres = check_overlap(temp_train_list, line)
            if pres == 0:
                train_list.remove(line)
                i -= 1
                test_list.append(line)
                tr_transfers += 1
        i += 1

    return test_list, train_list, n_transfers - tr_transfers


def transfer_common_indep_try(test_list, train_list):
    i = 0
    tr_transfers = 0
    train_rem = []
    while i < len(train_list):
        line = train_list[i]
        pres = check_overlap(test_list, line)
        if pres == 1:
            train_rem.append(train_rem)
            train_list.remove(line)
            i -= 1
            test_list.append(line)
            tr_transfers += 1
        i += 1

    return test_list, train_list, tr_transfers, train_rem


def transfer(test_list, train_list):  # Test to train
    train_list, test_list, com_comp1, test_rem = transfer_common_indep_try(train_list, test_list)
    # Transfer in opposite dir same number
    # Transfer bigger complexes to balance size distributions
    # test_list, train_list = transfer_same_dist(test_list, train_list, com_comp,  test_rem)
    # test_list, train_list, rem_transfers = transfer_final(test_list, train_list, com_comp)

    test_list, train_list, com_comp2, test_rem = transfer_common_indep_try(test_list, train_list)
    print("1", com_comp1)
    print("2", com_comp2)
    extra = com_comp1 - com_comp2
    train_rem = train_list[:extra]
    train_list = train_list[extra:]
    test_list = test_list + train_rem
    return test_list, train_list, com_comp1 - com_comp2


def split_ratio(perm_lines, ratio):
    split_pt_start = float(ratio[0]) / (ratio[0] + ratio[1])
    split_ind = int(round(len(perm_lines) * split_pt_start))
    train_list = [line for line in perm_lines[0:split_ind]]
    test_list = [line for line in perm_lines[split_ind:]]
    n_iters = 0
    print(len(train_list), len(test_list))
    while check_independence(test_list, train_list) != "Independent":
        # Do until train and test sets are completely separated
        test_list, train_list, com_comp = transfer(test_list, train_list)
        print(com_comp)

        print("Niter = ", n_iters)
        n_iters += 1
        if n_iters > 20:
            break
    # Handle case when empty / imbalanced
    print(check_independence(test_list, train_list))
    return train_list, test_list


def split_meth_orig(perm_lines, inputs):
    fact = inputs['fact']  # 0.99
    split_pt = int(round(len(perm_lines) * fact))
    train_list = [line for line in perm_lines[0:split_pt]]
    test_list = [line for line in perm_lines[split_pt:]]
    # Start with something that has a biased size distribution !!

    sizes = [len(line) for line in train_list]
    train_mean = np_mean(sizes)

    # Transferring some of the smaller complexes to the test list
    train_list_lower_mean = [line for line in train_list if len(line) < train_mean]
    perc_transfer = inputs['perc_transfer']  # 0.3 # You can optimize these parameters !
    to_transfer = train_list_lower_mean[:int(round(len(train_list_lower_mean) * perc_transfer))]
    test_list = test_list + to_transfer

    # Now remove from train set
    for line in to_transfer:
        train_list.remove(line)

    # Finding complexes in train that share an edge with a complex in test
    com_comp = 10
    while com_comp != 0:  # Do until train and test sets are completely separated

        # Removing super huge complexes also (nodes >30 ) from test set
        test_list = [line for line in test_list if len(line) < 30]

        # REMOVE OVERLAP B/W TRAIN AND TEST DATA
        # Remove complexes from train set sharing two proteins with test set
        train_rem = []
        train_rem_append = train_rem.append
        com_comp = 0
        for train_line in train_list:
            pres = 0
            for test_line in test_list:
                common = len(set(train_line.edges()).intersection(set(test_line.edges)))
                if common >= 1:
                    pres = 1
                    break
            if pres == 1:
                train_rem_append(train_line)
                com_comp += 1

        logging_info("No. of train complexes transferred = %s", str(com_comp))
        test_list = test_list + train_rem
        for t_line in train_rem:
            train_list.remove(t_line)
    return train_list, test_list


def merge_overlapped(list_comp):
    logging_info("Merging complexes...")

    fin_list = list(list_comp)

    n = len(fin_list)
    if n <= 1:
        return fin_list

    n_changes = 1
    while n_changes != 0:
        if len(fin_list) == 1:
            logging_debug("only one complex")
            break
        n_changes = 0
        ind = 0
        while ind < n:
            if len(fin_list) == 1:
                logging_debug("only one complex")
                break
            else:
                comp = fin_list[ind]
                temp_list = list(fin_list)
                temp_list.remove(comp)
                OS_comp = [jaccard_coeff(set(comp), set(comp2)) for comp2 in temp_list]

                OS_max_ind = int(np_argmax(OS_comp))
                OS_max = OS_comp[OS_max_ind]
                max_over_comp = temp_list[OS_max_ind]
                OS_max_ind_fin = fin_list.index(max_over_comp)

                if OS_max >= 0.6:
                    n_changes += 1
                    # Merge and find score. If score is higher than individual complexes
                    # Keep as new complex
                    merge_comp = comp + max_over_comp

                    fin_list.append(merge_comp)
                    fin_list.remove(comp)
                    fin_list.remove(max_over_comp)
                    n -= 1
                    if OS_max_ind_fin <= ind:
                        ind -= 1
                else:
                    ind += 1
        logging_info("No. of changes = %s", str(n_changes))

    logging_info("Finished filtering complexes.")

    return fin_list


def split_train_test_complexes(inputs, G):
    sep = inputs['sep']
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']

    all_complexes_path = inputs['dir_nm'] + inputs["comf_nm_all"]

    with open(all_complexes_path) as f:
        raw_lines = f.readlines()

    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("No. of raw complexes = ", len(raw_lines), file=fid)
    # Remove Nones and ensure min 3 nodes in complex after that
    lines1 = []
    lines1_append = lines1.append
    for line in raw_lines:
        line_noNone = [item for item in line.rstrip("\r\n").split(sep) if item != "None"]
        if len(line_noNone) > 2:
            lines1_append(line_noNone)

    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("No. of complexes after removing Nones and min 3 members = ", len(lines1), file=fid)

    # Min length and uniqueness checks and connectedness checks in graph
    lines_new = []
    lines_new_append = lines_new.append
    n_small = 0
    n_disconnected = 0
    for comp in lines1:
        Gsub = G.subgraph(comp)
        if not (len(Gsub.nodes()) > 2):
            n_small += 1
            continue
        if not nx_is_connected(Gsub):
            n_disconnected += 1
            continue
        if comp not in lines_new:
            lines_new_append(comp)

    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("No. of complexes after connected, uniqueness and min length checks= ", len(lines_new), file=fid)
        print("No. of small discarded complexes = ", n_small, file=fid)
        print("No. of disconnected complexes = ", n_disconnected, file=fid)
    # Removing redundancy from list:
    # merge complexes with jaccard similarity greater than 0.6 
    lines1 = merge_overlapped(rand_perm(lines_new))

    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("No. of complexes after merging overlapping = ", len(lines1), file=fid)

    complexes = [G.subgraph(comp) for comp in lines1]

    perm_lines = rand_perm(complexes)
    ratio = (70, 30)
    # train_list, test_list = split_meth_orig(perm_lines, inputs)
    train_list, test_list = split_ratio(perm_lines, ratio)
    plot_size_dists(train_list, test_list, out_comp_nm)
    with open(out_comp_nm + "_train_complexes_new.txt", "w") as f:
        for line in train_list:
            f.write(sep.join(line) + "\n")

    with open(out_comp_nm + "_test_complexes_new.txt", "w") as f:
        for line in test_list:
            f.write(sep.join(line) + "\n")
    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("Split ratio = %.3f" % str(float(len(train_list)) / len(test_list)), file=fid)
        # print("Initial train_test split = ", fact, file=fid)
        # print("Percentage of low sizes transferred from train to test = ", perc_transfer, file=fid)
    return train_list, test_list


def preprocess_complexes(complex_path, sep, G):
    with open(complex_path) as f:
        raw_lines = f.readlines()
    # Remove Nones and ensure min 3 nodes in complex after that
    complex_list = []
    lines1_append = complex_list.append
    for line in raw_lines:
        line_noNone = [item for item in line.rstrip("\r\n").split(sep) if item != "None"]
        if len(line_noNone) > 2:
            lines1_append(line_noNone)

    # Keeping only complexes with nodes greater than 2 FOR TOPOLOGICAL FEATURES ONLY CASE
    complex_graphs = []
    complex_graphs_append = complex_graphs.append
    for comp in complex_list:
        Gsub = G.subgraph(comp)
        if not (len(Gsub.nodes()) > 2):
            continue
        if not nx_is_connected(Gsub):
            continue
        if Gsub not in complex_graphs:
            complex_graphs_append(Gsub)

    return complex_graphs


def read_complexes(inputs, G):
    split_flag = inputs['split_flag']
    sep = inputs['sep']
    complex_path = inputs['dir_nm'] + inputs['comf_nm']
    test_complex_path = inputs['dir_nm'] + inputs['comf_test_nm']
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']

    logging_info("Reading complexes...")

    split_flag = split_flag
    if split_flag == 1:
        complex_graphs, test_complex_graphs = split_train_test_complexes(inputs, G)
    else:
        complex_graphs = preprocess_complexes(complex_path, sep, G)
        test_complex_graphs = preprocess_complexes(test_complex_path, sep, G)
        # Change later to extract features of complexes not present in network from CORUM edge predictor
    logging_info("Finished reading complexes...")

    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("No. of training complexes = ", len(complex_graphs), file=fid)
        print("No. of test complexes = ", len(test_complex_graphs), file=fid)
   
   # Plotting size distributions of complexes in the network 
    plot_size_dists(complex_graphs, test_complex_graphs, out_comp_nm)
    
    known_complexes = complex_graphs + test_complex_graphs

    known_complex_nodes_list = [comp.nodes() for comp in known_complexes]
    test_known_complex_nodes_list = [comp.nodes() for comp in test_complex_graphs]
    train_known_complex_nodes_list = [comp.nodes() for comp in complex_graphs]

    logging_info("Writing known complexes.")
    write_known_comps(known_complex_nodes_list,known_complexes,"", out_comp_nm, inputs['dir_nm'])
    write_known_comps(train_known_complex_nodes_list,complex_graphs,"_train", out_comp_nm, inputs['dir_nm'])
    write_known_comps(test_known_complex_nodes_list,test_complex_graphs,"_test", out_comp_nm, inputs['dir_nm'])
    logging_info("Finished writing known complexes.")
    
    known_complex_nodes = [item for sublist in known_complex_nodes_list for item in sublist]
    prot_list = set(known_complex_nodes)
    
    return known_complex_nodes_list, complex_graphs, test_complex_graphs, prot_list, test_known_complex_nodes_list, train_known_complex_nodes_list
