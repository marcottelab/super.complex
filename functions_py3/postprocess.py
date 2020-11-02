# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 17:22:47 2020

@author: Meghana
"""
from logging import info as logging_info, debug as logging_debug
from networkx import number_of_nodes as nx_number_of_nodes, write_weighted_edgelist as nx_write_weighted_edgelist, \
    Graph as nx_Graph
from jaccard_coeff import jaccard_coeff
from numpy import argmax as np_argmax, argsort as np_argsort
from sample import get_score
from pickle import load as pickle_load
from convert_humap_ids2names import convert2names_wscores

def filter_overlapped(list_comp, inputs):
    logging_info("Filtering complexes...")

    # Sort by size 
    sizes = [nx_number_of_nodes(comp) for comp in list_comp]

    sorted_ind = np_argsort(sizes)  # ascending order.

    list_comp = [list_comp[i] for i in sorted_ind]

    fin_list = list(list_comp)
    list_comp2 = list(list_comp)

    # print(len(list_comp))
    # Ensure ascending order 
    for comp in list_comp:
        OS_comp = []
        list_comp2.remove(comp)

        if len(list_comp2):
            for comp2 in list_comp2:
                Overlap = jaccard_coeff(comp.nodes(), comp2.nodes())

                OS_comp.append(Overlap)

            OS_max = max(OS_comp)
            # print(OS_max)

            if OS_max > inputs['over_t']:
                fin_list.remove(comp)

    logging_info("Finished filtering complexes.")

    return fin_list


def merge_filter_overlapped_score(list_comp, model, scaler, inputs, G):
    logging_info("Filtering complexes...")

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
                OS_comp = [jaccard_coeff(comp[0], comp2[0]) for comp2 in temp_list]

                OS_max_ind = int(np_argmax(OS_comp))
                OS_max = OS_comp[OS_max_ind]
                max_over_comp = temp_list[OS_max_ind]
                OS_max_ind_fin = fin_list.index(max_over_comp)

                if OS_max >= inputs['over_t']:
                    n_changes += 1
                    n -= 1
                    # Merge and find score. If score is higher than individual complexes 
                    # Keep as new complex                                        
                    merge_comp_nodes = comp[0].union(max_over_comp[0])
                    
                    # Rather than subgraph operation which requires the full graph, 
                    # you can add only additional edges from the node adjacency lists
                    merge_comp = nx_Graph(G.subgraph(merge_comp_nodes), comp_score=0)

                    (score_merge, comp_bool) = get_score(merge_comp, model, scaler, inputs['model_type'])
                    merge_comp.graph['comp_score'] = score_merge
                    sc1 = comp[1]
                    sc2 = max_over_comp[1]
                    if score_merge > sc1 and score_merge > sc2:
                        fin_list.append((frozenset(merge_comp.nodes()), merge_comp.graph['comp_score']))
                        fin_list.remove(comp)
                        fin_list.remove(max_over_comp)
                        if OS_max_ind_fin <= ind:
                            ind -= 1

                            # Otherwise: remove lower scoring complex
                    elif sc1 <= sc2:
                        fin_list.remove(comp)
                    else:
                        fin_list.remove(max_over_comp)
                        if OS_max_ind_fin > ind:
                            ind += 1
                else:
                    ind += 1
        logging_info("No. of changes = %s", str(n_changes))

    logging_info("Finished filtering complexes.")

    return fin_list


def postprocess(pred_comp_list, modelfname, scaler, inputs, G, prot_list, train_prot_list, test_prot_list):
    with open(modelfname, 'rb') as f:
        model = pickle_load(f)
    if len(pred_comp_list) == 0:
        return pred_comp_list

    # Removing complexes with only two nodes 
    # Finding unique complexes
    fin_list_graphs = set([(comp, score) for comp, score in pred_comp_list if len(comp) > 2])

    logging_info("Finished sampling complexes.")

    # Filtering complexes with high overlap with bigger complexes 
    fin_list_graphs = merge_filter_overlapped_score(fin_list_graphs, model, scaler, inputs, G)
    # Sort by scores
    fin_list_graphs = sorted(fin_list_graphs, key=lambda x: x[1], reverse=True)
    logging_info("Writing predicted complexes.")
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']

    if inputs['dir_nm'] == "humap":
        convert2names_wscores(fin_list_graphs, out_comp_nm + '_pred_names.out',
                              G, out_comp_nm + '_pred_edges_names.out')
    tot_pred_edges_unique_max_comp_prob = {}
    with open(out_comp_nm + '_pred.out', "w") as fn:
        with open(out_comp_nm + '_pred_edges.out', "wb") as f_edges:
            fn_write = fn.write
            f_edges_write = f_edges.write
            for index in range(len(fin_list_graphs)):
                tmp_graph_nodes = fin_list_graphs[index][0]
                tmp_score = fin_list_graphs[index][1]
                for node in tmp_graph_nodes:
                    fn_write("%s " % node)

                fn_write("%.3f" % tmp_score)
                tmp_graph = G.subgraph(tmp_graph_nodes)
                nx_write_weighted_edgelist(tmp_graph, f_edges)
                tmp_graph_edges = tmp_graph.edges()

                for edge in tmp_graph_edges:
                    edge_set = frozenset([edge[0], edge[1]])
                    if edge_set in tot_pred_edges_unique_max_comp_prob:
                        tot_pred_edges_unique_max_comp_prob[edge_set] = max(tot_pred_edges_unique_max_comp_prob[edge_set], tmp_score)
                    else:
                        tot_pred_edges_unique_max_comp_prob[edge_set] = tmp_score
                fn_write("\n")
                f_edges_write("\n".encode())

    with open(out_comp_nm + '_tot_pred_edges_unique_max_comp_prob.out', "w") as f:
        with open(out_comp_nm + '_tot_pred_edges_unique_max_comp_prob_inKnown.out', "w") as f_inKnown:
            with open(out_comp_nm + '_tot_pred_edges_unique_max_comp_prob_inKnown_train.out', "w") as f_inKnown_train:
                with open(out_comp_nm + '_tot_pred_edges_unique_max_comp_prob_inKnown_test.out', "w") as f_inKnown_test:
                    for edge_key in tot_pred_edges_unique_max_comp_prob:
                        edge = list(edge_key)
                        edge_score = tot_pred_edges_unique_max_comp_prob[edge_key]
                        f.write(edge[0] + "\t" + edge[1] + "\t" + "%.3f" % edge_score + "\n")

                        if edge[0] in prot_list and edge[1] in prot_list:
                            f_inKnown.write(edge[0] + "\t" + edge[1] + "\t" + "%.3f" % edge_score + "\n")
                        if edge[0] in train_prot_list and edge[1] in train_prot_list:
                            f_inKnown_train.write(edge[0] + "\t" + edge[1] + "\t" + "%.3f" % edge_score + "\n")
                        if edge[0] in test_prot_list and edge[1] in test_prot_list:
                            f_inKnown_test.write(edge[0] + "\t" + edge[1] + "\t" + "%.3f" % edge_score + "\n")

    logging_info("Finished writing predicted complexes.")

    return fin_list_graphs
