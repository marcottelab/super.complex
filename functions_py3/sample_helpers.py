from numpy import argmax as np_argmax, ceil as np_ceil
from numpy.random import uniform as rand_uniform
from random import choice as random_choice, sample as random_sample
from pickle import load as pickle_load
from create_feat_mat_1 import create_feat_mat_1
from logging import info as logging_info, debug as logging_debug
from networkx import Graph as nx_Graph


def update_neig_list(neig_list, node, folNm, g1_nodes):
    # remove newly added node from neighbor list
    del neig_list[node]
    # get its max neighbor and weight and store in dict
    g1_nodes = set(g1_nodes)
    with open(folNm + "/" + node, 'rb') as f:
        temp = pickle_load(f)
    # Remove neighbors already in graph - one small computation to save memory
    temp_fin = set(temp) - g1_nodes
    temp = dict([neig for neig in temp.items() if neig[0] in temp_fin])
    for neig in temp:
        if neig in neig_list:
            temp_wt_edg = temp[neig]['weight']
            if neig_list[neig]['weight'] < temp_wt_edg:  # Max weight of neigs
                neig_list[neig]['weight'] = temp_wt_edg
            neig_list[neig]['graph_neigs'].append((node, temp_wt_edg))
        else:
            neig_list[neig] = temp[neig]
            neig_list[neig]['graph_neigs'] = [(node, temp[neig]['weight'])]

    return neig_list


def read_neigs(g1_nodes, folNm):
    neig_list = {}
    for node in g1_nodes:
        # get its max neighbor and weight and store in dict
        with open(folNm + "/" + node, 'rb') as f:
            temp = pickle_load(f)
        for neig in temp:
            if neig in neig_list:
                temp_wt_edg = temp[neig]['weight']
                if neig_list[neig]['weight'] < temp_wt_edg:  # Max weight of neigs
                    neig_list[neig]['weight'] = temp_wt_edg
                neig_list[neig]['graph_neigs'].append((node, temp_wt_edg))
            else:
                neig_list[neig] = temp[neig]
                neig_list[neig]['graph_neigs'] = [(node, temp[neig]['weight'])]

    # Remove neighbors already in graph - one small computation to save memory
    neig_fin = set(neig_list) - set(g1_nodes)
    neig_list = dict([neig for neig in neig_list.items() if neig[0] in neig_fin])
    return neig_list


def add_top_neig(g1, thres_neig, folNm, inputs, model, scaler):
    g1_nodes = g1.nodes()
    neig_list = read_neigs(g1_nodes, folNm)
    rand_flag = 0
    if inputs["use_all_neigs"] == 0:
        # Don't check all neighbors - just a subset if number of neighbors is large
        if len(neig_list) > thres_neig:  # Make 500
            neig_list = dict(random_sample(neig_list.items(), thres_neig))
            rand_flag = 1
    if not neig_list:  # Checking if empty
        logging_debug("No more neighbors to add")
        return g1, 0, None, None, None, rand_flag

    node_to_add, score, compBool, rand_flag = find_imp_neig(neig_list, g1, inputs['perc'], model, scaler, inputs,
                                                            inputs['explore_prob'], rand_flag)

    g1 = add_newnode(g1, node_to_add, neig_list[node_to_add]['graph_neigs'])
    return g1, 1, node_to_add, score, compBool, rand_flag


def add_top_neig_update(g1, thres_neig, folNm, inputs, model, scaler, neig_list):
    neig_list_orig = None
    rand_flag = 0
    if not neig_list:  # Checking if empty
        logging_debug("No more neighbors to add")
        return g1, 0, None, None, None, rand_flag, neig_list
    if inputs["use_all_neigs"] == 0:
        # Don't check all neighbors - just a subset if number of neighbors is large
        if len(neig_list) > thres_neig:  # Make 500
            neig_list_orig = neig_list
            neig_list = dict(random_sample(neig_list.items(), thres_neig))
            rand_flag = 1

    node_to_add, score, compBool, rand_flag = find_imp_neig(neig_list, g1, inputs['perc'], model, scaler, inputs,
                                                            inputs['explore_prob'], rand_flag)

    g1 = add_newnode(g1, node_to_add, neig_list[node_to_add]['graph_neigs'])

    if neig_list_orig is None:
        neig_list = update_neig_list(neig_list, node_to_add, folNm, g1.nodes())
    else:
        neig_list = update_neig_list(neig_list_orig, node_to_add, folNm, g1.nodes())
    return g1, 1, node_to_add, score, compBool, rand_flag, neig_list


def add_newnode(g1, node_to_add, neig_wts):
    for node, wt_edge in neig_wts:
        g1.add_edge(node_to_add, node, weight=wt_edge)

    return g1


def get_score(g1, model, scaler, mod_type):
    feats = create_feat_mat_1(g1)

    if mod_type == "NN":
        feats = scaler.transform(feats)

        preds = model.predict(feats)
        pred = preds[0]
        comp_bool = np_argmax(pred)
        score_curr = pred[1]

    else:  # mod_type == "tpot"
        comp_bool = model.predict(feats)
        score_curr = model.predict_proba(feats)[:, 1]

    return float(score_curr), comp_bool


def starting_edge_update(folNm, seed_node):
    with open(folNm + "/" + seed_node, 'rb') as f:
        neig_list = pickle_load(f)

    for neighh in neig_list:
        neig_list[neighh]['graph_neigs'] = [(seed_node, neig_list[neighh]['weight'])]
    cd = 1
    if not neig_list:
        cd = 0
        return cd, None, None
    imp_neig = max(neig_list.items(), key=lambda elem: elem[1]['weight'])[0]
    # Largest weight neighbor - gives the most confident graphs
    wt_edge = neig_list[imp_neig]['weight']

    g1 = nx_Graph()
    g1.add_edge(seed_node, imp_neig, weight=wt_edge)
    neig_list = update_neig_list(neig_list, imp_neig, folNm, g1.nodes())
    return cd, g1, neig_list


def starting_edge(folNm, seed_node):
    with open(folNm + "/" + seed_node, 'rb') as f:
        neig_list = pickle_load(f)

    cd = 1
    if not neig_list:
        cd = 0
        return cd, None
    imp_neig = max(neig_list.items(), key=lambda elem: elem[1]['weight'])[0]
    # Largest weight neighbor - gives the most confident graphs
    wt_edge = neig_list[imp_neig]['weight']

    g1 = nx_Graph()
    g1.add_edge(seed_node, imp_neig, weight=wt_edge)
    return cd, g1


def pick_top_weight_neigs(neig_list, perc, inputs):
    n_maxs = len(neig_list)
    if n_maxs == 0:
        return None
    if n_maxs > inputs['min_thres_neig_sorted']:
        n_maxs = int(np_ceil(perc * len(neig_list)))

    imp_neigs = dict(sorted(neig_list.items(), key=lambda elem: elem[1]['weight'], reverse=True)[:n_maxs])
    return imp_neigs


def find_imp_neig(imp_neigs, g1, perc, model, scaler, inputs, explore_prob, rand_flag):
    score_fin = None
    comp_bool_fin = None
    if len(imp_neigs) == 1:
        imp_neig = list(imp_neigs.keys())[0]
    else:
        cur_trial = rand_uniform(low=0.0, high=1.0)
        if cur_trial <= explore_prob:
            logging_debug("Exploring with low probability")  # Move to top for efficiency and remove del max
            imp_neig = random_choice(list(imp_neigs.keys()))
            rand_flag = 1
        else:
            if inputs["use_all_neigs"] == 0:
                imp_neigs = pick_top_weight_neigs(imp_neigs, perc, inputs)

            for neig in imp_neigs:
                # Add to graph
                g1 = add_newnode(g1, neig, imp_neigs[neig]['graph_neigs'])
                # Check score
                (score_curr, comp_bool) = get_score(g1, model, scaler, inputs['model_type'])
                imp_neigs[neig]['compScore'] = score_curr
                imp_neigs[neig]['compBool'] = comp_bool
                g1.remove_node(neig)
            imp_neig_tup = max(imp_neigs.items(), key=lambda elem: elem[1]['compScore'])
            imp_neig = imp_neig_tup[0]
            score_fin = imp_neig_tup[1]['compScore']
            comp_bool_fin = imp_neig_tup[1]['compBool']
    return imp_neig, score_fin, comp_bool_fin, rand_flag
