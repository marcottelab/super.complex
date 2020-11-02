# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 16:16:35 2019

@author: Meghana
"""
# import os
# os.environ['JOBLIB_START_METHOD'] = "forkserver"
# export JOBLIB_START_METHOD="forkserver"
from numpy import exp as np_exp
from numpy.random import uniform as rand_uniform
from tqdm import tqdm
from joblib import Parallel, delayed
from pickle import load as pickle_load, dump as pickle_dump
from logging import info as logging_info, debug as logging_debug
from os import mkdir as os_mkdir, path as os_path,remove as os_remove
from distutils.dir_util import copy_tree
from shutil import copyfile,rmtree
from networkx import Graph as nx_Graph
from multiprocessing import cpu_count as mul_cpu_count
from sample_helpers import *


def search_max_neig(seed_node, scaler, par_inputs_fn):
    with open(par_inputs_fn, 'rb') as f:
        inputs = pickle_load(f)
    with open(inputs['modelfname'], 'rb') as f:
        model = pickle_load(f)  # Seed node
    logging_debug("Seed node is", seed_node)
    folNm = inputs['folNm']
    folNm_out = inputs['folNm_out']
    max_nodes = inputs["max_size"]
    score_curr = 0
    cd, g1 = starting_edge(folNm, seed_node)
    if cd == 0:
        return
    while len(g1) < max_nodes:
        # print(len(g1))
        logging_debug("Adding next node")

        neig_list = read_neigs(g1.nodes(), folNm)
        if not neig_list:  # Checking if empty
            logging_debug("No more neighbors to add")
            break
        node_to_add = max(neig_list.items(), key=lambda elem: elem[1]['weight'])[0]
        g1 = add_newnode(g1, node_to_add, neig_list[node_to_add]['graph_neigs'])

        score_prev = score_curr

        (score_curr, comp_bool) = get_score(g1, model, scaler, inputs['model_type'])

        if comp_bool == 0:
            logging_debug("Complex found")

            # Remove the node last added                
            g1.remove_node(node_to_add)
            score_curr = score_prev
            break
    with open(folNm_out + "/" + seed_node, 'wb') as f:
        pickle_dump((frozenset(g1.nodes()), score_curr), f)


def search_top_neigs(seed_node, scaler, par_inputs_fn):
    # Picks out of a subset of its neighbors and adds the best node
    # logging_debug("No. of nodes in g = ",len(G))
    # Assigning original graph to temporary variable
    with open(par_inputs_fn, 'rb') as f:
        inputs = pickle_load(f)
    with open(inputs['modelfname'], 'rb') as f:
        model = pickle_load(f)
    folNm = inputs['folNm']
    folNm_out = inputs['folNm_out']
    cd, g1 = starting_edge(folNm, seed_node)
    if cd == 0:
        return
    score_curr = 0
    max_nodes = inputs["max_size"]
    thres_neig = inputs["thres_neig"]  # Threshold on number of neighbors to consider
    while len(g1) < max_nodes:

        score_prev = score_curr
        logging_debug("Adding next node")
        g1, cc, node_to_add, score_curr, comp_bool, rand_flag = add_top_neig(g1, thres_neig, folNm, inputs,
                                                                                        model, scaler)
        if (score_curr is None) or (comp_bool is None):
            score_curr, comp_bool = get_score(g1, model, scaler, inputs['model_type'])
        if cc == 0:
            break
        if comp_bool == 0:
            logging_debug("Complex found")

            # Remove the node last added                
            g1.remove_node(node_to_add)
            score_curr = score_prev
            break

    with open(folNm_out + "/" + seed_node, 'wb') as f:
        pickle_dump((frozenset(g1.nodes()), score_curr), f)


def met(g1, model, scaler, inputs, score_prev):
    max_nodes = inputs["max_size"] - len(g1)

    num_iter = 1
    last_iter_imp = 0
    thres_neig = inputs["thres_neig"]
    prob_metropolis = inputs["prob_metropolis"]
    folNm = inputs['folNm']
    met_low_prob_acc = 0
    while num_iter < max_nodes:  # Limiting number of iteration rounds

        logging_debug("Adding next node")
        # neig_list_old = neig_list
        # g1, cc, node_to_add, score_curr, comp_bool, rand_flag, neig_list = add_top_neig(g1, thres_neig, folNm, inputs, model, scaler, neig_list)
        g1, cc, node_to_add, score_curr, comp_bool, rand_flag = add_top_neig(g1, thres_neig, folNm, inputs, model, scaler)
        if (score_curr is None) or (comp_bool is None):
            score_curr, comp_bool = get_score(g1, model, scaler, inputs['model_type'])

        if cc == 0:
            break
        if comp_bool == 0:
            logging_debug("Complex found")

            # Remove the node last added                
            g1.remove_node(node_to_add)
            break

        cur_trial = rand_uniform(low=0.0, high=1.0)
        if score_curr < score_prev:
            if cur_trial > prob_metropolis:
                # Remove the node last added
                g1.remove_node(node_to_add)
                # neig_list = neig_list_old
            else:
                logging_debug("Accepting with low probability")
                met_low_prob_acc += 1
                rand_flag = 1
        elif score_curr > score_prev:
            last_iter_imp = num_iter

        if (num_iter - last_iter_imp) > 10:  # Has been a long time since a score improvement
            logging_debug("Long time since score improvement")
            break

        score_prev = score_curr

        num_iter += 1
    logging_debug("No. of low probability acceptances = ")
    logging_debug(str(met_low_prob_acc))

    # print(g1.nodes())
    # print(g1.edges())
    return frozenset(g1.nodes()), score_prev


def search_metropolis_clique_start(scaler, par_inputs_fn, G_clique):
    # Picks out of a subset of its neighbors and adds the best node
    # print(seed_clique)
    with open(par_inputs_fn, 'rb') as f:
        inputs = pickle_load(f)
    with open(inputs['modelfname'], 'rb') as f:
        model = pickle_load(f)
    g1 = nx_Graph(G_clique)

    # Finding score
    score_prev, comp_bool = get_score(g1, model, scaler, inputs['model_type'])

    # Removing starting points which are not complexes    
    if comp_bool == 0:
        return
    a, b = met(g1, model, scaler, inputs, score_prev)
    name = " ".join([str(n) for n in g1.nodes()])
    with open(folNm_out + "/" + name, 'wb') as f:
        pickle_dump((a, b), f)


def search_metropolis(seed_node, scaler, par_inputs_fn):
    # Picks out of a subset of its neighbors and adds the best node
    with open(par_inputs_fn, 'rb') as f:
        inputs = pickle_load(f)
    with open(inputs['modelfname'], 'rb') as f:
        model = pickle_load(f)
    folNm = inputs['folNm']
    folNm_out = inputs['folNm_out']

    cd, g1 = starting_edge(folNm, seed_node)
    if cd == 0:
        return

    a, b = met(g1, model, scaler, inputs, 0)

    with open(folNm_out + "/" + seed_node, 'wb') as f:
        pickle_dump((a, b), f)


def search_isa(seed_node, scaler, par_inputs_fn):  # Picks out of a subset of its neighbors and adds the best node
    with open(par_inputs_fn, 'rb') as f:
        inputs = pickle_load(f)
    with open(inputs['modelfname'], 'rb') as f:
        model = pickle_load(f)

    folNm = inputs['folNm']
    folNm_out = inputs['folNm_out']

    score_prev = 0
    cd, g1 = starting_edge(folNm, seed_node)
    if cd == 0:
        return
    max_nodes = inputs["max_size"] - len(g1)

    num_iter = 1
    last_iter_imp = 0
    thres_neig = inputs["thres_neig"]
    T = inputs["T0"]  # T0 value
    alpha = inputs["alpha"]

    while num_iter < max_nodes:  # Limiting number of iteration rounds

        logging_debug("Adding next node")
        # neig_list_old = neig_list
        # g1, cc, node_to_add, score_curr, comp_bool, rand_flag, neig_list = add_top_neig(g1, thres_neig, folNm, inputs, model, scaler, neig_list)
        g1, cc, node_to_add, score_curr, comp_bool, rand_flag = add_top_neig(g1, thres_neig, folNm, inputs, model, scaler)
        if (score_curr is None) or (comp_bool is None):
            score_curr, comp_bool = get_score(g1, model, scaler, inputs['model_type'])

        if cc == 0:
            break

        if comp_bool == 0:
            logging_debug("Complex found")

            # Remove the node last added                
            g1.remove_node(node_to_add)
            break

        cur_trial = rand_uniform(low=0.0, high=1.0)
        if score_curr < score_prev:

            prob_isa = np_exp((score_curr - score_prev) / T)
            if cur_trial > prob_isa:
                # Remove the node last added
                g1.remove_node(node_to_add)
                # neig_list = neig_list_old
            else:
                logging_debug("Accepting with low probability")
                rand_flag = 1
        elif score_curr > score_prev:
            last_iter_imp = num_iter

        if (num_iter - last_iter_imp) > 10:  # Has been a long time since a score improvement
            logging_debug("Long time since score imporovement")
            break

        score_prev = score_curr
        num_iter += 1
        T = float(T) / alpha

    with open(folNm_out + "/" + seed_node, 'wb') as f:
        pickle_dump((frozenset(g1.nodes()), score_prev), f)


def sample(inputs, G, modelfname, scaler, seed_nodes, max_size,transfer2tmp):
    seed_mode = inputs['seed_mode']

    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    logging_info("Sampling complexes...")

    num_cores = mul_cpu_count()

    if 'num_cores' in inputs:
        num_cores = inputs['num_cores']

    # num_comp = 10

    num_comp = len(seed_nodes)

    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("No. of cores = ", num_cores, file=fid)
        print("No. of seeds for complex search = ", num_comp, file=fid)

    search_method = inputs["search_method"]
    folNm_out = "/tmp/" + out_comp_nm + "_orig_comps"
    folNm = inputs['dir_nm'] + "/neig_dicts"


    if not os_path.exists("/tmp/"):
        os_mkdir("/tmp/")
    if not os_path.exists("/tmp/" + inputs['dir_nm']):
        os_mkdir("/tmp/" + inputs['dir_nm'])
    if not os_path.exists("/tmp/" + out_comp_nm[:-3]):
        os_mkdir("/tmp/" + out_comp_nm[:-3])
        
    out_comp_nm_model = inputs['dir_nm'] + inputs['model_dir']
    
    if not os_path.exists("/tmp/" + out_comp_nm_model[:-3]):
        os_mkdir("/tmp/" + out_comp_nm_model[:-3])

    if not os_path.exists(folNm_out):
        os_mkdir(folNm_out)
  
    tmp_prefix_read = ""

    if transfer2tmp:
        tmp_prefix_read = "/tmp/"
        
    par_inputs = {"use_all_neigs": inputs["use_all_neigs"], "min_thres_neig_sorted": inputs["min_thres_neig_sorted"],
                  "modelfname": tmp_prefix_read + modelfname, "folNm_out": folNm_out, "folNm": tmp_prefix_read + folNm,
                  "model_type": inputs["model_type"],
                  "max_size": max_size, "perc": inputs["perc"], "thres_neig": inputs["thres_neig"],
                  "T0": inputs["T0"], "alpha": inputs["alpha"], "prob_metropolis": inputs["prob_metropolis"],
                  "explore_prob": inputs["explore_prob"], "feats": inputs["feats"]}
        
    par_inputs_fn = tmp_prefix_read + inputs['dir_nm'] + "/res_par_inputs"

    if transfer2tmp:
        if not os_path.exists("/tmp/" + modelfname):
            copyfile(modelfname, "/tmp/" + modelfname)
    
        if not os_path.exists("/tmp/" + folNm):
            os_mkdir("/tmp/" + folNm)
            copy_tree(folNm, "/tmp/" + folNm)

    with open(par_inputs_fn, 'wb') as f:
        pickle_dump(par_inputs, f)

    if inputs["run_mode"] == "parallel":
        # method = "threads" method = "processes" # better since we are not releasing the GIL ? prefer not required
        # when backend is specified. (prefer = method in parallel args)
        back = 'loky'  # loky and multiprocessing for processes and threading for threads

        if "backend" in inputs:
            back = inputs['backend']
        if seed_mode == "cliques":
            Parallel(n_jobs=num_cores, backend=back)(
                delayed(search_metropolis_clique_start)(scaler, par_inputs_fn, G.subgraph(clique)) for clique in
                tqdm(seed_nodes))
        elif search_method == "isa":
            Parallel(n_jobs=num_cores, backend=back)(
                delayed(search_isa)(node, scaler, par_inputs_fn) for node in tqdm(seed_nodes))

        elif search_method == "metropolis":
            Parallel(n_jobs=num_cores, backend=back)(
                delayed(search_metropolis)(node, scaler, par_inputs_fn) for node in tqdm(seed_nodes))

        elif search_method == "search_top_neigs":
            Parallel(n_jobs=num_cores, backend=back)(
                delayed(search_top_neigs)(node, scaler, par_inputs_fn) for node in tqdm(seed_nodes))
        else:
            Parallel(n_jobs=num_cores, backend=back)(
                delayed(search_max_neig)(node, scaler, par_inputs_fn) for node in tqdm(seed_nodes))

        '''# With multiprocessing pool = multiprocessing.Pool(processes=num_cores) pred_comp_list = list(tqdm(
        pool.imap_unordered(partial(search_isa,model=model,scaler=scaler,inputs=par_inputs),seed_nodes), 
        total=len(seed_nodes))) 
        
        # Asynchronous
        #r = pool.map_async(partial(search_isa,model=model,scaler=scaler,inputs=par_inputs),seed_nodes)
        #r.wait()
        #pred_comp_list = r.get()
        '''
    else:
        if seed_mode == "cliques":
            for clique in tqdm(seed_nodes):
                search_metropolis_clique_start(scaler, par_inputs_fn, G.subgraph(clique))
        elif search_method == "isa":
            for node in tqdm(seed_nodes):
                search_isa(node, scaler, par_inputs_fn)
        elif search_method == "metropolis":
            for node in tqdm(seed_nodes):
                search_metropolis(node, scaler, par_inputs_fn)
        elif search_method == "search_top_neigs":
            for node in tqdm(seed_nodes):
                search_top_neigs(node, scaler, par_inputs_fn)
        else:
            for node in tqdm(seed_nodes):
                search_max_neig(node, scaler, par_inputs_fn)

    if transfer2tmp:
        # Delete copied contents in tmp folder
        rmtree("/tmp/" + folNm)
        os_remove("/tmp/" + modelfname)
    
    return num_comp
