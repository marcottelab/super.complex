# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 16:16:35 2019

@author: Meghana
"""
#import os
#os.environ['JOBLIB_START_METHOD'] = "forkserver"
#export JOBLIB_START_METHOD="forkserver"
from numpy import argmax as np_argmax, ceil as np_ceil, argsort as np_argsort,exp as np_exp, sort as np_sort
from numpy.random import uniform as rand_uniform, permutation as rand_perm
from operator import itemgetter as operator_itemgetter
from tqdm import tqdm
from random import choice as random_choice, sample as random_sample
from joblib import Parallel, delayed
from pickle import load as pickle_load, dump as pickle_dump
from create_feat_mat_1 import create_feat_mat_1
from logging import info as logging_info,debug as logging_debug
#from functools import partial
from os import mkdir as os_mkdir, path as os_path
from networkx import Graph as nx_Graph, find_cliques as nx_find_cliques
from multiprocessing import cpu_count as mul_cpu_count
from distutils.dir_util import copy_tree
from shutil import copyfile


def get_score(g1,model,scaler,mod_type):
    feats = create_feat_mat_1(g1)
        
    if(mod_type == "tpot"):
        comp_bool = model.predict(feats)
        score_curr = model.predict_proba(feats)[:,1]

    elif(mod_type == "NN"):
        feats = scaler.transform(feats)  
            
        preds = model.predict(feats)
        pred = preds[0]
        comp_bool = np_argmax(pred)
        score_curr = pred[1]
    
    return (score_curr,comp_bool)
    
def search_max_neig(seed_node,scaler,par_inputs_fn):
    with open(par_inputs_fn,'rb') as f:                
        inputs = pickle_load(f)      
    with open(inputs['modelfname'],'rb') as f:
        model = pickle_load(f)    # Seed node
    logging_debug("Seed node is",seed_node)
    folNm = inputs['folNm']
    with open(folNm+"/"+seed_node,'rb') as f:
        neig_list = pickle_load(f)
    folNm_out = inputs['folNm_out']
    
    if not neig_list:
        return         
            
    imp_neig = max(neig_list) # Largest weight neighbor - gives the most confident graphs 
    wt = neig_list[imp_neig]
    wt_edge = wt['weight']
    
    score_curr = 0
    g1=nx_Graph()
    g1.add_edge(seed_node,imp_neig,weight=wt_edge)
    
    max_nodes = inputs["max_size"]
    
    
    while True:

        logging_debug("Adding next node")
            
        imp_neigs = dict()
        g1_nodes = g1.nodes()
        for node in g1_nodes:
            # get its max neighbor and weight and store in dict 
            with open(folNm+"/"+node,'rb') as f:
                neig_list = pickle_load(f)          
            
            # Remove neighbors already in graph - one small computation to save memory
            neig_fin = set(neig_list) -  set(g1_nodes)
            neig_list = dict([neig for neig in list(neig_list.items()) if neig[0] in neig_fin])
            
            if not neig_list: # Checking if empty
                break
            imp_neig = max(neig_list)
            wt = neig_list[imp_neig]
            wt_edge = wt['weight']
            imp_neigs[imp_neig] = wt_edge
                
        if not imp_neigs:
            logging_debug("No more neighbors to add")
            break
        
        node_to_add = max(imp_neigs) # Check again that this is the max 
        #ADD ALL EDGES OF NEW NODE TO ORIG GRAPH
        
        with open(folNm+"/"+node_to_add,'rb') as f:
            its_neig_list = pickle_load(f)              
            
        orig_nodes = g1.nodes()
        for node in orig_nodes:
            if node in its_neig_list:
                wt = its_neig_list[node]
                wt_edge = wt['weight']
                g1.add_edge(node_to_add,node,weight=wt_edge)
            
        if len(g1) > max_nodes:
            logging_debug("Max size exceeded")
            break
        
        score_prev = score_curr       
        
        (score_curr,comp_bool) = get_score(g1,model,scaler,inputs['model_type'])

        if comp_bool == 0:
            logging_debug("Complex found")
        
            # Remove the node last added                
            g1.remove_node(node_to_add)
            score_curr = score_prev
            break        
    with open(folNm_out+"/"+seed_node,'wb') as f:  
        pickle_dump((list(g1.nodes()),score_curr),f)

def find_imp_neig(neig_list,explore_prob):
    if len(neig_list) == 1:
        imp_neig = list(neig_list.keys())[0]
    else:    
        cur_trial = rand_uniform(low=0.0,high=1.0)
        if  cur_trial <= explore_prob: 
            logging_debug("Exploring with low probability") # Move to top for efficiency and remove del max
            imp_neig = random_choice(list(neig_list.keys()))
        else:    
            imp_neig = max(iter(neig_list.items()), key=lambda elem: elem[1]['compScore'])[0]
    return imp_neig
    

def find_max_neig(neig_list,g1,perc,model,scaler,inputs):

    n_maxs = len(neig_list)
    if n_maxs == 0:
        return None    
    if n_maxs > 10:
        # Ascending order         
        n_maxs = int(np_ceil(perc*len(neig_list)))    
    
    neig_key_list = [k for k in neig_list]
    neig_wt_list = [float(neig_list[k]['weight']) for k in neig_list]
    sorted_ind = np_argsort(neig_wt_list)
    sorted_wts = [{'weight':val} for val in np_sort(neig_wt_list)][-n_maxs:]
    sorted_neig_keys = [neig_key_list[i] for i in sorted_ind][-n_maxs:]
    imp_neigs = dict(zip(sorted_neig_keys,sorted_wts))
    
    folNm = inputs['folNm']
    
    if len(imp_neigs) == 1:
        imp_neig = list(imp_neigs.keys())[0]
        wt = imp_neigs[imp_neig]
        wt_edge = wt['weight']          
        node_to_add = imp_neig
        #ADD ALL EDGES OF NEW NODE TO ORIG GRAPH
        with open(folNm+"/"+node_to_add,'rb') as f:
            its_neig_list = pickle_load(f)      
                    
        orig_nodes = g1.nodes()
        all_nodesWedges = set(orig_nodes).intersection(its_neig_list)
        
        for node in all_nodesWedges:
            wt = its_neig_list[node]
            wt_edge = wt['weight']
            g1.add_edge(node_to_add,node,weight=wt_edge)        
        (score_imp_neig,comp_bool) = get_score(g1,model,scaler,inputs['model_type'])
        g1.remove_node(node_to_add)
    else:
        scores = {}
        for neig in imp_neigs:
            # Add to graph 
            wt = imp_neigs[neig]
            wt_edge = wt['weight']          
            node_to_add = neig
            #ADD ALL EDGES OF NEW NODE TO ORIG GRAPH
            with open(folNm+"/"+node_to_add,'rb') as f:
                its_neig_list = pickle_load(f)                              
            orig_nodes = g1.nodes()
            all_nodesWedges = set(orig_nodes).intersection(its_neig_list)
            
            for node in all_nodesWedges:
                wt = its_neig_list[node]
                wt_edge = wt['weight']
                g1.add_edge(node_to_add,node,weight=wt_edge)
            # Check score
            (score_curr,comp_bool) = get_score(g1,model,scaler,inputs['model_type'])         
            scores[neig] = score_curr
            g1.remove_node(node_to_add)
            
            
        imp_neig = max(iter(scores.items()), key=operator_itemgetter(1))[0]  
        score_imp_neig = scores[imp_neig]

    return(imp_neig,score_imp_neig)

def search_top_neigs(seed_node,scaler,par_inputs_fn): # Picks out of a subset of its neighbors and adds the best node 
    pass
    '''
    #logging_debug("No. of nodes in g = ",len(G))
    # Assigning original graph to temporary variable
    with open(par_inputs_fn,'rb') as f:                
        inputs = pickle_load(f)      
    with open(inputs['modelfname'],'rb') as f:
        model = pickle_load(f) 
    folNm = inputs['folNm']
    with open(folNm+"/"+seed_node,'rb') as f:
        neig_list = pickle_load(f)
    folNm_out = inputs['folNm_out']
        
    if not neig_list:
        return

    imp_neig = max(neig_list) # Largest weight neighbor - gives the most confident graphs 
    wt = neig_list[imp_neig]
    wt_edge = wt['weight']
    score_curr = 0
    g1=nx_Graph()
    g1.add_edge(seed_node,imp_neig,weight=wt_edge)
    
    max_nodes = inputs["max_size"]
    
    thres_neig = inputs["thres_neig"] # Threshold on number of neighbors to consider 
    while True:

        logging_debug("Adding next node")
            
        imp_neigs = dict()
        g1_nodes = g1.nodes()
        for node in g1_nodes:
            # get its max neighbor and weight and store in dict 
            with open(folNm+"/"+node,'rb') as f:
                neig_list = pickle_load(f)        
            
            # Remove neighbors already in graph - one small computation to save memory
            neig_fin = set(neig_list) -  set(g1_nodes)
            neig_list = dict([neig for neig in neig_list.items() if neig[0] in neig_fin])
            
            # Don't check all neighbors - just a subset if number of neighbors is large
            if len(neig_list) > thres_neig: # Make 500 
                neig_list = dict(random_sample(neig_list.items(),thres_neig))
            if not neig_list: # Checking if empty
                break
            imp_neig,max_score = find_max_neig(neig_list,g1,inputs['perc'],model,scaler,inputs)

            wt = neig_list[imp_neig]
            wt_edge = wt['weight']
            
            imp_neigs[imp_neig] = {'weight': wt_edge, 'compScore' : max_score}
                
        if not imp_neigs:
            logging_debug("No more neighbors to add")
            break
        
        node_to_add = find_imp_neig(imp_neigs,inputs['explore_prob'])
        #ADD ALL EDGES OF NEW NODE TO ORIG GRAPH
        with open(folNm+"/"+node_to_add,'rb') as f:
            its_neig_list = pickle_load(f)                  
        orig_nodes = g1.nodes()
        all_nodesWedges = set(orig_nodes).intersection(its_neig_list)
                
        for node in all_nodesWedges:
            wt = its_neig_list[node]
            wt_edge = wt['weight']
            g1.add_edge(node_to_add,node,weight=wt_edge)
            
        if len(g1) > max_nodes:
            logging_debug("Max size exceeded")
            break
        score_prev = score_curr 
        (score_curr,comp_bool) = get_score(g1,model,scaler,inputs['model_type']) 

        if comp_bool == 0:
            logging_debug("Complex found")
        
            # Remove the node last added                
            g1.remove_node(node_to_add)
            score_curr = score_prev
            break
       
    with open(folNm_out+"/"+seed_node,'wb') as f:  
        pickle_dump((g1.nodes(),score_curr),f)
    '''

def met(g1,model,scaler,inputs,score_prev):
    # Assigning original graph to temporary variable
    for edge in g1.edges():
        (node1,node2) = edge
    
    max_nodes = inputs["max_size"] - len(g1)
    
    num_iter = 1
    last_iter_imp = 0        
    thres_neig = inputs["thres_neig"]
    prob_metropolis = inputs["prob_metropolis"]
    rem_nodes = []
    folNm = inputs['folNm']
    
    while num_iter < max_nodes: # Limiting number of iteration rounds 

        logging_debug("Adding next node")
            
        imp_neigs = dict()
        g1_nodes = g1.nodes()
        for node in g1_nodes:
            # get its max neighbor and weight and store in dict 
            with open(folNm+"/"+node,'rb') as f:
                neig_list = pickle_load(f)    
            
            # Remove neighbors already in graph - one small computation to save memory
            neig_fin = set(neig_list) -  set(list(g1_nodes)+rem_nodes)
            neig_list = dict([neig for neig in neig_list.items() if neig[0] in neig_fin])
            
            # Don't check all neighbors - just a subset if number of neighbors is large
            if len(neig_list) > thres_neig: # Make 500 
                neig_list = dict(random_sample(neig_list.items(),thres_neig))
            if not neig_list: # Checking if empty
                break
            imp_neig,max_score = find_max_neig(neig_list,g1,inputs['perc'],model,scaler,inputs)

            wt = neig_list[imp_neig]
            wt_edge = wt['weight']
            
            imp_neigs[imp_neig] = {'weight': wt_edge, 'compScore' : max_score}
                
        if not imp_neigs:
            logging_debug("No more neighbors to add")
            break
        
        node_to_add = find_imp_neig(imp_neigs,inputs['explore_prob'])
        #ADD ALL EDGES OF NEW NODE TO ORIG GRAPH
        with open(folNm+"/"+node_to_add,'rb') as f:
            its_neig_list = pickle_load(f)                  
        orig_nodes = g1.nodes()
        all_nodesWedges = set(orig_nodes).intersection(its_neig_list)
                
        for node in all_nodesWedges:
            wt = its_neig_list[node]
            wt_edge = wt['weight']
            g1.add_edge(node_to_add,node,weight=wt_edge)
            
        (score_curr,comp_bool) = get_score(g1,model,scaler,inputs['model_type'])    

        if comp_bool == 0:
            logging_debug("Complex found")
        
            # Remove the node last added                
            g1.remove_node(node_to_add)
            break

        cur_trial = rand_uniform(low=0.0,high=1.0)
        if score_curr < score_prev:
            if  cur_trial > prob_metropolis:  
            # Remove the node last added                
                g1.remove_node(node_to_add)
                rem_nodes.append(node_to_add)
            # since edges from this node to complex have been removed from tempG it will not be revisited 
            else:
                logging_debug("Accepting with low probability")
        elif score_curr > score_prev:
            last_iter_imp = num_iter

        if (num_iter - last_iter_imp)> 10: # Has been a long time since a score improvement
            logging_debug("Long time since score imporovement") 
            break
        
        score_prev = score_curr

        num_iter += 1
    
    #print(g1.nodes())
    #print(g1.edges())
    return (g1.nodes(),score_prev)
    
def search_metropolis_clique_start(scaler,par_inputs_fn,G_clique): # Picks out of a subset of its neighbors and adds the best node 
    #print(seed_clique) 
    with open(par_inputs_fn,'rb') as f:                
        inputs = pickle_load(f)      
    with open(inputs['modelfname'],'rb') as f:
        model = pickle_load(f)       
    g1=nx_Graph(G_clique)
    
    # Finding score 
    (score_prev,comp_bool) = get_score(g1,model,scaler,inputs['model_type']) 
                
    # Removing starting points which are not complexes    
    if comp_bool == 0:
        return ([],0)
              
    g1 = met(g1,model,scaler,inputs,score_prev)
    
    return g1
    
def search_metropolis(seed_node,scaler,par_inputs_fn): # Picks out of a subset of its neighbors and adds the best node 
     
    
    with open(par_inputs_fn,'rb') as f:                
        inputs = pickle_load(f)      
    with open(inputs['modelfname'],'rb') as f:
        model = pickle_load(f)  
    folNm = inputs['folNm']
    with open(folNm+"/"+seed_node,'rb') as f:
        neig_list = pickle_load(f)
    folNm_out = inputs['folNm_out']
    '''    
    if not neig_list:
        return
            
    imp_neig = max(neig_list) # Largest weight neighbor - gives the most confident graphs 
    wt = neig_list[imp_neig]
    wt_edge = wt['weight']
    
    g1=nx_Graph() # Starting complex of 2 nodes is assigned score 0 arbitly
    g1.add_edge(seed_node,imp_neig,weight=wt_edge)
    
    a,b = met(g1,model,scaler,inputs,0)
    
    with open(folNm_out+"/"+seed_node,'wb') as f:  
        pickle_dump((a,b),f)
    '''
        
def search_isa(seed_node,scaler,par_inputs_fn): # Picks out of a subset of its neighbors and adds the best node
    pass
    '''
    with open(par_inputs_fn,'rb') as f:                
        inputs = pickle_load(f)  
    with open(inputs['modelfname'],'rb') as f:
        model = pickle_load(f)
      
    # Assigning original graph to temporary variable
    
    folNm = inputs['folNm']
    with open(folNm+"/"+seed_node,'rb') as f:
        neig_list = pickle_load(f)
    folNm_out = inputs['folNm_out']
        
    if not neig_list:
        return
    imp_neig = max(neig_list) # Largest weight neighbor - gives the most confident graphs 
    wt = neig_list[imp_neig]
    wt_edge = wt['weight']
    
    score_prev = 0
    g1=nx_Graph() # Starting complex of 2 nodes is assigned score 0 arbitly
    g1.add_edge(seed_node,imp_neig,weight=wt_edge)
    
    max_nodes = inputs["max_size"] - len(g1)
    
    num_iter = 1
    last_iter_imp = 0        
    thres_neig = inputs["thres_neig"]
    T = inputs["T0"] # T0 value 
    alpha = inputs["alpha"]
    rem_nodes = []
    
    while num_iter < max_nodes: # Limiting number of iteration rounds 

        logging_debug("Adding next node")
            
        imp_neigs = dict()
        g1_nodes = g1.nodes()
        for node in g1_nodes:
            # get its max neighbor and weight and store in dict 
            with open(folNm+"/"+node,'rb') as f:
                neig_list = pickle_load(f)    
            
            # Remove neighbors already in graph - one small computation to save memory
            neig_fin = set(neig_list) -  set(list(g1_nodes) + rem_nodes)
            neig_list = dict([neig for neig in neig_list.items() if neig[0] in neig_fin])
            
            # Don't check all neighbors - just a subset if number of neighbors is large
            if len(neig_list) > thres_neig: # Make 500 
                neig_list = dict(random_sample(neig_list.items(),thres_neig))
            if not neig_list: # Checking if empty
                break
            imp_neig,max_score = find_max_neig(neig_list,g1,inputs['perc'],model,scaler,inputs)

            wt = neig_list[imp_neig]
            wt_edge = wt['weight']
            
            imp_neigs[imp_neig] = {'weight': wt_edge, 'compScore' : max_score}
                
        if not imp_neigs:
            logging_debug("No more neighbors to add")
            break
        
        node_to_add = find_imp_neig(imp_neigs,inputs['explore_prob'])
        #ADD ALL EDGES OF NEW NODE TO ORIG GRAPH
        with open(folNm+"/"+node_to_add,'rb') as f:
            its_neig_list = pickle_load(f)                  
        orig_nodes = g1.nodes()
        all_nodesWedges = set(orig_nodes).intersection(its_neig_list)
                
        for node in all_nodesWedges:
            wt = its_neig_list[node]
            wt_edge = wt['weight']
            g1.add_edge(node_to_add,node,weight=wt_edge)
        
        (score_curr,comp_bool) = get_score(g1,model,scaler,inputs['model_type'])         

        if comp_bool == 0:
            logging_debug("Complex found")
        
            # Remove the node last added                
            g1.remove_node(node_to_add)
            break

        cur_trial = rand_uniform(low=0.0,high=1.0)
        if score_curr < score_prev:

            prob_isa = np_exp((score_curr - score_prev)/T)
            if  cur_trial > prob_isa:  
            # Remove the node last added                
                g1.remove_node(node_to_add)
                rem_nodes.append(node_to_add)
            # since edges from this node to complex have been removed from tempG it will not be revisited 
            else:
                logging_debug("Accepting with low probability")
        elif score_curr > score_prev:
            last_iter_imp = num_iter

        if (num_iter - last_iter_imp)> 10: # Has been a long time since a score improvement
            logging_debug("Long time since score imporovement") 
            break
        
        score_prev = score_curr
        num_iter += 1
        T = T/alpha
    
    #print(g1.nodes())
    #print(g1.edges())
    with open(folNm_out+"/"+seed_node,'wb') as f:
        pickle_dump((g1.nodes(),score_prev),f)            
    '''

def sample_dummy(inputs,G,modelfname,scaler,max_size_train,max_size_test,prot_list):
    num_comp = inputs['num_comp']
    seed_mode = inputs['seed_mode']
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    logging_info("Sampling complexes...")   
    max_size_thres = inputs['max_size_thres']
    max_size = max(max_size_test,max_size_train)
    if(max_size >= max_size_thres):
        max_size = max_size_thres

    with open(out_comp_nm+'_metrics.out',"a") as fid:        
        print("Max size of complexes = ",max_size,file = fid)  # NOT actual max size since you merge later
    num_cores = mul_cpu_count()

    if 'num_cores' in inputs:
        num_cores = inputs['num_cores']
    
    pred_comp_list = []
    #num_comp = 10

    if(seed_mode == "all_nodes"):
        seed_nodes = list(G.nodes())
    elif(seed_mode == "n_nodes"):
        seed_nodes = rand_perm(list(G.nodes()))[:num_comp]
    elif(seed_mode == "all_nodes_known_comp"):
        seed_nodes = prot_list
    elif(seed_mode == "cliques"):
        clique_list = list(nx_find_cliques(G))
        to_rem = []
        # Removing 2 node and big complexes
        for comp in clique_list:
            if (len(comp) <=2 or len(comp) >= max_size):
                to_rem.append(comp)
            
        for comp in to_rem:
            clique_list.remove(comp)

            
        seed_nodes = clique_list # Remove duplicates later.
    
    num_comp = len(seed_nodes)
 
    with open(out_comp_nm+'_metrics.out',"a") as fid: 
        print("No. of cores = ", num_cores,file = fid)               
        print("No. of seeds for complex search = ",num_comp,file = fid)        

    search_method = inputs["search_method"]
    folNm_out = "/tmp/" + out_comp_nm+"_orig_comps"
    folNm = inputs['dir_nm']+"/neig_dicts"
    
        
    par_inputs = {"modelfname": "/tmp/" + modelfname,"folNm_out":folNm_out,"folNm": "/tmp/" + folNm,"model_type":inputs["model_type"],
                  "max_size":max_size,"perc":inputs["perc"],"thres_neig":inputs["thres_neig"],
                  "T0":inputs["T0"],"alpha":inputs["alpha"],"prob_metropolis":inputs["prob_metropolis"],
                  "explore_prob":inputs["explore_prob"],"feats":inputs["feats"]}
    
    par_inputs_fn =  "/tmp/" + inputs['dir_nm']+"/res_par_inputs"
        
    if not os_path.exists("/tmp/"):
        os_mkdir("/tmp/")
    if not os_path.exists("/tmp/"+ inputs['dir_nm']):
        os_mkdir("/tmp/"+ inputs['dir_nm'])
    if not os_path.exists("/tmp/"+out_comp_nm[:-3]):
        os_mkdir("/tmp/" + out_comp_nm[:-3])
        
    if not os_path.exists(folNm_out):
        os_mkdir(folNm_out)    
        
    if not os_path.exists("/tmp/"+ modelfname):
        copyfile(modelfname, "/tmp/" + modelfname)

    if not os_path.exists("/tmp/" + folNm):
        copy_tree(folNm, "/tmp/" + folNm)
        
    with open(par_inputs_fn,'wb') as f:
        pickle_dump(par_inputs,f)
    if inputs["run_mode"] == "parallel":
        #method = "threads" 
        #method = "processes" # better since we are not releasing the GIL ? prefer not required when backend is specified. (prefer = method in parallel args)
        back = 'loky' # loky and multiprocessing for processes and threading for threads
        
        if "backend" in inputs:
            back = inputs['backend']
        if(seed_mode == "cliques"):
            pred_comp_list = Parallel(n_jobs=num_cores,backend=back)(delayed(search_metropolis_clique_start)(scaler,par_inputs_fn,G.subgraph(clique)) for clique in tqdm(clique_list)) 
        elif(search_method == "isa"):
            Parallel(n_jobs=num_cores,backend=back)(delayed(search_isa)(node,scaler,par_inputs_fn) for node in tqdm(seed_nodes))         
                
        elif(search_method == "metropolis"):     
            Parallel(n_jobs=num_cores,backend=back)(delayed(search_metropolis)(node,scaler,par_inputs_fn) for node in tqdm(seed_nodes))
            
        elif(search_method == "search_top_neigs"):         
            Parallel(n_jobs=num_cores,backend=back)(delayed(search_top_neigs)(node,scaler,par_inputs_fn) for node in tqdm(seed_nodes))           
        else:       
            Parallel(n_jobs=num_cores,backend=back)(delayed(search_max_neig)(node,scaler,par_inputs_fn) for node in tqdm(seed_nodes))           
      
        '''
        # With multiprocessing
        pool = multiprocessing.Pool(processes=num_cores)
        pred_comp_list = list(tqdm(pool.imap_unordered(partial(search_isa,model=model,scaler=scaler,inputs=par_inputs),seed_nodes), total=len(seed_nodes)))
        
        # Asynchronous
        #r = pool.map_async(partial(search_isa,model=model,scaler=scaler,inputs=par_inputs),seed_nodes)
        #r.wait()
        #pred_comp_list = r.get()
        '''
    else:
        if(seed_mode == "cliques"):
            for clique in tqdm(clique_list):
                search_metropolis_clique_start(scaler,par_inputs_fn,G.subgraph(clique)) 
        elif(search_method == "isa"):
            for node in tqdm(seed_nodes):
                search_isa(node,scaler,par_inputs_fn)                   
        elif(search_method == "metropolis"):     
            for node in tqdm(seed_nodes):
                search_metropolis(node,scaler,par_inputs_fn)                 
        elif(search_method == "search_top_neigs"):         
            for node in tqdm(seed_nodes):
                search_top_neigs(node,scaler,par_inputs_fn)                 
        else:       
            for node in tqdm(seed_nodes):
                search_max_neig(node,scaler,par_inputs_fn)                 
      
    pred_comp_list = [ pickle_load(open(folNm_out+"/"+seed_node,'rb'))  for seed_node in seed_nodes if os_path.exists(folNm_out+"/"+seed_node)]
    return pred_comp_list,num_comp
