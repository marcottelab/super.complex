# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 18:13:47 2021

@author: Meghana
"""
from os import path as os_path,chdir as os_chdir

os_chdir(os_path.dirname(os_path.abspath(__file__)))
from sys import path as sys_path

# insert at 1, 0 is the script path (or '' in REPL)
sys_path.insert(1, './functions_py3/')
from jaccard_coeff import jaccard_coeff
from numpy import percentile as np_percentile
from read_complexes import preprocess_complexes
#from pickle import load as pickle_load
import networkx

def NA_threshold(set1, set2):
    ls1 = len(set1)
    ls2 = len(set2)
    if ls1 == 0 and ls2 == 0:
        return 1
    inter = float(len(set1.intersection(set2)))
    
    a = inter/ls1
    b = inter/ls2
    return max(a,b)

def get_overlap_threshold(inputs,pp_flag,G=[]):
    if pp_flag:
        comps = preprocess_complexes(inputs['dir_nm'] + inputs['comf_nm'], ' ', G)
        comps = [set(list(g.nodes)) for g in comps]
    else:
        with open(inputs['dir_nm'] + inputs['comf_nm']) as f:
            comps = [set(line.rstrip().split()) for line in f.readlines()]   
        
    n_comps = len(comps)
    jcs = []
    for i in range(n_comps):
        for j in range(i+1,n_comps):
            jc= jaccard_coeff(comps[i], comps[j])
            jcs.append(jc)
                        
    jcs = [jc for jc in jcs if jc != 0]
    
    if len(jcs) == 0:
        return 0
    
    q1 = np_percentile(jcs,25)
    q2 = np_percentile(jcs,50)
        
    return float(q1 + 0.5*(q2 - q1))

def get_overlap_threshold_qi(inputs,pp_flag,G=[]):
    if pp_flag:
        comps = preprocess_complexes(inputs['dir_nm'] + inputs['comf_nm'], ' ', G)
        comps = [set(list(g.nodes)) for g in comps]
    else:    
        with open(inputs['dir_nm'] + inputs['comf_nm']) as f:
            comps = [set(line.rstrip().split()) for line in f.readlines()]   
        
    n_comps = len(comps)
    jcs = []
    for i in range(n_comps):
        for j in range(i+1,n_comps):
            jc= NA_threshold(comps[i], comps[j])
            jcs.append(jc)
            
    jcs = [jc for jc in jcs if jc != 0]
    
    if len(jcs) == 0:
        return 0    
    
    q1 = np_percentile(jcs,25)
    min_jc = np_percentile(jcs,2)
        
    return float(min_jc + (q1 - min_jc)/2.0)

#inputs = dict()
##inputs['dir_nm'] = 'humap'
##inputs['comf_nm'] = '/res_train_complexes_new_73_more.txt'
#
#inputs['dir_nm'] = 'yeast'
#inputs['comf_nm'] = '/TAP-MS.txt'
##inputs['comf_nm'] = '/mips.txt'
#
##inputs['dir_nm'] = 'toy_network'
##inputs['comf_nm'] = '/train_complexes.txt'
#
#pp_flag = 1
#inputs['graph_files_dir']='/graph_files'
#myGraphName = inputs['dir_nm'] + inputs['graph_files_dir']+ "/res_myGraph"
#with open(myGraphName, 'rb') as f:
#    myGraph = pickle_load(f)
#
#sol1 = get_overlap_threshold(inputs,pp_flag,myGraph)
#sol2=get_overlap_threshold_qi(inputs,pp_flag,myGraph)
#print(sol1)
#print(sol2)