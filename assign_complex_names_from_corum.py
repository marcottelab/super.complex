# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 11:42:23 2021

@author: Meghana
"""
import pickle
import re

def jaccard_coeff(set1, set2):
    ls1 = len(set1)
    ls2 = len(set2)
    if ls1 == 0 and ls2 == 0:
        return 1
    inter = len(set1.intersection(set2))
    return float(inter) / (ls1 + ls2 - inter)

import pandas as pd
 
df = pd.read_excel('./humap/CORUM_latest/corum_2019.xlsx')
# some names are different in df['subunits(Gene name)'], so get from converted names instead
with open('./humap/corum_namesoriginal.txt') as f:
    orig_corum_complexes_gene_names_converted = f.readlines()
comp2name = list(zip(df['ComplexName'][:2916],orig_corum_complexes_gene_names_converted))

comp2name = dict([(frozenset(elt_set[1].rstrip().split(' ')),elt_set[0]) for elt_set in comp2name])

results_folder = "./humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_qi_o0.375"        
results_file = results_folder + "/res_pred_names.out"

with open(results_file) as f:
    complexes = [frozenset(line.rstrip().split()[:-1]) for line in f.readlines()]
   
complexes_with_names = dict()
name_list = []
comp_num2name = dict()
perfect_matches_num = 0 
name_duplicate_nums = dict()
for i,comp in enumerate(complexes):
    max_jc = 0
    closest_match = None
    name = None
    # get closest match
    for corum_comp in comp2name:
        jc = jaccard_coeff(comp, corum_comp)
        if jc > max_jc:
            max_jc = jc
            closest_match = corum_comp
            
    if closest_match:
        if float("{:.1f}".format(max_jc)) == 1.0:
            name = comp2name[closest_match]
            perfect_matches_num += 1
        else:
            if float("{:.1f}".format(max_jc)) == 0.0:
                name = "{:.2f}".format(0.05) + ' ' + comp2name[closest_match]
            else:    
                name = "{:.1f}".format(max_jc) + ' ' + comp2name[closest_match]
            
    if not name:
        name = ' '.join(comp)
        
    if name+'\n' in name_list: # Ensuring unique names
        if name not in name_duplicate_nums:
            name_duplicate_nums[name] = 1
        else:
            name_duplicate_nums[name] += 1
        # get last number in list
        num = name_duplicate_nums[name]+1
        if num == '':
            num=1
        name = name + ' match '+ str(num)
    complexes_with_names[comp] = name
    comp_num2name['Complex'+str(i+1)]=name
    name_list.append(name+'\n')

with open(results_folder+'/res_pred_complex_num2name.pkl','wb') as f:
    pickle.dump(comp_num2name,f)
    
with open(results_folder+'/res_pred_complex_set2name.pkl','wb') as f:
    pickle.dump(complexes_with_names,f)
    
with open(results_folder+'/res_metrics.out','a') as f:  
    print('No. of perfectly recalled original CORUM matches = ',perfect_matches_num,file = f)     # 10
    print('No. of duplicate complexes (match with corum)= ',sum(name_duplicate_nums.values()),file = f)     # 10    

with open(results_folder+'/res_pred_complex_names.txt','w') as f:
    f.writelines(name_list)        
    

