# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 11:37:38 2021

@author: Meghana
"""
from publish_complex_list import write_comp2prot
import pickle
from argparse import ArgumentParser as argparse_ArgumentParser
import pandas as pd

parser = argparse_ArgumentParser("Input parameters")
parser.add_argument("--sars_cov2_map", default="../convert_ids/SARS_COV2_Map.xlsx", help="Input parameters file name")
parser.add_argument("--suffix_corum", default='', help="original")
parser.add_argument("--dir", default="../humap/", help="Input parameters file name")
parser.add_argument("--prot2url_file", default="../convert_ids/name2uniprotURL_full_humap_from_gene_names.pkl", help="Input parameters file name")


args = parser.parse_args()
    
##suffix_corum = ''
#suffix_corum = 'original'
suffix_corum=args.suffix_corum

df = pd.read_excel(args.sars_cov2_map,header=1)

edges = list(zip(df['PreyGene'],df['Bait']))

prot2sars = dict()
for edge in edges:
    prot = edge[0]
    sars = edge[1]
    if prot not in prot2sars:
        prot2sars[prot] = set([sars])
    else:
        prot2sars[prot].add(sars)

with open(args.dir +'corum_names' + suffix_corum + '.txt') as f:
    corum_names = [set(comp.rstrip().split()) for comp in f.readlines()]
    
if suffix_corum != 'original':     
    with open(args.dir + 'cleaned_corum_complex_num2name.pkl','rb') as f:
        comp_num2name_corum = pickle.load(f)           
else:
    with open(args.dir + 'original_corum_complex_num2name.pkl','rb') as f:
        comp_num2name_corum = pickle.load(f)       
    
complex_name_list = ['Complex' + str(i) for i in range(1,len(corum_names)+1)]

complex2prot = dict(zip(complex_name_list,corum_names))    
            
covid_comp2prot = dict()

with open("../convert_ids/covid_interactors.txt") as f:
    covid_prots = set([line.rstrip() for line in f.readlines()])

n_cov_prots = dict()   
for comp, comp_prots in complex2prot.items():
    num_cov_prots = len(covid_prots.intersection(comp_prots))
    if num_cov_prots:
        covid_comp2prot[comp] = comp_prots
        n_cov_prots[comp] = num_cov_prots        
        
write_comp2prot(complex2prot,comp_num2name_corum,suffix='',write_html_file = "./Complexes/" + suffix_corum + "CORUM_Complex2proteins", folder_prefix=suffix_corum,n_cov_prots=dict(),prot2sars=dict(),logfile='logs.txt',prot2url_file=args.prot2url_file)
write_comp2prot(covid_comp2prot,comp_num2name_corum,suffix='_covid',write_html_file = "./Complexes/"+ suffix_corum + "CORUM_Complex2proteins", folder_prefix=suffix_corum,n_cov_prots=n_cov_prots,prot2sars=prot2sars,logfile='logs.txt',prot2url_file=args.prot2url_file)