# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 11:37:38 2021

@author: Meghana
"""
from publish_complex_list import write_comp2prot
import pickle
from argparse import ArgumentParser as argparse_ArgumentParser

parser = argparse_ArgumentParser("Input parameters")
parser.add_argument("--suffix_corum", default='', help="original")
args = parser.parse_args()
    
##suffix_corum = ''
#suffix_corum = 'original'
suffix_corum=args.suffix_corum

with open('../humap/corum_names' + suffix_corum + '.txt') as f:
    corum_names = [set(comp.rstrip().split()) for comp in f.readlines()]
    
if suffix_corum != 'original':     
    with open('../humap/cleaned_corum_complex_num2name.pkl','rb') as f:
        comp_num2name_corum = pickle.load(f)           
else:
    with open('../humap/original_corum_complex_num2name.pkl','rb') as f:
        comp_num2name_corum = pickle.load(f)       
    
complex_name_list = ['Complex' + str(i) for i in range(1,len(corum_names)+1)]

complex2prot = dict(zip(complex_name_list,corum_names))    
            
covid_comp2prot = dict()

with open("../convert_ids/covid_interactors.txt") as f:
    covid_prots = set([line.rstrip() for line in f.readlines()])
    
for comp, comp_prots in complex2prot.items():
    if len(covid_prots.intersection(comp_prots)):
        covid_comp2prot[comp] = comp_prots
        
write_comp2prot(complex2prot,comp_num2name_corum,suffix='',write_html_file = "./Complexes/" + suffix_corum + "CORUM_Complex2proteins", folder_prefix=suffix_corum)
write_comp2prot(covid_comp2prot,comp_num2name_corum,suffix='_covid',write_html_file = "./Complexes/"+ suffix_corum + "CORUM_Complex2proteins", folder_prefix=suffix_corum)