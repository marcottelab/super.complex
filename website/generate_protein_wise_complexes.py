# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 10:14:01 2021

@author: Meghana
"""

from os import path as os_path, mkdir as os_mkdir, chdir as os_chdir
import pickle
import argparse


def get_uniprot_urls(prot2complex,prot2url_file="../convert_ids/name2uniprotURL_full_humap_from_gene_names.pkl"):
    with open(prot2url_file,'rb') as f:
        prot2url = pickle.load(f)    
        links = dict()
    
    for key in prot2complex:
        if key in prot2url:
            links[key] = '<a href="' + prot2url[key] + '">' + key + '</a>'
        else:
            links[key] = key
            
    return links


def write_prot2comp(prot2complex,comp_num2name,write_path,suffix='',write_html_file = "./Complexes/Protein2complex", folder_pref='Predicted/',prot2url_file="../convert_ids/name2uniprotURL_full_humap_from_gene_names.pkl"):
    if not os_path.exists(write_path):
        with open(write_path ,'w') as f:
            f.writelines([key + ' : ' + ' '.join(val)+'\n' for key,val in prot2complex.items()])
    
    prefix = ''
    suffix2 = 'predictions'

    link = 'All proteins with corresponding predicted complexes are listed below. Only SARS-COV2 interacting proteins and their complexes can be found <a href="Protein2complex_covid.html"> here</a>. <br>'
    if 'CORUM' in write_html_file:
        link = 'All proteins with corresponding CORUM complexes are listed below. Only SARS-COV2 interacting proteins and their complexes can be found <a href="CORUM_Protein2complex_covid.html"> here</a>. <br>'
        if folder_pref == 'original':
            folder_pref = 'originalCORUM/'
        else:
            folder_pref = 'CORUM/'
        suffix2 = 'CORUM'
        
    if suffix == '_covid':
        prefix = 'SARS-COV2 interacting '
        link = ''
    pre = '''
    
    <!doctype html>
    
    <html>
    
    <head>
    <title> ''' + prefix + '''Protein-wise complex list - ''' + suffix2 + '''</title>
    
        </head>
        
    
        
        <body>''' + link +  '''
        <section>
    	
     <table style="width:100%">
      <tr>
        <th>''' + prefix + '''Protein</th>
        <th>Complex list</th>
      </tr>
    '''    
    
    post = '''     </section>
    
    
        </body>
        </html>
        '''
        
    links = get_uniprot_urls(prot2complex,prot2url_file)
    
    # sort alphabetical
    prot2complex_items = sorted(prot2complex.items())
    
    with open(write_html_file  + suffix + ".html",'w') as f:
        f.write(pre)
        f.writelines([ '<tr>' + '<td>' + links[key]  + ' </td> <td>' + ' | '.join(['<a href="' +folder_pref+ node + '.html">' + comp_num2name[node] + '</a>' for node in val])+'</td> </tr>' for key,val in prot2complex_items])
        f.write(post)    
    

def get_prot2complex(complexes):
    prot2complex = dict()
    
    for i,comp in enumerate(complexes):
        for node in comp:
            if node not in prot2complex:
                prot2complex[node] = ['Complex' + str(i+1)]
            else:
                prot2complex[node].append('Complex' + str(i+1))
    
    return prot2complex
    

def main():
    parser = argparse.ArgumentParser("Input parameters")
    parser.add_argument("--results_folder", default="../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_qi_o0.375", help="Input parameters file name")
    parser.add_argument("--id2name_file_pkl", default="../convert_ids/humap_gene_id_name_map_updated.pkl", help="Input parameters file name")
    parser.add_argument("--prot2url_file", default="../convert_ids/name2uniprotURL_full_humap_from_gene_names.pkl", help="Input parameters file name")
    
    
    args = parser.parse_args()
    # results_folder = "../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_qi_o0.375"        
    # #results_folder = "../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_best"
    results_folder = args.results_folder
    results_file = results_folder + "/res_pred_names.out"
    
    with open(results_file) as f:
        complexes = [line.rstrip().split()[:-1] for line in f.readlines()]
        # recall order is in decreasing score size
        
    with open(args.id2name_file_pkl,'rb') as f:
        id2name=pickle.load(f)
    # Update gene names if any
    complexes = [[prot if (prot not in id2name) or (id2name[prot] == '-') else id2name[prot] for prot in comp] for comp in complexes ]
    
    
    prot2complex = get_prot2complex(complexes)
                
    covid_prot2complex = dict()
    
    with open("../convert_ids/covid_interactors.txt") as f:
        covid_prots = set([line.rstrip() for line in f.readlines()])
        
    for key in prot2complex:
        if key in covid_prots:
            covid_prot2complex[key] = prot2complex[key]
            
    with open(results_folder+'/res_pred_complex_num2name.pkl','rb') as f:
        comp_num2name = pickle.load(f)
            
    write_prot2comp(prot2complex,comp_num2name,write_path = results_folder + '/res_pred_prot2comp' +  ".out",suffix='',write_html_file = "./Complexes/Protein2complex", folder_pref='Predicted/', prot2url_file=args.prot2url_file)
    write_prot2comp(covid_prot2complex,comp_num2name, write_path = results_folder + '/res_pred_prot2comp' + '_covid' + ".out",suffix='_covid',write_html_file = "./Complexes/Protein2complex", folder_pref='Predicted/', prot2url_file=args.prot2url_file)
    
if __name__ == "__main__":
    main()