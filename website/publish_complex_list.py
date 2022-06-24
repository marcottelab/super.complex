# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 10:14:01 2021

@author: Meghana
"""
import pickle
import pandas as pd
import argparse


def get_uniprot_urls(complex2prot,prot2url_file="../convert_ids/name2uniprotURL_full_humap_from_gene_names.pkl"):
    with open(prot2url_file,'rb') as f:
        prot2url = pickle.load(f)    
        links = dict()
    
    for key in set.union(*list(complex2prot.values())):
        
        if key in prot2url:
            links[key] = '<a href="' + prot2url[key] + '">' + key + '</a>'
        else:
            links[key] = key
            
    return links


def write_comp2prot(complex2prot,comp2names,suffix='',write_html_file = "./Complexes/Complex2proteins", folder_prefix='',n_cov_prots=dict(),prot2sars=dict(),logfile='logs.txt',prot2url_file="../convert_ids/name2uniprotURL_full_humap_from_gene_names.pkl"):
    
    prefix = ''
    folder_pref = 'Predicted/'
    suffix2 = 'predictions'

    link = 'Protein-wise complex membership list can be viewed <a href="Protein2complex.html"> here</a>. <br> All predicted complexes are listed below with their corresponding proteins. Only SARS-COV2 interacting protein complexes from the predictions can be found <a href="Complex2proteins_covid.html"> here</a>. <br>'
    if 'CORUM' in write_html_file:
        link = 'Protein-wise complex membership list can be viewed <a href="CORUM_Protein2complex.html"> here</a>. <br> All CORUM complexes are listed below with their corresponding proteins. Only SARS-COV2 interacting protein complexes from CORUM can be found <a href="CORUM_Complex2proteins_covid.html"> here</a>. <br>'
        folder_pref = folder_prefix+'CORUM/'
        suffix2 = 'CORUM'
        
    pre = '''
    
    <!doctype html>
    
    <html>
    
    <head>
    <title> ''' + prefix + '''Protein complex list - ''' + suffix2 + '''</title>
    
        </head>
        
    
        
        <body>''' + link +  '''
        <section>
    	
     <table style="width:100%">
      <tr>
      <th>''' + prefix + '''Complex </th>
        <th>Protein membership in complex</th>
      '''
      
    if suffix == '_covid':
        prefix = 'SARS-COV2 interacting '
        link = ''  
        pre = pre + '''  <th>No. of proteins interacting with SARS-COV2</th>  
        <th>No. of interacting SARS-COV2 proteins</th>
        '''
        complex2prot_items = sorted(complex2prot.items(),key = lambda x: -n_cov_prots[x[0]])
        comp2nsars = dict([(comp,len(set([lt_item for lt in [list(prot2sars[prot]) for prot in prots if prot in prot2sars] for lt_item in lt ])))for comp,prots in complex2prot_items])
   
        counts_nsars = dict()
        for comp,val in comp2nsars.items():
            if val not in counts_nsars:
                counts_nsars[val] = 1
            else:
                counts_nsars[val] += 1
                
        with open(logfile,'a') as f:  
            print('Counts of the no. of SARS COV2 proteins an individual complex interacts with in ' + folder_pref+write_html_file  + suffix,file=f)                
            print(counts_nsars,file=f)    
            
        counts_nsars = dict()
        for comp,val in n_cov_prots.items():
            if val not in counts_nsars:
                counts_nsars[val] = 1
            else:
                counts_nsars[val] += 1
                
        with open(logfile,'a') as f:  
            print('Counts of the no. of proteins in an individual complex interacting with SARS COV2 in ' + folder_pref+write_html_file  + suffix,file=f)                
            print(sorted(counts_nsars.items(),key=lambda x:x[0]),file=f)                
                
    pre = pre+'''
      </tr>
    '''    
    
    post = '''     </section>
    
    
        </body>
        </html>
        '''
        
    links = get_uniprot_urls(complex2prot,prot2url_file)
    
    with open(write_html_file  + suffix + ".html",'w') as f:
        f.write(pre)
        if suffix == '_covid':
            f.writelines([ '<tr>' + '<td>' + '<a href="' +folder_pref+ key + '.html">' + comp2names[key] + '</a>'  + ' </td> <td>' + ' '.join([links[prot] for prot in val])+'</td> <td>' + str(n_cov_prots[key])+'</td><td>' + str(comp2nsars[key])+'</td></tr>' for key,val in complex2prot_items])
        else:
            f.writelines([ '<tr>' + '<td>' + '<a href="' +folder_pref+ key + '.html">' + comp2names[key] + '</a>'  + ' </td> <td>' + ' '.join([links[prot] for prot in val])+'</td> </tr>' for key,val in complex2prot.items()])
        f.write(post)    

def main():    
    parser = argparse.ArgumentParser("Input parameters")
    parser.add_argument("--sars_cov2_map", default="../convert_ids/SARS_COV2_Map.xlsx", help="Input parameters file name")
    parser.add_argument("--results_folder", default="../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_qi_o0.375", help="Input parameters file name")
    parser.add_argument("--id2name_file_pkl", default="../convert_ids/humap_gene_id_name_map_updated.pkl", help="Input parameters file name")
    parser.add_argument("--prot2url_file", default="../convert_ids/name2uniprotURL_full_humap_from_gene_names.pkl", help="Input parameters file name")
    parser.add_argument("--comp_num2name", default="/res_pred_complex_num2name.pkl", help="Input parameters file name")
    
    args = parser.parse_args()
    
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
            
    # prot2sars_cts = dict([(prot,len(sarss)) for prot,sarss in prot2sars.items()])
    
    # Check - should it be names or gene ids          
    results_folder = args.results_folder     
    results_file = results_folder + "/res_pred_names.out"
    #results_file = "../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_best/res_pred_names.out"
    #results_file = "../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_o0.25/res_pred_names.out"
    
    with open(results_folder+args.comp_num2name,'rb') as f:
        comp_num2name = pickle.load(f)
        
    with open(results_file) as f:
        complexes = [set(line.rstrip().split()[:-1]) for line in f.readlines()]
        # recall order is in decreasing score size
        
    with open(args.id2name_file_pkl,'rb') as f:
        id2name=pickle.load(f)
    # Update gene names if any
    complexes = [set([prot if (prot not in id2name) or (id2name[prot] == '-') else id2name[prot] for prot in comp]) for comp in complexes ]
        
        
    complex_name_list = ['Complex' + str(i) for i in range(1,len(complexes)+1)]
    
    complex2prot = dict(zip(complex_name_list,complexes))
                
    covid_comp2prot = dict()
    
    with open("../convert_ids/covid_interactors.txt") as f:
        covid_prots = set([line.rstrip() for line in f.readlines()])
        
    n_cov_prots = dict()
    for comp, comp_prots in complex2prot.items():
        num_cov_prots = len(covid_prots.intersection(comp_prots))
        if num_cov_prots:
            covid_comp2prot[comp] = comp_prots
            n_cov_prots[comp] = num_cov_prots
            
    write_comp2prot(complex2prot,comp_num2name,suffix='',write_html_file = "./Complexes/Complex2proteins", folder_prefix='',n_cov_prots=dict(),prot2sars=dict(),logfile='logs.txt',prot2url_file=args.prot2url_file)
    write_comp2prot(covid_comp2prot,comp_num2name,suffix='_covid',write_html_file = "./Complexes/Complex2proteins", folder_prefix='',n_cov_prots=n_cov_prots,prot2sars=prot2sars,logfile=results_folder+'/res_metrics_more.out',prot2url_file=args.prot2url_file)
    
if __name__ == "__main__":
    main()