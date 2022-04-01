# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 10:14:01 2021

@author: Meghana
"""
import pickle
from os import path as os_path, mkdir as os_mkdir, chdir as os_chdir


def get_uniprot_urls(prot2complex):
    with open("../convert_ids/name2uniprotURL_full_humap_from_gene_names.pkl",'rb') as f:
        prot2url = pickle.load(f)    
        links = dict()
    
    for key in prot2complex:
        if key in prot2url:
            links[key] = '<a href="' + prot2url[key] + '">' + key + '</a>'
        else:
            links[key] = key
            
    return links

def write_prot2comp(prot2complex,comp_num2name,write_path,suffix='',write_html_file = "./Complexes/Protein2complex", folder_pref='Predicted/',cov_link="Protein2complex_covid.html",logfile='logs.txt'):
    if not os_path.exists(write_path):
        with open(write_path ,'w') as f:
            f.writelines([key + ' : ' + ' '.join(val)+'\n' for key,val in prot2complex.items()])
    
    prefix = ''
    suffix2 = 'predictions'

    link = 'All proteins with corresponding predicted complexes are listed below. Only SARS-COV2 interacting proteins and their complexes can be found <a href=' + cov_link + '> here</a>. <br>'
    if 'CORUM' in write_html_file:
        cov_link="CORUM_"+cov_link
        link = 'All proteins with corresponding CORUM complexes are listed below. Only SARS-COV2 interacting proteins and their complexes can be found <a href=' + cov_link + '> here</a>. <br>'
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
    '''
    pre2='''  
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
        
    links = get_uniprot_urls(prot2complex)
    
    score_wise_dicts = dict()
    for prot in prot2complex:
        if prot in name2annot:
            score = name2annot[prot]
        else:
            score = 'Unknown'
        
        if score not in score_wise_dicts:
            score_wise_dicts[score] = dict([(prot,prot2complex[prot])])
        else:
            score_wise_dicts[score][prot] = prot2complex[prot]
    
    with open(write_html_file  + suffix + ".html",'w') as f:
        f.write(pre)  
        
    scores_list = ['Unknown',0,1,2,3,4,5]
    with open(logfile,'a') as f:  
        print('Protein score and number of proteins in ' + folder_pref+write_html_file  + suffix,file=f)                
    
    for score in scores_list:
        if score in score_wise_dicts:
            # sort alphabetical
            prot2complex_items = sorted(score_wise_dicts[score].items())
            
            with open(logfile,'a') as f:  
                print(score,len(prot2complex_items),file=f)
            pre_fin = "\n<b>Score: " + str(score) + ' </b> <br>\n'+ pre2
            with open(write_html_file  + suffix + ".html",'a') as f:
                f.write(pre_fin)  
                f.writelines([ '<tr>' + '<td>' + links[key]  + ' </td> <td>' + ' | '.join(['<a href="' +folder_pref+ node + '.html">' + comp_num2name[node] + '</a>' for node in val])+'</td> </tr>' for key,val in prot2complex_items])
                f.write("</table> <br>\n")
    with open(write_html_file  + suffix + ".html",'a') as f:
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
    
results_folder = "../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_qi_o0.375"        
#results_folder = "../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_best"
results_file = results_folder + "/res_pred_names.out"

with open(results_file) as f:
    complexes = [line.rstrip().split()[:-1] for line in f.readlines()]
    # recall order is in decreasing score size

with open('../convert_ids/humap_gene_id_name_map_updated.pkl','rb') as f:
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
    
with open('../convert_ids/name2annot_full_humap_from_gene_names.pkl','rb') as f:
    name2annot = pickle.load(f)    
        
write_prot2comp(prot2complex,comp_num2name,write_path = results_folder + '/res_pred_prot2comp' +  ".out",suffix='',write_html_file = "./Complexes/Protein2complex_annotated", folder_pref='Predicted/',cov_link="Protein2complex_annotated_covid.html",logfile=results_folder+'/res_metrics_more.out')
write_prot2comp(covid_prot2complex,comp_num2name, write_path = results_folder + '/res_pred_prot2comp' + '_covid' + ".out",suffix='_covid',write_html_file = "./Complexes/Protein2complex_annotated", folder_pref='Predicted/',cov_link="Protein2complex_annotated_covid.html",logfile=results_folder+'/res_metrics_more.out')