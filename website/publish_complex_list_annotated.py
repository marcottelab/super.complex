# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 10:14:01 2021

@author: Meghana
"""
import pickle
import numpy as np
from publish_complex_list import get_uniprot_urls
import argparse


def get_complex_annot_score(comp,name2annot):
    prot_annot_scores = [name2annot[prot] if prot in name2annot else 0 for prot in comp]
    
    n_sum = len(prot_annot_scores)
    
    min_annot = min(prot_annot_scores)
    median_annot = np.median(prot_annot_scores)
    counts_dict = dict([(i,0) for i in range(0,7)])
    for score in prot_annot_scores:
        counts_dict[score] += 1
    
    char_score = cand_score = 0
    if n_sum:
        char_score = 0.2*sum([key*val for key, val in counts_dict.items()])/n_sum 
    
    n_low = counts_dict[0] + counts_dict[1] + counts_dict[2] + counts_dict[3]
    n_high = counts_dict[4] + counts_dict[5]
    if n_low and n_high:
        cand_score = (counts_dict[5]+0.9*counts_dict[4])*(counts_dict[0] + counts_dict[1] + 0.9*counts_dict[2] + 0.8*counts_dict[3])/(n_low*n_high)
    
    return char_score,cand_score, min_annot, median_annot
    

def write_comp2prot(complex2prot,comp2names,suffix='',write_html_file = "./Complexes/Complex2proteins", folder_prefix='',cov_link="Complex2proteins_covid.html",logfile='logs.txt',prot2url_file="../convert_ids/name2uniprotURL_full_humap_from_gene_names.pkl",name2annot=None):
    
    prefix = ''
    folder_pref = 'Predicted/'
    suffix2 = 'predictions'

    link = 'Protein-wise complex membership list can be viewed <a href="Protein2complex.html"> here</a>. <br> All predicted complexes are listed below with their corresponding proteins. Only SARS-COV2 interacting protein complexes from the predictions can be found <a href=' + cov_link + '> here</a>. <br> Click on a column header to sort. <br>'
    if 'CORUM' in write_html_file:
        link = 'Protein-wise complex membership list can be viewed <a href="CORUM_Protein2complex.html"> here</a>. <br> All CORUM complexes are listed below with their corresponding proteins. Only SARS-COV2 interacting protein complexes from CORUM can be found <a href=' + cov_link + '> here</a>. <br> Click on a column header to sort. <br>'
        folder_pref = folder_prefix+'CORUM/'
        suffix2 = 'CORUM'
        
    if suffix == '_covid':
        prefix = 'SARS-COV2 interacting '
        link = ''
    pre = '''
    
    <!doctype html>
    
    <html>
    
    <head>
    <style>
th{
  border: 1px solid black;
}
</style>
    <title> ''' + prefix + '''Protein complex list - ''' + suffix2 + '''</title>
    
        </head>
        
    
        
        <body>''' + link +  '''
        <section>
    	
     <table style="width:100%">
     <script src="https://www.kryogenix.org/code/browser/sorttable/sorttable.js"></script>

     <table class="sortable">
     <thead>
      <tr>
        <th>''' + prefix + '''Complex </th>
        <th>Protein membership in complex</th>
        <th>Characterization/Annotation score of complex</th>
        <th>Candidate score (for annotated uncharacterized proteins)</th>
        <th>Min Annotation score</th>
        <th>Median Annotation score</th>
      </tr>
      </thead>
       <tbody>
    '''    
    
    post = ''' </tbody>  </table>  </section>
    
    
        </body>
        </html>
        '''
        
    links = get_uniprot_urls(complex2prot,prot2url_file)
    
    for key,val in complex2prot.items():
        complex2prot[key] = dict()
        complex2prot[key]['prots'] = val
        complex2prot[key]['char_score'],complex2prot[key]['cand_score'], complex2prot[key]['min_annot'], complex2prot[key]['median_annot']=get_complex_annot_score(val,name2annot)
    
    complex2prot_items = sorted(complex2prot.items(), key = lambda x: x[1]['min_annot'])   
    with open(write_html_file  + suffix + ".html",'w') as f:
        f.write(pre)
        f.writelines([ '<tr>' + '<td>' + '<a href="' +folder_pref+ key + '.html">' + comp2names[key] + '</a>'  + ' </td> <td>' + ' '.join([links[prot] for prot in val['prots']])+'</td><td>' + "{:0.2f}".format(complex2prot[key]['char_score']) +'</td><td>' + "{:0.2f}".format(complex2prot[key]['cand_score']) + '</td><td>'+"{:0.0f}".format(complex2prot[key]['min_annot'])+ '</td><td>'+"{:0.0f}".format(complex2prot[key]['median_annot']) + '</td></tr>' for key,val in complex2prot_items])
        f.write(post)    


def main():
    parser = argparse.ArgumentParser("Input parameters")
    parser.add_argument("--results_folder", default="../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_qi_o0.375", help="Input parameters file name")
    parser.add_argument("--name2annot_file", default="../convert_ids/name2annot_full_humap_from_gene_names.pkl", help="Input parameters file name")
    parser.add_argument("--prot2url_file", default="../convert_ids/name2uniprotURL_full_humap_from_gene_names.pkl", help="Input parameters file name")
    parser.add_argument("--comp_num2name", default="/res_pred_complex_num2name.pkl", help="Input parameters file name")
    parser.add_argument("--id2name_file_pkl", default="../convert_ids/humap_gene_id_name_map_updated.pkl", help="Input parameters file name")
    
    
    args = parser.parse_args()
      
    # Check - should it be names or gene ids          
    results_folder = args.results_folder
    # results_folder = "../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_qi_o0.375"        
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
        
    for comp, comp_prots in complex2prot.items():
        n_cov_prots = len(covid_prots.intersection(comp_prots))
        if n_cov_prots:
            covid_comp2prot[comp] = comp_prots
    
    with open(args.name2annot_file,'rb') as f:
        name2annot = pickle.load(f)    
            
    write_comp2prot(complex2prot,comp_num2name,suffix='',write_html_file = "./Complexes/Complex2proteins_annotated", folder_prefix='',cov_link="Complex2proteins_annotated_covid.html",logfile='logs.txt',prot2url_file = args.prot2url_file,name2annot=name2annot)
    write_comp2prot(covid_comp2prot,comp_num2name,suffix='_covid',write_html_file = "./Complexes/Complex2proteins_annotated", folder_prefix='',cov_link="Complex2proteins_annotated_covid.html",logfile='logs.txt',prot2url_file = args.prot2url_file,name2annot=name2annot)
    
    
if __name__ == '__main__':
    main()