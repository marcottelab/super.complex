# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 11:37:38 2021

@author: Meghana
"""
from generate_protein_wise_complexes import get_prot2complex, write_prot2comp
from sys import path as sys_path
# insert at 1, 0 is the script path (or '' in REPL)
sys_path.insert(1, '../functions_py3/')

from convert_humap_ids2names import read_gene_id_map
from networkx import read_weighted_edgelist,write_weighted_edgelist, relabel_nodes
import pickle
from argparse import ArgumentParser as argparse_ArgumentParser

parser = argparse_ArgumentParser("Input parameters")
parser.add_argument("--suffix_corum", default='original', help="original")
args = parser.parse_args()
    
##suffix_corum = ''
#suffix_corum = 'original'
suffix_corum=args.suffix_corum

def write_corum_complex_edges(corum,id2name,suffix): 
    # overlayed using gene ids on humap with gene ids and then converted to gene names using the map
    G = read_weighted_edgelist('../humap/test_graph.txt')    

    with open('../humap/corum_edges_humap' + suffix + '.txt', "wb") as f_edges:
        f_edges_write = f_edges.write
        
        with open('../humap/corum_edges_humap_names' + suffix + '.txt', "wb") as f_edges_names:
            f_edges_names_write = f_edges_names.write

            for tmp_graph_nodes in corum:
                tmp_graph = G.subgraph(tmp_graph_nodes)
                tmp_graph_names = relabel_nodes(tmp_graph,id2name)
                if len(tmp_graph.edges()) == 0:
                    if tmp_graph_nodes[0] == 'None':
                        id2name[tmp_graph_nodes[0]]='None'                   
                    if id2name[tmp_graph_nodes[0]] == '-':
                        id2name[tmp_graph_nodes[0]]=tmp_graph_nodes[0]
                    dummy_edge_name = str(id2name[tmp_graph_nodes[0]]) + ' ' + str(id2name[tmp_graph_nodes[0]]) + ' 0\n'                    
                    dummy_edge = str(tmp_graph_nodes[0]) + ' ' + str(tmp_graph_nodes[0]) + ' 0\n'
                    f_edges_write(dummy_edge.encode())
                    f_edges_names_write(dummy_edge_name.encode())                    
                else:
                    write_weighted_edgelist(tmp_graph, f_edges)
                    write_weighted_edgelist(tmp_graph_names, f_edges_names)
                
                f_edges_write("\n".encode())
                f_edges_names_write("\n".encode())
    
# ./humap/humap_gene_id_name_map.txt has error entries - check generation of this file - ex:  '55739', 'P)HX', 'dehydratase(NAXD'   

#with open('./convert_ids/humap_gene_id_name_map.txt') as f:
#    id2name2 = dict([line.rstrip().split() for line in f.readlines()[1:]])

# Replace Nones with original names 
id2name = read_gene_id_map(id_name_map_path="../convert_ids/humap_gene_id_name_map.txt")

if suffix_corum != 'original':
    with open('../humap/res_known_complex_nodes_list','rb') as f:
        corum = pickle.load(f)
        
    with open('../humap/cleaned_corum_complex_num2name.pkl','rb') as f:
        comp_num2name_corum = pickle.load(f)           
else:
    with open('../humap/all_complexes.txt') as f:
        corum = [comp.rstrip().split() for comp in f.readlines()]
    with open('../humap/original_corum_complex_num2name.pkl','rb') as f:
        comp_num2name_corum = pickle.load(f)           
        
write_corum_complex_edges(corum,id2name,suffix_corum)
       
corum_names = [[id2name[node] if node in id2name else node for node in comp if node.strip() != 'None'] for comp in corum]
    
with open('../humap/corum_names' + suffix_corum + '.txt','w') as f:
    f.writelines([' '.join(comp) + '\n' for comp in corum_names])

prot2complex = get_prot2complex(corum_names)
            
covid_prot2complex = dict()

with open("../convert_ids/covid_interactors.txt") as f:
    covid_prots = set([line.rstrip() for line in f.readlines()])
    
for key in prot2complex:
    if key in covid_prots:
        covid_prot2complex[key] = prot2complex[key]
              
write_prot2comp(prot2complex,comp_num2name_corum,write_path="../humap/CORUM_latest/res_prot2comp" + suffix_corum,suffix='',write_html_file = "./Complexes/" + suffix_corum + 'CORUM_Protein2complex', folder_pref=suffix_corum)
write_prot2comp(covid_prot2complex,comp_num2name_corum, write_path="../humap/CORUM_latest/res_prot2comp"+ suffix_corum,suffix='_covid',write_html_file = "./Complexes/" + suffix_corum + 'CORUM_Protein2complex',folder_pref=suffix_corum)