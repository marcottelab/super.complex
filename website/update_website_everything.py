# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 02:04:38 2021

@author: Meghana
"""
from os import system as os_system
from sys import path as sys_path
import networkx as nx

sys_path.insert(1, '../functions_py3/')
from convert_humap_ids2names import convert2names_wscores

# generate_names_flag = 1
generate_names_flag = 0

# Input variables
sars_cov2_map = "../convert_ids/SARS_COV2_Map.xlsx"
comp_num2name = "/res_pred_complex_num2name.pkl"
name2annot_covid_file = "../convert_ids/name2annot_covid_interactors_from_gene_names.pkl"
prot2url_covid_file = "../convert_ids/name2uniprotURL_covid_interactors_from_gene_names.pkl"
prot2comp_cov_file_name = "/res_pred_prot2comp_covid.out"

# Input variables - hu.MAP
# results_folder = "../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_qi_o0.375"
# name2annot_file = "../convert_ids/name2annot_full_humap_from_gene_names.pkl"
# prot2url_file = "../convert_ids/name2uniprotURL_full_humap_from_gene_names.pkl"
# id2name_file = "../convert_ids/humap_gene_id_name_map.txt"
# id2name_file_pkl = "../convert_ids/humap_gene_id_name_map_updated.pkl"
# dir_ = "../humap/"

# Input variables - hu.MAP2
results_folder = "../humap2/results_allfeats_train_diff_cutoffs_merged_o0.1"
name2annot_file = "../convert_ids/name2annot_full_humap2_from_gene_names.pkl"
prot2url_file = "../convert_ids/name2uniprotURL_full_humap2_from_gene_names.pkl"
id2name_file = "../convert_ids/humap2_gene_id_name_map_updated.txt"
id2name_file_pkl = "../convert_ids/humap2_gene_id_name_map_updated.pkl"
dir_ = "../humap2/"

if generate_names_flag:
    G = nx.read_weighted_edgelist(dir_+'test_graph.txt')
    
    with open(results_folder+ '/res_pred.out', "r") as fn:
        fin_list_graphs = [(set(line.split()[:-1]), float(line.split()[-1])) for line in fn.readlines()]
    
    convert2names_wscores(fin_list_graphs, results_folder + '/res_pred_names.out', G, results_folder + '/res_pred_edges_names.out',id2name_file)
            
    os_system("python" + " ../assign_complex_names_from_corum.py --corum_files_folder ../humap"+ " --results_folder " + results_folder)


os_system("python" + " ./publish_complex_list.py"+ " --sars_cov2_map " + sars_cov2_map+ " --results_folder " + results_folder+ " --id2name_file_pkl " + id2name_file_pkl+ " --prot2url_file " + prot2url_file+ " --comp_num2name " + comp_num2name)
os_system("python" + " ./publish_complex_list_annotated.py"+ " --results_folder " + results_folder+ " --name2annot_file " + name2annot_file+ " --prot2url_file " + prot2url_file+ " --comp_num2name " + comp_num2name+ " --id2name_file_pkl " + id2name_file_pkl)
os_system("python" + " ./generate_protein_wise_complexes.py"+ " --results_folder " + results_folder+ " --id2name_file_pkl " + id2name_file_pkl+ " --prot2url_file " + prot2url_file)
os_system("python" + " ./generate_protein_wise_complexes_sorted_annot_scores.py" + " --results_folder " + results_folder + " --name2annot_file " + name2annot_file + " --prot2url_file " + prot2url_file + " --id2name_file_pkl " + id2name_file_pkl + " --comp_num2name " + comp_num2name)
os_system("python" + " ./generate_complexes_html_predicted_with_sars_cov_links.py"+ " --sars_cov2_map " + sars_cov2_map+ " --name2annot_file " + name2annot_file+ " --results_folder " + results_folder+ " --id2name_file_pkl " + id2name_file_pkl+ " --prot2url_file " + prot2url_file+ " --comp_num2name " + comp_num2name)
os_system("python" + " ./generate_sars2predicted_map_complexes_names.py"+ " --sars_cov2_map " + sars_cov2_map+ " --results_folder " + results_folder+ " --id2name_file " + id2name_file+ " --prot2url_file " + prot2url_file)
os_system("python" + " ./generate_sars_human_predicted_map_with_complexes_only_mapped_prots_sars_cov_protein_wise.py"+ " --sars_cov2_map " + sars_cov2_map + " --results_folder " + results_folder + " --name2annot_file " + name2annot_file + " --name2annot_covid_file " + name2annot_covid_file + " --prot2url_file " + prot2url_file + " --prot2url_covid_file " +  prot2url_covid_file + " --prot2comp_cov_file_name " + prot2comp_cov_file_name + " --comp_num2name " + comp_num2name)

os_system("python" + " ./publish_complex_list_CORUM.py" + " --dir "+ dir_ + " --prot2url_file " + prot2url_file+ " --sars_cov2_map " + sars_cov2_map)
os_system("python" + " ./generate_protein_wise_complexes_CORUM.py"+ " --id2name_file " + id2name_file + " --dir "+ dir_ + " --prot2url_file " + prot2url_file)
os_system("python" + " ./generate_complexes_html_CORUM.py"+ " --name2annot_file " + name2annot_file+ " --dir "+ dir_ + " --prot2url_file " + prot2url_file)

os_system("python" + " ./publish_complex_list_CORUM.py --suffix_corum original"+ " --dir "+ dir_ + " --prot2url_file " + prot2url_file+ " --sars_cov2_map " + sars_cov2_map)
os_system("python" + " ./generate_protein_wise_complexes_CORUM.py --suffix_corum original"+ " --id2name_file " + id2name_file + " --dir "+ dir_ + " --prot2url_file " + prot2url_file)
os_system("python" + " ./generate_complexes_html_CORUM.py --suffix_corum original"+ " --name2annot_file " + name2annot_file+ " --dir "+ dir_ + " --prot2url_file " + prot2url_file)