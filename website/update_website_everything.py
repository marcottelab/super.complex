# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 02:04:38 2021

@author: Meghana
"""
from os import system as os_system

os_system("python" + " ./generate_sars2predicted_map_complexes_names.py")
os_system("python" + " ./publish_complex_list.py")
os_system("python" + " ./publish_complex_list_annotated.py")
os_system("python" + " ./generate_sars_human_predicted_map_with_complexes_only_mapped_prots_sars_cov_protein_wise.py")
os_system("python" + " ./generate_protein_wise_complexes.py")
os_system("python" + " ./generate_protein_wise_complexes_sorted_annot_scores.py")
os_system("python" + " ./generate_complexes_html_predicted_with_sars_cov_links.py")
os_system("python" + " ./publish_complex_list_CORUM.py")
os_system("python" + " ./generate_complexes_html_CORUM.py")
os_system("python" + " ./generate_protein_wise_complexes_CORUM.py")
os_system("python" + " ./publish_complex_list_CORUM.py --suffix_corum original")
os_system("python" + " ./generate_complexes_html_CORUM.py --suffix_corum original")
os_system("python" + " ./generate_protein_wise_complexes_CORUM.py --suffix_corum original")