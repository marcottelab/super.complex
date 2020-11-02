# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 13:45:21 2020

@author: Meghana
"""

import argparse

'''
def alphabetize(edges):
    alphEdg = []
    for edge in edges:
        words = edge.split()
        ...
        
    return alphEdg
'''


def main():
    parser = argparse.ArgumentParser("Input parameters")
    parser.add_argument("--known_edge_list",
                        default="../humap/results_50step/results_50step_isa_best/res_known_edges.out",
                        help="Known edge list")
    parser.add_argument("--pred_comps", default="../humap/results_50step/results_50step_isa_best/res_pred.out",
                        help="Predicted complex list")
    parser.add_argument("--pred_edge_list",
                        default="../humap/results_50step/results_50step_isa_best/res_pred_edges.out",
                        help="Predicted edge list")
    parser.add_argument("--pred_edge_list_out",
                        default="../humap/results_50step/results_50step_isa_best/res_pred_edges_list_inKnown.out",
                        help="Predicted edge list")
    parser.add_argument("--known_edge_list_out",
                        default="../humap/results_50step/results_50step_isa_best/res_known_edges_list.out",
                        help="Known edge list")

    args = parser.parse_args()
    # separator should be " "
    with open(args.known_edge_list, 'r') as f:
        known_edges = f.readlines()
        known_edges = [edge for edge in known_edges if edge != "\n"]

    known_prot_list = []

    with open(args.known_edge_list_out, 'w') as f:
        for edge in known_edges:
            words = edge.split(" ")
            f.write(words[0] + "\t" + words[1] + "\n")
            known_prot_list.append(words[0])
            known_prot_list.append(words[1])

    known_prot_list = list(set(known_prot_list))

    with open("../known_prot_list.txt", 'w') as f:
        for prot in known_prot_list:
            f.write(prot + "\n")

    with open(args.pred_comps, 'r') as f:
        pred_comps = f.readlines()
        comp_scores = [line.split(" ")[-1] for line in pred_comps if line != "\n"]

    with open(args.pred_edge_list, 'r') as f:
        pred_edges = f.readlines()
        fin = []
        compNo = 0
        for edge in pred_edges:
            if edge == "\n":
                compNo += 1
                continue
            fin.append(edge + " " + comp_scores[compNo])

    '''    
    alphKnown = alphabetize(known_edges)
    alphPred = alphabetize(pred_edges)
    # No need - just sort !
    '''
    '''
    with open(args.pred_edge_list_out, 'w') as f:
        f.write("ID1"+"\t"+"ID2"+"\t"+"P_1"+"\n")
        for edge in fin:
            words = edge.split(" ")
            f.write(words[0]+"\t"+words[1]+"\t"+words[3])
    
    with open(args.pred_edge_list_out, 'w') as f:
        f.write("ID1"+"\t"+"ID2"+"\t"+"P_1"+"\n")
        for edge in fin:
            words = edge.split(" ")
            if words[0] in known_prot_list and words[1] in known_prot_list:
                f.write(words[0]+"\t"+words[1]+"\t"+words[3])    
    '''


main()
