# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 15:37:24 2019

@author: Meghana
"""

from seaborn import distplot as sns_distplot
from numpy import zeros as np_zeros, count_nonzero as np_count_nonzero, sum as np_sum, argmax as np_argmax
from logging import info as logging_info
from matplotlib.pyplot import figure as plt_figure, savefig as plt_savefig, close as plt_close, xlabel as plt_xlabel, title as plt_title
from convert_humap_ids2names import convert2names_wscores_matches
from collections import Counter


def write_best_matches(best_matches_for_known,out_comp_nm,dir_nm,suffix):
       
    sorted_matches = sorted(best_matches_for_known,key=lambda x: x[2],reverse=True)
    if dir_nm == "humap":
        convert2names_wscores_matches(sorted_matches, out_comp_nm + suffix + '_known_pred_matches_names.out')
 
    with open(out_comp_nm + suffix + '_known_pred_matches.out', "w") as fn:
        fn_write = fn.write
        fn_write("Known complex nodes ||| Predicted complex nodes ||| Match F1 score ||| Complex score \n")

        for index in range(len(sorted_matches)):
            known_graph_nodes = sorted_matches[index][0]
            pred_graph_nodes = sorted_matches[index][1]
            match_score = sorted_matches[index][2]
            complex_score = sorted_matches[index][3]
            for node in known_graph_nodes:
                fn_write("%s " % node)
            fn_write(" ||| ")
            for node in pred_graph_nodes:
                fn_write("%s " % node)
            
            fn_write(" ||| ")
            fn_write("%.3f" % match_score)
            fn_write(" ||| ")
            fn_write("%.3f" % float(complex_score))            
            fn_write("\n")

def plot_f1_scores(best_matches,out_comp_nm,suffix,prefix):
    # plot histogram of F1 scores
    max_f1_scores = [match[2] for match in best_matches]
    
    avged_f1_score = sum(max_f1_scores)/float(len(max_f1_scores))
    
    f1_score_counts = Counter()
    
    for score in max_f1_scores:
        f1_score_counts[score] += 1
        
    n_perfect_matches = 0
    if 1 in f1_score_counts:
        n_perfect_matches = f1_score_counts[1]
        
    n_no_matches = 0
    if 0 in f1_score_counts:
        n_no_matches = f1_score_counts[0]    
                
    if len(set(max_f1_scores)) > 1:
        fig = plt_figure()
        sns_distplot(max_f1_scores, hist=True)
        plt_xlabel("F1 score")
        plt_title(prefix + "F1 scores distributions")
        plt_savefig(out_comp_nm +suffix+ '_f1_scores_histogram.png')
        plt_close(fig)    
        
    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print(prefix, file=fid)
        print("Averaged F1 score = %.3f" % avged_f1_score, file=fid)
        print("No. of perfectly recalled matches = %d" % n_perfect_matches, file=fid)
        print("No. of matches not recalled at all = %d" % n_no_matches, file=fid)     
    return avged_f1_score
        
    
def one2one_matches(known_complex_nodes_list, fin_list_graphs, N_pred_comp, N_test_comp,out_comp_nm,suffix,dir_nm):

    Metric = np_zeros((N_test_comp, N_pred_comp))

    for i, test_complex in enumerate(known_complex_nodes_list):
        for j, pred_complex in enumerate(fin_list_graphs):
            T = set(test_complex)
            P = pred_complex[0]
            C = len(T.intersection(P))
            A = len(P.difference(T))
            B = len(T.difference(P))
            
            Precision = float(C) / (A + C)
            Recall = float(C) / (B + C)
            
            if Precision == Recall == 0:
                F1_score = 0
            else:
                F1_score = 2 * Precision * Recall / (Precision + Recall)            
            
            Metric[i, j] = F1_score

    max_indices_i = np_argmax(Metric, axis=0)
    best_matches_4predicted = [(fin_list_graphs[j][0],known_complex_nodes_list[i],Metric[i,j],fin_list_graphs[j][1]) for j,i in enumerate(max_indices_i)]

    max_indices_j = np_argmax(Metric, axis=1)
    best_matches_4known = [(fin_list_graphs[j][0],known_complex_nodes_list[i],Metric[i,j],fin_list_graphs[j][1]) for i,j in enumerate(max_indices_j)]
    
    avged_f1_score4known = plot_f1_scores(best_matches_4known,out_comp_nm,'_best4known'+suffix,'Best predicted match for known complexes - ')
    avged_f1_score4pred = plot_f1_scores(best_matches_4predicted,out_comp_nm,'_best4predicted'+suffix,'Best known match for predicted complexes - ')
    
    avg_f1_score = (avged_f1_score4known + avged_f1_score4pred)/2
    
    write_best_matches(best_matches_4known,out_comp_nm,dir_nm,'_best4known' + suffix)
    write_best_matches(best_matches_4predicted,out_comp_nm,dir_nm,'_best4predicted' + suffix)

    return avg_f1_score

def node_comparison_prec_recall(known_complex_nodes_list, fin_list_graphs, N_pred_comp, N_test_comp, p):
    N_matches_test = 0

    Metric = np_zeros((N_test_comp, N_pred_comp))

    for i, test_complex in enumerate(known_complex_nodes_list):
        N_match_pred = 0
        for j, pred_complex in enumerate(fin_list_graphs):
            T = set(test_complex)
            P = pred_complex[0]
            C = len(T.intersection(P))
            A = len(P.difference(T))
            B = len(T.difference(P))

            if float(C) / (A + C) > p and float(C) / (B + C) > p:
                Metric[i, j] = 1
                N_match_pred = N_match_pred + 1

        if N_match_pred > 0:
            N_matches_test = N_matches_test + 1

    Recall = float(N_matches_test) / N_test_comp

    N_matches_pred = np_count_nonzero(np_sum(Metric, axis=0))
    Precision = float(N_matches_pred) / N_pred_comp

    if Precision == Recall == 0:
        F1_score = 0
    else:
        F1_score = 2 * Precision * Recall / (Precision + Recall)

    return Precision, Recall, F1_score

def plot_size_dists(known_complex_nodes_list, fin_list_graphs, sizes_orig, out_comp_nm):
    sizes_known = [len(comp) for comp in known_complex_nodes_list]
    # Size distributions
    sizes_new = [len(comp[0]) for comp in fin_list_graphs]
    fig = plt_figure()
    if len(set(sizes_known)) <= 1:
        return
    sns_distplot(sizes_known, hist=False, label="known")
    if len(set(sizes_orig)) <= 1:
        return
    sns_distplot(sizes_orig, hist=False, label="predicted")
    if len(set(sizes_new)) <= 1:
        return
    sns_distplot(sizes_new, hist=False, label="predicted_known_prots")
    plt_xlabel("Complex Size")
    plt_title("Complex size distributions")
    plt_savefig(out_comp_nm + '_size_dists_known_pred.png')
    plt_close(fig)


def remove_unknown_prots(fin_list_graphs_orig, prot_list):
    # Remove all proteins in predicted complexes that are not present in known complex protein list
    fin_list_graphs = []
    for comp in fin_list_graphs_orig:
        comp = (comp[0].intersection(prot_list), comp[1])

        if len(comp[0]) > 2:  # Removing complexes with only one,two or no nodes
            fin_list_graphs.append(comp)
    return fin_list_graphs


def compute_metrics(known_complex_nodes_list, fin_list_graphs,out_comp_nm,N_test_comp,N_pred_comp,inputs,suffix):

    if N_test_comp != 0 and N_pred_comp != 0:
        Precision, Recall, F1_score = node_comparison_prec_recall(known_complex_nodes_list,fin_list_graphs, N_pred_comp, N_test_comp, inputs["eval_p"])
        
        avg_f1_score = one2one_matches(known_complex_nodes_list, fin_list_graphs, N_pred_comp, N_test_comp,out_comp_nm,suffix,inputs['dir_nm'])
        
        with open(out_comp_nm + '_metrics.out', "a") as fid:     
            print("Net Averaged F1 score = %.3f" % avg_f1_score, file=fid)
            print("Prediction Precision = %.3f" % Precision, file=fid)
            print("Prediction Recall = %.3f" % Recall, file=fid)
            print("Prediction F1 score = %.3f" % F1_score, file=fid)    
    
def eval_complex(rf=0, rf_nm=0, inputs={}, known_complex_nodes_list=[], prot_list=[], fin_list_graphs=[],suffix="both"):
    # rf - read flag to read complexes from file
    logging_info("Evaluating complexes...")
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']

    if rf == 1:
        if rf_nm == 0:
            rf_nm = out_comp_nm + '_pred.out'
        with open(rf_nm) as fn:
            fin_list_graphs = [(set(line.rstrip('\n').split()),1) for line in fn]  # Space separated text only
            # Just list of list of nodes

    sizes_orig = [len(comp[0]) for comp in fin_list_graphs]

    N_pred_comp = len(fin_list_graphs)
    if N_pred_comp == 0:
        return
    N_test_comp = len(known_complex_nodes_list)

    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("No. of known complexes = ", N_test_comp, file=fid) 
        print("No. of predicted complexes = ", N_pred_comp, file=fid)       
        print("\n -- Metrics on complexes with all proteins -- ", file=fid)       
    
    compute_metrics(known_complex_nodes_list, fin_list_graphs, out_comp_nm,N_test_comp,N_pred_comp,inputs,suffix+'_all_prots')            
    
    fin_list_graphs = remove_unknown_prots(fin_list_graphs, prot_list)
    plot_size_dists(known_complex_nodes_list, fin_list_graphs, sizes_orig, out_comp_nm)
    
    N_pred_comp = len(fin_list_graphs)
    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("No. of predicted complexes after removing non-gold std proteins = ", N_pred_comp, file=fid)      
        print("\n -- Metrics on complexes with only gold std proteins -- ", file=fid)   
    
    compute_metrics(known_complex_nodes_list, fin_list_graphs, out_comp_nm,N_test_comp,N_pred_comp,inputs,suffix+'_gold_std_prots')            
    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("-- Finished writing main metrics -- \n", file=fid)   

    logging_info("Finished Evaluating complexes.")
