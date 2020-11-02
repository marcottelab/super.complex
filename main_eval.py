# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 17:36:07 2020

@author: Meg_94
"""
from time import time as time_time

start_time = time_time()

from matplotlib import use as mpl_use

mpl_use('Agg')  # Issues warning on spyder - don't worry abt it
from os import path as os_path, mkdir as os_mkdir, chdir as os_chdir, system as os_system
os_chdir(os_path.dirname(os_path.abspath(__file__)))
from sys import path as sys_path

# insert at 1, 0 is the script path (or '' in REPL)
sys_path.insert(1, './functions_py3/')
from yaml import load as yaml_load, dump as yaml_dump, Loader as yaml_Loader
from argparse import ArgumentParser as argparse_ArgumentParser
from eval_complex import eval_complex
# from random_walk_control import control

from logging import basicConfig as logging_basicConfig, INFO as logging_INFO, DEBUG as logging_DEBUG
from pickle import load as pickle_load
from complex_comparison import *
from main_postprocess import get_prot_list
from sizewise_scores import sizewise_scores
from matplotlib.pyplot import subplots as plt_subplots, savefig as plt_savefig


def run_metrics(gold_standard_complexes, predicted_clusters, out_comp_nm, pref):
    predicted_clusters = [comp[0] for comp in predicted_clusters]
    excluded_complexes = []
    cplx_compare = ComplexComparison(gold_standard_complexes, predicted_clusters, remove_non_gold_standard_proteins=True, exclusion_complexes=excluded_complexes, normalize_by_combinations=True, pseudocount=0.00001)

    #print "Sensitivity: %s" % cplx_compare.sensitivity()
    #print "PPV: %s" % cplx_compare.ppv()
    #print "ACC: %s" % cplx_compare.acc()
    #print "MMR: %s" % cplx_compare.mmr()
    #print "PWMMR: %s" % cplx_compare.pwmmr()
    #print "MMR_PWMMR_hmean: %s" % cplx_compare.mmr_pwmmr_hmean()
    #print "Precision measure: %s" % cplx_compare.precision_measure()
    #print "Recall measure: %s" % cplx_compare.recall_measure()
    #print "Precision Recall product: %s" % cplx_compare.precision_recall_product()
    #ccmm = cplx_compare.clique_comparison_metric_mean()
    #print "Clique Precision Mean: %s Recall Mean: %s" % (ccmm['precision_mean'],ccmm['recall_mean'])
    #ccmm = cplx_compare.clique_comparison_metric_mean(weighted=True)
    #print "Clique Weighted Precision Mean: %s Weighted Recall Mean: %s" % (ccmm['precision_mean'],ccmm['recall_mean'])
    #print "Clique Weighted hmean (F-weighted K-Clique): %s" % (hmean([ccmm['precision_mean'],ccmm['recall_mean']]))

    sensitivity = cplx_compare.sensitivity()
    ppv = cplx_compare.ppv()
    acc = cplx_compare.acc()
    mmr = cplx_compare.mmr()
    pwmmr = cplx_compare.pwmmr()
    pwmmr_hmean = cplx_compare.mmr_pwmmr_hmean()
    precision_measure = cplx_compare.precision_measure()
    recall_measure = cplx_compare.recall_measure()
    precision_recall_product = cplx_compare.precision_recall_product()
    ccmm = cplx_compare.clique_comparison_metric_mean()
    clique_pr_mean = ccmm['precision_mean']
    clique_re_mean = ccmm['recall_mean']
    clique_f1grand = cplx_compare.clique_comparison_metric_grandf1score(mean_func=np.mean)
    wccmm = cplx_compare.clique_comparison_metric_mean(weighted=True)
    clique_weighted_pr_mean = wccmm['precision_mean']
    clique_weighted_re_mean = wccmm['recall_mean']
    clique_weighted_hmean = hmean([wccmm['precision_mean'],wccmm['recall_mean']])
    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("Sensitivity\tPPV\tACC\tMMR\tPWMMR\tMMR_PWMMR_hmean", file=fid)
        print("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (sensitivity, ppv, acc, mmr, pwmmr, pwmmr_hmean), file = fid )
        print("Precision\tRecall\tPrecision Recall product", file=fid)
        print("%.3f\t%.3f\t%.3f\n" % (precision_measure, recall_measure, precision_recall_product), file = fid )
        print("Clique Precision Mean\tRecall Mean\tF-Grand K-Clique\tClique Weighted Precision Mean\tWeighted Recall Mean\tClique Weighted hmean (F-weighted K-Clique)\n", file=fid)
        print("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (clique_pr_mean, clique_re_mean, clique_f1grand, clique_weighted_pr_mean, clique_weighted_re_mean, clique_weighted_hmean), file = fid )

    plot_filename = out_comp_nm + "_eval_metrics_hist" + pref + ".png"

    if plot_filename != None:
        f, subplots = plt_subplots(3)
        subplots[0].hist(cplx_compare.max_matching_ratio_distribution())
        subplots[0].set_title('MMR')
        subplots[1].hist(cplx_compare.sensitivity_distribution())
        subplots[1].set_title('Sensitivity')
        subplots[2].hist(cplx_compare.ppv_distribution())
        subplots[2].set_title('PPV')
        plt_savefig(plot_filename)

def main():
    parser = argparse_ArgumentParser("Input parameters")
    parser.add_argument("--input_file_name", default="input_humap.yaml", help="Input parameters file name")
    parser.add_argument("--out_dir_name", default="/results", help="Output directory name, by default - /results")
    parser.add_argument("--seed_mode", help="Seed mode - specify 'cliques' for the cliques algo")
    parser.add_argument("--search_method", help="Sampling algorithm")
    parser.add_argument("--model_dir", help="Directory containing model")
    parser.add_argument("--python_command", default="python", help="python / python3")
    parser.add_argument("--read_flag", default=0, help="1 when you want to read from file for evaluation")
    parser.add_argument("--complex_file_name", default="humap/results_2stageclustering_comparison/humap_2stage_clustering_res.txt", help="complexes file name")

    args = parser.parse_args()
    rf = args.read_flag
    rf_nm = args.complex_file_name
    
    with open(args.input_file_name, 'r') as f:
        inputs = yaml_load(f, yaml_Loader)

    if args.model_dir:
        inputs['model_dir'] = args.model_dir
    # Override output directory name if same as gen
    if inputs['out_comp_nm'] == "/results/res":
        if not os_path.exists(inputs['dir_nm'] + args.out_dir_name):
            os_mkdir(inputs['dir_nm'] + args.out_dir_name)
        inputs['out_comp_nm'] = args.out_dir_name + "/res"

    with open(inputs['dir_nm'] + inputs['out_comp_nm'] + "_input_eval.yaml", 'w') as outfile:
        yaml_dump(inputs, outfile, default_flow_style=False)

    logging_basicConfig(filename=inputs['dir_nm'] + inputs['out_comp_nm'] + "_logs.yaml", level=logging_INFO)
    # fin_list_graphs = control(myGraph,inputs,n=50)



    # eval_complex(rf,rf_nm,inputs,known_complex_nodes_list,prot_list,myGraph,fin_list_graphs)

    known_complex_nodes_listfname = inputs['dir_nm'] + "/res_known_complex_nodes_list"

    protlistfname = inputs['dir_nm'] + "/res_protlist"
    with open(protlistfname, 'rb') as f:
        prot_list = pickle_load(f)
    with open(known_complex_nodes_listfname, 'rb') as f:
        known_complex_nodes_list = pickle_load(f)
    
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']

    if not rf:
        with open(inputs['dir_nm'] + inputs["out_comp_nm"] + '_pred.out', "r") as fn:
            lines = fn.readlines()
    
        fin_list_graphs = []
        for line in lines:
            words = line.split()
            fin_list_graphs.append((set(words[:-1]), words[-1]))
            N_pred_complexes = len(fin_list_graphs)
            
        with open(out_comp_nm + '_metrics.out', "a") as fid:
            print("No. of predicted complexes = ", N_pred_complexes, file=fid)
        if N_pred_complexes == 0:
            print("0 predicted complexes")
            return            

    pythonCommand = args.python_command
    if rf == 1:
        eval_complex(rf, rf_nm, inputs, known_complex_nodes_list, prot_list)
    else:
        start_time_eval = time_time()
        with open(out_comp_nm + '_metrics.out', "a") as fid:
            print("\n --- On training set ---", file=fid)

        train_complex_path = inputs['dir_nm'] + "/res_train_known_complex_nodes_list"
        try:
            with open(train_complex_path + "_prot_list",'rb') as f:
                train_prot_list = pickle_load(f)
        except:
            train_prot_list = get_prot_list(train_complex_path)
    
        with open(train_complex_path, 'rb') as f:
            train_complex_list = pickle_load(f)
        
        eval_complex(rf, rf_nm, inputs, train_complex_list, train_prot_list, fin_list_graphs,"_train")
        try:
            run_metrics(train_complex_list, fin_list_graphs, out_comp_nm, "_train")
        except:
            print("Error in running additional metrics for train")
            
        test_complex_path = inputs['dir_nm'] + "/res_test_known_complex_nodes_list"
        try:
            with open(test_complex_path + "_prot_list",'rb') as f:
                test_prot_list = pickle_load(f)
        except:
            test_prot_list = get_prot_list(test_complex_path)
    
        with open(test_complex_path, 'rb') as f:
            test_complex_list = pickle_load(f)
            
        with open(out_comp_nm + '_metrics.out', "a") as fid:
            print("\n --- On test set ---", file=fid)
        eval_complex(rf, rf_nm, inputs, test_complex_list, test_prot_list, fin_list_graphs,"_test")
        try:
            run_metrics(test_complex_list, fin_list_graphs, out_comp_nm, "_test")
        except:
            print("Error in running additional metrics for test")
        with open(out_comp_nm + '_metrics.out', "a") as fid:
            print("\n --- On both sets ---", file=fid)
        eval_complex(rf, rf_nm, inputs, known_complex_nodes_list, prot_list, fin_list_graphs,"_both")
        try:
            run_metrics(known_complex_nodes_list, fin_list_graphs, out_comp_nm, "")
        except:
            print("Error in running additional metrics for both")
        if not os_path.exists(out_comp_nm + "_edge_pr_files"):
            os_mkdir(out_comp_nm + "_edge_pr_files")
        for pref in ["", "_train", "_test"]:
	    # model dir not outcompnm
    	    out_comp_nm_model = inputs['dir_nm'] + inputs['model_dir']
    	    results_wprob = out_comp_nm + '_tot_pred_edges_unique_max_comp_prob_inKnown' + pref + '.out'
    	    input_pos = out_comp_nm_model + pref + '_tot_known_edges_unique.out'
    	    outfname = out_comp_nm + "_edge_pr_files/res"+ '_edge_pr_curve' + pref
    	    os_system(pythonCommand + " functions_py3/prcurve_overlay_noneg.py --labels All " + "--results_wprob " + results_wprob +
                      " --input_positives " + input_pos +
                      " --output_file " + outfname)

        fname = out_comp_nm + "_pred.out"
        figname = out_comp_nm + "_sizewise_scores_pred.png"
        sizewise_scores(fname, figname)
        eval_time = time_time() - start_time_eval

        tot_time = time_time() - start_time

        # Write to yaml file instead
        with open(out_comp_nm + '_runtime_performance.out', "a") as fid:
            print("--- Runtime performance ---", file=fid)
            print("Evaluate complex time (s) = ", eval_time, "[", round(100 * float(eval_time) / tot_time, 2), "%]",
                  file=fid)
            print("Total time (s) = ", tot_time, file=fid)


if __name__ == '__main__':
    main()
