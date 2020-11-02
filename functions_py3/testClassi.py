# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 18:05:09 2019

@author: Meghana
"""
from numpy import savetxt as np_savetxt, argmax as np_argmax, array as np_array
from sklearn.metrics import average_precision_score as sklearn_metrics_average_precision_score, \
    precision_recall_curve as sklearn_metrics_precision_recall_curve
from matplotlib import use as mpl_use

mpl_use('Agg')  # Issues warning on spyder - don't worry abt it
import matplotlib.pyplot as plt
from logging import info as logging_info
import pandas as pd


def generate_freq_dict(X_pos_test, res_pos, lbl):
    # dictionary with size as key and values as total and correctly classified
    sizes_pos = list(X_pos_test['nodes'])
    freqs = {}
    for i, val in enumerate(sizes_pos):
        val = int(val)
        correct_flag = int(res_pos[i] == lbl)
        if val not in freqs:
            freqs[val] = [1, correct_flag]
        else:
            freqs[val][0] += 1
            freqs[val][1] += correct_flag

    pos_sizes = sorted(freqs.keys())
    pos_accs = [float(freqs[sz][1])/freqs[sz][0] for sz in pos_sizes]
    pos_tots = [freqs[sz][0] for sz in pos_sizes]
    return pos_sizes, pos_accs, pos_tots


def analyze_sizewise_accuracies(X_pos_test, res_pos, X_neg_test, res, out_filename):
    pos_sizes, pos_accs, pos_tots = generate_freq_dict(X_pos_test, res_pos, 1)
    neg_sizes, neg_accs, neg_tots = generate_freq_dict(X_neg_test, res, 0)

    fig = plt.figure()
    plt.plot(pos_sizes, pos_accs, 'yo' )
    for i, txt in enumerate(pos_tots):
        plt.annotate(txt, (pos_sizes[i], pos_accs[i]))
    plt.plot(neg_sizes, neg_accs, 'b+')
    for i, txt in enumerate(neg_tots):
        plt.annotate(txt, (neg_sizes[i], neg_accs[i]))
    plt.legend(['positives', 'negatives'])
    #plt.xticks(range(min(min(pos_sizes), min(neg_sizes))-1,max(max(pos_sizes), max(neg_sizes))+1, 1))
    plt.xlabel('Sizes')
    plt.ylabel('Accuracies')
    plt.title('Size wise accuracies')
    plt.savefig(out_filename + "_annotated.png")
    plt.close(fig)

    fig = plt.figure()
    plt.plot(pos_sizes, pos_accs, 'yo' )
    plt.plot(neg_sizes, neg_accs, 'b+')
    plt.legend(['positives', 'negatives'])
    #plt.xticks(range(min(min(pos_sizes), min(neg_sizes))-1,max(max(pos_sizes), max(neg_sizes))+1, 1))
    plt.xlabel('Sizes')
    plt.ylabel('Accuracies')
    plt.title('Size wise accuracies')
    plt.savefig(out_filename)
    plt.close(fig)


def plot_pr_curve(test_p, test_r, test_aps, out_comp_nm):
    fig = plt.figure()
    plt.plot(test_r, test_p)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title('Precision-Recall curve: AP={0:0.2f}'.format(test_aps))
    plt.savefig(out_comp_nm + '_pr_curve.png')
    # plt.show()            # Does not work on pod
    plt.close(fig)


def calc_metrics(res, res_pos, n_neg, n_pos):
    TP = sum(res_pos)  # assuming negatives are 0s
    FN = n_pos - TP
    if n_pos != 0:
        acc = TP / float(n_pos)
    else:
        acc = 'NA'

    TN = sum([res[ind] == 0 for ind in range(len(res))])
    FP = n_neg - TN
    if n_neg != 0:
        acc_neg = TN / float(n_neg)  # assuming negatives are 0s
    else:
        acc_neg = 'NA'

    Recall = float(TP) / (TP + FN)  # Just accuracy of test positives
    Precision = float(TP) / (TP + FP)
    if Recall == Precision == 0:
        F1_score = 0
    else:
        F1_score = 2 * Precision * Recall / (Precision + Recall)
    return acc, acc_neg, Recall, Precision, F1_score

def test_classi(model, scaler, inputs, X_pos_test, X_neg_test, test_complex_graphs, X_test, y_test):
    logging_info("Evaluating test complexes...")
    model_type = inputs['model_type']
    out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
    res = None
    if model_type == "tpot":
        res_pos = model.predict(X_pos_test)
        res = model.predict(X_neg_test)

        if hasattr(model, 'decision_function'):
            score = model.decision_function(X_pos_test)
            np_savetxt(out_comp_nm + '_test_pos_score.out', score)
            # print("Scores for positive complexes are",score)
            score = model.decision_function(X_neg_test)
            np_savetxt(out_comp_nm + '_test_neg_score.out', score)

            # Write the else case

    elif model_type == "NN":

        X_pos_test = scaler.transform(X_pos_test)

        preds = model.predict(X_pos_test)
        res_pos = [np_argmax(pred) for pred in preds]
        score = np_array([pred[1] for pred in preds])
        np_savetxt(out_comp_nm + '_test_pos_score.out', score)

        X_neg_test = scaler.transform(X_neg_test)
        preds = model.predict(X_neg_test)
        res = [np_argmax(pred) for pred in preds]
        # Score of being negative !!
        score = np_array([pred[0] for pred in preds])
        np_savetxt(out_comp_nm + '_test_pos_score.out', score)
        # print("Scores for negative complexes are",score)

    n_pos = len(test_complex_graphs)
    n_neg = len(X_neg_test)
    
    analyze_sizewise_accuracies(X_pos_test, res_pos, X_neg_test, res, out_comp_nm + '_size_wise_accuracies_test.png')

    acc, acc_neg, Recall, Precision, F1_score = calc_metrics(res,res_pos,n_neg, n_pos)


    if model_type == "tpot":
        test_fit_probs = model.predict_proba(X_test)[:, 1]

    elif model_type == "NN":
        X_test = scaler.transform(X_test)
        preds = model.predict(X_test)
        test_fit_probs = np_array([pred[1] for pred in preds])

    test_aps = sklearn_metrics_average_precision_score(y_test, test_fit_probs)
    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("Test set average precision score = %.3f" % test_aps, file=fid)

    test_p, test_r, _ = sklearn_metrics_precision_recall_curve(y_test, test_fit_probs)
    plot_pr_curve(test_p, test_r, test_aps, out_comp_nm)

    with open(out_comp_nm + '_metrics.out', "a") as fid:
        print("Accuracy for test positive complexes = %.3f" % acc, file=fid)
        print("Accuracy for test negative complexes = %.3f" % acc_neg,
              file=fid)  # Really just tells you complex or not for random graphs
        print("Test Precision = %.3f" % Precision, file=fid)
        print("Test Recall = %.3f" % Recall, file=fid)
        print("Test F1 score = %.3f" % F1_score, file=fid)

    logging_info("Finished evaluating test complexes.")
