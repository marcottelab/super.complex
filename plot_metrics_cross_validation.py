# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 13:13:10 2021

@author: Meghana
"""
from matplotlib import use as mpl_use

mpl_use('Agg')  # Issues warning on spyder - don't worry abt it
from glob import glob
from argparse import ArgumentParser as argparse_ArgumentParser
import re
import pandas as pd
from matplotlib.pyplot import figure as plt_figure, savefig as plt_savefig, close as plt_close, plot as plt_plot, legend as plt_legend,xlabel as plt_xlabel,ylabel as plt_ylabel,title as plt_title


parser = argparse_ArgumentParser("Input parameters")
parser.add_argument("--direct", default="humap", help="Input parameters file name")
parser.add_argument("--main_folder", default="/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_qi_o", help="Input parameters file name")
parser.add_argument("--out_file_suffix", default="", help="out files suffix")
parser.add_argument("--parameter", default="overlap_threshold_qi", help="Input parameters file name")
args = parser.parse_args()


def get_metrics(lines):
    metrics = dict()
    metric_keys = ["No. of predicted complexes after removing non-gold std proteins ","MMR F1 score ","Net F1 score ","Unbiased accuracy","Prediction F1 score ","MMR Precision ","MMR Recall ","Prediction Precision ", "Prediction Recall "]
    for line in lines:  # Reverse order
        words = line.strip().split("=")
        for key in metric_keys:
            if words[0] == key:
                if float(words[1]) > 0:
                    metrics[key] = float(words[1])  
                break
    Sn_PPV_line = lines[-9].split()
    metrics['Sn-PPV acc'] = Sn_PPV_line[2]
    metrics['MMR_PWMMR_hmean'] = Sn_PPV_line[-1]
    
    PR_line = lines[-6].split()
    metrics['Precision Recall product'] = PR_line[-1]   
    
    Clique_line = lines[-2].split()
    metrics['F-Grand K-Clique'] = Clique_line[2]
    metrics['F-weighted K-Clique'] = Clique_line[-1]
        
    return metrics


# level1subd = './humap/*/res_metrics*'
allsubd = './' + args.direct + args.main_folder + '*/res_metrics*'
# fname = "./humap/results_73_neg_same_size_distmetropolis/res_metrics.out"
max_f1_score = 0
max_fname = ""
all_sets = []

precisions_qi = []
recalls_qi = []
precisions_MMR = []
recalls_MMR = []

if args.parameter:
    par = args.parameter
for fname in glob(allsubd, recursive=True):
    with open(fname) as f:
        lines = f.readlines()    
    metrics = get_metrics(lines[-35:])
    param_vals = re.findall('\d*\.?\d+',fname)
    over_t = param_vals[-1]
    metrics[par] = over_t
    all_sets.append(metrics)    
    precisions_MMR.append(metrics["MMR Precision "])    
    recalls_MMR.append(metrics["MMR Recall "])    
    precisions_qi.append(metrics["Prediction Precision "])    
    recalls_qi.append(metrics["Prediction Recall "])    
    
# Read humap 2 stage clustering results 
if args.direct == 'humap':
    with open('./humap/results_2stageclustering_comparison/res_metrics.out') as f:
        lines = f.readlines()    
    metrics = get_metrics(lines[-35:])
    metrics[par] = '2stage clustering'
    all_sets.append(metrics)      

o_f = args.out_file_suffix
df = pd.DataFrame(all_sets)
df = df.set_index(par)
df = df.sort_index()
df.to_csv('./' + args.direct + '/' + par + o_f + '_metrics.csv')

if args.direct == 'humap':
    df = df.drop('2stage clustering')
df[list(df.columns)] = df[list(df.columns)].astype(float)
print(df.idxmax())
#print(df)
fig = plt_figure()
df.plot(subplots=True, style='-.',layout=(7,2),figsize=(15,15))  
plt_savefig('./' + args.direct + '/' + par + o_f + '_metrics.png')
plt_close(fig) 

# Plot PR curves for MMR and Qi et al based metrics 
lists = sorted(zip(*[recalls_MMR, precisions_MMR]))
recalls_MMR, precisions_MMR = list(zip(*lists))

lists = sorted(zip(*[recalls_qi, precisions_qi]))
recalls_qi, precisions_qi = list(zip(*lists))
fig = plt_figure()
plt_plot(recalls_MMR,precisions_MMR,'k.-') 
plt_plot(recalls_qi,precisions_qi,'b.-') 
plt_xlabel('Recall')
plt_ylabel('Precision')
plt_title('PR curve over ' + par)
plt_legend(['MMR','Qi'])
plt_savefig('./' + args.direct + '/pr_' + par + o_f + '_MMR_qi.png')
plt_close(fig)  
