# sup_graph
--NEW--
Available experiments:
1. for hu.MAP - use input file input_humap.yaml
2. for toy network use input_toy.yaml
3. for hu.MAP2 use input_humap2.yaml, 
4. for yeast, use input_yeast_comb.yaml, 
5. for yeast, for training on TAP complexes and testing on MIPS, use input_yeast_tap.yaml
6. for yeast, for training on MIPS complexes and testing on TAP, use input_yeast_mips.yaml

Installation:
pip install -r requirements_py3.txt
Example Bash script to run Super.Complex pipeline: 

#!/bin/bash
mtype=tap
input_file_name=input_yeast_$mtype.yaml
out_dir_name=/results_conn_$mtype
out_dir_name_full=yeast$out_dir_name
tpot_tmp_dir=TPOT/tpot_tmp_yeast_$mtype
classifier_file=$out_dir_name_full/res_classifiers.txt
train_dat=$out_dir_name_full/res_train_dat.csv
test_dat=$out_dir_name_full/res_test_dat.csv
python3 main_read.py --input_file_name $input_file_name --out_dir_name $out_dir_name
python3 train.py --input_file_name $input_file_name --out_dir_name $out_dir_name --mode gen
mkdir $tpot_tmp_dir
python3 TPOT/train_TPOT3.py --training_data $train_dat --testing_data $test_dat --outfile $out_dir_name_full/res_tpot_best_pipeline.py --outfile_classifiers $classifier_file --outfile_fig $out_dir_name_full/res_classifiers_pr.png --generations 50 --population_size 50 --n_jobs 20 --temp_dir $tpot_tmp_dir
python3 train.py --input_file_name $input_file_name --out_dir_name $out_dir_name --train_feat_mat $train_dat --test_feat_mat $test_dat --classifier_file $classifier_file
meths=( search_top_neigs metropolis isa cliques )
for meth in "${meths[@]}"
do
out_dir_name_meth=$out_dir_name$meth
python3 partition_search_seeds.py --input_file_name $input_file_name --out_dir_name $out_dir_name_meth --search_method $meth
python3 main_sample.py --input_file_name $input_file_name --out_dir_name $out_dir_name_meth --search_method $meth
python3 main_postprocess.py --input_file_name $input_file_name --out_dir_name $out_dir_name_meth
python3 main_eval.py --input_file_name $input_file_name --out_dir_name $out_dir_name_meth
done



--OLD VERSION - pls ask author if you want this version --
For each step except step 3, run: 
python main.py --input_file_name input_humap.yaml --out_dir_name /results_your_name
The input file input_humap.yaml needs to be modified at each step according to the instruction before running the above command.

1. Split complex list into train and test:
Set option split_flag = 1
Play around with parameters fact and perc_transfer until train test size distributions in figure are the same. Also check that number of training complexes is not too low by looking at the res_metrics.out file.
Set options comf_nm and comf_test_nm with these two files. All the above paths are set relative to the directory specified in the option:dir_nm

2. Generate feature matrices:
Set option split_flag = 0
Set option mode = gen
Option use_full = 1 if you want random walks generated on the full network, else set it to 0

3. Find classifier:
Run train_TPOT.py with the train and test 
python train_TPOT.py --training_data ../humap/results_your_name/res_train_dat.csv --testing_data ../humap/results_your_name/res_test_dat.csv --outfile tpot_output.py --outfile_fig Classifiers.png --outfile_classifiers tpot_classifiers.txt --temp_dir tpot_temp/ --generations 50 --population_size 50 --n_jobs 10 --warm_start

Look at Classifiers.png and tpot_classifiers.txt to decide best classifier and set it in the option: best_classifier
Set classifier_file = tpot_classifiers.txt

4. Search:
Set mode = non_gen 
model_type = tpot
feature matrices paths

For the cliques algorithm search, set seed_mode = cliques
For other algos, it can be all_nodes_known_comp, all_nodes or n_nodes and algo needs to be set in search_method

Requirements installation:
python -m pip install -r requirements.txt --user
Works for Python version < 3.8
