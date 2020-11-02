# sup_graph
Installation:
Requires python3 
python3 -m pip install -r requirements_py3.txt --user

Available experiments:
1. for toy network use input_toy.yaml
2. for hu.MAP - use input file input_humap.yaml
3. for hu.MAP2 use input_humap2.yaml, 
4. for yeast, use input_yeast_comb.yaml, 
5. for yeast, for training on TAP complexes and testing on MIPS, use input_yeast_tap.yaml
6. for yeast, for training on MIPS complexes and testing on TAP, use input_yeast_mips.yaml

To run the pipeline on a new network, construct a similar input file (as input_toy.yaml) containing the required inputs.

1. Specify input options relating to the network: 
Set options dir_nm - directory containing the network, netf_nm - file name of the network

2. Specify input options relating to the known communities:  
If you already have separated known communities into train and test communities, specify their paths in the options comf_nm and comf_test_nm (relative to the directory specified in the option:dir_nm) 

Otherwise, 
Split complex list into train and test:
Set option split_flag = 1
Play around with parameters fact and perc_transfer until train test size distributions in figure are the same. Also check that number of training complexes is not too low by looking at the res_metrics.out file.
Set options comf_nm and comf_test_nm with these two files. All the above paths are set relative to the directory specified in the option:dir_nm
Make sure to change the option split_flag back to 0 after this step

Example Bash script to run Super.Complex pipeline after the above 2 steps for yeast, training on TAP complexes and testing on MIPS : 

#!/bin/bash  
mtype=tap  
input_file_name=input_yeast_$mtype.yaml  
out_dir_name=/results_conn_$mtype  
out_dir_name_full=yeast$out_dir_name  
tpot_tmp_dir=TPOT/tpot_tmp_yeast_$mtype  
classifier_file=$out_dir_name_full/res_classifiers.txt  
train_dat=$out_dir_name_full/res_train_dat.csv  
test_dat=$out_dir_name_full/res_test_dat.csv  

echo Reading network...
python3 main_read.py --input_file_name $input_file_name --out_dir_name $out_dir_name  

echo Generating feature matrices for known communities...
python3 train.py --input_file_name $input_file_name --out_dir_name $out_dir_name --mode gen  

echo Finding the best community fitness function...
mkdir $tpot_tmp_dir  
python3 TPOT/train_TPOT3.py --training_data $train_dat --testing_data $test_dat --outfile $out_dir_name_full/res_tpot_best_pipeline.py --outfile_classifiers $classifier_file --outfile_fig $out_dir_name_full/res_classifiers_pr.png --generations 50 --population_size 50 --n_jobs 20 --temp_dir $tpot_tmp_dir  

echo Training the best community fitness function...
python3 train.py --input_file_name $input_file_name --out_dir_name $out_dir_name --train_feat_mat $train_dat --test_feat_mat $test_dat --classifier_file $classifier_file  

meths=( search_top_neigs metropolis isa cliques )  
for meth in "${meths[@]}"  
do  
out_dir_name_meth=$out_dir_name$meth  

echo Partitioning graph nodes across compute nodes...
python3 partition_search_seeds.py --input_file_name $input_file_name --out_dir_name $out_dir_name_meth --search_method $meth  

echo Growing communities...
python3 main_sample.py --input_file_name $input_file_name --out_dir_name $out_dir_name_meth --search_method $meth  

echo Merging very similar communities...
python3 main_postprocess.py --input_file_name $input_file_name --out_dir_name $out_dir_name_meth  

echo Comparing predicted and known communities...
python3 main_eval.py --input_file_name $input_file_name --out_dir_name $out_dir_name_meth  
done  

---xxx---

For each of the scripts, optional arguments can be viewed by running:
python3 script_name.py --help

Add the desired arguments with each of the commands directly on the terminal. If the same argument is present in the input options file, the command line argument value will prevail and be considered in the program.  

The best models for each of the experiments can be made available on request for transferring to different applications. 
