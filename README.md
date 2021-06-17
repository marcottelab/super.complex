# Super.Complex

## Installation:
Requires python3  
Requirements installation:
python3 -m pip install -r requirements_py3.txt --user
Works for Python version < 3.8

Note: For faster installation, if you don't plan to use neural networks, you can skip tensorflow installation, by removing the line with tensorflow from the file requirements_py3.txt (None of the experiments provided below need neural networks)

## Available experiments:
1. for toy network use input_toy.yaml
2. for hu.MAP - use input file input_humap.yaml
3. for yeast, use input_yeast_comb.yaml, 
4. for yeast, for training on TAP complexes and testing on MIPS, use input_yeast_tap.yaml
5. for yeast, for training on MIPS complexes and testing on TAP, use input_yeast_mips.yaml

## Tests:
pip install -U pytest

cd /path_to_super.complex_directory

pytest

## Instructions:
To run the pipeline on a new network, construct a similar input file (as input_toy.yaml) containing the required inputs.

1. Specify input options relating to the network: 
Set options dir_nm - directory containing the network, netf_nm - file name of the network

2. Specify input options relating to the known communities:  
If you already have separated known communities into train and test communities, specify their paths in the options comf_nm and comf_test_nm (relative to the directory specified in the option:dir_nm)  
Otherwise, 
Split complex list into train and test:
Set option split_flag = 1
Verify that train test size distributions in figure are the similar. Also check that number of training complexes is not too low by looking at the res_metrics.out file.
Set options comf_nm and comf_test_nm with these two files. All the above paths are set relative to the directory specified in the option:dir_nm
Make sure to change the option split_flag back to 0 after this step

Example Bash script to run Super.Complex pipeline after the above 2 steps:
This is for yeast, training on TAP complexes and testing on MIPS. To run on your network, simply modify the input_file_name and out_dir_name parameters

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
python3 main1_read.py --input_file_name $input_file_name --out_dir_name $out_dir_name  

echo Generating feature matrices for known communities...  
python3 main2_train.py --input_file_name $input_file_name --out_dir_name $out_dir_name --mode gen  

echo Finding the best community fitness function...  
mkdir $tpot_tmp_dir  
python3 TPOT/train_TPOT3.py --training_data $train_dat --testing_data $test_dat --outfile $out_dir_name_full/res_tpot_best_pipeline.py --outfile_classifiers $classifier_file --outfile_fig $out_dir_name_full/res_classifiers_pr.png --generations 50 --population_size 50 --n_jobs 20 --temp_dir $tpot_tmp_dir  

echo Training the best community fitness function...  
python3 main2_train.py --input_file_name $input_file_name --out_dir_name $out_dir_name --train_feat_mat $train_dat --test_feat_mat $test_dat --classifier_file $classifier_file  

meths=( search_top_neigs metropolis isa cliques )  
for meth in "${meths[@]}"  
do  
out_dir_name_meth=$out_dir_name$meth  

echo Partitioning graph nodes across compute nodes...  
python3 main3_partition_search_seeds.py --input_file_name $input_file_name --out_dir_name $out_dir_name_meth --search_method $meth  

echo Growing communities...  
python3 main4_sample.py --input_file_name $input_file_name --out_dir_name $out_dir_name_meth --search_method $meth  

echo Merging very similar communities...  
python3 main5_postprocess.py --input_file_name $input_file_name --out_dir_name $out_dir_name_meth  

echo Comparing predicted and known communities...  
python3 main6_eval.py --input_file_name $input_file_name --out_dir_name $out_dir_name_meth  
done  

python3 get_best_f1_score.py

### Additional tips:
For each of the scripts, optional arguments can be viewed by running:
python3 script_name.py --help

Add the desired arguments with each of the commands directly on the terminal. If the same argument is present in the input options file, the command line argument value will prevail and be considered in the program.  

For further cross validation in each of the methods, in the bash script, loop through different values for explore_prob and prob_metropolis if using metropolis, and T0 and alpha if using iterative simulated annealing. We recommend setting max_size_thres based on domain knowledge of the network, to specify the maximum number of steps the algorithm is allowed to take while growing the community. Alternately, performing a box plot of the sizes of known communities and choosing a large size that is not an outlier works well for this. We recommend setting the use_all_neigs to 1, so that all neighbors of a node are explored before deciding which to choose to grow the community. If this is not efficient enough for your compute system, you can set the parameters thres_neig, min_thres_neig_sorted and perc to limit the number of neighbors explored.


The best models for each of the experiments are available on zenodo at 

https://doi.org/10.5281/zenodo.4814944

These can be used for transferring learning to different applications. The input data for human and yeast data is also available on zenodo. Parameters yielding the best results are specified in the paper and in the input files.

### Website:
Interactive visualizations of results can be constructed by running the file update_everything_webiste.py in the websites folder after specifying paths to learned community results. 
