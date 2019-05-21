# sup_graph
For each step except step 3, run: 
python main.py --input_file_name input_humap.yaml --out_dir_name /results_your_name
or for python 3 (if you want to use a neural network):
python3 main3.py --input_file_name input_humap.yaml --out_dir_name /results_your_name 

The input file input_humap.yaml needs to be modified at each step according to the instructions before running the above command.
More information on each of the input options is given in all of the example input files. 

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