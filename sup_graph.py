# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 21:42:41 2018
Copyright (C) - Unpublished work of Meghana Venkata Palukuri
All Rights Reserved
Do not redistribute without the express permission of the author 
@author: Meg_94
"""
from __future__ import print_function
from joblib import Parallel, delayed
import multiprocessing 
import networkx as nx
import numpy as np
import sklearn.metrics
import matplotlib
matplotlib.use('Agg')     # Issues warning on spyder - don't worry abt it
import matplotlib.pyplot as plt
import seaborn as sns
import operator
from tqdm import tqdm
import random
from scipy.stats import ttest_ind

#from profilehooks import profile
def unwrap_self_isa(arg, **kwarg):
    return Sup_graph.search_isa(*arg, **kwarg)

def unwrap_self_metropolis(arg, **kwarg):
    return Sup_graph.search_metropolis(*arg, **kwarg)

def unwrap_self_top_neigs(arg, **kwarg):
    return Sup_graph.search_top_neigs(*arg, **kwarg)

def unwrap_self_max_neig(arg, **kwarg):
    return Sup_graph.search_max_neig(*arg, **kwarg)

class Sup_graph:
    
    def __init__(self,inputs):
        print("Initializing...")
        self.inputs = inputs
        self.network_path = inputs['dir_nm'] + inputs['netf_nm']
        self.complex_path = inputs['dir_nm'] + inputs['comf_nm']
        self.test_complex_path = inputs['dir_nm'] + inputs['comf_test_nm']
        self.out_comp_nm = inputs['dir_nm'] + inputs['out_comp_nm']
        self.model_type = inputs['model_type']
        self.split_flag = inputs['split_flag']
        self.use_full = inputs['use_full']
        self.model_name = inputs['model_name']
        self.train_feat_mat = inputs['train_feat_mat']
        self.test_feat_mat = inputs['test_feat_mat']
        if self.model_name == "tpot_select":
            self.classifier_file = inputs['classifier_file']
            self.best_classifier = inputs['best_classifier']
        self.sep = inputs['sep']
        self.scale_factor = inputs['scale_factor']
        self.over_t = inputs['over_t']
        self.perc = inputs['perc']
        self.explore_prob = inputs['explore_prob']
        self.max_size = inputs['max_size'] # Will get overwritten by actual max size from data if starting from the start 
        self.all_complexes_path = inputs['dir_nm'] + inputs["comf_nm_all"]
        with open(self.out_comp_nm+'_metrics.out',"w") as fid:
            print("Network:",self.network_path,file = fid)
        print("Finished Initializing")

     
    def random_walk(self,seed_node,n): # in graph G 
        neg_comp = nx.Graph()
        neg_comp.add_node(seed_node)
        node_num = 1
        pres_node = seed_node
        while node_num<n:
            neig_list = self.G[pres_node]
            if not neig_list:
                print("No neig")
                break
            if(len(neig_list) != 1):
            	new_node = np.random.choice(neig_list)
            else:
                temp = neig_list.keys()
                new_node = list(temp)[0]
            wt = neig_list[new_node]
            wt_edge = wt['weight']
            neg_comp.add_edge(pres_node,new_node,weight=wt_edge)
            pres_node = new_node
            node_num = node_num + 1
	    #print(pres_node) 
        return neg_comp
    
    def random_graph(self,n): # random seed 
        seed_node_rand = np.random.choice(self.G.nodes())
        rand_graph = self.random_walk(seed_node_rand,n)
        return rand_graph

    def jaccard_coeff(self,list1, list2):
        inter = len(list(set(list1).intersection(set(list2))))
        union = (len(list1) + len(list2)) - inter
        return float(inter)/union

    def split_train_test_complexes(self):
        with open(self.all_complexes_path) as f:
            raw_lines = f.readlines()
            
        lines1_orig = [line.rstrip("\r\n").split(self.sep) for line in raw_lines if len(line.rstrip("\r\n").split(self.sep)) > 2]
            
        # Remove Nones 
        lines1 = []
        for line in lines1_orig:
            line_noNone = [item for item in line if item != "None"]
            if line_noNone:
                lines1.append(line_noNone)
            
        lines1=np.random.permutation(lines1)   
        
        # Removing redundancy from list: 
        # merge complexes with jaccard similarity greater than 0.6 
        # Loop since it is a forward merging algorithm 
        n_merge = 1
        while (n_merge != 0):
            n_merge = 0
            complex_list = []
            
            lines2 = list(lines1) # Assign as new variable to avoid call by ref
            for line1 in lines1:
                lines2_rem = list(lines2)
                lines2_rem.remove(line1)
                for line2 in lines2_rem:
                    jc = self.jaccard_coeff(line1,line2)
                    
                    if jc >= 0.6:
                        # merge and store in line1
                        n_merge +=1 
                        line1 = line1 + line2
                
                complex_list.append(line1)
            lines1 = list(complex_list)
            print("No. of merges=",n_merge)  
        
        complexes = [self.G.subgraph(complex) for complex in complex_list if len(self.G.subgraph(complex).nodes())>2 and len(self.G.subgraph(complex).edges())>=2 ]                    
        perm_lines=np.random.permutation(complexes)        
        fact = self.inputs['fact'] #0.99
        train_list = [line for line in perm_lines[0:int(round(len(perm_lines)*fact))]]
        test_list = [line for line in perm_lines[int(round(len(perm_lines)*fact)):]]            

        # Start with something that has a biased size distribution !! 
        
        sizes = [len(line) for line in train_list]
        train_mean = np.mean(sizes)
        
        # Transferring some of the smaller complexes to the test list 
        train_list_lower_mean = [line for line in train_list if len(line) < train_mean]
        perc_transfer = self.inputs['perc_transfer'] #0.3 # You can optimize these parameters !
        to_transfer = train_list_lower_mean[:int(round(len(train_list_lower_mean)*perc_transfer))]
        test_list = test_list + to_transfer
        
        # Now remove from train set 
        for line in to_transfer:
            train_list.remove(line)        
            
        # Finding complexes in train that share an edge with a complex in test 
        com_comp = 10
        while (com_comp !=0): # Do until train and test sets are completely separated
        
            # Removing super huge complexes also (nodes >30 ) from test set         
            test_list = [line for line in test_list if len(line)<30]    
            
            # REMOVE OVERLAP B/W TRAIN AND TEST DATA 
            # Remove complexes from train set sharing two proteins with test set 
            train_rem = []
            for train_line in train_list:
                pres = 0
                for test_line in test_list:
                    common = len(list(set(train_line.edges()).intersection(set(test_line.edges))))
                    if common >= 1:
                        pres += 1
                if pres > 0:
                    train_rem.append(train_line)
                        
            com_comp = len(train_rem)
            print("No. of train complexes transferred = ",com_comp)
                    
            for t_line in train_rem:        
                test_list.append(t_line)
        
            for line in train_rem:
                train_list.remove(line)      
                
        self.complex_graphs =  train_list  
        self.test_complex_graphs = test_list

        
        with open(self.out_comp_nm+"_train_complexes_wo2s.txt","w") as f:
            for line in train_list:
                f.write(self.sep.join(line)+"\n")
            
        with open(self.out_comp_nm+"_test_complexes_wo2s.txt","w") as f:
            for line in test_list:
                f.write(self.sep.join(line)+"\n")   
        with open(self.out_comp_nm+'_metrics.out',"a") as fid:
            print("No. of training complexes = ",len(train_list),file = fid)
            print("No. of test complexes = ",len(test_list),file = fid)            
            print("Initial train_test split = ",fact,file = fid)    
            print("Percentage of low sizes transferred from train to test = ",perc_transfer,file = fid)                

 #   @profile
    def read_graphs(self):
        
        import re
        # Extracting subnetwork of the original PPIN containing proteins in the known complexes
        with open(self.complex_path) as f:
            data = f.read()
            id_list_train = re.findall(r"[\w']+", data)
            
        with open(self.test_complex_path) as f:
            data = f.read()
            id_list_test = re.findall(r"[\w']+", data)
        
        prot_list_temp = id_list_train + id_list_test
        
        prot_list = []
        for prot in prot_list_temp:
            if prot not in prot_list:
                prot_list.append(prot)
     
        print("Reading network...")        
        G_full = nx.read_weighted_edgelist(self.network_path,nodetype=str);

        print("Finished Reading network")  
        
        # Removing self loops 
        G_full.remove_edges_from(G_full.selfloop_edges())

        print("Extracting sub network with proteins from complexes...")          

        if self.use_full == 1:
            self.G = G_full
        else:
            self.G =   G_full.subgraph(prot_list)      

        n_nodes_net = nx.number_of_nodes(self.G)
        n_edges_net = nx.number_of_edges(self.G)
        
        with open(self.out_comp_nm+'_metrics.out',"a") as fid:
            print("No. of nodes in network = ",n_nodes_net,file = fid)
            print("No. of edges in network = ",n_edges_net,file = fid)
        
        print("Reading complexes...")
        
        split_flag = self.split_flag
        if (split_flag == 1):
            self.split_train_test_complexes()
        else:
            with open(self.complex_path) as f:
                complex_list = [line.rstrip('\n').split(self.sep) for line in f ]# Space separated text only
                
                #self.complex_graphs = [self.G.subgraph(complex) for complex in complex_list if self.G.subgraph(complex).nodes() ]                    
                
                # Keeping only complexes with nodes greater than 2 FOR TOPOLOGICAL FEATURES ONLY CASE 
    
                self.complex_graphs = [self.G.subgraph(complex) for complex in complex_list if len(self.G.subgraph(complex).nodes())>2 and len(self.G.subgraph(complex).edges())>=2 ]                    
            
            with open(self.test_complex_path) as f:
                complex_list = [line.rstrip('\n').split(self.sep) for line in f ]# Space separated text only
                
                #self.test_complex_graphs = [self.G.subgraph(complex) for complex in complex_list if self.G.subgraph(complex).nodes() ]                    
                
                # Keeping only complexes with nodes greater than 2 FOR TOPOLOGICAL FEATURES ONLY CASE 
                self.test_complex_graphs = [self.G.subgraph(complex) for complex in complex_list if len(self.G.subgraph(complex).nodes())>2 and len(self.G.subgraph(complex).edges())>=2]                    
            # Change later to extract features of complexes not present in network from CORUM edge predictor 
        print("Finished Reading complexes...")  
        
        known_complexes = self.complex_graphs + self.test_complex_graphs
        
        self.known_complex_nodes_list = [list(complex.nodes()) for complex in known_complexes]
        
        known_complex_nodes = [item for sublist in self.known_complex_nodes_list for item in sublist]
        self.prot_list=[]
        for protein in known_complex_nodes:
            if protein not in self.prot_list:
                self.prot_list.append(protein)         
        # Plotting size distributions of complexes in the network 
        train_sizes = [len(comp) for comp in self.complex_graphs]
        test_sizes = [len(comp) for comp in self.test_complex_graphs]        
        fig = plt.figure()
        sns.distplot(train_sizes,hist=False,label="train")
        sns.distplot(test_sizes,hist=False,label="test")
        plt.savefig(self.out_comp_nm+'_size_dists_train_test.png')        
        plt.close(fig)        
        
        # T-test for sizes 
        #tval, pval = ttest_ind(train_sizes, test_sizes)
        #print("pval",pval)
        #print("tval",tval)

      
            
    def create_feat_mat(self,graph_list,n_feats):
        dens_pos = [nx.density(graph) for graph in graph_list]
        nodes_pos = [nx.number_of_nodes(graph) for graph in graph_list]
        
        #CC statistics - mean and max  - faster to use a big loop mostly
        CC_max = [max(nx.clustering(graph).values()) for graph in graph_list]
        CC_mean = [np.mean(list(nx.clustering(graph).values())) for graph in graph_list]
        CC_var = [np.var(list(nx.clustering(graph).values())) for graph in graph_list]
        
        
        # Degree correlation - avg degree of the neighborhood 
        DC_mean= [np.mean(list(nx.average_neighbor_degree(graph).values())) for graph in graph_list]
        DC_var= [np.var(list(nx.average_neighbor_degree(graph).values())) for graph in graph_list]
        DC_max= [max(nx.average_neighbor_degree(graph).values()) for graph in graph_list]
        
        
        #Degree statistics 
        degree_mean=[]
        degree_max = []
        degree_median = []
        degree_var = []
        for graph in graph_list:
            degrees = [val for (node, val) in graph.degree()]
            degree_mean.append(np.mean(degrees))
            degree_median.append(np.median(degrees))
            degree_max.append(max(degrees))
            degree_var.append(np.var(degrees))
        
        # Edge weight statistics 
        edge_wt_mean=[]
        edge_wt_max = []
        edge_wt_var=[]
        for graph in graph_list:
            edge_wts = [val for (node1,node2, val) in graph.edges.data('weight')]
            edge_wt_mean.append(np.mean(edge_wts))
            edge_wt_var.append(np.var(edge_wts))            
            edge_wt_max.append(max(edge_wts))        
            
        # First 3 singular values 
        sv1=[]
        sv2=[]
        sv3=[]
        
        for graph in graph_list:
            A_mat = nx.to_numpy_matrix(graph)
            svs = np.linalg.svd(A_mat,full_matrices=False,compute_uv=False)
            
            if(len(svs) >=3):
                sv1.append(svs[0])
                sv2.append(svs[1])
                sv3.append(svs[2])
            elif(len(svs) >=2):
                sv1.append(svs[0])
                sv2.append(svs[1])
                sv3.append(0)                
            else:
                sv1.append(svs[0])
                sv2.append(0)
                sv3.append(0)                
        
        #feat_mat = np.vstack((dens_pos,nodes_pos)).T
        feat_mat = np.vstack((dens_pos,nodes_pos,degree_max,degree_mean,degree_median,degree_var,CC_max,CC_mean,CC_var,edge_wt_mean,edge_wt_max,edge_wt_var,DC_mean,DC_var,DC_max,sv1,sv2,sv3)).T
        
        #feat_mat = np.vstack((dens_pos,nodes_pos,degree_max,degree_mean,degree_var,CC_max,CC_mean,CC_var,edge_wt_mean,edge_wt_max,edge_wt_var,DC_mean,DC_var,sv2,sv3)).T
        if(n_feats == 1):
            feat_mat = np.array(dens_pos).reshape(-1,1)
            
        return feat_mat
    
  #  @profile
    def feature_extract(self,feats):
        self.n_feats = feats
        
        import pandas as pd
        mode = self.inputs['mode']
        #mode = "non_gen" # Change to gen if you want to generate matrices 
        if self.model_type == "tpot" and mode == "non_gen": # CHANGE X_POS, Y_POS later !!!!
            # Read X,y from file
            df_full = pd.read_csv(self.train_feat_mat)
            self.y = df_full.pop('complex')
            self.X = df_full
            
            neg_start_ind = self.y[self.y == 0].index[0]
            self.X_pos = self.X.iloc[0:neg_start_ind]
            self.y_pos = self.y[0:neg_start_ind]
            self.X_neg = self.X.iloc[neg_start_ind:]
            self.y_neg = self.y[neg_start_ind:]
            
            df_full = pd.read_csv(self.test_feat_mat)
            self.y_test = df_full.pop('complex')
            self.X_test = df_full
            
            neg_start_ind = self.y_test[self.y_test == 0].index[0]
            self.X_pos_test = self.X_test.iloc[0:neg_start_ind]
            self.y_pos_test = self.y_test[0:neg_start_ind]
            self.X_neg_test = self.X_test.iloc[neg_start_ind:]
            self.y_neg_test = self.y_test[neg_start_ind:]
        else:               
    
            print("Feature extraction...")        
            
            self.X_pos = self.create_feat_mat(self.complex_graphs,self.n_feats)
           
            n_pos = len(self.X_pos)
            
            dims = self.X_pos.shape
            n_feats = dims[1]
            with open(self.out_comp_nm+'_metrics.out',"a") as fid:        
                print("No. of features = ",n_feats,file = fid)          
                print("No. of train positive complexes = ",n_pos,file = fid) 
                
            print("Constructing negative complexes...") 
            
            n_pos = len(self.complex_graphs)
            sizes = [len(complex) for complex in self.complex_graphs]
            min_size = min(sizes)
            max_size = max(sizes)
            self.max_size_train = max_size
            n_each = int(np.ceil( n_pos/(max_size - min_size + 1))) # Uniform distribution 
            n_each=int(np.ceil(n_each*self.scale_factor))    
            
            # Negative complexes - random walks from random seeds with incremental steps from min_size to max_size
            neg_comp_list_init = []
            for n in range(min_size,max_size+1): # Generate with same distribution as sizes later
                for k in range(n_each):
                    # Pick a random seed node 
                    #seed_node_rand = np.random.choice(self.G.nodes())
                    
                    #neg_comp = self.random_walk(seed_node_rand,n)
                    neg_comp = self.random_graph(n)
                    #print(neg_comp.nodes())
                    #print(neg_comp.edges())
                    neg_comp_list_init.append(neg_comp)  
            
            # Remove random walks with 1 or 2 nodes 
            self.neg_comp_list = []
            for walk in neg_comp_list_init:
                if len(walk) >= 3:
                    self.neg_comp_list.append(walk)
                    
            print("Finished constructing negative complexes")                  
    
            self.X_neg = self.create_feat_mat(self.neg_comp_list,self.n_feats)        
             
            #Removing negative feature rows that exactly match any row in positives
            for ind in range(n_pos):
                matching_inds = np.where((self.X_neg == self.X_pos[ind]).all(axis=1))
                self.X_neg = np.delete(self.X_neg,matching_inds,axis=0)
                for index in sorted(list(matching_inds[0]), reverse=True):
                    del self.neg_comp_list[index]             
    
            n_neg = len(self.X_neg)
            #print(n_neg)
            # HHANDLE CASE WHEN n_neg = 0 !!!!!
            with open(self.out_comp_nm+'_metrics.out',"a") as fid:                
                print("No. of train negative complexes = ",n_neg,file = fid)
            
            with open(self.out_comp_nm+'_neg_train.out',"w") as fn:
                with open(self.out_comp_nm+'_neg_train_edges.out',"wb") as f_edges:
                    for index in range(len(self.neg_comp_list)):
                        for node in self.neg_comp_list[index].nodes():
                            fn.write("%s " % node)            
                        nx.write_weighted_edgelist(self.neg_comp_list[index],f_edges)
                        f_edges.write("\n".encode())
                        fn.write("\n")
                 
             
            self.X=np.vstack((self.X_pos,self.X_neg))
            
            y_pos = [1 for i in range(n_pos)]
            y_neg = [0 for i in range(n_neg)]            
            y = y_pos + y_neg 
            self.y = np.array(y) 
            self.y_pos = np.array(y_pos) 
            self.y_neg = np.array(y_neg)   
            
            # Writing raw training data to csv in tpot format 
            feat_list = ["dens","nodes","degree_max","degree_mean","degree_median","degree_var","CC_max","CC_mean","CC_var","edge_wt_mean","edge_wt_max","edge_wt_var","DC_mean","DC_var","DC_max","sv1","sv2","sv3","complex"]
    
            dat = np.hstack((self.X,self.y[:,None]))
            import pandas as pd 
            df = pd.DataFrame(dat)
            df.to_csv(path_or_buf=self.out_comp_nm+"_train_dat.csv",index=False,header=feat_list)
            
            X_pos_test = self.create_feat_mat(self.test_complex_graphs,self.n_feats)
                    
            # Negative tests
            
            #N_neg = len(self.test_complex_graphs)
            #neg_test_comp = [self.random_graph(n) for n in range(N_neg)]
            n_pos = len(self.test_complex_graphs)
            sizes = [len(complex) for complex in self.test_complex_graphs]
            min_size = min(sizes)
            max_size = max(sizes)
            self.max_size_test = max_size        
            
            self.max_size = max(self.max_size_test,self.max_size_train)
            with open(self.out_comp_nm+'_metrics.out',"a") as fid:        
                print("Max size of complexes = ",self.max_size,file = fid)  
                
            n_each = int(np.ceil( n_pos/(max_size - min_size + 1)))
            n_each=int(np.ceil(n_each*self.scale_factor))
            print("Constructing test negative complexes...")        
    
            # Negative complexes - random walks from random seeds with incremental steps from min_size to max_size
            neg_test_comp_init = []
            for n in range(min_size,max_size+1): # Generate with same distribution as sizes later
                for k in range(n_each):
                    # Pick a random seed node 
                    #seed_node_rand = np.random.choice(self.G.nodes())
                    
                    #neg_comp = self.random_walk(seed_node_rand,n)
                    neg_comp = self.random_graph(n)
                    #print(neg_comp.nodes())
                    #print(neg_comp.edges())
                    neg_test_comp_init.append(neg_comp) 
                    
            # Remove random walks with 1 or 2 nodes 
            neg_test_comp = []
            for walk in neg_test_comp_init:
                if len(walk) >= 3:
                    neg_test_comp.append(walk)                
                    
            print("Finished constructing test negative complexes.")     
        
            
            X_neg_test = self.create_feat_mat(neg_test_comp,self.n_feats)
            
            self.X_allpos=np.vstack((self.X_pos,X_pos_test))
            n_allpos = len(self.X_allpos)
            
            #Removing negative feature rows that exactly match any row in positives
            for ind in range(n_allpos):
                matching_inds = np.where((X_neg_test == self.X_allpos[ind]).all(axis=1))
                X_neg_test = np.delete(X_neg_test,matching_inds,axis=0) 
                for index in sorted(list(matching_inds[0]), reverse=True):
                    del neg_test_comp[index] 
                    
            with open(self.out_comp_nm+'_neg_test.out',"w") as fn:
                with open(self.out_comp_nm+'_neg_test_edges.out',"wb") as f_edges:
                    for index in range(len(neg_test_comp)):
                        for node in neg_test_comp[index].nodes():
                            fn.write("%s " % node)            
                        nx.write_weighted_edgelist(neg_test_comp[index],f_edges)
                        fn.write("\n")
                        f_edges.write("\n".encode())                    
    
    
            n_neg = len(X_neg_test)
            
            y_pos = [1 for i in range(n_pos)]
            y_neg = [0 for i in range(n_neg)]            
            y = y_pos + y_neg 
            self.y_test = np.array(y)         
            
            self.X_test=np.vstack((X_pos_test,X_neg_test))    
            self.X_pos_test = X_pos_test
            self.X_neg_test = X_neg_test
            
            # Writing raw test data to csv in tpot format 
            dat = np.hstack((self.X_test,self.y_test[:,None]))
            df = pd.DataFrame(dat)
            df.to_csv(path_or_buf=self.out_comp_nm+"_test_dat.csv",index=False,header=feat_list)
    
            with open(self.out_comp_nm+'_metrics.out',"a") as fid:        
                print("No. of test positive complexes = ",n_pos,file = fid)     
                print("No. of test negative complexes = ",n_neg,file = fid)  
                
        print("Finished Feature extraction")        
        

   # @profile
    def train_model(self,model_name):
        
        if(self.model_type == "tpot"):
            print("Training model...",self.model_type)        
            
            from sklearn.pipeline import make_pipeline
            	
            if(model_name == "tpot_select"):
                best_classifier = self.best_classifier 
                from tpot import TPOTClassifier
                from deap import creator
                import re
                tpot = TPOTClassifier()
                with open(self.classifier_file) as f:
                    raw_lines = f.readlines()
                    words = [line.rstrip("\n").split("\t") for line in raw_lines[1:]]
                    classifiers = [word[0] for word in words]
                    pipelines = [word[2] for word in words]
                    
                for i,classi in enumerate(classifiers):
                    if classi == best_classifier:
                        # if model is linearsvc then convert to svc
                        pipeline_string = pipelines[i]
                        #print(pipeline_string)
                        # convert pipeline string to scikit-learn pipeline object
                        deap_pipeline = creator.Individual.from_string(pipeline_string, tpot._pset)
                        clf = tpot._toolbox.compile(expr=deap_pipeline)		
                        if classi == "LinearSVC":
                            #print(clf)
                            n = len(clf.steps)
                            linsvc = str(clf.steps.pop(n-1))
                            #print(linsvc)
                            #print(clf)
                            match = re.search(r"C=(\d*.\d*)",linsvc)
                            C_val = float(match.group(1))
                            #print(C_val)
                            from sklearn.svm import SVC
                            clf.steps.append(('svc',SVC(kernel='linear',probability=True,C=C_val,tol=1e-05)))
                            #print(clf)   
            
            if(model_name == "SVM"):
                print("Training model...",model_name)                 
                # Imports from tpot output             
                from sklearn.preprocessing import StandardScaler
                #from sklearn.svm import LinearSVC
                from sklearn.svm import SVC
                
                # Pipeline from tpot 
                #clf = make_pipeline(StandardScaler(), LinearSVC(random_state=0, tol=1e-5))
                # Cross validate with C vals - default is 1
                # LinearSVC does not have a predict_proba function 
                clf = make_pipeline(StandardScaler(), SVC(kernel='linear',probability=True,random_state=0, tol=1e-5))
            elif(model_name == "estimator_SVM"):            
                '''
                from sklearn.ensemble import GradientBoostingClassifier
                from sklearn.pipeline import make_pipeline, make_union
                #from sklearn.svm import LinearSVC
                from sklearn.svm import SVC
                from sklearn.tree import DecisionTreeClassifier
                from tpot.builtins import StackingEstimator
                from sklearn.preprocessing import FunctionTransformer
                from copy import copy

                # Score on the training set was:0.972240704784
                #clf = make_pipeline(make_union(StackingEstimator(estimator=GradientBoostingClassifier(learning_rate=0.001, max_depth=1, max_features=0.65, min_samples_leaf=16, min_samples_split=11, n_estimators=100, subsample=0.5)),FunctionTransformer(copy)),StackingEstimator(estimator=GradientBoostingClassifier(learning_rate=0.001, max_depth=1, max_features=0.65, min_samples_leaf=16, min_samples_split=11, n_estimators=100, subsample=0.5)),StackingEstimator(estimator=DecisionTreeClassifier(criterion="gini", max_depth=4, min_samples_leaf=16, min_samples_split=14)),LinearSVC(C=5.0, dual=False, loss="squared_hinge", penalty="l2", tol=1e-05))

                clf = make_pipeline(make_union(StackingEstimator(estimator=GradientBoostingClassifier(learning_rate=0.001, max_depth=1, max_features=0.65, min_samples_leaf=16, min_samples_split=11, n_estimators=100, subsample=0.5)),FunctionTransformer(copy)),StackingEstimator(estimator=GradientBoostingClassifier(learning_rate=0.001, max_depth=1, max_features=0.65, min_samples_leaf=16, min_samples_split=11, n_estimators=100, subsample=0.5)),StackingEstimator(estimator=DecisionTreeClassifier(criterion="gini", max_depth=4, min_samples_leaf=16, min_samples_split=14)),SVC(kernel='linear',probability=True,C=5.0,tol=1e-05))
                '''
                from sklearn.ensemble import GradientBoostingClassifier
                from sklearn.feature_selection import SelectFwe, f_classif
                from sklearn.linear_model import LogisticRegression
                from sklearn.pipeline import make_pipeline, make_union
                #from sklearn.svm import LinearSVC
                from tpot.builtins import StackingEstimator
                from xgboost import XGBClassifier

                # Score on the training set was:0.968003998605
                #clf = make_pipeline(StackingEstimator(estimator=GradientBoostingClassifier(learning_rate=0.1, max_depth=9, max_features=0.05, min_samples_leaf=2, min_samples_split=17, n_estimators=100, subsample=1.0)),SelectFwe(score_func=f_classif, alpha=0.02),StackingEstimator(estimator=LogisticRegression(C=1.0, dual=True, penalty="l2")),StackingEstimator(estimator=XGBClassifier(learning_rate=0.001, max_depth=7, min_child_weight=16, n_estimators=100, nthread=1, subsample=0.65)),LinearSVC(C=1.0, dual=True, loss="squared_hinge", penalty="l2", tol=0.001))

                clf = make_pipeline(StackingEstimator(estimator=GradientBoostingClassifier(learning_rate=0.1, max_depth=9, max_features=0.05, min_samples_leaf=2, min_samples_split=17, n_estimators=100, subsample=1.0)),SelectFwe(score_func=f_classif, alpha=0.02),StackingEstimator(estimator=LogisticRegression(C=1.0, dual=True, penalty="l2")),StackingEstimator(estimator=XGBClassifier(learning_rate=0.001, max_depth=7, min_child_weight=16, n_estimators=100, nthread=1, subsample=0.65)),SVC(kernel='linear',probability=True,C=1.0,tol=0.001))
            elif(model_name == "log_reg"):  
                print("Training model...",model_name)                                 
                # Imports from tpot output             
                from sklearn.ensemble import ExtraTreesClassifier
                from sklearn.linear_model import LogisticRegression
                from tpot.builtins import StackingEstimator, ZeroCount            
    
                # Pipeline from tpot 
                # Score on humap was:0.986160063433
                clf = make_pipeline(ZeroCount(),StackingEstimator(estimator=ExtraTreesClassifier(bootstrap=False, criterion="entropy", max_features=0.6, min_samples_leaf=4, min_samples_split=6, n_estimators=100)),    LogisticRegression(C=15.0, dual=False, penalty="l2"))            

            elif(model_name == "extra_trees"):
                from sklearn.ensemble import ExtraTreesClassifier
                from tpot.builtins import StackingEstimator
                '''
                from sklearn.naive_bayes import GaussianNB
                from sklearn.svm import LinearSVC
                
                # Score on the training set was:0.980618836533
                clf = make_pipeline(StackingEstimator(estimator=GaussianNB()),StackingEstimator(estimator=LinearSVC(C=1.0, dual=False, loss="squared_hinge", penalty="l1", tol=0.001)),ExtraTreesClassifier(bootstrap=False, criterion="entropy", max_features=0.7, min_samples_leaf=1, min_samples_split=16, n_estimators=100))
                '''
                from sklearn.pipeline import make_pipeline, make_union
                from sklearn.preprocessing import Normalizer
                from sklearn.preprocessing import FunctionTransformer
                from copy import copy

                # Score on the training set was:0.948305771055
                clf = make_pipeline(make_union(FunctionTransformer(copy),make_pipeline(StackingEstimator(estimator=ExtraTreesClassifier(bootstrap=False, criterion="gini", max_features=0.25, min_samples_leaf=8, min_samples_split=11, n_estimators=100)),Normalizer(norm="l1"))),StackingEstimator(estimator=ExtraTreesClassifier(bootstrap=False, criterion="entropy", max_features=0.75, min_samples_leaf=15, min_samples_split=18, n_estimators=100)),ExtraTreesClassifier(bootstrap=True, criterion="entropy", max_features=0.85, min_samples_leaf=5, min_samples_split=4, n_estimators=100))

            elif(model_name == "rand_forest"):  
                print("Training model...",model_name)                                 
                # Imports from tpot output             
                from sklearn.ensemble import RandomForestClassifier
                from sklearn.feature_selection import VarianceThreshold
                from sklearn.preprocessing import PolynomialFeatures           
    
                # Pipeline from tpot 
                # Score on humap was:0.986160063433
                clf = make_pipeline(VarianceThreshold(threshold=0.05),PolynomialFeatures(degree=2, include_bias=False, interaction_only=False),RandomForestClassifier(bootstrap=False, criterion="entropy", max_features=0.35, min_samples_leaf=1, min_samples_split=11, n_estimators=100)
)            

            clf.fit(self.X, self.y)
            
            print("Finished Training model")        
            print("Evaluating training accuracy...")        
            #Training accuracy 
            
            
            acc_overall_train = clf.score(self.X,self.y)
            acc_pos_train = clf.score(self.X_pos,self.y_pos)
            acc_neg_train = clf.score(self.X_neg,self.y_neg)
            
            train_fit_probs = clf.predict_proba(self.X)[:,1]
            train_aps = sklearn.metrics.average_precision_score(self.y,train_fit_probs)
            with open(self.out_comp_nm+'_metrics.out',"a") as fid:          
                print("Training set average precision score = ",train_aps,file = fid)                        
            
            #res = clf.predict(self.X_neg)
            #print(res)
            #TN = sum([res[ind] == 0 for ind in range(len(res))])
            #acc_neg = TN /float(len(self.X_neg)) # assuming negatives are 0s
            #print("Accuracy for train negative complexes = ",acc_neg) #Really just tells you complex or not for random graphs
         
            # Checking results 
            #X_new = np.array([1,0.9,0.8,0.75,0.7,0.6,0.5]) # Rescale this now !
            # array([1, 1, 1, 0, 0, 0, 0])
            #print(clf.predict(X_new.reshape(-1,1)))
            
            self.model = clf
            
            if hasattr(self.model, 'decision_function'):
                score = self.model.decision_function(self.X_neg)
                np.savetxt(self.out_comp_nm+'_train_neg_score.out',score)
                score = self.model.decision_function(self.X_pos)
                np.savetxt(self.out_comp_nm+'_train_pos_score.out',score)
            
        elif(self.model_type == "NN"):
            
            # Standardizing the feature matrix 
            from sklearn import preprocessing 
            self.scaler = preprocessing.StandardScaler().fit(self.X)
            
            self.X = self.scaler.transform(self.X)
            
            # Scaling X_pos and X_neg as well now for testing with them later
            self.X_pos = self.scaler.transform(self.X_pos)
            self.X_neg = self.scaler.transform(self.X_neg)   
            
            import tensorflow as tf
            from tensorflow import keras
            
            #tf.enable_eager_execution() # Fix ensuing errors 

            
            print("Training model...",self.model_type)

            # multi-layer perceptron
            #for most problems, one could probably get decent performance (even without a second optimization step) by setting the hidden layer configuration using just two rules: (i) number of hidden layers equals one; and (ii) the number of neurons in that layer is the mean of the neurons in the input and output layers. 
            print()
            dims = self.X.shape
            n_feats = dims[1]
            n_classes = 2
            print("No. of nodes in input layer = ", n_feats)
            print("No. of nodes in output layer (since softmax) = ", n_classes)
            hidden_nodes = int((n_feats + n_classes)/2)
            print("No. of nodes in the one hidden layer = ", hidden_nodes)            
            model = keras.Sequential([keras.layers.Dense(n_feats, activation = tf.nn.relu),keras.layers.Dense(hidden_nodes, activation = tf.nn.relu), keras.layers.Dense(n_classes, activation = tf.nn.softmax)])                  
            #model = keras.Sequential([keras.layers.Dense(n_feats, activation = tf.nn.relu), keras.layers.Dense(n_classes, activation = tf.nn.softmax)])                  
            model.compile(optimizer='adam',loss = 'sparse_categorical_crossentropy',metrics=['accuracy'])
            N_epochs = 1000
            model.fit(self.X, self.y, epochs = N_epochs)
            with open(self.out_comp_nm+'_metrics.out',"a") as fid:          
                print("No. of epochs = ",N_epochs,file = fid)
            
            print("Finished Training model")        
            print("Evaluating training accuracy...")
            loss_overall, acc_overall_train = model.evaluate(self.X,self.y)
            loss_pos, acc_pos_train = model.evaluate(self.X_pos,self.y_pos)
            loss_neg, acc_neg_train = model.evaluate(self.X_neg,self.y_neg) 

            self.model = model  

        print("Finished Evaluating training accuracy.")  
        with open(self.out_comp_nm+'_metrics.out',"a") as fid:          
            print("Accuracy overall train = ",acc_overall_train,file = fid)
            print("Accuracy positive train = ",acc_pos_train,file = fid)
            print("Accuracy negative train = ",acc_neg_train,file = fid)                    
    
    #@profile
    def test(self):
        print("Evaluating test complexes...")        
        
        if(self.model_type == "tpot"):
            res_pos = self.model.predict(self.X_pos_test)
            res = self.model.predict(self.X_neg_test)            
            
            if hasattr(self.model, 'decision_function'):            
                score = self.model.decision_function(self.X_pos_test)
                np.savetxt(self.out_comp_nm+'_test_pos_score.out',score)
                #print("Scores for positive complexes are",score)
                score = self.model.decision_function(self.X_neg_test)
                np.savetxt(self.out_comp_nm+'_test_neg_score.out',score)                
                
            # Write the else case 
            
        elif(self.model_type == "NN"):

            X_pos_test = self.scaler.transform(self.X_pos_test)        
            
            preds = self.model.predict(X_pos_test)
            res_pos = [np.argmax(pred) for pred in preds]
            score = np.array([pred[1] for pred in preds])
            np.savetxt(self.out_comp_nm+'_test_pos_score.out',score)
            
            X_neg_test = self.scaler.transform(self.X_neg_test)             
            preds = self.model.predict(X_neg_test)
            res = [np.argmax(pred) for pred in preds]            
            # Score of being negative !!
            score = np.array([pred[0] for pred in preds])
            np.savetxt(self.out_comp_nm+'_test_pos_score.out',score)  
        #print("Scores for negative complexes are",score)
            
        n_pos = len(self.test_complex_graphs)
        n_neg = len(self.X_neg_test)
                   
        TP = sum(res_pos) # assuming negatives are 0s
        FN = n_pos - TP
        acc = TP/float(n_pos) 
        
        TN = sum([res[ind] == 0 for ind in range(len(res))])
        FP = n_neg - TN
        acc_neg = TN /float(n_neg) # assuming negatives are 0s
        
        Recall = float(TP)/(TP+FN) # Just accuracy of test positives
        Precision = float(TP)/(TP + FP)
        F1_score = 2*Precision*Recall/(Precision + Recall)
        
        if(self.model_type == "tpot"):
            test_fit_probs = self.model.predict_proba(self.X_test)[:,1]
            
        elif(self.model_type == "NN"):
            self.X_test = self.scaler.transform(self.X_test)             
            preds = self.model.predict(self.X_test)
            test_fit_probs = np.array([pred[1] for pred in preds])           
            
        test_aps = sklearn.metrics.average_precision_score(self.y_test,test_fit_probs)
        with open(self.out_comp_nm+'_metrics.out',"a") as fid:          
            print("Training set average precision score = ",test_aps,file = fid)    
        
        test_p, test_r, _ = sklearn.metrics.precision_recall_curve(self.y_test, test_fit_probs)
            

        fig = plt.figure()        
        plt.plot(test_r,test_p)
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0.0, 1.05])
        plt.xlim([0.0, 1.0])
        plt.title('Precision-Recall curve: AP={0:0.2f}'.format(test_aps))
        plt.savefig(self.out_comp_nm+'_pr_curve.png')
        #plt.show()            # Does not work on pod 
        plt.close(fig)        

        with open(self.out_comp_nm+'_metrics.out',"a") as fid:        
            print("Accuracy for test positive complexes = ",acc,file = fid)                   
            print("Accuracy for test negative complexes = ",acc_neg,file = fid) #Really just tells you complex or not for random graphs    
            print("Test Precision = ",Precision,file = fid)
            print("Test Recall = ",Recall,file = fid)
            print("Test F1 score = ",F1_score,file = fid)
        

        print("Finished evaluating test complexes.")        
        
    def get_score(self,g1):
        g_list=[g1]
        feats = self.create_feat_mat(g_list,self.n_feats)
            
        if(self.model_type == "tpot"):
            comp_bool = self.model.predict(feats)
            score_curr = self.model.predict_proba(feats)[:,1]

        elif(self.model_type == "NN"):
            feats = self.scaler.transform(feats)  
                
            preds = self.model.predict(feats)
            pred = preds[0]
            comp_bool = np.argmax(pred)
            score_curr = pred[1]
            print(score_curr)
        
        return (score_curr,comp_bool)
        
    def search_max_neig(self,seed_node):
        
        # Seed node
        #print("Seed node is",seed_node)
        
        # Assigning original graph to temporary variable
        tempG = self.G.copy()
        
        neig_list= tempG[seed_node]
        if not neig_list:
            return []
        imp_neig = max(neig_list) # Largest weight neighbor - gives the most confident graphs 
        wt = neig_list[imp_neig]
        wt_edge = wt['weight']
        
        score_curr = 0
        g1=nx.Graph(comp_score=score_curr)
        g1.add_edge(seed_node,imp_neig,weight=wt_edge)
        tempG.remove_edge(seed_node,imp_neig)
        
        max_nodes = self.max_size
        
        
        while True:

            #print("Adding next node")
                
            imp_neigs = dict()
            for node in g1.nodes():
                # get its max neighbor and weight and store in dict 
                neig_list= tempG[node]
                if not neig_list: # Checking if empty
                    break
                imp_neig = max(neig_list)
                wt = neig_list[imp_neig]
                wt_edge = wt['weight']
                imp_neigs[imp_neig] = wt_edge
                    
            if not imp_neigs:
                #print("No more neighbors to add")
                break
            
            node_to_add = max(imp_neigs) # Check again that this is the max 
            #ADD ALL EDGES OF NEW NODE TO ORIG GRAPH
            its_neig_list= tempG[node_to_add]
                
            orig_nodes = g1.nodes()
            for node in list(orig_nodes):
                if node in its_neig_list:
                    wt = its_neig_list[node]
                    wt_edge = wt['weight']
                    g1.add_edge(node_to_add,node,weight=wt_edge)
                    tempG.remove_edge(node_to_add,node)
                
            if len(g1) > max_nodes:
                print("Max size exceeded")
                break
            
            score_prev = score_curr       
            
            (score_curr,comp_bool) = self.get_score(g1)

            if comp_bool == 0:
            #print("Complex found")
            
                # Remove the node last added                
                g1.remove_node(node_to_add)
                score_curr = score_prev
                break
            
        g1.graph['comp_score'] = score_curr
        
        #print(g1.nodes())
        #print(g1.edges())
        
        return g1

    def find_imp_neig(self,neig_list,g1,tempG):
        sorted_neig_list = sorted(neig_list.items(), key=operator.itemgetter(1))
        #perc = 0.5
        n_maxs = int(np.ceil(self.perc*len(neig_list)))
        imp_neigs = dict(sorted_neig_list[-n_maxs:])
        
        if len(imp_neigs) == 1:
            imp_neig = imp_neigs.keys()[0]
        else:
            explore_prob = self.explore_prob  
            cur_trial = np.random.uniform(low=0.0,high=1.0)
            if  cur_trial <= explore_prob: 
                print("Exploring with low probability") # Move to top for efficiency and remove del max
                imp_neig = random.choice(list(imp_neigs.keys()))
            else:
                scores = {}
                for neig in imp_neigs:
                    # Add to graph 
                    wt = imp_neigs[neig]
                    wt_edge = wt['weight']          
                    temp_g1 = g1.copy()
                    node_to_add = neig
                    #ADD ALL EDGES OF NEW NODE TO ORIG GRAPH
                    its_neig_list= tempG[node_to_add]
                                
                    orig_nodes = temp_g1.nodes()
                    for node in list(orig_nodes):
                        if node in its_neig_list:
                            wt = its_neig_list[node]
                            wt_edge = wt['weight']
                            temp_g1.add_edge(node_to_add,node,weight=wt_edge)
                    # Check score
                    (score_curr,comp_bool) = self.get_score(temp_g1)         
                    scores[neig] = score_curr
                    
                max_neig = max(scores.iteritems(), key=operator.itemgetter(1))[0]  
                imp_neig = max_neig          

        return(imp_neig)
    
    def search_top_neigs(self,seed_node): # Picks out of a subset of its neighbors and adds the best node 
        
        # Assigning original graph to temporary variable
        tempG = self.G.copy()
        
        neig_list= tempG[seed_node]
        if not neig_list:
            return []
        imp_neig = max(neig_list) # Largest weight neighbor - gives the most confident graphs 
        wt = neig_list[imp_neig]
        wt_edge = wt['weight']
        score_curr = 0
        g1=nx.Graph(comp_score=score_curr)
        g1.add_edge(seed_node,imp_neig,weight=wt_edge)
        tempG.remove_edge(seed_node,imp_neig)
        
        max_nodes = self.max_size
        
        thres_neig = self.inputs["thres_neig"] # Threshold on number of neighbors to consider 
        while True:

            #print("Adding next node")
                
            imp_neigs = dict()
            for node in g1.nodes():
                # get its max neighbor and weight and store in dict 
                neig_list= tempG[node]
                # Don't check all neighbors - just a subset if number of neighbors is large
                if len(neig_list) > thres_neig: # Make 500 
                    neig_list = dict(random.sample(list(neig_list.items()),thres_neig))
                if not neig_list: # Checking if empty
                    break
                imp_neig = self.find_imp_neig(neig_list,g1,tempG)
    
                wt = neig_list[imp_neig]
                wt_edge = wt['weight']
                
                imp_neigs[imp_neig] = {'weight': wt_edge}
                    
            if not imp_neigs:
                #print("No more neighbors to add")
                break
            
            node_to_add = self.find_imp_neig(imp_neigs,g1,tempG)
            #ADD ALL EDGES OF NEW NODE TO ORIG GRAPH
            its_neig_list= tempG[node_to_add]
                
            orig_nodes = g1.nodes()
            for node in list(orig_nodes):
                if node in its_neig_list:
                    wt = its_neig_list[node]
                    wt_edge = wt['weight']
                    g1.add_edge(node_to_add,node,weight=wt_edge)
                    tempG.remove_edge(node_to_add,node)
                
            if len(g1) > max_nodes:
                print("Max size exceeded")
                break
            score_prev = score_curr 
            (score_curr,comp_bool) = self.get_score(g1) 

            if comp_bool == 0:
            #print("Complex found")
            
                # Remove the node last added                
                g1.remove_node(node_to_add)
                score_curr = score_prev
                break
        g1.graph['comp_score'] = score_curr            
        
        #print(g1.nodes())
        #print(g1.edges())
        
        return g1
    
    def met(self,g1):
        # Assigning original graph to temporary variable
        score_prev = g1.graph['comp_score']
        tempG = self.G.copy()
        for edge in g1.edges():
            (node1,node2) = edge
            tempG.remove_edge(node1,node2)
        
        max_nodes = self.max_size - len(g1)
        
        num_iter = 1
        last_iter_imp = 0        
        thres_neig = self.inputs["thres_neig"]
        prob_metropolis = self.inputs["prob_metropolis"]
        while num_iter < max_nodes: # Limiting number of iteration rounds 

            #print("Adding next node")
                
            imp_neigs = dict()
            for node in g1.nodes():
                # get its max neighbor and weight and store in dict 
                neig_list= tempG[node]
                # Don't check all neighbors - just a subset if number of neighbors is large
                if len(neig_list) > thres_neig: # Make 500 
                    neig_list = dict(random.sample(list(neig_list.items()),thres_neig))
                if not neig_list: # Checking if empty
                    break
                imp_neig = self.find_imp_neig(neig_list,g1,tempG)
    
                wt = neig_list[imp_neig]
                wt_edge = wt['weight']
                
                imp_neigs[imp_neig] = {'weight': wt_edge}
                    
            if not imp_neigs:
                #print("No more neighbors to add")
                break
            
            node_to_add = self.find_imp_neig(imp_neigs,g1,tempG)
            #ADD ALL EDGES OF NEW NODE TO ORIG GRAPH
            its_neig_list= tempG[node_to_add]
                
            orig_nodes = g1.nodes()
            for node in list(orig_nodes):
                if node in its_neig_list:
                    wt = its_neig_list[node]
                    wt_edge = wt['weight']
                    g1.add_edge(node_to_add,node,weight=wt_edge)
                    tempG.remove_edge(node_to_add,node)
                
            (score_curr,comp_bool) = self.get_score(g1)    

            if comp_bool == 0:
            #print("Complex found")
            
                # Remove the node last added                
                g1.remove_node(node_to_add)
                break

            cur_trial = np.random.uniform(low=0.0,high=1.0)
            if score_curr < score_prev:
                if  cur_trial > prob_metropolis:  
                # Remove the node last added                
                    g1.remove_node(node_to_add)
                # since edges from this node to complex have been removed from tempG it will not be revisited 
                else:
                    print("Accepting with low probability")
            elif score_curr > score_prev:
                last_iter_imp = num_iter

            if (num_iter - last_iter_imp)> 10: # Has been a long time since a score improvement
                print("Long time since score imporovement") 
                break
            
            score_prev = score_curr

            num_iter += 1
        
        #print(g1.nodes())
        #print(g1.edges())
        g1.graph['comp_score'] = score_prev
        return g1
        
    def search_metropolis_clique_start(self,seed_clique): # Picks out of a subset of its neighbors and adds the best node 
        #print(seed_clique) 
        g1=nx.Graph(self.G.subgraph(seed_clique),comp_score=0)
        
        # Finding score 
        (score_prev,comp_bool) = self.get_score(g1) 
                    
        # Removing starting points which are not complexes    
        if comp_bool == 0:
            return []
            
        g1.graph['comp_score'] = score_prev
          
        g1 = self.met(g1)
        
        return g1
        
    def search_metropolis(self,seed_node): # Picks out of a subset of its neighbors and adds the best node 
        
        neig_list= self.G[seed_node]
        if not neig_list:
            return []
        imp_neig = max(neig_list) # Largest weight neighbor - gives the most confident graphs 
        wt = neig_list[imp_neig]
        wt_edge = wt['weight']
        
        g1=nx.Graph(comp_score=0) # Starting complex of 2 nodes is assigned score 0 arbitly
        g1.add_edge(seed_node,imp_neig,weight=wt_edge)
        
        g1 = self.met(g1)
        
        return g1
    
    def search_isa(self,seed_node): # Picks out of a subset of its neighbors and adds the best node
        
        # Assigning original graph to temporary variable
        tempG = self.G.copy()
        
        neig_list= tempG[seed_node]
        if not neig_list:
            return []
        imp_neig = max(neig_list) # Largest weight neighbor - gives the most confident graphs 
        wt = neig_list[imp_neig]
        wt_edge = wt['weight']
        
        score_prev = 0
        g1=nx.Graph(comp_score=score_prev) # Starting complex of 2 nodes is assigned score 0 arbitly
        g1.add_edge(seed_node,imp_neig,weight=wt_edge)
        tempG.remove_edge(seed_node,imp_neig)
        
        max_nodes = self.max_size - len(g1)
        
        num_iter = 1
        last_iter_imp = 0        
        thres_neig = self.inputs["thres_neig"]
        T = self.inputs["T0"] # T0 value 
        alpha = self.inputs["alpha"]
        while num_iter < max_nodes: # Limiting number of iteration rounds 

            #print("Adding next node")
                
            imp_neigs = dict()
            for node in g1.nodes():
                # get its max neighbor and weight and store in dict 
                neig_list= tempG[node]
                # Don't check all neighbors - just a subset if number of neighbors is large
                if len(neig_list) > thres_neig: # Make 500 
                    neig_list = dict(random.sample(list(neig_list.items()),thres_neig))
                if not neig_list: # Checking if empty
                    break
                imp_neig = self.find_imp_neig(neig_list,g1,tempG)
    
                wt = neig_list[imp_neig]
                wt_edge = wt['weight']
                
                imp_neigs[imp_neig] = {'weight': wt_edge}
                    
            if not imp_neigs:
                #print("No more neighbors to add")
                break
            
            node_to_add = self.find_imp_neig(imp_neigs,g1,tempG)
            #ADD ALL EDGES OF NEW NODE TO ORIG GRAPH
            its_neig_list= tempG[node_to_add]
                
            orig_nodes = g1.nodes()
            for node in list(orig_nodes):
                if node in its_neig_list:
                    wt = its_neig_list[node]
                    wt_edge = wt['weight']
                    g1.add_edge(node_to_add,node,weight=wt_edge)
                    tempG.remove_edge(node_to_add,node)     
            
            (score_curr,comp_bool) = self.get_score(g1)         

            if comp_bool == 0:
            #print("Complex found")
            
                # Remove the node last added                
                g1.remove_node(node_to_add)
                break

            cur_trial = np.random.uniform(low=0.0,high=1.0)
            if score_curr < score_prev:

                prob_isa = np.exp((score_curr - score_prev)/T)
                if  cur_trial > prob_isa:  
                # Remove the node last added                
                    g1.remove_node(node_to_add)
                # since edges from this node to complex have been removed from tempG it will not be revisited 
                else:
                    print("Accepting with low probability")
            elif score_curr > score_prev:
                last_iter_imp = num_iter

            if (num_iter - last_iter_imp)> 10: # Has been a long time since a score improvement
                print("Long time since score imporovement") 
                break
            
            score_prev = score_curr
            num_iter += 1
            T = T/alpha
        
        #print(g1.nodes())
        #print(g1.edges())
        g1.graph['comp_score'] = score_prev
        return g1
    def filter_overlapped(self,list_comp):
        print("Filtering complexes...")        

        # Sort by size 
        sizes = [nx.number_of_nodes(comp) for comp in list_comp]
        
        sorted_ind = np.argsort(sizes) # ascending order.
        
        list_comp = [list_comp[i] for i in sorted_ind]

        fin_list = list(list_comp)
        list_comp2 = list(list_comp)
        
        #print(len(list_comp))
        # Ensure ascending order 
        for comp in list_comp:
            OS_comp = []            
            list_comp2.remove(comp)
            
            if(len(list_comp2)):
                for comp2 in list_comp2:
                    
                    Overlap = self.jaccard_coeff(comp.nodes(),comp2.nodes())
                    
                    OS_comp.append(Overlap)
                
                OS_max = max(OS_comp)
                #print(OS_max)
                
                if( OS_max > self.over_t):
                    fin_list.remove(comp)
        
        #print(len(list_comp))
        print("Finished filtering complexes.")        

        return fin_list
    
    def merge_filter_overlapped_score(self,list_comp):
        print("Filtering complexes...")        

        fin_list = list(list_comp)
         
        n = len(fin_list)
        while True:
            n_changes = 0
            ind = 0
            while ind < n:
                comp = fin_list[ind]
                temp_list = list(fin_list)
                OS_comp = []            
                temp_list.remove(comp)
                
                if(len(temp_list)):
                    for comp2 in temp_list:
                        
                        #print("1:",comp.nodes())
                        #print("2:",comp2.nodes())
                        Overlap = self.jaccard_coeff(comp.nodes(),comp2.nodes())
                        #print("Overlap:",Overlap)
                        OS_comp.append(Overlap)
                    
                    OS_max = max(OS_comp)
                    OS_max_ind = np.argmax(OS_comp)
                    max_over_comp = fin_list[OS_max_ind]
                    
                    #print("OSmax=",OS_max)
    
                    if( OS_max >= self.over_t):
                        n_changes += 1                    
                        # Merge and find score. If score is higher than individual complexes 
                        # Keep as new complex                                        
                        merge_comp_nodes = list(set(comp.nodes()).union(set(max_over_comp.nodes())))
                        merge_comp = nx.Graph(self.G.subgraph(merge_comp_nodes),comp_score=0)
                        
                        (score_merge,comp_bool) = self.get_score(merge_comp)
                        merge_comp.graph['comp_score'] = score_merge
                        sc1 = comp.graph['comp_score']
                        sc2 = max_over_comp.graph['comp_score']                    
                        if(score_merge > sc1 and score_merge > sc2):
                            fin_list.append(merge_comp)
                            fin_list.remove(comp)
                            del fin_list[OS_max_ind]
                            if OS_max_ind <= ind:
                                ind -= 1                                
    
                        # Otherwise: remove lower scoring complex
                        elif sc1 <= sc2:
                            fin_list.remove(comp)
                        else:
                            del fin_list[OS_max_ind]
                            if OS_max_ind > ind:
                                ind += 1
                    else:
                        ind += 1
                        
                else:
                    print("temp_list is empty")
                    ind += 1
                n = len(fin_list)

            print(n_changes)
            if n_changes == 0:
                break
        
        #print(len(list_comp))
        print("Finished filtering complexes.")          

        return fin_list

    '''    
    def processInput(self,i):
        return i * i                
    '''        
    #@profile
    def sample(self,num_comp,run_mode,seed_mode):
        print("Sampling complexes...")   
        
        if(self.max_size >= 30):
            self.max_size = 30
        num_cores = multiprocessing.cpu_count()

        #import time 
        #start_time_sample = time.time()  
        
        pred_comp_list = []
        pred_comp_node_list = []
        #num_comp = 10

        if(seed_mode == "all_nodes"):
            seed_nodes = list(self.G.nodes())
        elif(seed_mode == "n_nodes"):
            seed_nodes = np.random.permutation(list(self.G.nodes()))[:num_comp]
        elif(seed_mode == "all_nodes_known_comp"):
            seed_nodes = self.prot_list
        elif(seed_mode == "cliques"):
            clique_list = list(nx.find_cliques(self.G))
            to_rem = []
            # Removing 2 node and big complexes
            for comp in clique_list:
                if (len(comp) <=2 or len(comp) >= self.max_size):
                    to_rem.append(comp)
                
            for comp in to_rem:
                clique_list.remove(comp)

                
            seed_nodes = clique_list # Remove duplicates later.
        
        num_comp = len(seed_nodes)
 
        with open(self.out_comp_nm+'_metrics.out',"a") as fid: 
            print("No. of cores = ", num_cores,file = fid)               
            print("No. of seeds for complex search = ",num_comp,file = fid)        

        search_method = self.inputs["search_method"]
        
        if(seed_mode == "cliques"):
            pred_comp_list = [self.search_metropolis_clique_start(clique) for clique in tqdm(clique_list)] 
        elif(search_method == "isa"):
            if(run_mode == "parallel"):
                pred_comp_list = Parallel(n_jobs=num_cores,prefer="threads")(delayed(unwrap_self_isa)(node) for node in tqdm(zip([self]*num_comp,seed_nodes)))
            else:
                pred_comp_list = [self.search_isa(node) for node in tqdm(seed_nodes)]              
                
        elif(search_method == "metropolis"):
            if(run_mode == "parallel"):            
                pred_comp_list = Parallel(n_jobs=num_cores,prefer="threads")(delayed(unwrap_self_metropolis)(node) for node in tqdm(zip([self]*num_comp,seed_nodes)))
            else:
                pred_comp_list = [self.search_metropolis(node) for node in tqdm(seed_nodes)]              
        elif(search_method == "search_top_neigs"):
            if(run_mode == "parallel"):            
                pred_comp_list = Parallel(n_jobs=num_cores,prefer="threads")(delayed(unwrap_self_top_neigs)(node) for node in tqdm(zip([self]*num_comp,seed_nodes)))
            else:
                pred_comp_list = [self.search_top_neigs(node) for node in tqdm(seed_nodes)]              
        else:
            if(run_mode == "parallel"):            
                pred_comp_list = Parallel(n_jobs=num_cores,prefer="threads")(delayed(unwrap_self_max_neig)(node) for node in tqdm(zip([self]*num_comp,seed_nodes)))
            else:
                pred_comp_list = [self.search_max_neig(node) for node in tqdm(seed_nodes)]              
      
                    
        '''
        # Does'nt work for some reason I think 
        pool = multiprocessing.Pool(4)
        seeds = [all_nodes[num] for num in range(num_comp)]
        pred_comp_list = pool.map(self.search,seeds)
        '''
        
        '''
        # Example parallel process
        inputs = range(10)
        
        res1 = [self.processInput(i) for i in inputs]
        print(res1)

        results = Parallel(n_jobs=num_cores)(delayed(self.processInput)(i) for i in inputs)
        print(results)
        '''
       
        # Removing complexes with only two nodes 
        pred_comp_list = [comp for comp in pred_comp_list if len(comp)>2]
        pred_comp_node_list = [list(comp.nodes()) for comp in pred_comp_list]
        # Finding unique ccomplexes     
        #pred_comp_node_list = np.array(pred_comp_node_list)
        #print("all",pred_comp_node_list)
        
        #(fin_list,indices) = np.unique(pred_comp_node_list,axis = 0, return_index=True)
        
        fin_list=[]
        indices=[]
        for j,comp in enumerate(pred_comp_node_list):
            if comp not in fin_list:
                fin_list.append(comp)        
                indices.append(j)
        #print("fin",fin_list)
        self.fin_list_graphs = [pred_comp_list[i] for i in indices]
        
        #self.fin_list_graphs = pred_comp_list

        print("Finished sampling complexes.")  
        #with open(self.out_comp_nm+'_metrics.out',"a") as fid:                
        #    print("--- Sampling Runtime  = ",(time.time() - start_time_sample)," seconds ---",file = fid)

 
        # Filtering complexes with high overlap with bigger complexes 
        self.fin_list_graphs = self.merge_filter_overlapped_score(self.fin_list_graphs)
        
        print("Writing predicted complexes.")        
        
        with open(self.out_comp_nm+'_pred.out',"w") as fn:
            with open(self.out_comp_nm+'_pred_edges.out',"wb") as f_edges:
                for index in range(len(self.fin_list_graphs)):
                    for node in self.fin_list_graphs[index].nodes():
                        fn.write("%s " % node)     
                    
                    fn.write("%.3f" % self.fin_list_graphs[index].graph['comp_score'])                        
                    nx.write_weighted_edgelist(self.fin_list_graphs[index],f_edges)
                    fn.write("\n")
                    f_edges.write("\n".encode())

        print("Finished writing predicted complexes.")    
        
        print("Writing known complexes.")        
        known_complexes = self.complex_graphs + self.test_complex_graphs
        with open(self.out_comp_nm+'_known.out',"w") as fn:
            with open(self.out_comp_nm+'_known_edges.out',"wb") as f_edges:
                for index in range(len(self.known_complex_nodes_list)):
                    for node in self.known_complex_nodes_list[index]:
                        fn.write("%s " % node)
                    nx.write_weighted_edgelist(known_complexes[index],f_edges)
                    fn.write("\n")
                    f_edges.write("\n".encode())

        print("Finished writing known complexes.")       

    #@profile
    def eval_complex(self,rf=0,rf_nm = 0):
        print("Evaluating complexes...")     
        
        if(rf == 1):
            if(rf_nm == 0):
                rf_nm = self.out_comp_nm+'_pred.out'
            with open(rf_nm) as fn:
                complex_list = [line.rstrip('\n').split() for line in fn ]# Space separated text only
                self.fin_list_graphs = [nx.Graph(self.G.subgraph(complex)) for complex in complex_list]                    
                # Just list of list of nodes 
                #self.fin_list_graphs = [nx.Graph(complex) for complex in complex_list]                    
        
        sizes_orig = [len(comp) for comp in self.fin_list_graphs]  
        
        p = self.inputs["eval_p"]
        
        N_pred_comp = len(self.fin_list_graphs)
        with open(self.out_comp_nm+'_metrics.out',"a") as fid:                
            print("No. of predicted complexes = ",N_pred_comp,file = fid)
        # Remove all proteins in predicted complexes that are not present in known complex protein list
        
        comp_remove=[]
        for comp in self.fin_list_graphs:
            to_remove=[]
            for node in comp.nodes():
                if(node not in self.prot_list):
                    to_remove.append(node)
            comp.remove_nodes_from(to_remove)
            
            if(len(comp) <= 1): # Removing complexes with only one node or none 
                comp_remove.append(comp)
        
        self.fin_list_graphs = [graph for graph in self.fin_list_graphs if graph not in comp_remove]

        sizes_known = [len(comp) for comp in self.known_complex_nodes_list]
        # Size distributions 
        sizes_new = [len(comp) for comp in self.fin_list_graphs]        
        fig = plt.figure()
        sns.distplot(sizes_known,hist=False,label="known")
        sns.distplot(sizes_orig,hist=False,label="predicted")     
        sns.distplot(sizes_new,hist=False,label="predicted_known_prots")  
        
        plt.savefig(self.out_comp_nm+'_size_dists_known_pred.png')        
        #plt.show()
        plt.close(fig)            
        

        N_test_comp = len(self.known_complex_nodes_list)
        N_pred_comp = len(self.fin_list_graphs)
        with open(self.out_comp_nm+'_metrics.out',"a") as fid:        
            print("No. of predicted complexes after removing non-gold std proteins = ",N_pred_comp,file = fid)        

        N_matches_test = 0
        
        Metric = np.zeros((N_test_comp,N_pred_comp))
        
        for i,test_complex in enumerate(self.known_complex_nodes_list):
            N_match_pred=0
            for j,pred_complex in enumerate(self.fin_list_graphs):
                T = set(test_complex)
                P = set(pred_complex.nodes())
                C = len(T.intersection(P))
                A = len(P.difference(T))
                B = len(T.difference(P))
                
                if(float(C)/(A+C) > p and float(C)/(B+C) > p):
                    Metric[i,j] = 1; 
                    N_match_pred=N_match_pred+1
                    
            if(N_match_pred > 0):
                N_matches_test=N_matches_test+1
                
       
        Recall = float(N_matches_test)/N_test_comp
        
        N_matches_pred = np.count_nonzero(np.sum(Metric,axis=0))
        Precision = float(N_matches_pred)/N_pred_comp
        
        with open(self.out_comp_nm+'_metrics.out',"a") as fid:                
            print("No. of known complexes = ",N_test_comp,file = fid)         
            print("Prediction Precision = ",Precision,file = fid)
            print("Prediction Recall = ",Recall,file = fid)        
            try:
                F1_score = 2*Precision*Recall/(Precision + Recall)
                print("Prediction F1 score = ",F1_score,file = fid)            
            except:
                print("Error in calculating F1 score - likely divide by 0")
                
                        
        print("Finished Evaluating complexes.")        
        
