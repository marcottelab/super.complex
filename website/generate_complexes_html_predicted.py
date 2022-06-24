# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 11:04:51 2021

@author: Meghana
"""

import pickle
import argparse
from os import mkdir as os_mkdir, path as os_path,remove as os_remove
from distutils.dir_util import copy_tree
from shutil import copyfile,rmtree

parser = argparse.ArgumentParser("Input parameters")
parser.add_argument("--sars_cov2_map", default="../convert_ids/SARS_COV2_Map.xlsx", help="Input parameters file name")
parser.add_argument("--results_folder", default="../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_qi_o0.375", help="Input parameters file name")
parser.add_argument("--id2name_file", default="../convert_ids/humap_gene_id_name_map_updated.pkl", help="Input parameters file name")
parser.add_argument("--prot2url_file", default="../convert_ids/name2uniprotURL_full_humap_from_gene_names.pkl", help="Input parameters file name")
parser.add_argument("--name2annot_file", default="../convert_ids/name2annot_full_humap_from_gene_names.pkl", help="Input parameters file name")


args = parser.parse_args()

with open(args.name2annot_file,'rb') as f:
    name2annot = pickle.load(f)
    
    
pref = '''
<!doctype html>

<html>

<head>
'''

circle_post = '''
style: [
    {
        selector: 'node',
        style: {
		'font-size': 30,
        'font-weight': 'bold',
        'shape': function( ele ){ 
                                if ( ele.data('type') == 'cov_link' ){return 'star'}                                
  return 'rectangle';
  },            
            'background-color': 'mapData(score, 0, 5, yellow, green)',
            'label': 'data(id)'        }
    },
							{
							selector: 'edge',
							style: {
								'width': 'mapData(weight, 0.002, 1, 0.5, 2.5)',
								'line-color': 'blue',
								'opacity': 0.5
							}
						}
	]
    
});

cy.layout({
    name: 'circle'
}).run();
    
cy.on('tap', 'node', function(){
  try { // your browser may block popups
    window.open( this.data('href') );
  } catch(e){ // fall back on url change
    window.location.href = this.data('href');
  }
});      


    </script>
'''

fcose_post = '''
style: [
    {
        selector: 'node',
        style: {
		'font-size': 10,
		'width':'10',
		'height': '10',	
        'font-weight': 'bold',
        'shape': function( ele ){ 
                                if ( ele.data('type') == 'cov_link' ){return 'star'}                                
  return 'rectangle';
  },            
            'background-color': 'mapData(score, 0, 5, yellow, green)',
            'label': 'data(id)'        }
    },
							{
							selector: 'edge',
							style: {
								'width': 'mapData(weight, 0.002, 1, 0.3, 1)',
								'line-color': 'blue',
								'opacity': 0.5
							}
						}
	]
    
});

cy.layout({
    'name': 'fcose',
		'animate':false,
	'nodeDimensionsIncludeLabels':true,
	'nodeSeparation': 1000
}).run();
    
cy.on('tap', 'node', function(){
  try { // your browser may block popups
    window.open( this.data('href') );
  } catch(e){ // fall back on url change
    window.location.href = this.data('href');
  }
});      


    </script>
'''
results_folder = args.results_folder
#results_folder = "../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_qi_o0.375"        
results_file = results_folder + "/res_pred_names.out"
results_file_edges = results_folder + "/res_pred_edges_names.out"
#results_file = "../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_best/res_pred_names.out"

#results_file = "../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_o0.25/res_pred_names.out"

with open(results_file) as f:
    complexes_nodes = f.read()
    
with open(results_file_edges) as f:
    complexes_edgesF = f.read()
    
with open("../convert_ids/covid_interactors.txt") as f:
    covid_interactors = set([elt.rstrip() for elt in f.readlines()])
    
with open(args.prot2url_file,'rb') as f:
    prot2url = pickle.load(f)
    
complexes = complexes_nodes.split("\n")
complexes_edges = complexes_edgesF.split("\n\n")[:-1]

edge_nm = ""
num = 0
covid_complex_ids = []
comp_num = 0
for comp,comp_edges in zip(complexes,complexes_edges):
    comp_num += 1
    pre = pref + '<title> Complex ' + str(comp_num) + '</title>\n'
    
    pre = pre + '''
		<script src="https://unpkg.com/cytoscape/dist/cytoscape.min.js"></script>

		<script src="https://unpkg.com/layout-base/layout-base.js"></script>
<script src="https://unpkg.com/cose-base/cose-base.js"></script>
<script src="https://unpkg.com/cytoscape-fcose/cytoscape-fcose.js"></script>
    </head>
    
    <style>
        #cy {
            width: 90%;
            height: 90%;
            position: absolute;
            middle: 0px;
            left: 0px;
        }
    </style>
    
    <body>
    <section>
    '''
    
    pre = pre + '<a href="Complex' + str(comp_num+1) + '.html">Next</a>\n'
    
    pre = pre + '''
     </section>
        <div id="cy"></div>
        <script>
          var cy = cytoscape({
      container: document.getElementById('cy'),
      '''

    res = pre + "\n" + "elements:{nodes: ["
    num +=1
    words = comp.split(" ")
    nodes = words[:-1]
    #score = words[-1]
    nNodes = len(nodes)
    n_covid_interactors = 0
    complex_covid_interactors = []
    for i,node in enumerate(nodes):
            
        annot_score = 0
        if node in name2annot:
            annot_score = name2annot[node]         
            
        # covid check 
        if node in covid_interactors:
            n_covid_interactors += 1
            complex_covid_interactors.append(node)
            type_cov = "'cov_link'"
        else:
            type_cov = "''"
            
        res += "{data: { id: '" + node + "', score: " + str(annot_score) + ", type: " + type_cov
        
        if node in prot2url:
            res += ', href: "' + prot2url[node] + '"'
            
        res += " }}"
        if i != nNodes-1:
            res += ",\n"
        else:
            res += "], \n edges: ["
            
    edges = comp_edges.split("\n")
    complex_score = float(edges[-1].split()[-1])
    edges = edges[:-1]
    nEdges = len(edges)
    for j,edge in enumerate(edges):
        words = edge.split(" ")
        nodes = words[:-1]
        score = words[-1]
        res += "{ data: { source: '" + nodes[0] +"', target: '" + nodes[1] + "', directed: 'false', weight: " + score + "} }"
        if j != nEdges-1:
            res += ",\n"
        else:
            res += "]}" 
    
    density = float(2*nEdges)/(nNodes*(nNodes-1))
    if (density > 0.8):
        post = circle_post
    else:
        post = fcose_post    
    res += ",\n" + post
    
    res = res + '''
    	<section>
	<br>
	<br>
	<br> 
    ''' +    'Complex Score (from Super.Complex) = ' + "{:.5f}".format(complex_score)

    
    if n_covid_interactors > 0:
        covid_complex_ids.append(comp_num)
        res = res + '''
        <br> 
        <br>
        Connection to SARS-COV2: 
        <br>
        ''' + str(n_covid_interactors) + ' interacting protein(s): <br>' + ' '.join(complex_covid_interactors)
    
    res = res + '''
    </section>
    </body>
    </html>
    '''
    
    with open("./Complexes/Predicted/Complex"+ str(num) + ".html","w") as f:
        f.write(res)