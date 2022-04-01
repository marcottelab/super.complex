# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 11:18:26 2021

@author: Meghana
"""

import pandas as pd

df = pd.read_excel('../convert_ids/SARS_COV2_Map.xlsx',header=1)

edges = list(zip(df['Bait'],df['PreyGene'],df['MIST']))

nodes_sars = set(list(df['Bait']))

nodes_human_prots = set(list(df['PreyGene']))

results_folder = "../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_qi_o0.375_reproduce"        
prot2comp_cov_file = results_folder + '/res_pred_prot2comp_covid.out'
# generate protein to complexes edges
with open(prot2comp_cov_file) as f:
    lines = [line.rstrip().split() for line in f.readlines()]
    
edges_comp = []
nodes_comp = []
prots_in_comps = set()
for line in lines:
    prot = line[0]
    prots_in_comps.add(prot)
    for comp in line[2:]:
        edges_comp.append((prot,comp,0.99))
        nodes_comp.append(comp)

nodes_comp = set(nodes_comp)

nodes_human_prots = nodes_human_prots.intersection(prots_in_comps)

edges = [edge for edge in edges if edge[1] in nodes_human_prots]    

import pickle
with open('../convert_ids/name2annot_full_humap_from_gene_names.pkl','rb') as f:
    name2annot = pickle.load(f)
    
    
pref = '''
<!doctype html>

<html>

<head>
'''

post = '''
style: [
    {
        selector: 'node',
        style: {
            shape: function( ele ){ 
                                if ( ele.data('type') == 'complex' ){return 'rectangle'}
                                if ( ele.data('type') == 'sars' ){return 'triangle'}
                                
  return 'circle';
  },
            'background-color': function( ele ){ 
                                if ( ele.data('type') == 'complex' ){return 'blue'}
                                if ( ele.data('type') == 'sars' ){return 'red'}             
  return 'yellow';
  },
            label: 'data(id)'
        }
    },
							{
							selector: 'edge',
							style: {
								'width': 'mapData(weight, 0.5, 1, 0.1, 2.5)',
								'line-color': function( ele ){ 
                                if ( ele.data('type') == 'prot2complex' ){return 'blue'}
  return 'red';
  },
								'opacity': 0.5
							}
						}
	]
    
});

cy.layout({
    name: 'concentric',
  concentric: function( ele ){ // returns numeric value for each node, placing higher nodes in levels towards the centre
  return ele.data('score')+ ele.degree()/1000;
  },
  equidistant: true, 
  levelWidth: function( nodes ){ // the variation of concentric values in each level
  return 1;
  },
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
    
with open("../convert_ids/name2uniprotURL_full_humap_from_gene_names.pkl",'rb') as f:
    prot2url = pickle.load(f)
    

edge_nm = ""
num = 0
covid_complex_ids = []

pre = pref + '<title> SARS-COV2-human map ' +  '</title>\n'

pre = pre + '''
    <script src="./cytoscape.js"></script>
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

pre = pre + '\n'

pre = pre + '''
 </section>
    <div id="cy"></div>
    <script>
      var cy = cytoscape({
  container: document.getElementById('cy'),
  '''

res = pre + "\n" + "elements:{nodes: ["
num +=1
n_covid_interactors = 0
complex_covid_interactors = []
#edges = edges + edges_comp

for i,node in enumerate(nodes_sars):
    res += "{data: { id: '" + node + "', score: " + str(2) + ", type: " + " 'sars' "
    
    res += ', href: "' + "./SARS_COV2prots/" + node + '.html' '"'
    
    
    if node in prot2url:
        res += ', href: "' + prot2url[node] + '"'
        
    res += " }}"
    res += ",\n"
    
for i,node in enumerate(nodes_human_prots):
    res += "{data: { id: '" + node + "', score: " + str(1) + ", type: " + " 'human' "
    
    if node in prot2url:
        res += ', href: "' + prot2url[node] + '"'
        
    res += " }}"
    res += ",\n"    


nNodes = len(nodes_comp)
for i,node in enumerate(nodes_comp):
    res += "{data: { id: '" + node + "', score: " + str(0) + ", type: " + "'complex'"
    
    res += ', href: "' + "./Predicted/" + node + '.html' '"'
        
    res += " }}"
    if i != nNodes-1:
        res += ",\n"
    else:
        res += "], \n edges: ["

for j,elts in enumerate(edges_comp):
    res += "{ data: { source: '" + elts[0] +"', target: '" + elts[1] + "', directed: 'false'," + "type: 'prot2complex'," + " weight: " + str(elts[2]) + "} }"
    res += ",\n"
        
edges = edges[:-1]
nEdges = len(edges)
for j,elts in enumerate(edges):
    res += "{ data: { source: '" + elts[0] +"', target: '" + elts[1] + "', directed: 'false', " + "type: 'sars2prot'," + "weight: " + str(elts[2]) + "} }"
    if j != nEdges-1:
        res += ",\n"
    else:
        res += "]}" 
res += ",\n" + post

res = res + '''
	<section>
	<br>
	<br>
	<br> 
''' 

    
res = res + '''
</section>
</body>
</html>
'''

with open("./Complexes/SARS_COV2_Map_only_mapped_prots"+ ".html","w") as f:
    f.write(res)