# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 11:18:26 2021

@author: Meghana
"""

import pandas as pd

df = pd.read_excel('../convert_ids/SARS_COV2_Map.xlsx',header=1)

edges = list(zip(df['Bait'],df['PreyGene'],df['MIST']))

nodes = set(list(df['Bait'])+list(df['PreyGene']))

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
            shape: 'rectangle',
            'background-color': 'mapData(score, 0, 5, yellow, green)',
            label: 'data(id)'
        }
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
    name: 'concentric'
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

pre = pref + '<title> Complex ' +  '</title>\n'

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
nNodes = len(nodes)
n_covid_interactors = 0
complex_covid_interactors = []
for i,node in enumerate(nodes):
     
    res += "{data: { id: '" + node + "', score: " + str(0) 
    
    if node in prot2url:
        res += ', href: "' + prot2url[node] + '"'
        
    res += " }}"
    if i != nNodes-1:
        res += ",\n"
    else:
        res += "], \n edges: ["
        
edges = edges[:-1]
nEdges = len(edges)
for j,elts in enumerate(edges):
    res += "{ data: { source: '" + elts[0] +"', target: '" + elts[1] + "', directed: 'false', weight: " + str(elts[2]) + "} }"
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

with open("./Complexes/SARS_COV2_Map"+ str(num) + ".html","w") as f:
    f.write(res)