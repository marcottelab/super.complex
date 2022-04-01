# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 11:04:51 2021

@author: Meghana
"""

import pickle
from argparse import ArgumentParser as argparse_ArgumentParser

parser = argparse_ArgumentParser("Input parameters")
parser.add_argument("--suffix_corum", default='', help="original")
args = parser.parse_args()
    
##suffix_corum = ''
#suffix_corum = 'original'
suffix_corum=args.suffix_corum
with open('../convert_ids/name2annot_full_humap_from_gene_names.pkl','rb') as f:
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
            shape: 'rectangle',
            'background-color': 'mapData(score, 0, 5, yellow, green)',
            label: 'data(id)'
        }
    },
							{
							selector: 'edge',
							style: {
								'width': 'mapData(weight, 0, 1, 0, 2.5)',
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
            shape: 'rectangle',
            'background-color': 'mapData(score, 0, 5, yellow, green)',
            label: 'data(id)'
        }
    },
							{
							selector: 'edge',
							style: {
								'width': 'mapData(weight, 0, 1, 0, 1)',
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

with open('../humap/corum_names' + suffix_corum + '.txt') as f:
    complexes_nodes = f.read()
    
with open('../humap/corum_edges_humap_names' + suffix_corum + '.txt') as f:
    complexes_edgesF = f.read()
    
with open("../convert_ids/covid_interactors.txt") as f:
    covid_interactors = set([elt.rstrip() for elt in f.readlines()])
    
with open("../convert_ids/name2uniprotURL_full_humap_from_gene_names.pkl",'rb') as f:
    prot2url = pickle.load(f)
    
if suffix_corum != 'original':
    with open('../humap/cleaned_corum_complex_num2name.pkl','rb') as f:
        comp_num2name_corum = pickle.load(f)           
else:
    with open('../humap/original_corum_complex_num2name.pkl','rb') as f:
        comp_num2name_corum = pickle.load(f)     
    
complexes = complexes_nodes.split("\n")
complexes_edges = complexes_edgesF.split("\n\n")[:-1]

edge_nm = ""
num = 0
covid_complex_ids = []
comp_num = 0
for comp,comp_edges in zip(complexes,complexes_edges):
    comp_num += 1
    pre = pref + '<title> ' + comp_num2name_corum['Complex' + str(comp_num)] + '</title>\n'
    
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
    nodes = comp.split(" ")
    #score = words[-1]
    nNodes = len(nodes)
    n_covid_interactors = 0
    complex_covid_interactors = []
    for i,node in enumerate(nodes):
        
        # covid check 
        if node in covid_interactors:
            n_covid_interactors += 1
            complex_covid_interactors.append(node)
            
        annot_score = 0
        if node in name2annot:
            annot_score = name2annot[node]            
        res += "{data: { id: '" + node + "', score: " + str(annot_score) 
        
        if node in prot2url:
            res += ', href: "' + prot2url[node] + '"'
            
        res += " }}"
        if i != nNodes-1:
            res += ",\n"
        else:
            res += "], \n edges: ["
            
    edges = comp_edges.split("\n")

    nEdges = len(edges)
    nEdges_tot = nEdges
    for j,edge in enumerate(edges):
        words = edge.split(" ")
        nodes = words[:-1]
        score = words[-1]
        res += "{ data: { source: '" + nodes[0] +"', target: '" + nodes[1] + "', directed: 'false', weight: " + score + "} }"
        if float(score) == 0:
            nEdges_tot -= 1
        if j != nEdges-1:
            res += ",\n"
        else:
            res += "]}"
    density = float(2*nEdges_tot)/(nNodes*(nNodes-1))
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
    ''' + comp_num2name_corum['Complex' + str(comp_num)]  + '	<br> Unweighted Density = ' + "{:.2f}".format(density)

    
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
    
    with open("./Complexes/"+suffix_corum+"CORUM/Complex"+ str(num) + ".html","w") as f:
        f.write(res)