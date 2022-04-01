import pickle
from os import mkdir as os_mkdir, path as os_path,remove as os_remove
from distutils.dir_util import copy_tree
from shutil import copyfile,rmtree
import pandas as pd

df = pd.read_excel('../convert_ids/SARS_COV2_Map.xlsx',header=1)

edges = list(zip(df['Bait'],df['PreyGene'],df['MIST']))
nodes_sars = set(list(df['Bait']))


# SARS to prot
sars2prot_edges = dict()
sars2prot_nodes = dict()

for edge in edges:
    prot = edge[1]
    sars = edge[0]
    wt = edge[2]
    if sars not in sars2prot_edges:
        sars2prot_nodes[sars] = [prot]
        sars2prot_edges[sars] = dict()
    else:
        sars2prot_nodes[sars].append(prot)
    sars2prot_edges[sars][prot] = wt         # if duplicates exist, last value will be used
        
        
# prot to comp
results_folder = "../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_qi_o0.375"        
results_file = results_folder + "/res_pred_prot2comp_covid.out"

with open(results_file) as f:
    prot2comp = dict([(line2[0].rstrip(' '), line2[1].split()) for line2 in (line.rstrip().split(':') for line in f.readlines()) if line2])
        
# sars to comp 
sars_graph_edges = []
sars2comp = dict()
sars2comp_edges = dict()
for sars in sars2prot_nodes:
    prots = sars2prot_nodes[sars]
    
    if sars not in sars2comp_edges:
        sars2comp_edges[sars] = dict()   
        
    for prot in prots:
        if prot in prot2comp:
            if sars not in sars2comp:
                sars2comp[sars] = prot2comp[prot][:]
            else:
                sars2comp[sars] = sars2comp[sars]  + prot2comp[prot][:]
    # only mapped proteins
 
            for comp in prot2comp[prot]:
                if comp not in sars2comp_edges[sars]:
                    sars2comp_edges[sars][comp] = sars2prot_edges[sars][prot]
                else:
                    sars2comp_edges[sars][comp] += sars2prot_edges[sars][prot]

# Normalize weights between 0 and 1
max_wt = max([sars2comp_edges[sarsss][comppp] for sarsss in sars2comp_edges for comppp in sars2comp_edges[sarsss]])   
for sars in sars2comp_edges:
    for comp in sars2comp_edges[sars]:
        sars2comp_edges[sars][comp] = sars2comp_edges[sars][comp]/max_wt
        sars_graph_edges.append((sars,comp,sars2comp_edges[sars][comp]))
            

nodes_comp = set() # Only complexes !!

for comp_list in sars2comp.values():
    for comp in comp_list:
        nodes_comp.add(comp)
        

    
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
                                if ( ele.data('type') == 'complex' ){return 'circle'}
                                if ( ele.data('type') == 'sars' ){return 'triangle'}
                                
  return 'rectangle';
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
								'width': 'mapData(weight, 0.1, 1, 0.1, 2.5)',
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
    name: 'cose',
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

        
nEdges = len(sars_graph_edges)
for j,elts in enumerate(sars_graph_edges):
    res += "{ data: { source: '" + elts[0] +"', target: '" + elts[1] + "', directed: 'false'," + "type: 'prot2complex'," + " weight: " + str(elts[2]) + "} }"
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

with open("./Complexes/SARS_COV2_Map_only_mapped_complexes"+ ".html","w") as f:
    f.write(res)
        
    
                
            

                
