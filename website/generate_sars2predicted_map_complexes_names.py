import pickle
from os import mkdir as os_mkdir, path as os_path,remove as os_remove
from distutils.dir_util import copy_tree
from shutil import copyfile,rmtree
import pandas as pd
import argparse

def get_sars2url():
    df = pd.read_excel('../convert_ids/SARS2URL.xlsx',header=0)
    return dict(zip(df['SARS'],df['URL']))

def jaccard_coeff(set1, set2):
    ls1 = len(set1)
    ls2 = len(set2)
    if ls1 == 0 and ls2 == 0:
        return 1
    inter = len(set1.intersection(set2))
    return float(inter) / (ls1 + ls2 - inter)

def get_comp2comp_edges(nodes_comp, results_folder): # get complex to complex edges as jc b/w nodes
    with open(results_folder+'/res_pred_names.out') as f:
        comp2prots = dict([('Complex'+str(i+1),set(line.rstrip().split()[:-1])) for i,line in enumerate(f.readlines()) if 'Complex'+str(i+1) in nodes_comp])
    
    comp_edges=[]
    n_comps = len(nodes_comp)
    nodes_comp = list(nodes_comp)
    for i in range(n_comps):
        for j in range(i+1,n_comps):
            jc = jaccard_coeff(comp2prots[nodes_comp[i]],comp2prots[nodes_comp[j]])
            if jc != 0:
                comp_edges.append((nodes_comp[i],nodes_comp[j],jc))
                
    return comp_edges

def main():
    parser = argparse.ArgumentParser("Input parameters")
    parser.add_argument("--sars_cov2_map", default="../convert_ids/SARS_COV2_Map.xlsx", help="Input parameters file name")
    parser.add_argument("--results_folder", default="../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_qi_o0.375", help="Input parameters file name")
    parser.add_argument("--id2name_file", default="../convert_ids/humap_gene_id_name_map_updated.pkl", help="Input parameters file name")
    parser.add_argument("--prot2url_file", default="../convert_ids/name2uniprotURL_full_humap_from_gene_names.pkl", help="Input parameters file name")
    
    
    args = parser.parse_args()
        
    #sars2url = get_sars2url()
    
    df = pd.read_excel(args.sars_cov2_map,header=1)
    
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
    results_folder = args.results_folder
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
                
    
    with open(results_folder+'/res_pred_complex_num2name.pkl','rb') as f:
        comp_num2name = pickle.load(f)
        
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
            				'font-size': 20,
            'font-weight': 'bold',
                'shape': function( ele ){ 
                                    if ( ele.data('type') == 'complex' ){return 'circle'}
                                    if ( ele.data('type') == 'sars' ){return 'triangle'}
                                    
      return 'rectangle';
      },
                'background-color': function( ele ){ 
                                    if ( ele.data('type') == 'complex' ){return 'blue'}
                                    if ( ele.data('type') == 'sars' ){return 'red'}             
      return 'yellow';
      },
                  'text-wrap': 'wrap',
        'text-max-width': '250',
                'label': 'data(id)'
            }
        },
    							{
    							selector: 'edge',
    							style: {
    								'width': 'mapData(weight, 0, 1, 0, 5)',
    								'line-color': function( ele ){ 
                                    if ( ele.data('type') == 'prot2complex' ){return 'red'}
      return 'blue';
      },
    								'opacity': 0.5
    							}
    						}
    	]
        
    });
    
    cy.layout({
        'name': 'cose',
    	'animate':false,
    	'nodeDimensionsIncludeLabels':true,
    	'padding': 1,
    	'componentSpacing': 20,
    	 'idealEdgeLength': function( edge ){ return 10; },
    	 // Divisor to compute edge forces
      'edgeElasticity': function( edge ){ return 20; },
       'fit': true,
      'gravity' : 100   
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
        
    
        
    
    # edge_nm = ""
    num = 0
    # covid_complex_ids = []
    
    pre = pref + '<title> SARS-COV2-human map ' +  '</title>\n'
    
    pre = pre + '''
    		<script src="https://unpkg.com/cytoscape/dist/cytoscape.min.js"></script>
    
    		<script src="https://unpkg.com/layout-base/layout-base.js"></script>
    <script src="https://unpkg.com/cose-base/cose-base.js"></script>
    <script src="https://unpkg.com/cytoscape-fcose/cytoscape-fcose.js"></script>
    </head>
    
    <style>
        #cy {
            width: 1920px;
            height: 1080px;
            position: absolute;
                top: 100px;
                right: -300px;	
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
    # n_covid_interactors = 0
    # complex_covid_interactors = []
    #edges = edges + edges_comp
    
    comp_edges=get_comp2comp_edges(nodes_comp, results_folder)
    
    for i,node in enumerate(nodes_sars):
        res += "{data: { id: '" + node + "', score: " + str(2) + ", type: " + " 'sars' "
        
        res += ', href: "' + "./SARS_COV2prots/" + node + '.html' '"'
            
        res += " }}"
        res += ",\n"  
    
    nNodes = len(nodes_comp)
    for i,node in enumerate(nodes_comp):
        res += "{data: { id:"+'"' + comp_num2name[node] + '"' + ", score: " + str(0) + ", type: " + "'complex'"
        
        res += ', href: "' + "./Predicted/" + node + '.html' '"'
            
        res += " }}"
        if i != nNodes-1:
            res += ",\n"
        else:
            res += "], \n edges: ["
    
    for j,elts in enumerate(comp_edges):
        res += '{ data: { source: "' + comp_num2name[elts[0]] +'", target: "' + comp_num2name[elts[1]] +'"'+ ", directed:" + " 'false'," + "type: 'complex2complex'," + " weight: " + str(elts[2]) + "} }"
        
        res += ",\n"
            
    nEdges = len(sars_graph_edges)
    for j,elts in enumerate(sars_graph_edges):
        res += '{ data: { source: "' + elts[0] +'", target: "' + comp_num2name[elts[1]] +'"'+ ", directed:" + " 'false'," + "type: 'prot2complex'," + " weight: " + str(elts[2]) + "} }"
        
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
    
    with open("./Complexes/SARS_COV2_Map_only_mapped_complexes_names"+ ".html","w") as f:
        f.write(res)
        
if __name__ == "__main__":
    main()
            
        
                    
                
    
                    
