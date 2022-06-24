# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 11:18:26 2021

@author: Meghana
"""
import pickle
import pandas as pd
import argparse


def jaccard_coeff(set1, set2):
    ls1 = len(set1)
    ls2 = len(set2)
    if ls1 == 0 and ls2 == 0:
        return 1
    inter = len(set1.intersection(set2))
    return float(inter) / (ls1 + ls2 - inter)


def get_comp2comp_edges(nodes_comp,comp2prots): # get complex to complex edges as jc b/w nodes
    
    comp_edges=[]
    n_comps = len(nodes_comp)
    nodes_comp = list(nodes_comp)
    for i in range(n_comps):
        for j in range(i+1,n_comps):
            jc = jaccard_coeff(comp2prots[nodes_comp[i]],comp2prots[nodes_comp[j]])
            if jc != 0:
                comp_edges.append((nodes_comp[i],nodes_comp[j],jc))
                
    return comp_edges


def get_sars2url():
    df = pd.read_excel('../convert_ids/SARS2URL.xlsx',header=0)
    return dict(zip(df['SARS'],df['URL']))

def main():
    parser = argparse.ArgumentParser("Input parameters")
    parser.add_argument("--sars_cov2_map", default="../convert_ids/SARS_COV2_Map.xlsx", help="Input parameters file name")
    parser.add_argument("--results_folder", default="../humap/results_73_neg_unif_10xisa_e0.01_T01.75_a0.005_qi_o0.375", help="Input parameters file name")
    parser.add_argument("--name2annot_file", default="../convert_ids/name2annot_full_humap_from_gene_names.pkl", help="Input parameters file name")
    parser.add_argument("--name2annot_covid_file", default='../convert_ids/name2annot_covid_interactors_from_gene_names.pkl', help="Input parameters file name")
    parser.add_argument("--prot2url_file", default="../convert_ids/name2uniprotURL_full_humap_from_gene_names.pkl", help="Input parameters file name")
    parser.add_argument("--prot2url_covid_file", default="../convert_ids/name2uniprotURL_covid_interactors_from_gene_names.pkl", help="Input parameters file name")
    parser.add_argument("--prot2comp_cov_file_name", default="/res_pred_prot2comp_covid.out", help="Input parameters file name")
    parser.add_argument("--comp_num2name", default="/res_pred_complex_num2name.pkl", help="Input parameters file name")
    
    args = parser.parse_args()
    
        
    sars2url = get_sars2url()
        
    df = pd.read_excel(args.sars_cov2_map,header=1)
    
    edges = list(zip(df['Bait'],df['PreyGene'],df['MIST']))
    
    nodes_sars = set(list(df['Bait']))
    
    nodes_human_prots = set(list(df['PreyGene']))
    
    sars_prot_wise_edges = dict()
    for edge in edges:
        sars_prot = edge[0]
        if sars_prot not in sars_prot_wise_edges:
            sars_prot_wise_edges[sars_prot] = [edge]  
        else:
            sars_prot_wise_edges[sars_prot].append(edge)
    
    results_folder = args.results_folder      
    prot2comp_cov_file = results_folder + args.prot2comp_cov_file_name
    
    with open(results_folder+args.comp_num2name,'rb') as f:
        comp_num2name = pickle.load(f)
        
    # generate protein to complexes edges
    with open(prot2comp_cov_file) as f:
        lines = [line.rstrip().split() for line in f.readlines()]
    
    edges_comp = []
    nodes_comp = []
    
    for line in lines:
        prot = line[0]
        for comp in line[2:]:
            edges_comp.append((prot,comp,0.99))
            nodes_comp.append(comp)
            
    prot2comps = dict([(line[0],line[2:]) for line in lines])
    
    prot2comp_edges = dict()
    for prot in prot2comps:
        prot2comp_edges[prot] = [(prot,comp,0.99) for comp in prot2comps[prot]]
        
    prots_in_comps = set(list(prot2comps.keys()))
    
    
    nodes_comp = set(nodes_comp)
    
    nodes_human_prots_all = nodes_human_prots.intersection(prots_in_comps)
    
    edges = [edge for edge in edges if edge[1] in nodes_human_prots_all]    
    
    with open(args.name2annot_file,'rb') as f:
        name2annot = pickle.load(f)
        
    with open(args.name2annot_covid_file,'rb') as f:
        name2annot2 = pickle.load(f)
        
    with open(args.prot2url_file,'rb') as f:
        prot2url = pickle.load(f)    
            
    with open(args.prot2url_covid_file,'rb') as f:
        prot2url2 = pickle.load(f)       
        
    for key, val in name2annot2.items():
        name2annot[key] = val
        
    for key, val in prot2url2.items():
        prot2url[key] = val
        
        
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
                          'text-wrap': 'wrap',
        'text-max-width': '300',
                shape: function( ele ){ 
                                    if ( ele.data('type') == 'complex' ){return 'circle'}
                                    if ( ele.data('type') == 'sars' ){return 'triangle'}
                                    
      return 'star';
      },
                'background-color': function( ele ){ 
                                    if ( ele.data('type') == 'complex' ){return 'blue'}
                                    if ( ele.data('type') == 'sars' ){return 'red'}             
                                    if ( ele.data('score_annot') == 1 ){return 'rgb(201, 160, 54)'} 
    								if ( ele.data('score_annot') == 2 ){return 'rgb(211, 186, 63)'}  
                                     if ( ele.data('score_annot') == 3 ){return 'rgb(160, 171, 105)'}   
                                     if ( ele.data('score_annot') == 4 ){return 'rgb(78, 148, 70)'}   
                                     if ( ele.data('score_annot') == 5 ){return 'rgb(52, 136, 74)'}   
                                     return 'rgb(180, 131, 40)';
      },
                label: 'data(id)'
            }
        },
    							{
    							selector: 'edge',
    							style: {
    								'width': 'mapData(weight, 0, 1, 0, 5)',
    								'line-color': function( ele ){ 
                                    if ( ele.data('type') == 'prot2complex' ){return 'orange'}
                                      if ( ele.data('type') == 'complex2complex' ){return 'blue'}
      return 'red';
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
    	'componentSpacing': 1,
    	 'idealEdgeLength': function( edge ){ return 1; },
    	 // Divisor to compute edge forces
      'edgeElasticity': function( edge ){ return 20; }
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
    
    with open(results_folder+'/res_pred_names.out') as f:
        comp2prots = dict([('Complex'+str(i+1),set(line.rstrip().split()[:-1])) for i,line in enumerate(f.readlines()) if 'Complex'+str(i+1) in nodes_comp])
    
    
    for sars_prot in nodes_sars:
    
        pre = pref + '<title>' + sars_prot +  '</title>\n'
        
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
        
        #res += "{data: { id: '" + sars_prot + "', score: " + str(2) + ", type: " + " 'sars' "
    
        res += '{data: { id: "' + sars_prot + '", score: ' + str(2) + ', type: ' + ' "sars" '
        if sars_prot in sars2url:
            res += ', href: "' + sars2url[sars_prot] + '"'
        res += " }}"
        res += ",\n"
        
        edges_sars_cov_prot = sars_prot_wise_edges[sars_prot]
        
        nodes_human_prots = [edgee[1] for edgee in edges_sars_cov_prot]
        
        nodes_comp = []
        edges_comp = []
        
        for prot in nodes_human_prots:
            if prot in prot2comps:
                nodes_comp = nodes_comp + prot2comps[prot]
                edges_comp = edges_comp + prot2comp_edges[prot]
            
        nodes_comp = set(nodes_comp)
        nNodes_comp = len(nodes_comp)
        
        nNodes = len(nodes_human_prots)
        
        for i,node in enumerate(nodes_human_prots):
            annot_score = 0
            if node in name2annot:
                annot_score = name2annot[node]          
            #res += "{data: { id: '" + node + "', score: " + str(1)+ ", score_annot: " + str(annot_score) + ", type: " + " 'human' "
            res += '{data: { id: "' + node + '", score: ' + str(1)+ ', score_annot: ' + str(annot_score) + ', type: ' + ' "human" '
            
            if node in prot2url:
                res += ', href: "' + prot2url[node] + '"'
                
            res += " }}"
            if i != nNodes-1:
                res += ",\n"
            else:
                if nNodes_comp == 0:
                    res += "], \n edges: ["  
                else:
                    res += ",\n"
    
        for i,node in enumerate(nodes_comp):
            #res += "{data: { id: '" + comp_num2name[node] + "', score: " + str(0) + ", type: " + "'complex'"
            res += '{data: { id: "' + comp_num2name[node] + '", score: ' + str(0) + ', type: ' + '"complex"'
            
            res += ', href: "' + "../Predicted/" + node + '.html' '"'
                
            res += " }}"
            if i != nNodes_comp-1:
                res += ",\n"
            else:
                res += "], \n edges: ["
        
        comp_edges=get_comp2comp_edges(nodes_comp,comp2prots)
                
        for j,elts in enumerate(comp_edges):
            #res += "{ data: { source: '" + comp_num2name[elts[0]] +"', target: '"  + comp_num2name[elts[1]] +"', directed: 'false'," + "type: 'complex2complex'," + " weight: " + str(elts[2]) + "} }"
            res += '{ data: { source: "' + comp_num2name[elts[0]] +'", target: "'  + comp_num2name[elts[1]] +'", directed: "false",' + 'type: "complex2complex",' + ' weight: ' + str(elts[2]) + '} }'
            
            res += ",\n"            
        
        for j,elts in enumerate(edges_comp):
            #res += "{ data: { source: '" + elts[0] +"', target: '" + comp_num2name[elts[1]] + "', directed: 'false'," + "type: 'prot2complex'," + " weight: " + str(elts[2]) + "} }"
            res += '{ data: { source: "' + elts[0] +'", target: "' + comp_num2name[elts[1]] + '", directed: "false",' + 'type: "prot2complex",' + ' weight: ' + str(elts[2]) + '} }'
            
            res += ",\n"
                
        nEdges = len(edges_sars_cov_prot)
        for j,elts in enumerate(edges_sars_cov_prot):
            #res += "{ data: { source: '" + elts[0] +"', target: '" + elts[1] + "', directed: 'false', " + "type: 'sars2prot'," + "weight: " + str(elts[2]) + "} }"
            res += '{ data: { source: "' + elts[0] +'", target: "' + elts[1] + '", directed: "false", ' + 'type: "sars2prot",' + 'weight: ' + str(elts[2]) + '} }'
            
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
        
        with open("./Complexes/SARS_COV2prots/" + sars_prot + ".html","w") as f:
            f.write(res)

if __name__ == "__main__":
    main()