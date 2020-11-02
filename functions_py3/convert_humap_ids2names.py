from networkx.relabel import relabel_nodes as nx_relabel_nodes
from networkx import write_weighted_edgelist as nx_write_weighted_edgelist


def get_word_names(words, id_name_map):
    word_names = []
    for node in words:
        if node == "None":
            continue
        if node not in id_name_map or id_name_map[node] == '-':
            print(node)
            node_nm = node
        else:
            node_nm = id_name_map[node]
        word_names.append(node_nm)
    return frozenset(word_names)


def read_gene_id_map():
    ids = []
    names = []
    with open("./humap/humap_gene_id_name_map.txt") as f:
        for line in f:
            entries = line.split("\t")
            ids.append(entries[0])
            names.append(entries[1])
    del ids[0]
    del names[0] # Since header

    id_name_map = dict(zip(ids,names))
    return id_name_map


def convert_nodes(complexes, filename, id_name_map):
    unique_lines = set()
    for comp in complexes:
        word_names = get_word_names(comp, id_name_map)
        unique_lines.add(word_names)
    lines = []
    for myset in unique_lines:
        lines.append(" ".join(myset))
    edge_nm = "\n".join(lines)

    # Remember order is not maintained and duplicate names and complexes have been removed
    with open(filename, "w") as f:
        f.write(edge_nm)


def convert_edges(complex_graphs, filename, id_name_map):
    with open(filename, "wb") as f_edges:
        f_edges_write = f_edges.write
        for comp in complex_graphs:
            # Convert graph names
            comp = nx_relabel_nodes(comp, id_name_map)
            nx_write_weighted_edgelist(comp, f_edges)
            f_edges_write("\n".encode())


def convert2names(complexes,filename, complex_graphs, filename_edges):
    id_name_map = read_gene_id_map()
    convert_edges(complex_graphs, filename_edges, id_name_map)
    convert_nodes(complexes, filename, id_name_map)


def convert_nodes_wscore_unordered(complexes, filename, id_name_map):
    unique_lines = set()
    for comp in complexes:
        word_names = set(get_word_names(comp[0], id_name_map))
        word_names.add(comp[1])
        unique_lines.add(frozenset(word_names))
    lines = []
    for myset in unique_lines:
        score = 1
        myset = set(myset)
        for word in myset:
            if isinstance(word, float):
                score = word
                myset.remove(word)
                break

        lines.append(" ".join(myset) + " " + str(score))
    edge_nm = "\n".join(lines)

    # Remember order is not maintained and duplicate names and complexes have been removed
    with open(filename, "w") as f:
        f.write(edge_nm)

def convert_nodes_wscore(complexes, filename, id_name_map):
    # Order by score
    unique_lines = set()
    for comp in complexes:
        word_names = set(get_word_names(comp[0], id_name_map))
        word_names.add(comp[1])
        unique_lines.add(frozenset(word_names))
    lines = []
    for myset in unique_lines:
        score = 1
        myset = set(myset)
        for word in myset:
            if isinstance(word, float):
                score = word
                myset.remove(word)
                break
        lines.append((myset, score))
        
    lines = sorted(lines, key=lambda x: x[1], reverse=True)

    fin_lines = []
    for line in lines:
        fin_lines.append(" ".join(line[0]) + " " + str(line[1]))
    edge_nm = "\n".join(fin_lines)

    # Remember order is not maintained and duplicate names and complexes have been removed
    with open(filename, "w") as f:
        f.write(edge_nm)
        
def convert_nodes_matches_wscore(complex_matches, filename, id_name_map):
    # Order by score
    unique_lines = set()
    
    fin_lines = ["Known complex nodes ||| Predicted complex nodes ||| Match F1 score ||| Complex score \n"]    
    for comp in complex_matches:
        known_graph_nodes = comp[0]
        pred_graph_nodes = comp[1]
        match_score = comp[2]
        complex_score = comp[3]        
        
        known_graph_nodes = list(set(get_word_names(known_graph_nodes, id_name_map)))
        
        if frozenset(known_graph_nodes) in unique_lines:
            continue
        
        unique_lines.add(frozenset(known_graph_nodes))
        
        pred_graph_nodes = list(set(get_word_names(pred_graph_nodes, id_name_map)))
        
        fin_lines.append(" ".join(known_graph_nodes) + " ||| " + " ".join(pred_graph_nodes) + " ||| " + str(match_score) + " ||| " + str(complex_score))
    edge_nm = "\n".join(fin_lines)

    # Remember order is not maintained and duplicate names and complexes have been removed
    with open(filename, "w") as f:
        f.write(edge_nm)        
        
def convert_edges_wscore(complexes, G, filename, id_name_map):
    with open(filename, "wb") as f_edges:
        f_edges_write = f_edges.write
        for comp in complexes:
            # Convert graph names
            comp_nam = nx_relabel_nodes(G.subgraph(comp[0]), id_name_map)
            nx_write_weighted_edgelist(comp_nam, f_edges)
            eos = "Score = " + str(comp[1]) + "\n\n"
            f_edges_write(eos.encode())


def convert2names_wscores(complexes, filename, G, filename_edges):
    id_name_map = read_gene_id_map()
    convert_edges_wscore(complexes, G, filename_edges, id_name_map)
    convert_nodes_wscore(complexes, filename, id_name_map)
    
def convert2names_wscores_matches(complex_matches, filename):
    id_name_map = read_gene_id_map()
    convert_nodes_matches_wscore(complex_matches, filename, id_name_map)    
