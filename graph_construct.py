from prody import *
from pylab import *
import numpy as np
from os.path import basename
import fnmatch
import os
import networkx as nx
import matplotlib.pyplot as plt


def graph_construct(pdb_file):
    pdb = parsePDB(str(pdb_file))
    file_name_wh_ex = str(os.path.splitext(pdb_file)[0])
    #carbon_betas = pdb.select("name CB").copy()
    #gly_calphas = pdb.select("resname GLY and calpha").copy()
    #nodes = carbon_betas + gly_calphas
    betas_gly_alphas = pdb.select("(name CB and protein) or (name CA and resname GLY)").copy()
    nodes = betas_gly_alphas
    nodes = sortAtoms(nodes, "resnum")
    nodes_range = len(nodes)
    nodes_list = nodes.getResindices()
    ia_list = []
    for i in range(nodes_range-1):
        for j in range(i+1, nodes_range):
            dist = calcDistance(nodes[i], nodes[j])
            if dist >= 6.7:
                continue
            ia_list.append((nodes[i].getIndex(), nodes[j].getIndex()))
    #print(sorted(ia_list))
    print(len(ia_list))
    #print(nodes.getResindices())
    print(len(nodes))
    #print(nodes.getResnums())
    #print(nodes.getDataLabels())
    #print(nodes[0].getFlagLabels(which="user"))
    protein_graph = nx.Graph()
    dj_path_matrix = np.zeros((int(nodes_range),int(nodes_range)))
    protein_graph.add_nodes_from(nodes_list)
    protein_graph.add_edges_from(ia_list)
    for i in range(nodes_range):
        for j in range(nodes_range):
            dj_path_matrix[i,j] = nx.dijkstra_path_length(protein_graph,i,j)

    print((dj_path_matrix/(nodes_range*(nodes_range-1))).sum())

    #for i in nodes_list:
    #    protein_graph.node[i]["residue number"] = int(nodes[i].getResnum())
    #print(nx.dijkstra_path(protein_graph,1,54))
    print(nx.average_shortest_path_length(protein_graph))
    print(nodes[dj_path_matrix.max()].getResnum())
    #all_pairs = nx.all_pairs_shortest_path_length(protein_graph)
    #print(len(all_pairs))
    #print(nx.betweenness_centrality(protein_graph, normalized=False))
    #print(nodes[dijkstra_index].getResnums())
    #print(nodes.getResnums()[0])

    betwen_cent = nx.betweenness_centrality(protein_graph)
    for i in range(nodes_range):
        betwen_cent[nodes[i].getResnum()] = betwen_cent.pop(i)
    plt.plot(betwen_cent.keys(),betwen_cent.values())
    plt.title(str(file_name_wh_ex)+"_BC.png")
    plt.savefig(str(file_name_wh_ex)+"_BC.png",dpi=300, bbox_inches='tight')
    plt.close()

    avg_len_per_node_list = []
    for i in range(nodes_range):
        avg_len_per_node = float((dj_path_matrix[i,:].sum())/nodes_range-1)
        avg_len_per_node_list.append(float(avg_len_per_node))
    plt.plot(betwen_cent.keys(),avg_len_per_node_list)
    plt.title(str(file_name_wh_ex)+"_L.png")
    plt.savefig(str(file_name_wh_ex)+"_L.png",dpi=300, bbox_inches='tight')
    plt.close()

    plt.close()
    nx.draw(protein_graph, with_labels=True, node_size=900)
    plt.savefig(str(file_name_wh_ex)+"_path.png")
    plt.close()
    return

for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*.pdb'):
        pdb = file
        graph_construct(str(pdb))
