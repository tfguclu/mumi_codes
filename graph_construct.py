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
    carbon_betas = pdb.select("name CB").copy()
    gly_calphas = pdb.select("resname GLY and calpha").copy()
    nodes = carbon_betas + gly_calphas
    nodes = nodes.copy()
    nodes = sortAtoms(nodes, "resnum")
    nodes_range = len(nodes)
    nodes_list = nodes.getResnums()
    ia_list = []
    for i in range(nodes_range-1):
        for j in range(i+1, nodes_range):
            dist = calcDistance(nodes[i], nodes[j])
            if dist >= 6.7:
                continue
            ia_list.append((nodes[i].getResnum(), nodes[j].getResnum()))
    #print(sorted(ia_list))
    print(len(ia_list))
    #print(nodes.getResindices())
    print(len(nodes))
    #print(nodes.getResnums())
    #print(nodes.getDataLabels())
    #print(nodes[0].getFlagLabels(which="user"))
    protein_graph = nx.Graph()
    protein_graph.add_nodes_from(nodes_list)
    protein_graph.add_edges_from(ia_list)
    #for i in nodes_list:
    #    protein_graph.node[i]["residue number"] = int(nodes[i].getResnum())
    print(nx.dijkstra_path(protein_graph,1,54))
    #print(nodes[dijkstra_index].getResnums())
    #print(nodes.getResnums()[0])
    nx.draw(protein_graph, with_labels=True)
    plt.savefig("5rxn_wt_path.png")
    plt.close()
    return

graph_construct("5rxn_lf_prot.pdb")
