from prody import *
from pylab import *
import numpy as np
from os.path import basename
import fnmatch
import os
import networkx as nx
import matplotlib.pyplot as plt
from natsort import natsorted



def nat_graph_construct(pdb_file):


    pdb = parsePDB(str(pdb_file))
    file_name_wh_ex = str(os.path.splitext(pdb_file)[0])
    betas_gly_alphas_ligands = pdb.select("(name CB and protein) or (name CA and resname GLY)").copy()
    nodes = betas_gly_alphas_ligands
    nodes = sortAtoms(nodes, "resnum")
    nodes_range = len(nodes)
    nodes_list = range(0,int(len(nodes)))
    nodes_list_graph = [x+1 for x in nodes_list]

##########################################

    ia_list = []
    for i in range(nodes_range-1):
        for j in range(i+1, nodes_range):
            dist = calcDistance(nodes[i], nodes[j])
            if dist >= 6.7:
                continue
            ia_list.append((nodes[i].getIndex(), nodes[j].getIndex()))

    protein_graph = nx.Graph()
    protein_graph.add_nodes_from(nodes_list)
    protein_graph.add_edges_from(ia_list)

    nx.write_gml(protein_graph,str(file_name_wh_ex)+"_graph.gml")
    nx.write_graphml(protein_graph,str(file_name_wh_ex)+"_graph.graphml")

###########################################

    dj_path_matrix = np.zeros((int(nodes_range),int(nodes_range)))
    for i in range(nodes_range):
        for j in range(nodes_range):
            dj_path_matrix[i,j] = nx.dijkstra_path_length(protein_graph,i,j)
    np.savetxt(str(file_name_wh_ex)+"_dj_path_matrix.dat", dj_path_matrix)
    avg_len_per_node_list = (np.sum(dj_path_matrix, axis=0))/(int(nodes_range)-1)


    betwen_cent = nx.edge_betweenness_centrality(protein_graph, normalized=True, weight="0")
    betwen_cent = np.asarray(list(betwen_cent.values()))

############################################

    plt.plot(nodes_list_graph, betwen_cent)
    plt.title(str(file_name_wh_ex)+"_BC", fontsize=18)
    plt.xlabel('Node Indices', fontsize=16)
    plt.ylabel('BC', fontsize=16)
    plt.savefig(str(file_name_wh_ex)+"_BC.png",dpi=300, bbox_inches='tight')
    plt.close()


    plt.plot(nodes_list_graph, avg_len_per_node_list)
    plt.title(str(file_name_wh_ex)+"_L", fontsize=18)
    plt.xlabel('Node Indices', fontsize=16)
    plt.ylabel('L', fontsize=16)
    plt.savefig(str(file_name_wh_ex)+"_L.png",dpi=300, bbox_inches='tight')
    plt.close()

############################################

    avg_len_per_node_list = avg_len_per_node_list.reshape(1,int(nodes_range))
    betwen_cent = betwen_cent.reshape(1,int(nodes_range))
    np.savetxt(str(file_name_wh_ex)+"_avg_lenght.dat", avg_len_per_node_list)
    np.savetxt(str(file_name_wh_ex)+"_bc.dat", betwen_cent)
    return



# file_names_sorted = []
# for file in os.listdir('.'):
#     if fnmatch.fnmatch(file, "*.pdb"):
#         file_name = file
#         file_names_sorted.append(file_name)
# file_names_sorted = natsorted(file_names_sorted)
# print(file_names_sorted)
# file_range = len(file_names_sorted)
# for main_looper in range(int(file_range)):
#     nat_graph_construct(str(file_names_sorted[main_looper]))

nat_graph_construct("1be9_chainA_minim_laststep.pdb")