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
    betas_gly_alphas = pdb.select("(name CB and protein) or (name CA and resname GLY)").copy()
    nodes = betas_gly_alphas
    nodes = sortAtoms(nodes, "resnum")
    nodes_range = len(nodes)
    nodes_list = nodes.getResindices()
    nodes_resnum_list = list(nodes.getResnums())


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


    dj_path_matrix = np.zeros((int(nodes_range),int(nodes_range)))
    for i in range(nodes_range):
        for j in range(nodes_range):
            dj_path_matrix[i,j] = nx.dijkstra_path_length(protein_graph,i,j)
    np.savetxt(str(file_name_wh_ex)+"_dj_path_matrix.dat", dj_path_matrix)
    avg_len_per_node_list = []
    for i in range(nodes_range):
        avg_len_per_node = float((dj_path_matrix[i,:].sum())/nodes_range-1)
        avg_len_per_node_list.append(float(avg_len_per_node))


    betwen_cent = nx.betweenness_centrality(protein_graph, normalized=False)
    betwen_cent = list(betwen_cent.values())
    betwen_cent_norm_fac = ((int(nodes_range)-1)*int(nodes_range))/2
    betwen_cent = [x / betwen_cent_norm_fac for x in betwen_cent]


    plt.plot(nodes_resnum_list,betwen_cent)
    plt.title(str(file_name_wh_ex)+"_BC", fontsize=18)
    plt.xlabel('Residue Numbers', fontsize=16)
    plt.ylabel('BC', fontsize=16)
    plt.savefig(str(file_name_wh_ex)+"_BC.png",dpi=300, bbox_inches='tight')
    plt.close()



    plt.plot(nodes_resnum_list, avg_len_per_node_list)
    plt.title(str(file_name_wh_ex)+"_L", fontsize=18)
    plt.xlabel('Residue Numbers', fontsize=16)
    plt.ylabel('L', fontsize=16)
    plt.savefig(str(file_name_wh_ex)+"_L.png",dpi=300, bbox_inches='tight')
    plt.close()


    nx.draw(protein_graph, with_labels=True, node_size=900)
    plt.savefig(str(file_name_wh_ex)+"_path.png")
    plt.close()

    nat_file_name_wh_ex = file_name_wh_ex
    nat_len_per_node_list = np.asarray(avg_len_per_node_list)
    nat_bc_list = np.asarray(betwen_cent)
    nodes_length = nodes_range

    return (nodes_length, nat_file_name_wh_ex, nat_bc_list, nat_len_per_node_list)


def graph_construct(nat_len_per_node_list, pdb_pattern):


    delta_lenght_per_node_dictionary = {}

    file_names_sorted = []
    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, str(pdb_pattern)):
            file_name = file
            file_names_sorted.append(file_name)
    file_names_sorted = natsorted(file_names_sorted)
    print(file_names_sorted)
    file_range = len(file_names_sorted)


    for i in range(int(file_range)):

        pdb = parsePDB(str(file_names_sorted[i]))
        file_name_wh_ex = str(os.path.splitext(file_names_sorted[i])[0])
        betas_gly_alphas = pdb.select("(name CB and protein) or (name CA and resname GLY)").copy()
        nodes = betas_gly_alphas
        nodes = sortAtoms(nodes, "resnum")
        nodes_range = len(nodes)
        nodes_list = nodes.getResindices()
        nodes_resnum_list = list(nodes.getResnums())


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


        dj_path_matrix = np.zeros((int(nodes_range),int(nodes_range)))
        for i in range(nodes_range):
            for j in range(nodes_range):
                dj_path_matrix[i,j] = nx.dijkstra_path_length(protein_graph,i,j)
        np.savetxt(str(file_name_wh_ex)+"_dj_path_matrix.dat", dj_path_matrix)
        avg_len_per_node_list = []
        for i in range(nodes_range):
            avg_len_per_node = float((dj_path_matrix[i,:].sum())/nodes_range-1)
            avg_len_per_node_list.append(float(avg_len_per_node))

        nat_len_per_node_list = np.asarray(nat_len_per_node_list)
        avg_len_per_node_list = np.asarray(avg_len_per_node_list)
        delta_lenght_per_node = np.absolute(np.subtract(nat_len_per_node_list, avg_len_per_node))


        betwen_cent = nx.betweenness_centrality(protein_graph, normalized=False)
        betwen_cent = list(betwen_cent.values())
        betwen_cent_norm_fac = ((int(nodes_range)-1)*int(nodes_range))/2
        betwen_cent = [x / betwen_cent_norm_fac for x in betwen_cent]


        plt.plot(nodes_resnum_list,betwen_cent)
        plt.title(str(file_name_wh_ex)+"_BC", fontsize=18)
        plt.xlabel('Residue Numbers', fontsize=16)
        plt.ylabel('BC', fontsize=16)
        plt.savefig(str(file_name_wh_ex)+"_BC.png",dpi=300, bbox_inches='tight')
        plt.close()



        plt.plot(nodes_resnum_list, avg_len_per_node_list)
        plt.title(str(file_name_wh_ex)+"_L", fontsize=18)
        plt.xlabel('Residue Numbers', fontsize=16)
        plt.ylabel('L', fontsize=16)
        plt.savefig(str(file_name_wh_ex)+"_L.png",dpi=300, bbox_inches='tight')
        plt.close()


        nx.draw(protein_graph, with_labels=True, node_size=900)
        plt.savefig(str(file_name_wh_ex)+"_path.png")
        plt.close()

        delta_lenght_per_node_dictionary[str(file_name_wh_ex)] = delta_lenght_per_node

    return delta_lenght_per_node_dictionary, nodes_resnum_list


nodes_length, nat_file_name_wh_ex, nat_bc_list, nat_len_per_node_list = nat_graph_construct("1be9_wt.pdb")

delta_lenght_per_node_dictionary, nodes_resnum_list = graph_construct(nat_len_per_node_list, "protein_res*_mutated_autopsf_wb_ionized_lf.pdb")

delta_lenght_per_node_array = np.array(delta_lenght_per_node_dictionary.values()) 
delta_lenght_average_per_node = np.mean(delta_lenght_per_node_array, axis=0)


plt.plot(nodes_resnum_list, delta_lenght_average_per_node)
plt.title("Delta L", fontsize=18)
plt.xlabel('Residue Numbers', fontsize=16)
plt.ylabel('$\Delta$ L', fontsize=16)
plt.savefig("average_L.png",dpi=300, bbox_inches='tight')
plt.close()

np.savetxt("delta_lenght_per_node", delta_lenght_per_node_array)