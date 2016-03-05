from prody import *
from pylab import *
import numpy as np
from os.path import basename
import fnmatch
import os
import networkx as nx
import matplotlib.pyplot as plt



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

    return (nat_file_name_wh_ex, nat_bc_list, nat_len_per_node_list)


def graph_construct(pdb_pattern,):


    average_lenght_list = np.zeros(nodes_length,nodes_length)


    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, str(pdb_pattern)):
            pdb_file = file


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

    return


list_of_nat, file_name_of_nat = nat_graph_construct("1BE9.pdb")
