from prody import *
from pylab import *
import numpy as np
from os.path import basename
import fnmatch
import os
import networkx as nx
import matplotlib.pyplot as plt
from natsort import natsorted




def graph_construct(pdb_pattern):

    nodes_number = 115
    delta_lenght_average_whole_array = np.zeros((int(nodes_number),int(nodes_number)))
    L_mutant_matrix = np.zeros((int(nodes_number),int(nodes_number)))

    nat = np.loadtxt("1be9_chainA_minim_laststep_avg_lenght.dat")

    file_names_sorted = []
    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, str(pdb_pattern)):
            file_name = file
            file_names_sorted.append(file_name)
    file_names_sorted = natsorted(file_names_sorted)
    print(file_names_sorted)
    file_range = len(file_names_sorted)


    for main_looper in range(int(file_range)):
        mut = np.loadtxt(str(file_names_sorted[main_looper]))
        diference = mut - nat
        diference = diference.reshape(1,int(nodes_number))


        L_mutant_matrix[main_looper,:] = mut
        delta_lenght_average_whole_array[main_looper,:] = diference


    nodes_resnum_list = (list(range(int(nodes_number))))
    np.savetxt("delta_lenght_average_whole_array.dat", delta_lenght_average_whole_array)
    delta_lenght_full_mean = np.absolute(np.mean(delta_lenght_average_whole_array, axis=0))
    np.savetxt("delta_lenght_full_mean.dat", delta_lenght_full_mean)
    np.savetxt("L_mutant_matrix.dat", L_mutant_matrix)
    plt.plot(nodes_resnum_list, delta_lenght_full_mean)
    plt.title("Delta L", fontsize=18)
    plt.xlabel('Node Indices', fontsize=16)
    plt.ylabel('$\Delta$ L', fontsize=16)
#    plt.ylim([0,1])
    plt.savefig("delta_L_all.png",dpi=300, bbox_inches='tight')
    plt.close()
    return

graph_construct("*ALA_laststep_avg_lenght.dat")