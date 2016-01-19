#!/usr/bin/env python

from prody import *
from pylab import *
import numpy as np
from os.path import basename
import fnmatch
import os


def lf_extract(pdb_name):
    file_name_wh_ex = str(os.path.splitext(pdb_name)[0])
    structure = parsePDB(str(pdb_name))
    traj = Trajectory(str(file_name_wh_ex)+".dcd")
    traj.link(structure)
    last_frame = traj[-1]
    writePDB(str(file_name_wh_ex)+"_lf.pdb", last_frame)
    return

for file in os.listdir('.'):
    if fnmatch.fnmatch(file, 'protein_res*_mutated_wb_autopsf.pdb'):
        pdb = file
        lf_extract(str(pdb))
