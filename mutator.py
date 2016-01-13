#!/usr/bin/env python

from prody import *
from pylab import *
import numpy as np
from os.path import basename
import fnmatch
import os


def mutator(pdb_file):

    structure = parsePDB(str(pdb_file))
    water = structure.select("water")
    ions = structure.select("ion")

    global water_box, protein

    water_box = water + ions
    water_box = water_box.copy()

    protein = structure.select("protein")
    protein = protein.copy()

    res_nums = protein.getResnums()
    res_nums = sorted(set(res_nums))
    res_nums = array(res_nums)

    for i in res_nums:

        protein_wh_res_sd = protein.select("protein and not resnum %s sidechain" % str(i))
        protein_wh_res_sd = protein_wh_res_sd.copy()

        res_id_mut_slc = protein_wh_res_sd.select("resnum %s" % str(i))
        res_id_mut_slc = res_id_mut_slc.setResnames("ALA")

        mutated_protein_structure = protein_wh_res_sd.copy()
        mutated_protein_structure_wb = mutated_protein_structure + water_box
        mutated_protein_structure_wb = mutated_protein_structure_wb.copy()
        writePDB("protein_res%s_mutated_wb.pdb" % str(i), mutated_protein_structure_wb)

mutator("1lit_min_last_frame.pdb")
