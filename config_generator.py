#!/usr/bin/env python

from prody import *
from pylab import *
import numpy as np
from os.path import basename
import fnmatch
import os



def config_gen(pdb_name,min_step_count):
    file_name_wh_ex = str(os.path.splitext(pdb_name)[0])
    structure = parsePDB(str(pdb_name))
    f = open(str(file_name_wh_ex)+"_config.conf", 'w')
    f.write("#############################################################\n")
    f.write("## JOB DESCRIPTION                                         ##\n")
    f.write("#############################################################\n")
    f.write("\n")
    f.write("#%s minimization" % file_name_wh_ex)
    f.write("\n")
    f.write("\n")
    f.write("#############################################################\n")
    f.write("## ADJUSTABLE PARAMETERS                                   ##\n")
    f.write("#############################################################\n")
    f.write("\n")
    f.write("%-10s\t\t\t%-10s\n" % ("structure",str(file_name_wh_ex)+".psf"))
    f.write("%-10s\t\t\t%-10s\n" % ("coordinates",str(file_name_wh_ex)+".pdb"))
    f.write("\n")
    f.write("\n")
    f.write("\n")
    f.write("set temperature    310\n")
    f.write("set outputname     %s\n" % (str(file_name_wh_ex)))
    f.write("\n")
    f.write("\n")
    f.write("firsttimestep      0\n")
    f.write("#############################################################\n")
    f.write("## SIMULATION PARAMETERS                                   ##\n")
    f.write("#############################################################\n")
    f.write("# Input\n")
    f.write("paraTypeCharmm      on\n")
    f.write("parameters          par_all27_prot_lipid.inp\n")
    f.write("temperature         $temperature\n")
    f.write("\n")
    f.write("\n")
    f.write("# Force-Field Parameters\n")
    f.write("""
exclude             scaled1-4
1-4scaling          1.0
cutoff              12.
switching           on
switchdist          10.
pairlistdist        13.5
""")
    f.write("\n")
    f.write("\n")
    f.write("""
# Integrator Parameters
numsteps            1
timestep            2.0  ;# 2fs/step
rigidBonds          all  ;# needed for 2fs steps
nonbondedFreq       1
fullElectFrequency  2  
stepspercycle       10
""")
    f.write("\n")
    f.write("\n")
    f.write("""
# Constant Temperature Control
langevin            on  ;# do langevin dynamics
langevinDamping     5     ;# damping coefficient (gamma) of 5/ps
langevinTemp        $temperature
langevinHydrogen    off    ;# don't couple langevin bath to hydrogens
""")
    f.write("\n")
    f.write("\n")

    # TODO: Create a class, please :ç

    ###################
    pdb_atom_coords = structure.getCoords()
    xmax = np.max(pdb_atom_coords[:, 0])
    ymax = np.max(pdb_atom_coords[:, 1])
    zmax = np.max(pdb_atom_coords[:, 2])
    xmin = np.min(pdb_atom_coords[:, 0])
    ymin = np.min(pdb_atom_coords[:, 1])
    zmin = np.min(pdb_atom_coords[:, 2])
    cb_vec1 = abs(xmax-xmin)+0.1
    cb_vec2 = abs(ymax-ymin)+0.1
    cb_vec3 = abs(zmax-zmin)+0.1
    cb_vec1_f = "{0:.1f}".format(cb_vec1)
    cb_vec2_f = "{0:.1f}".format(cb_vec2)
    cb_vec3_f = "{0:.1f}".format(cb_vec3)
    co_x = (xmax+xmin)/2
    co_y = (ymax+ymin)/2
    co_z = (zmax+zmin)/2
    co_x_f = "{0:.1f}".format(co_x)
    co_y_f = "{0:.1f}".format(co_y)
    co_z_f = "{0:.1f}".format(co_z)
    ###################

    f.write("# Periodic Boundary Conditions\n")
    f.write("cellBasisVector1     %s      0.0    0.0\n" % cb_vec1_f)
    f.write("cellBasisVector2     0.0       %s   0.0\n" % cb_vec2_f)
    f.write("cellBasisVector3     0.0       0.0    %s\n"% cb_vec3_f)
    f.write("cellOrigin           %s     %s     %s\n" % (co_x_f, co_y_f, co_z_f))
    f.write("\n")
    f.write("\n")
    f.write("wrapall              on\n")
    f.write("\n")
    f.write("\n")
    f.write("""
# PME (for full-system periodic electrostatics)
PME                 yes
pmeGridSpacing      1.0
#PMEGridSizeX        90
#PMEGridSizeY       108
#PMEGridSizeZ       108


# Constant Pressure Control (variable volume)
useGroupPressure      yes ;# needed for rigidBonds
useFlexibleCell       no
useConstantArea       no

langevinPiston        yes
langevinPistonTarget  1.01325 ;#  in bar -> 1 atm
langevinPistonPeriod  100.
langevinPistonDecay   50.
langevinPistonTemp    $temperature


# Output
outputName          $outputname

restartfreq         1000     ;# 
dcdfreq             1000
DCDfile             $outputname.dcd
xstFreq             1000
outputEnergies      1000
outputPressure      1000


#############################################################
## EXTRA PARAMETERS                                        ##
#############################################################


#############################################################
## EXECUTION SCRIPT                                        ##
#########################################################
""")
    f.write("\n")
    f.write("\n")
    f.write("minimize %s" % (str(min_step_count)))

    f.close()

    return

global step_size
step_size = 5000

for file in os.listdir('.'):
    if fnmatch.fnmatch(file, 'protein_res*_mutated_wb_autopsf.pdb'):
        pdb = file
        config_gen(str(pdb),str(step_size))

