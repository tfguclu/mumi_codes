#!/usr/bin/env python

from prody import *
from pylab import *
import numpy as np
from os.path import basename
import fnmatch
import os

f = open("config_joblist", 'w')

f.write("#!/bin/bash\n")
f.write("echo config jobs started\n")
file_names_sorted = []
for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*_config.conf'):
        file_names_sorted.append(str(os.path.splitext(file)[0])) # no need to create variable to append

file_names_sorted = sorted(file_names_sorted, key=str.lower) # sorting should executed once

range_list = len(file_names_sorted)

for i in range(range_list):
    conffile = str(file_names_sorted[i] +".conf")
    logfile = str(file_names_sorted[i] +".log")
    f.write("namd2 +p8 %s > %s\n" % (str(conffile), str(logfile)))
    f.write("echo job for %s is finished\n" % str(conffile))

f.write(" echo all jobs are finished")
f.close()
