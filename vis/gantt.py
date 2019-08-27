'''

Gantt chart

21 Aug 2019
updated: 27 Aug 2019

usage:

python gantt.py /Users/mukundraj/Desktop/work/projects/diy2-dlb-kd-new-copy/archive/bebop_190822/gantt constrained/48 8

'''

import numpy as np
from numpy import genfromtxt
import os
import matplotlib.pyplot as plt
import shutil
import sys
from matplotlib.ticker import MaxNLocator
import csv


ipath =  sys.argv[1] # /Users/mukundraj/Desktop/work/projects/diy2-dlb-kd-new/archive/bebop_190822/gantt/
prefix = sys.argv[2] # e.g. constrained/48
nprocs = sys.argv[3] # e.g. 8
ip_prefix = ipath+"/"+prefix+"/"
# ip_prefix = "/Users/mukundraj/Desktop/work/projects/diy2-dlb-kd-new/archive/bebop_190822/gantt/"+prefix+"/"
# nprocs = '8'


files = os.listdir(ip_prefix)

files = [x for x in files if x.find('gantt_'+nprocs)==0]

# Declaring a figure "gnt" 
fig, gnt = plt.subplots() 
gnt.grid(True) 



total_time = 0

for file in files:
	fullpath = ip_prefix+file
	print fullpath
	processed_pairs = []
	with open(fullpath) as csvfile:
		readCSV = csv.reader(csvfile)
		rank = int(file.split('_')[-1].split('.')[0])
		for row in readCSV:
			pairs = row

		for pair in pairs:
			pair = pair.split(':')
			pair = [int(pair[0]), float(pair[1])]
			total_time += pair[1]
			processed_pairs.append(pair)

	# plotting
	tim = 0
	for p in processed_pairs:

		if p[0] == 0:
			
			gnt.broken_barh([(0, p[1])], (10+10*rank, 10), facecolors =('purple'))
			tim += p[1]
		elif p[0] == 1:
			
			gnt.broken_barh([(tim, p[1])], (10+10*rank, 10), facecolors =('green'))
			tim += p[1]
		elif p[0] == 2:
			gnt.broken_barh([(tim, p[1])], (10+10*rank, 10), facecolors =('pink'))
			tim += p[1]
		elif p[0] == 3:
			gnt.broken_barh([(tim, p[1])], (10+10*rank, 10), facecolors =('orange'))
			tim += p[1]

gnt.set_ylim(0, 20 + int(nprocs)*10) 
  
# Setting X-axis limits 
gnt.set_xlim(0, 1.1*total_time/int(nprocs)) 

plt.show()