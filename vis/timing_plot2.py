'''
This version plots the timing plot for bebop output files

Started: 21 Aug 2019
Updated: 27 Aug 2019 : read input folder from command line


usage:
python timing_plot2.py /Users/mukundraj/Desktop/work/projects/diy2-dlb-kd-new-copy/archive/bebop_190823/nogantt

python timing_plot2.py /Users/mukundraj/Desktop/work/results/diy2-dlb-kd-new/bebop_190827/nogantt

'''

import numpy as np
import os
import matplotlib.pyplot as plt
import shutil
import sys
from matplotlib.ticker import MaxNLocator

# in_folder = "/Users/mukundraj/Desktop/work/projects/diy2-dlb-kd-new/archive/bebop_190822/nogantt/"
# in_folder = "/Users/mukundraj/Desktop/work/projects/diy2-dlb-kd-new-copy/archive/bebop_190823/nogantt/"

in_folder =  sys.argv[1] + "/"


# subfolders = ['prediction0/', 'prediction5/', 'prediction10/', 'prediction20/']
subfolders = ['prediction0/', 'prediction1/']
# subfolders = ['baseline/', 'constrained/']
mems = ['48/', '96/', '384/']

op_image_name = "timing.png"

def kfunc(val):
	return int(val.split('_')[1])

def get_timings(fpath):

	files = os.listdir(fpath)
	files = [file for file in files if file != ".DS_Store"]
	files.sort(key=kfunc)
	times = []
	time_trace = 0
	time_kdtree = 0
	time_prediction = 0
	load_balance_indicat = 0
	maxbyavg = 0
	for file in files:
		with open(fpath+file) as fp:
		   line = fp.readline()
		   while line:
		   		
		   		sline = line.split('\t')
		   		print sline
		   		if sline[0]=='time_init=' :
		   			time_init = float(sline[2].rstrip())
		   		if sline[0]=='time_purerun=' :
		   			time_purerun = float(sline[2].rstrip())
		   		if sline[0]=='time_all=' :
		   			time_all = float(sline[2].rstrip())
		   		if sline[0]=='time_io=' :
		   			time_io = float(sline[2].rstrip())
		   		if sline[0]=='time_run=' :
		   			time_run = float(sline[2].rstrip())
		   		if sline[0]=='_time_trace=' :
		   			time_trace = float(sline[2].rstrip())
		   		if sline[0]=='_time_kdtree=' :
		   			time_kdtree = float(sline[2].rstrip())
		   		if sline[0]=='_time_prediction=' :
		   			time_prediction = float(sline[2].rstrip())
		   		if sline[0]=='load_balance_indicat=' :
		   			load_balance_indicat = float(sline[1].rstrip())
		   		if sline[0]=='max/avg_integr_steps=' :
		   			maxbyavg = float(sline[1].rstrip())
		   		line = fp.readline()
		   times.append(time_trace)
	return times



plt.title(op_image_name)
plt.ylabel('time (s)')
plt.xlabel('nprocs')



for sf in subfolders:
	for mem in mems:
		folpath = in_folder+sf+mem
		tim = get_timings(folpath)
		if len(tim) == 7:
			xs = [1, 2, 4, 8, 16, 32, 64]
		elif len(tim) == 6:
			xs = [16, 32, 64, 128, 256, 512]	
		elif len(tim) == 10:
			xs = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
		elif len(tim) == 9:
			xs = [2, 4, 8, 16, 32, 64, 128, 256, 512]
		elif len(tim) == 5:
			xs = [32, 64, 128, 256, 512]

		if sf=='prediction0/':
			styl = '-'
		if sf=='prediction1/':
			styl = ':'
		elif sf=='prediction5/':
			styl = '--'
		elif sf=='prediction10/':
			styl = '-.'
		elif sf=='prediction20/':
			styl = ':'
		else:
			styl = "-"

		if mem == '48/':
			col = 'red'
		elif mem == '96/':
			col = 'blue'
		elif mem == '384/':
			col = 'green'
		elif mem == 'unlim/':
			col = 'orange'


		lab = sf+mem
		print lab, tim
		plt.plot(xs, tim, label=lab, linestyle=styl, linewidth=2, color=col)

# plt.gca().set_ylim(0, 400)

plt.gca().set_xscale('log')

plt.legend(loc='best', fontsize=12)
plt.gca().minorticks_off()
plt.gca().set_xticks(xs)
plt.gca().set_xticklabels(xs)
plt.show()
plt.savefig(in_folder + op_image_name)


