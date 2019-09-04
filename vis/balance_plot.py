'''
This version plots the balance plot for bebop output files

Started: 30 Aug 2019


usage:

python balance_plot.py /Users/mukundraj/Desktop/work/results/diy2-dlb-kd-new/bebop_190829c/nogantt 512

'''

import numpy as np
import os
import matplotlib.pyplot as plt
import shutil
import sys
from matplotlib.ticker import MaxNLocator
import ntpath
from shutil import copyfile

in_folder =  sys.argv[1] + "/"
nproc = sys.argv[2]


def get_time_and_balance(ffullpath):

	with open(ffullpath) as fp:
		line = fp.readline()
		while line:
		   	line = line.split(',')
		   	zero = float(line[0].split(':')[0])
		   	time = [float(x.split(':')[0])-zero for x in line[:-1]]
		   	balance = [float(x.split(':')[1]) for x in line[:-1]]
		   	line = fp.readline()

	return time, balance

subfolders = ['prediction0/', 'prediction1/']
mems = ['48/', '96/', '384/', 'unlim/']
op_image_name = 'balance_'+nproc

plt.title(op_image_name)
plt.ylabel(op_image_name)
plt.xlabel('round')

for sf in subfolders:
	for mem in mems:
		folpath = in_folder+sf+mem

		files = os.listdir(folpath)
		files = [file for file in files if 'balance' in file and nproc in file]


	
		time, balance = get_time_and_balance(folpath+files[0])

		if sf=='prediction0/':
			styl = '-'
		if sf=='prediction1/':
			styl = ':'
		elif sf=='prediction5/':
			styl = '--'
		elif sf=='prediction10/':
			styl = '-.'
		# elif sf=='prediction20/':
		# 	styl = ':'
		# else:
		# 	styl = "-"

		if mem == '48/':
			col = 'red'
		elif mem == '96/':
			col = 'blue'
		elif mem == '384/':
			col = 'green'
		elif mem == 'unlim/':
			col = 'orange'

		lab = sf+mem
		x = range(len(time))
		plt.plot(x, balance, label=lab, linestyle=styl, linewidth=2, color=col)


# plt.gca().set_xscale('log')

plt.legend(loc='best', fontsize=12)
plt.gca().minorticks_off()
# plt.gca().set_xticks(xs)
# plt.gca().set_xticklabels(xs)

plt.savefig(in_folder + op_image_name + ".png")

cur_fname = os.path.abspath(__file__)
scrpt_fname = in_folder+op_image_name+".py"
if cur_fname != scrpt_fname:
	copyfile(os.path.abspath(__file__), scrpt_fname)

plt.show()