'''
This version plots the balance plot for epochs and rounds.

Started: 23 Oct 2019


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
import csv

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


def parse_round_file(ffullpath):
	
	xs=[]
	ys=[]
	with open(ffullpath) as csvfile:
		readCSV = csv.reader(csvfile, delimiter=',')
		for row in readCSV:
			for col in row:
				col = col.split(':')
				if (len(col)==2):
					xs.append(col[0])
					ys.append(col[1])

	return xs, ys



def parse_epoch_file(ffullpath):
	xs=[]
	ys=[]
	with open(ffullpath) as csvfile:
		readCSV = csv.reader(csvfile, delimiter=',')
		for row in readCSV:
			for col in row:
				col = col.split(':')
				if (len(col)==3):
					xs.append(col[2])
					ys.append(col[1])
	return xs, ys


subfolders = ['prediction0/']#, 'prediction1/']
mems = ['96/']#, '96/', '384/', 'unlim/']
nprocs = ['8', '16', '32', '64']
op_image_name = 'balance_'

plt.title(op_image_name)
plt.ylabel(op_image_name)
plt.xlabel('round')
for sf in subfolders:
	for mem in mems:

		folpath = in_folder+sf+mem
		files = os.listdir(folpath)
		for nproc in nprocs:
			print nproc
			
			rfiles = [file for file in files if 'round' in file and nproc in file]
			efiles = [file for file in files if 'balance' in file and nproc in file and 'round' not in file]


			round_xs, round_ys = parse_round_file(folpath+rfiles[0])
			epoch_xs, epoch_ys = parse_epoch_file(folpath+efiles[0])
			plt.plot(round_xs, round_ys)
			plt.plot(epoch_xs, epoch_ys)


plt.show()