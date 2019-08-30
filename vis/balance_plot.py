'''
This version plots the balance plot for bebop output files

Started: 30 Aug 2019


usage:

python balance_plot.py /Users/mukundraj/Desktop/work/results/diy2-dlb-kd-new/bebop_190829c/nogantt

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

def get_time_and_balance():

	pass


subfolders = ['prediction0/', 'prediction1/']
mems = ['48/', '96/', '384/', 'unlim/']

for sf in subfolders:
	for mem in mems:
		folpath = in_folder+sf+mem
		files = os.listdir(folpath)
		files = [file for file in files if 'balance' in file]

		print files
