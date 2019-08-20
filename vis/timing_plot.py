# timing plot for baseline, 512, and 1024

import numpy as np
import os
import matplotlib.pyplot as plt
import shutil
import sys
from matplotlib.ticker import MaxNLocator


path_base = "/Users/mukundraj/Desktop/work/projects/diy2-dlb-kd-new/results/baseline/results/"
path_512 = "/Users/mukundraj/Desktop/work/projects/diy2-dlb-kd-new/results/constrained/512/"
path_1024 = "/Users/mukundraj/Desktop/work/projects/diy2-dlb-kd-new/results/constrained/1024/"

op_path = "/Users/mukundraj/Desktop/work/projects/diy2-dlb-kd-new/results/16m/"
op_image_name = "timing_16m.png"

def get_timings(fpath):

	files = os.listdir(fpath)
	files.sort()
	times = []
	for file in files:
		with open(fpath+file) as fp:
		   line = fp.readline()
		   while line:
		   		
		   		sline = line.split('\t')
		   		if sline[0]=='time_run=' :
		   			times.append(float(sline[2].rstrip()))
		   		line = fp.readline()

	return times


times_base = get_timings(path_base)
times_512 = get_timings(path_512)
times_1024 = get_timings(path_1024)

plt.title(op_image_name)
plt.ylabel('time (s)')
plt.xlabel('nprocs')


xs = [2, 4, 8, 16, 32]
plt.plot(xs, times_base, label = 'base')
plt.plot(xs, times_512, label = '512')
plt.plot(xs, times_1024, label = '1024')

plt.legend(loc='best', fontsize=12)


plt.savefig(op_path+op_image_name)





