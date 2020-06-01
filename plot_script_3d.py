#!/usr/bin/python

import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
import csv
import os
import subprocess
import math
from itertools import cycle

plt.rcParams.update({'font.size': 14})
#cmd = "./main 0.00 0.00 1.0 1.0 0 gridfile.dat"
#returned_value = os.system(cmd)  # returns the exit code in unix
returned_value = os.system("seq 0.0 0.001 3.0 > gridfile.dat")

ymin = -0.1; ymax = 6.0;
xmin = 0.0; xmax = 2.5;

qmid = [0.0, 0.0, 0.6];


#alpha_range = range(0,91,30);

t1 = math.degrees(math.acos(math.sqrt(2.0/3)));
t2 = math.degrees(math.acos(math.sqrt(1.0/3)));
#angles = [[0,0],[0,30],[0,60],[0,90],[30,90],[60,90],[90,90],[60,0],[30,0],[30,t1],[30,t2],[45,45]];

angles = [[90,0],[90,45],[90,90],[45,90],[0,0],[45,0],[55,45]];

colors = pl.cm.brg(np.linspace(0,1,len(angles)))

lines = ["-","--","-.",":"]
linecycler = cycle(lines)
	
for i in range(0,len(angles)):
	
	theta = angles[i][0]; phi = angles[i][1];
	cmd = "./main "+str(theta)+" "+str(phi)+" "+str(qmid[0])+" "+str(qmid[1])+" "+str(qmid[2])+" gridfile.dat | grep Energy"
	print(cmd)
	returned_value = os.system(cmd)
	#print(returned_value)
	
	hosts = subprocess.check_output(cmd, shell=True)
	error = float(hosts.split(' ')[1])-1
	error = abs(error)
	error2 = float(hosts.split(' ')[3])
	#print(error2)
	thetalabel = r'$(\theta,\varphi)=$'+'('+"{0:.0f}".format(theta)+r'$\degree$'+', '+"{:.0f}".format(phi)+r'$\degree$'+')'+r', $\varepsilon_1$ = '+"{:.1e}".format(error)
	thetatitle = r'Translation of parallel fibers'
	
	x = []
	y = []
	
	with open('outfile.gr','r') as csvfile:
	    next(csvfile); 
	    next(csvfile);
	    plots = csv.reader(csvfile, delimiter=' ')
	    for row in plots:
	        x.append(float(row[0]))
	        y.append(float(row[1]))
	
	plt.plot(x,y, next(linecycler), label=thetalabel, color=colors[i])
	#plt.plot(x,y, "-", label=thetalabel, color=colors[i])
	

plt.xlabel('r', fontsize=20)
plt.ylabel('g(r)', fontsize=20)
plt.ylim(ymin, ymax) 
plt.xlim(xmin, xmax) 
# ax[ip].set_title(thetatitle, fontsize=12)
plt.legend(loc=1, fontsize=12)


#plt.show()
#plt.plot(x,x, label='Loaded from file!')
#ip+=1

qmid = [0.0, 0.5, 0.0];
	
maintitle = "Distance distribution g(r) with "+r'$0.0 \leq \Delta x \leq 1.4$'
#plt.savefig('outfile.png')
#plt.suptitle(maintitle);
fileid = qmid[0]*9+qmid[1]*3+qmid[2];
filename = "plot_3d.eps"
plt.savefig(filename, format="eps", bbox_inches='tight')
#os.system("inkscape plot_pt.svg --export-pdf=plot_pt.pdf")
