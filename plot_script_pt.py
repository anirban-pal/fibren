#!/usr/bin/python

import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
import csv
import os
import subprocess
from itertools import cycle

plt.rcParams.update({'font.size': 14})
#cmd = "./main 0.00 0.00 1.0 1.0 0 gridfile.dat"
#returned_value = os.system(cmd)  # returns the exit code in unix
returned_value = os.system("seq 0.0 0.001 3.0 > gridfile.dat")

ymin = -0.1; ymax = 5.0;
xmin = 0.25; xmax = 4.0;

#qmid = [1.1, 0.0, 0.0]; ymax = 2.5; xmax = 2.5; xmin = 0.0;
#qmid = [1.1, 1.1, 0.0]; ymax = 4.0; xmax = 2.5; xmin = 0.5;
#qmid = [0.0, 1.1, 0.0]; ymax = 25.0; xmax = 2.0; xmin = 0.5;

#qmid = [1.1, 0.0, 1.1]; ymax = 3.5; xmax = 3.0; xmin = 0.5;
#qmid = [1.1, 1.1, 1.1]; ymax = 5.0; xmax = 3.0; xmin = 1.0;
#qmid = [0.0, 1.1, 1.1]; ymax = 10.0;

q0_range = [x / 100.0 for x in range(0, 151, 30)]
print(q0_range)
colors = pl.cm.brg(np.linspace(0,1,len(q0_range)))

#alpha_range = range(0,91,30);

alpha = 0.0; theta = 0.0;
i = 0;
	
lines = ["-","--","-.",":"]
linecycler = cycle(lines)
	
for q0x in q0_range:
	cmd = "./main "+str(theta)+" "+str(alpha)+" "+str(q0x)+" 0.5 0.0 gridfile.dat | grep Energy"
	print(cmd)
	returned_value = os.system(cmd)
	#print(returned_value)
	
	hosts = subprocess.check_output(cmd, shell=True)
	error = float(hosts.split(' ')[1])-1
	error = abs(error)
	error2 = float(hosts.split(' ')[3])
	#print(error2)
	thetalabel = r'$\Delta x$ = '+str(q0x)+r', $\varepsilon_1$ = '+"{:.1e}".format(error)+r', $\varepsilon_2$ = '+"{:.1e}".format(error2)
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
	i+=1

plt.xlabel('r', fontsize=20)
plt.ylabel('g(r)', fontsize=20)
plt.ylim(ymin, ymax) 
plt.xlim(xmin, xmax) 
# ax[ip].set_title(thetatitle, fontsize=12)
plt.legend(loc=1)


#plt.show()
#plt.plot(x,x, label='Loaded from file!')
#ip+=1

qmid = [0.0, 0.5, 0.0];
	
maintitle = "Distance distribution g(r) with "+r'$0.0 \leq \Delta x \leq 1.4$'
#plt.savefig('outfile.png')
#plt.suptitle(maintitle);
fileid = qmid[0]*9+qmid[1]*3+qmid[2];
filename = "plot_pt.eps"
plt.savefig(filename, format="eps", bbox_inches='tight')
#os.system("inkscape plot_pt.svg --export-pdf=plot_pt.pdf")
