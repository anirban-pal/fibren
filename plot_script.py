#!/usr/bin/python

import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
import csv
import os
import subprocess
from itertools import cycle

plt.rcParams.update({'font.size': 12})
#cmd = "./main 0.00 0.00 1.0 1.0 0 gridfile.dat"
#returned_value = os.system(cmd)  # returns the exit code in unix
returned_value = os.system("seq 0.0 0.001 3.0 > gridfile.dat")

ymin = -0.1; xmax = 3.0;

#qmid = [1.1, 0.0, 0.0]; ymax = 2.5; xmax = 2.5; xmin = 0.0;
#qmid = [1.1, 1.1, 0.0]; ymax = 4.0; xmax = 2.5; xmin = 0.5;
#qmid = [0.0, 1.1, 0.0]; ymax = 25.0; xmax = 2.0; xmin = 0.5;

#qmid = [1.1, 0.0, 1.1]; ymax = 3.5; xmax = 3.0; xmin = 0.5;
qmid = [1.1, 1.1, 1.1]; ymax = 5.0; xmax = 3.0; xmin = 1.0;

tmax = 180
tmin = 0
interval = 30
theta_range = range(tmin,tmax,interval)
colors = pl.cm.Dark2(np.linspace(0.5,1,(tmax-tmin)/interval))

fig, ax = plt.subplots(1, 4, figsize=(15,5))
ip = 0;

alpha_range = range(0,91,30);

for alpha in alpha_range:
	i = 0;
		
	lines = ["-","--","-.",":"]
	linecycler = cycle(lines)
		
	for theta in theta_range:
		cmd = "./main "+str(theta)+" "+str(alpha)+" "+str(qmid[0])+" "+str(qmid[1])+" "+str(qmid[2])+" gridfile.dat | grep Energy"
		print(cmd)
		#returned_value = os.system(cmd)
		#print(returned_value)
		
		hosts = subprocess.check_output(cmd, shell=True)
		error = float(hosts.split(' ')[1])-1
		error = abs(error)
		print(error)
		thetalabel = r'$\theta$ = '+str(theta)+r'$\degree$ $\varepsilon$ = '+"{:.1e}".format(error)
		thetatitle = r'$\alpha$ = '+str(alpha)+r'$\degree$'
		
		x = []
		y = []
		
		with open('outfile.gr','r') as csvfile:
		    next(csvfile); 
		    next(csvfile);
		    plots = csv.reader(csvfile, delimiter=' ')
		    for row in plots:
		        x.append(float(row[0]))
		        y.append(float(row[1]))
		
		ax[ip].plot(x,y, next(linecycler), label=thetalabel, color=colors[i])
		i+=1
	
	ax[ip].set_xlabel('r', fontsize=20)
	#ax[ip].set_ylabel('g(r)', fontsize=18)
	ax[ip].set_ylim(ymin, ymax) 
	ax[ip].set_xlim(xmin, xmax) 
	ax[ip].set_title(thetatitle, fontsize=12)
	ax[ip].legend(loc=9)
	
	
	#plt.show()
	#plt.plot(x,x, label='Loaded from file!')
	ip+=1
	
maintitle = "Distance distribution g(r) with "+r"$q_m$ = "+str(qmid)
#plt.savefig('outfile.png')
plt.suptitle(maintitle);
fileid = (qmid[0]*9+qmid[1]*3+qmid[2])*10;
filename = "outfile_"+str(fileid)+".pdf"
plt.savefig("outfile.svg", format="svg", bbox_inches='tight')
cmd2 = "inkscape outfile.svg --export-pdf="+filename
os.system(cmd2)
