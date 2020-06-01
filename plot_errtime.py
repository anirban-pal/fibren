#!/usr/bin/python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
from math import log

plt.rcParams.update({'font.size': 14})

filename = sys.argv[1]

file1 = open(filename, 'r') 
Lines = file1.readlines() 

for line in Lines: 
    line_curr = line.strip()
    #print(line_curr)

fib_ref=Lines[1:2]
fib_adapt=Lines[3:12]
fib_discrete=Lines[13:23]

#print(fib_ref)
#print(fib_adapt)
#print(fib_discrete)

Uref = fib_ref[0].strip().split()[0];
#print(Uref)

U_adapt = np.array([])
t_adapt = np.array([])
t_adapt2 = np.array([])

fac = 1.0/log(10);

for line in fib_adapt:
		U = line.strip().split()[0]
		t = line.strip().split()[2]
		ts = line.strip().split()[3]
		U_adapt = np.append(U_adapt,[abs(float(U)-float(Uref))])
		t_adapt = np.append(t_adapt,[float(t)])
		t_adapt2 = np.append(t_adapt2,[float(ts)])

U_discrete = np.array([])
t_discrete = np.array([])
t_discrete2 = np.array([])

for line in fib_discrete:
		U = line.strip().split()[0]
		t = line.strip().split()[2]
		ts = line.strip().split()[3]
		U_discrete = np.append(U_discrete,[abs(float(U)-float(Uref))])
		t_discrete = np.append(t_discrete,[float(t)])
		t_discrete2 = np.append(t_discrete2,[float(ts)])


#######################

filename = sys.argv[2]

file1 = open(filename, 'r') 
Lines = file1.readlines() 

for line in Lines: 
    line_curr = line.strip()
    #print(line_curr)

fib2_ref=Lines[1:2]
fib2_adapt=Lines[3:12]
fib2_discrete=Lines[13:23]

#print(fib_ref)
#print(fib_adapt)
#print(fib_discrete)

Uref = fib2_ref[0].strip().split()[0];
#print(Uref)

U2_adapt = np.array([])
t2_adapt = np.array([])
t2_adapt2 = np.array([])

for line in fib2_adapt:
		U = line.strip().split()[0]
		t = line.strip().split()[2]
		ts = line.strip().split()[3]
		U2_adapt = np.append(U2_adapt,[abs(float(U)-float(Uref))])
		t2_adapt = np.append(t2_adapt,[float(t)])
		t2_adapt2 = np.append(t2_adapt2,[float(ts)])

U2_discrete = np.array([])
t2_discrete = np.array([])
t2_discrete2 = np.array([])

for line in fib2_discrete:
		U = line.strip().split()[0]
		t = line.strip().split()[2]
		ts = line.strip().split()[3]
		U2_discrete = np.append(U2_discrete,[abs(float(U)-float(Uref))])
		t2_discrete = np.append(t2_discrete,[float(t)])
		t2_discrete2 = np.append(t2_discrete2,[float(ts)])



ef = r'$\epsilon_f=|(U/U_{ref})-1|$'
#######################
fig = plt.figure()

plt.xlabel("Error fraction "+ef, fontsize=20)
plt.ylabel("Computation time (sec)", fontsize=20)
plt.xscale("log", nonposx='clip')
plt.yscale("log", nonposy='clip')
plt.xlim([10**(-19),10**9])
plt.ylim([10**(-6),10**3])

plt.errorbar(U_adapt, t_adapt, fmt = 'rs', markersize=5, markerfacecolor='none', yerr=t_adapt2, label=r'$\theta = 0 \degree$ (fibren)')
plt.errorbar(U_discrete, t_discrete, fmt = 'bo', markersize=5, markerfacecolor='none', yerr=t_discrete2, label=r'$\theta = 0 \degree$ (discrete)')
plt.errorbar(U2_adapt, t2_adapt, fmt = 'rs', markersize=5, yerr=t2_adapt2, label=r'$\theta = 30 \degree$ (fibren)')
plt.errorbar(U2_discrete, t2_discrete, fmt = 'bo', markersize=5, yerr=t2_discrete2, label=r'$\theta = 30 \degree$ (discrete)')
plt.grid(linestyle='-', linewidth='0.1', color='black')
plt.legend(loc=3)

plt.savefig("plot_discrete.eps", format="eps", bbox_inches='tight')
#cmd2 = "inkscape outfile.svg --export-pdf=errtime.pdf"
#os.system(cmd2)
