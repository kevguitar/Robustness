"""
*** plot_comparison.py ***
Author: Kevin Wolz, date: 09/2018

Reads several data files (.rana or .rnum) separated by spaces as a command line
parser, creates a common plot of the result, and stores it inside a file 
"plots/comparison/<timestamp>.png". 

"""

import matplotlib.pyplot as plt
import argparse
from datetime import datetime

# Calculates the numerical derivative of y_list w.r.t. x_list and outputs it in
# another list.
def PlotDerivative(x_list, y_list):
	dydx = []
	for i in range(len(x_list) - 1):
		dy = y_list[i + 1] - y_list[i]
		dx = x_list[i + 1] - x_list[i]
		if dx != 0.0:
			dydx.append(dy/dx)
		else:
			print "Derivative could not be performed at point " + str(x_list[i])
			dydx.append(dydx[i-1])
	dydx.append(dydx[len(x_list) - 2])
			
	return dydx			

	
"""
________________________________________________________________________________

	MAIN
________________________________________________________________________________

"""
# Reads the filepaths as a space-separated parser arguments passed in the shell.
parser = argparse.ArgumentParser(description='Plot several data for \
                                 comparison. Accepts a space-separated list \
                                 containing the paths of the data files to be \
                                 analyzed, e.g. 15-Aug-18_22-21-27_30.rnum')
parsed, filenames = parser.parse_known_args()

for filename in filenames:
	datafile = open('data/' + filename, 'r')
	print filename
	robustness = []
	dist = []
	for line in datafile.readlines():
		robustness.append(float(line.split()[0]))
		dist.append(float(line.split()[1]))
	datafile.close()

	# Plots the distributions.
	bins = 200

	if '20.rnum' in filename:
		text = 'num WPL' #+ '(' + filename.split('/')[-1] + ')'
		plt.hist(robustness, bins, density=True, cumulative=False,
		 		 histtype='step', label=text, color='orange', alpha=0.8, lw=2)
	elif '30.rnum' in filename:
		text = 'num' #+ '(' + filename.split('/')[-1] + ')'
		plt.hist(robustness, bins, density=True, cumulative=False,
				 histtype='step', label=text, color='red', alpha=0.8, lw=2)
	elif '21.rana' in filename:
		text = 'ana' #+ 'WPL (' + filename.split('/')[-1] + ')'
		plt.plot(robustness, dist, label=text, color='cyan', lw=2)		
	else:
		print "Data file name not supported. " +
			  "Supports only suffixes .rnum and .rana."
		quit()
		
plt.title(r'Probability Distribution of $R$', fontsize=16)
plt.xlabel(r'Log-robustness $R$', fontsize=14)
plt.ylabel('Probability Density', fontsize=14)
plt.grid()
plt.legend(loc='upper left', fontsize=12)
plt.autoscale()
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

now = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
plt.savefig('plots/comparison/' + str(now) + '.png',bbox_inches='tight',
        	dpi=300)
plt.show()
plt.clf()
