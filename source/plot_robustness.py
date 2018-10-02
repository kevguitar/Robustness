"""
*** plot_robustness.py ***
Author: Kevin Wolz, date: 09/2018

Reads a data file (.rana or .rnum) as a command line parser, creates a plot of
the result, and stores it inside files "plots/<data_filename>.png" and 
"plots/<data_filename>.svg". 

"""

import matplotlib.pyplot as plt
import argparse

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

# Reads the filepath as a parser argument passed in the shell.
parser = argparse.ArgumentParser(description='Plot the data.')
parser.add_argument('filepath', type=str, help='The path of the data file to \
					be analyzed, e.g. ./data/15-Aug-18_22-21-27_30.rnum')
args = parser.parse_args()
filename = args.filepath


# Reads data from the indicated file.
datafile = open(filename, 'r')
robustness = []
dist = []
cumulative = []
summed = 0.0
count = 0
for line in datafile.readlines():
	count += 1
	term = float(line.split()[1])
	summed += term
	robustness.append(float(line.split()[0]))
	dist.append(term)
	cumulative.append(summed)
datafile.close()


# Plots the distributions.
bins = 150
text = filename.split('/')[-1]
plot = plt.axes()
subplot = plt.axes([0.15, 0.45, .3, .3])
R_lo, R_hi = 0.0, 0.0
y_lo, y_hi = 0.0, 0.0

if '20.rnum' in filename:
	n_plot, b_plot, p_plot = plot.hist(robustness, bins, density=True, 
	                                   cumulative=False, histtype='bar', 
	                                   ec='black', fc='orange', label=text, 
	                                   lw=0.5)
	plot.set_title(r'Numerical Distribution of $R$ (Weak Prior Limit)',
	 			   fontsize=16)
	n_splot, b_splot, p_splot = subplot.hist(robustness, bins, density=True, 
	                                         cumulative=True, histtype='step',
	                                         ec='orange', label=text, lw=1)
	for i in range(len(n_splot)):
		if (n_splot[i] > 0.05) and (R_lo == 0.0):
			R_lo = b_splot[i]
		if (n_splot[i] >= 0.95) and (R_hi == 0.0):
			R_hi = b_splot[i]
		if (R_lo == 0.0):
			plt.setp(p_plot[i], 'facecolor', 'blue')
		if (R_hi != 0.0):
			plt.setp(p_plot[i], 'facecolor', 'blue')

elif '30.rnum' in filename:
	n_plot, b_plot, p_plot = plot.hist(robustness, bins, density=True, 
	                                   cumulative=False, histtype='bar', 
	                                   ec='black', fc='red', label=text,
	         			   			   lw=0.5)
	plot.set_title(r'Numerical Distribution of $R$', fontsize=16)
	n_splot, b_splot, p_splot = subplot.hist(robustness, bins, density=True, 
	                                         cumulative=True, histtype='step', 
	                                         ec='red', label=text, lw=1)
	
	for i in range(len(n_splot)):
		if (n_splot[i] > 0.05) and (R_lo == 0.0):
			R_lo = b_splot[i]
		if (n_splot[i] >= 0.95) and (R_hi == 0.0):
			R_hi = b_splot[i]
		if (R_lo == 0.0):
			plt.setp(p_plot[i], 'facecolor', 'cyan')
		if (R_hi != 0.0):
			plt.setp(p_plot[i], 'facecolor', 'cyan')

elif '21.rana' in filename:
	# The data input for mode 21 is the probability distribution.
	cumulative = [cum/summed for cum in cumulative]
	plot.plot(robustness, dist, color='cyan', label=text, lw=2)
	plot.set_title(r'Analytical Distribution of $R$ (Weak Prior Limit)', 
	               fontsize=16)
	subplot.plot(robustness, cumulative, color='cyan', label=text, lw=1)

	for i in range(len(cumulative)):
		if (cumulative[i] > 0.05) and (R_lo == 0.0):
			R_lo = robustness[i]
			y_lo = dist[i]
		if (cumulative[i] >= 0.95) and (R_hi == 0.0):
			R_hi = robustness[i]
			y_hi = dist[i]
	plot.vlines(R_lo, 0.0, y_lo, colors='red')
	plot.vlines(R_hi, 0.0, y_hi, colors='red')

else:
	print "Data file name not supported. Supports only suffixes .rnum and .rana."
	quit()

subplot.text(0.1, 0.3, 
             r'$R_{hi} = %1.3f$'
             '\n'
             r'$R_{lo} = %1.3f$' %(R_hi, R_lo), 
             transform=subplot.transAxes,
             fontsize=12)
subplot.axhline(y=0.05, c='grey')
subplot.axhline(y=0.95, c='grey')
#subplot.set_ylim(0,1)
#subplot.set_xticks([])
subplot.set_yticks([])

plot.set_xlabel(r'Log-robustness $R$', fontsize=14)
plot.set_ylabel('Probability Density', fontsize=14)
plot.grid()
plot.text(0.0, 1.2, text, transform=subplot.transAxes, fontsize=12, 
          bbox={'facecolor':'white', 'edgecolor':'none'})
plot.autoscale()


#Saves plot in the local filesystem.	
plotpath = 'plots/' \
			+ filename.split('/')[-1]	
plt.savefig(plotpath + '.png', format='png')
plt.savefig(plotpath + '.svg', format='svg')	
plt.clf()

