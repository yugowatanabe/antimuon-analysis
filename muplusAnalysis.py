# Script for analysis of data for muplus experiment, using data generated from
# the ./muplus data collection program
# Author: Yugo Hamada Watanabe, 03/16/2018
# Use --help option for details on how to use script

# Un-comment two lines below if using script over SSH
#import matplotlib
#matplotlib.use('Agg')

import numpy as np							# general
import matplotlib.pyplot as plt				# for plotting
from matplotlib.dates import DateFormatter	# for using dates in plots
from scipy.optimize import curve_fit		# for fitting
from scipy.stats import chi2, chisquare				# for chi-square hypothesis test
import pandas as pd							# for dataframes
import sys, getopt							# for command line arguments
from warnings import simplefilter			# to manage warnings
from os.path import isfile					# to check for existence of files
from datetime import datetime				# to convert unix time to readable time

# Ignore FutureWarning type warnings
simplefilter(action='ignore', category=FutureWarning)

# Default options (global vars) that may be changed with user input
filename = 'none'						# Name of file containing data
datapath = '/home/p180f/Desktop/'		# Directory where data files are assumed to be stored
bgsub = 'none'							# File name for data for background subtraction
K = 2197			                    # Mu+ lifetime to be used in fitting process for magnet ON data
tdc_offset = 80
max_tdc = 8000							# Maximum TDC reading we want to consider (ns)
num_bins = 40							# Number of bins to partition data into
mag = 1;								# Whether we are working with magnet ON or magnet OFF data
save = 0								# For generated graphs, whether to save (1) or display only (0)
bin_size = round(max_tdc/num_bins)

# MAIN ROUTINE #
def main():

	# BAD BINS: bins to be removed from data (add post-hoc)
	badBins = [0,3]

	# Read and verify command line arguments, customizes global vars as necessary
	cl_args()

	# Print confirmation message detailing all of the variables to be used in analysis
	print(confirm_str())

	# Get dataframe including ALL, UP, and DOWN events, partitioned into defined bins
	# eventsdata is array of individual TDC times and their associated unix timestamps
	[df, bins, eventsdata] = readData(datapath, filename, max_tdc, num_bins)

	# Get midpoints of bins
	t = ((bins[1:]+bins[:-1])/2).astype(int)
	t += tdc_offset

	# Add uncertainty in each bin to dataframe
	df = getUnc(df)

	# Background subtraction if requested by user
	bg = 0
	if (bgsub != "none"): [df, bg] = bgsubtraction(df, bgsub)

	# Removing bad bins from data
	[df, bg, t] = removeBadBins(badBins, df, bg, t)

	# Set up plotting environment (6 subplots)
	[fig, ax] = setupPlot(mag);

	# Fitting curve to data and plotting
	[params, perrs] = plotData(ax, eventsdata, t, df, bg)

	# Get chi-square of the fit
	[chisquared, pval] = chisq(t, df, params, perrs)
	print("\nChi-squared value of fit of last plot: %0.02f +/- %0.02f" % (chisquared[0], chisquared[1]))
	print("P-value: %0.04f +/- %0.05f\n" % (pval[0], pval[1]))

	# Saving or showing plot
	if (save):
		print('Saving plot as %s.pdf in current directory\n' % (sys.argv[0])[0:-3])
		plt.savefig('%s.pdf' % (sys.argv[0])[0:-3])
		plt.close()
	else:
		plt.show()


#############
# FUNCTIONS #
#############

# Get and validate command line arguments
def cl_args():

	# Use global variables
	global filename, datapath, bgsub, K, max_tdc, num_bins, mag, save

	# Guarantee that a data file is passed in (which must end with .txt)
	if (len(sys.argv) <= 1):
		print('\nError: missing argument.')
		usage()
		sys.exit(2)
	else:
		filename = str(sys.argv[1])
		if (filename == "--help" or filename == "-h"):
			usage()
			sys.exit(1)

	# Get command line arguments
	try:
		opts, args = getopt.getopt(sys.argv[2:], "D:B:T:C:N:M:sh", \
		["datapath=", "bgsub=", "lifetime=", "maxtdc=", "numbins=", "mag=", "save", "help"])
	except getopt.GetoptError as err:
		print("\n%s" % err)  	# Print error if something bad about arguments
		usage()					# Print usage message
		sys.exit(2)

	# Do appropriate assignments if options passed in, and check for validity
	for o, a in opts:
		if o in ("-D", "--datapath"):
			datapath = a
			if (isfile(datapath+filename) == 0):
				print('\nError: data file %s does not exist. Try using --help.\n' % (datapath+filename))
				sys.exit(2)
		elif o in ("-B", "--bgsub"):
			bgsub = a
			if (isfile(datapath+bgsub) == 0):
				print('\nError: data file %s does not exist. Try using --help.\n' % (datapath+bgsub))
				sys.exit(2)
		elif o in ("-T", "--lifetime"):
			try:
				if (int(a) > 0): K = int(a)
				else: raise Exception('Negative value')
			except:
				print("\nError: mu+ lifetime must be positive integer.\n")
				sys.exit(2)
		elif o in ("-C", "--maxtdc"):
			try:
				if (int(a) > 0): max_tdc = int(a)
				else: raise Exception('Negative value')
			except:
				print("\nError: maxtdc value must be positive integer.\n")
				sys.exit(2)
		elif o in ("-N", "--numbins"):
			try:
				if (int(a) > 0): num_bins = int(a)
				else: raise Exception('Negative value')
			except:
				print("\nError: number of bins must be positive integer.\n")
				sys.exit(2)
		elif o in ("-M", "--mag"):
			try:
				if (int(a) != 0 and int(a) != 1): raise Exception('Not 1 or 0')
				else: mag = int(a)
			except:
				print("\nError: magnet option must be either 1 or 0.\n")
				sys.exit(2)
		elif o in ("-s", "--save"):
			save = 1
		elif o in ("-h", "--help"):
			usage()
			sys.exit(1)

	# Check if data file actually exists
	if (isfile(datapath+filename) == 0):
		print('\nError: data file %s does not exist. Try using --help.\n' % (datapath+filename))
		sys.exit(2)

# Build string giving details of analysis round
def confirm_str():
	print("\nSUMMARY OF OPTIONS FOR ANALYSIS RUN:\n")

	if (mag): magStr = "magnet ON data\n"
	else: magStr = "magnet OFF data\n"

	if (bgsub != "none"): bgsubStr = "-Using background subtraction from file %s\n" % bgsub
	else: bgsubStr = "-Without background subtraction\n"

	if (mag): KStr = "-Assuming theoretical mu+ lifetime of %0.0fns for scaling of data\n" % K
	else: KStr = ""

	helpStr = "(use --help to learn to modify any of these options)\n"

	return("-Reading %s-From file %s\n%s%s-Analyzing data points with lifetimes up to %0.0fns\n-Binning with %0.0f bins\n\n%s"\
	% (magStr, filename, bgsubStr, KStr, max_tdc, num_bins, helpStr))

# Read in data
def readData(datapath, filename, max_tdc, num_bins):
	# Description of data indices in a data file generated by the ./muplus program:
	# 0              1       2      3         4
	# CAMAC_SLOT_ID  NUMBER  NEVENT TDC_VALUE UNIX_TIME

	# Reading in data from text file
	data = np.loadtxt(datapath+filename)							# Load data as numpy array
	data[:,3] *= 20													# Convert TDC reading to ns
	dataf = pd.DataFrame(data[:,[0,3]], columns=('U+D','life_ns')) 	# Turn array into dataframe (relevant data only)
	dataf = dataf.loc[dataf['life_ns'] <= max_tdc]					# Keep only data with TDC readings <= max_tdc
	u_dataf = dataf.loc[dataf['U+D'] == 6]							# Dataframe of only UP events
	d_dataf = dataf.loc[dataf['U+D'] == 7]							# Dataframe of only DOWN events

	# Partitioning data into bins
	bins = np.linspace(0, max_tdc, num_bins)										# Define bins
	df = dataf.groupby(pd.cut(dataf.life_ns, bins)).count()['U+D'].to_frame()		# ALL events
	df['U'] = u_dataf.groupby(pd.cut(u_dataf.life_ns, bins)).count()['U+D'].values	# UP events
	df['D'] = d_dataf.groupby(pd.cut(d_dataf.life_ns, bins)).count()['U+D'].values	# DOWN events

	# Get unix time and tdc value data
	eventsdata = data[:,(3,4)]
	eventsdata = eventsdata[eventsdata[:,0] <= max_tdc]

	return [df, bins, eventsdata]

# Get uncertainty in number of events in each bin
def getUnc(df):
	df['err_U+D'] = np.ceil(df['U+D'].values**(1/2.0)).astype(int)		# Uncertainty given as sqrt(n), where n is # of events in bin
	df['err_U'] = np.ceil(df['U'].values**(1/2.0)).astype(int)
	df['err_D'] = np.ceil(df['D'].values**(1/2.0)).astype(int)

	return df

# Background subtraction functions
def bgsubtraction(df, bgsub):
	data_hrs = float(raw_input("\nINPUT the # of hrs of collected data: "))		# Get hours of data collected as input from user
	bg_hrs = float(raw_input("INPUT the # of hrs of background data: "))		# Get hours of background data as input
	mult = round(data_hrs/bg_hrs, 3)
	[bg, _, _] = readData(datapath, bgsub, max_tdc, num_bins)					# Read in background data
	bg.loc[:,:] *= mult															# Scale the background subtraction data to real data
	bg = getUnc(bg)																# Get uncertainty in background data

	df.loc[:,('U+D','U','D')] -= np.around(bg.loc[:,('U+D','U','D')])			# Subtract scaled background events

	# Can't have < 0 events (avoid 0 as well, causes issue when taking log)
	df.loc[df['U'] < 1, 'U'] = 1
	df.loc[df['D'] < 1, 'D'] = 1
	df['U+D'] = df['U'].values+df['D'].values

	return [df, bg]

# Remove specific bins from data
def removeBadBins(badBins, df, bg, t):
	df = df.drop(df.index[badBins])
	t = np.delete(t, badBins)

	if (isinstance(bg, pd.DataFrame)):
		bg = bg.drop(bg.index[badBins])

	return [df, bg, t]

# Build the plot environment
def setupPlot(mag):
	# Set some matplotlib options (for plots)
	plt.rc('axes', titlesize=15)     # fontsize of the axes title
	plt.rc('xtick', labelsize=8)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=8)    # fontsize of the tick labels

	fig = plt.figure()
	fig.subplots_adjust(left=0.1, bottom=0.05, right=0.97, top=0.9, wspace=0.25, hspace=None)

	if (mag):
		fig.suptitle('Magnet ON Data', fontsize=16)
	else:
		fig.suptitle('Magnet OFF Data', fontsize=16)

	# All event points plotted over time
	ax1 = fig.add_subplot(3,2,1)
	xfmt = DateFormatter('%m/%d')
	ax1.xaxis.set_major_formatter(xfmt)
	ax1.set_ylabel('TDC time (ns)')

	# Target out data
	ax2 = fig.add_subplot(3,2,2)
	ax2.set_ylabel('Targ out events')
	ax2.grid(linestyle='dashed')

	ax3 = fig.add_subplot(3,2,3)
	# ax3.set_ylim(-100, 800)
	ax4 = fig.add_subplot(3,2,4)
	# ax4.set_ylim(-100, 800)
	ax5 = fig.add_subplot(3,2,5)
	# ax5.set_ylim(-100, 1400)
	ax6 = fig.add_subplot(3,2,6)
	# ax6.set_ylim(-500, 650)

	ax = [ax1, ax2, ax3, ax4, ax5, ax6]

	return [fig, ax]

# Plot magnet ON data
def plotData(ax, eventsdata, t, df, bg):

	global K, bin_size, max_tdc
	X = np.linspace(0, max_tdc, 500, endpoint=True)			# x-axis points to plot on

	# Plot events over time
	date = [0]*len(eventsdata[:,1])
	date = [datetime.fromtimestamp(ts) for ts in eventsdata[:,1]]
	ax[0].plot(date,eventsdata[:,0],'k.',alpha=0.3)

	# Plot target out data if given
	if (isinstance(bg, pd.DataFrame)):
		ax[1].errorbar(t, bg['U'].values, xerr=bin_size/2, \
		fmt='r.', ms=3, ecolor="Red", linestyle="None", capsize=0, label='U events')
		ax[1].errorbar(t, bg['D'].values, xerr=bin_size/2, \
		fmt='b.', ms=3, ecolor="Blue", linestyle="None", capsize=0, label='D events')
		ax[1].legend(fontsize='x-small', loc='upper right', bbox_to_anchor=(1, 1), fancybox=False, shadow=False, numpoints=1)
	else:
		ax[1].text(0.2, 0.5, '(No target out data)')

	# Fit simple decay formulas to ALL, U, and D data
	[A, T, Aerr, Terr] = fitDecay(df, t)

	# Plot UP events
	if (isinstance(bg, pd.DataFrame)):
		ax[2].set_ylabel('U minus bkgrnd')
	else:
		ax[2].set_ylabel('U events')
	ax[2].grid(linestyle='dashed')
	ax[2].errorbar(t, df['U'].values, xerr=bin_size/2, yerr=df['err_U'].values, fmt='', ecolor="Red", linestyle="None", capsize=0)
	ufit_y = np.exp(-X/T[1])*A[1]
	ax[2].plot(X, ufit_y, 'k--', label='$y =  A_U e^{-t/T_U}$')
	ax[2].legend(fontsize='x-small', loc='lower left', bbox_to_anchor=(0, 0), fancybox=False, shadow=False)
	ax2str = "Entries:    %0.0f\nRMS:         %0.1f\n$A_U$:   %0.0f$\pm$%0.1f\n$T_U$:  %0.0f$\pm$%0.0f" %\
				(sum(df['U'].values), rms(df['U'].values, np.exp(-t/T[1])*A[1]), A[1], Aerr[1], T[1], Terr[1])
	ax[2].text(0.98, 0.96, ax2str,
        bbox={'facecolor':'white'}, style='oblique', horizontalalignment='right', verticalalignment='top', fontsize=7, transform = ax[2].transAxes)

	# Plot DOWN events
	if (isinstance(bg, pd.DataFrame)):
		ax[3].set_ylabel('D minus bkgrnd')
	else:
		ax[3].set_ylabel('D events')
	ax[3].grid(linestyle='dashed')
	ax[3].errorbar(t, df['D'].values, xerr=bin_size/2, yerr=df['err_D'].values, fmt='', ecolor="Blue", linestyle="None", capsize=0)
	dfit_y = np.exp(-X/T[2])*A[2]
	ax[3].plot(X, dfit_y, 'k--', label='$y =  A_D e^{-t/T_D}$')
	ax[3].legend(fontsize='x-small', loc='lower left', bbox_to_anchor=(0, 0), fancybox=False, shadow=False)
	ax3str = "Entries:    %0.0f\nRMS:         %0.1f\n$A_D$:   %0.0f$\pm$%0.1f\n$T_D$:  %0.0f$\pm$%0.0f" %\
				(sum(df['D'].values), rms(df['D'].values, np.exp(-t/T[2])*A[2]), A[2], Aerr[2], T[2], Terr[2])
	ax[3].text(0.98, 0.96, ax3str,
        bbox={'facecolor':'white'}, style='oblique', horizontalalignment='right', verticalalignment='top', fontsize=7, transform = ax[3].transAxes)

	# Plot U+D data
	UmD = df['U'].values - df['D'].values								# Get U-D values
	err_UmD = (df['err_U'].values + df['err_D'].values)					# Uncertainty associated with each U-D bin
	numUD = sum(df['U+D'].values)
	numUmD = sum(UmD)
	if (isinstance(bg, pd.DataFrame)):
		ax[4].set_ylabel('U+D minus bkgrnd')
	else:
		ax[4].set_ylabel('U+D events')
	ax[4].grid(linestyle='dashed')
	ax[4].errorbar(t, df['U+D'].values, xerr=bin_size/2, yerr=df['err_U+D'].values, fmt='', ecolor="Black", linestyle="None", capsize=0)
	fit_y = np.exp(-X/T[0])*A[0]
	ax[4].plot(X, fit_y, 'k--', label='$y = A e^{-t/T}$')
	ax[4].legend(fontsize='x-small', loc='lower left', bbox_to_anchor=(0, 0), fancybox=False, shadow=False)
	ax4str = "Entries:    %0.0f\nRMS:         %0.1f\n$A$:    %0.0f$\pm$%0.1f\n$T$:    %0.0f$\pm$%0.0f" %\
				(numUD, rms(df['U+D'].values, np.exp(-t/T[0])*A[0]), A[0], Aerr[0], T[0], Terr[0])
	ax[4].text(0.98, 0.96, ax4str,
        bbox={'facecolor':'white'}, style='oblique', horizontalalignment='right', verticalalignment='top', fontsize=7, transform = ax[4].transAxes)

	# Stop plots here if magnet off data
	if (mag == 0):
		ax[5].text(0.2, 0.5, '(Purposefully empty)')
		return [(A[0], T[0]), (Aerr[0], Terr[0])]

	# Plot scaled U-D data
	[C, B, w, err, linUmD, linerr_UmD] = fitDecayOsc(UmD, err_UmD, t)	# Fit curve to data, get curve parameters

	ax[5].set_ylabel('$e^{t/%0.0f}(U-D)$' % K)
	ax[5].grid(linestyle='dashed')
	linfit_y = model_func(X, C, B, w)
	ax[5].errorbar(t, linUmD, xerr=bin_size/2, yerr=linerr_UmD, fmt='', ecolor="Black", linestyle="None", capsize=0)
	ax[5].plot(X, linfit_y, 'g-', label='$y = B\cos(\omega t) + D$')
	ax[5].legend(fontsize='x-small', loc='lower left', bbox_to_anchor=(0, 0), fancybox=False, shadow=False)
	ax5str = "RMS:           %0.1f\n$D$:          %0.1f$\pm$%0.1f\n$B$:          %0.1f$\pm$%0.1f\n$\omega$:      %0.1E$\pm$%0.1E" %\
				(rms(model_func(t, C, B, w), linUmD), C, err[0], B, err[1], w, err[2])
	ax[5].text(0.02, 0.96, ax5str,
        bbox={'facecolor':'white'}, style='oblique', horizontalalignment='left', verticalalignment='top', fontsize=7, transform = ax[5].transAxes)

	return [(C, B, w), err]

# Get chi-squared value for fit to fitted distribution
def chisq(t, df, params, perrs):

	if (mag):
		[C, B, w] = params
		[Cerr, Berr, werr] = perrs
		O = np.exp(t/K)*np.array(df['U'].values - df['D'].values)	# Observed scaled U-D data
		E = model_func(t, C, B, w)									# Expected scaled U-D data
		Oerr = np.exp(t/K)*(df['err_U'].values + df['err_D'].values)
		Eerr = np.sqrt(Cerr**2 + ( (B*np.cos(w*t)*Berr/B)**2) + ((B*np.cos(w*t)*t*np.tan(w*t)*werr)**2) )	# Error in expected data due to parameter uncertainties
		dof = len(O) - 4
		chisquared = sum((O-E)**2/Oerr**2)
		errs = 1/Oerr**2 * 2*abs(O-E) * Eerr					# Contribution to error of each uncertainty in each expected data point
		chisquarederr = np.sqrt(sum(errs**2))
		pvalmin = 1. - chi2.cdf(chisquared+chisquarederr, df=dof)
		pvalmax = 1. - chi2.cdf(chisquared-chisquarederr, df=dof)
		pval = pvalmin + (pvalmax-pvalmin)/2.
		pvalerr = pvalmax-pval
	else:
		[A, T] = params
		[Aerr, Terr] = perrs
		O = np.log(df['U+D'].values)
		E = -t/T + np.log(A)
		Oerr = abs(df['err_U+D'].values/df['U+D'].values)
		Eerr = abs(t/T**2/A)*Terr*Aerr
		dof = len(O) - 3
		[chisquared, pval] = chisquare(O, E, ddof=dof)
		errs = 2*chisquared*Eerr
		errs = 1/Oerr**2 * 2*abs(O-E) * Eerr
		pvalerr = 0.01
		chisquarederr = np.sqrt(sum(errs**2))


	return [(chisquared, chisquarederr), (pval, pvalerr)]

# Fit equation to magnet OFF data
def fitDecay(df, t):
	linUD = np.log(df['U+D'].values)
	linU = np.log(df['U'].values)
	linD = np.log(df['D'].values)

	R = np.array([0., 0., 0.])		# Slope of log of data for ALL, UP, and DOWN
	C = np.array([0., 0., 0.])		# y-intercept of log of data for ALL, UP, and DOWN
	v = np.zeros((3,2,2))			# covariance matrix of each parameter

	[(R[0], C[0]), v[0]]  = np.polyfit(t, linUD, 1, cov=True)	# Fit data to first-order poly
	[(R[1], C[1]), v[1]] = np.polyfit(t, linU, 1, cov=True)
	[(R[2], C[2]), v[2]] = np.polyfit(t, linD, 1, cov=True)

	Rstdev = v[:,0,0]**(0.5)		# Std dev of parameters
	Cstdev = v[:,1,1]**(0.5)

	T = -1./R						# Get mean lifetime
	A = np.exp(C)					# Remove scaling constant
	Tstdev = Rstdev/R**2			# Propagation of error
	Astdev = np.exp(Cstdev)

	return [A, T, Astdev, Tstdev]

# Fit equation to magnet ON data (decay with oscillation)
def fitDecayOsc(UmD, err_UmD, t):
	global K
	linUmD = np.exp(t/K)*UmD				# Remove decay component to data by scaling
	linerr_UmD = np.exp(t/K)*err_UmD		# Scale errors as well
	p0 = [0, UmD[0], 2*np.pi/K]				# Our rough guess for parameters
	[C, B, w], v = curve_fit(model_func, t, linUmD, p0, sigma=linerr_UmD)	# Fit the data

	err = np.zeros(3)						# Uncertainty (stdev) of each parameter
	err[0] = v[0,0]**(0.5)
	err[1] = v[1,1]**(0.5)
	err[2] = v[2,2]**(0.5)

	return [C, B, w, err, linUmD, linerr_UmD]

# Model function we expect our data to follow
def model_func(t, A, B, w):
	return A + B*np.cos(w*t)

# Calculate root mean square
def rms(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

# Message detailing usage of script
def usage():
	print('\nUsage: python %s [name of data file] [OPTIONS]\n' % sys.argv[0])
	print('OPTIONS:\n')
	print('-D, --datapath=FILENAME     path to directory where data is stored. If not provided, assumes \n\t\t\t    default: %s\n' % datapath)
	print('-B, --bgsub=FILENAME        data for background subtraction. If not provided, no background \n\t\t\t    subtraction. Note: background data must be in same directory as data.\n')
	print('-T, --lifetime=#            mu+ mean lifetime to use in fitting (used only for magnet ON data). \n\t\t\t    If not provided, assumes T=%sns.\n' % K)
	print('-C, --maxtdc=#              upperbound of TDC count (in ns) to be considered in analysis. If \n\t\t\t    not provided, assume maxtdc=%sns.\n' % max_tdc)
	print('-N, --numbins=#             Number of bins to divide data into. If not provided, assumes \n\t\t\t    numbins=%s.\n' % num_bins)
	print('-M, --mag=[0 or 1]          Whether given data is for magnet ON or OFF. If not provided, assumes \n\t\t\t    mag=%s (ON).\n' % mag)
	print('-s, --save                  save the generated plots instead of displaying them.\n')


if __name__ == '__main__':
    main()
