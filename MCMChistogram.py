import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.colors as colors

"""
def histogram2d(x,y,nx,ny):
	hist,xedges,yedges= np.histogram2d(x,y,[nx,ny],normed=False)
	xe                = 0.5*(xedges[0:xedges.size-1]+xedges[1:xedges.size])
	ye                = 0.5*(yedges[0:yedges.size-1]+yedges[1:yedges.size])
	nbinx             = xe.size
	nbiny             = ye.size
	tothist           = np.sum(hist.astype(float))
	hist              = hist.astype(float)/tothist
	return hist,xe,ye
"""
	
def main():

	parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
	parser.add_argument("i0", type=int, help="index of first model param")
	parser.add_argument("i1", type=int, help="index of second model param")
	#parser.add_argument("a0fit", type=float, help="best fit param 1")
	#parser.add_argument("a1fit", type=float, help="best fit param 2")

	args = parser.parse_args()
	i0 = args.i0
	i1 = args.i1
	#a0fit = args.a0fit
	#a1fit = args.a1fit


	df = pd.read_csv('TRKMCMC_data.txt', sep=" ", header=None)
	ar = df.values
	ar = np.transpose(ar)

	a0 = ar[i0]
	a1 = ar[i1]
	
	R = a0.size
	bincount = np.int(np.sqrt(float(R)))
	hist, xe, ye = np.histogram2d(a0, a1, [bincount, bincount], normed=False)
	
	plt.figure(num=1,figsize=(12,6),dpi=100,facecolor='white')
	plt.subplot(223)
	plt.imshow(hist, cmap="binary", norm=colors.LogNorm(), interpolation='nearest', origin='low', aspect="auto", extent=[xe[0], xe[-1], ye[0], ye[-1]])
	#plt.plot(a0fit, a1fit, color="red")
	plt.xlabel("a0")
	plt.ylabel("a1")
	plt.title("Parameter Space plot")
	
	plt.subplot(221)
	n, bins, patches = plt.hist(a0, bincount, density=True, facecolor='b', alpha=0.75)
	plt.xlabel("a0")
	
	plt.subplot(222)
	n, bins, patches = plt.hist(a1, bincount, density=True, facecolor='r', alpha=0.75)
	plt.xlabel("a1")
	
	plt.show()
	
		
		
		
main()