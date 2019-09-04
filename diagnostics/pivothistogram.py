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
	parser.add_argument("name", type=str, help="filename")
	parser.add_argument("bound", type=str, help="+- bound")

	args = parser.parse_args()
	filename = args.name
	bound = float(args.bound)

	df = pd.read_csv(filename, sep=" ", header=None)
	ar = df.values
	ar = np.transpose(ar)

	data = ar[0]
	w = ar[1]
	
	R = data.size
	bincount = np.int(np.sqrt(float(R)))
	
	plt.figure(num=1,figsize=(12,6),dpi=100,facecolor='white')

	plt.hist(data, bins=bincount, range = (-bound, bound), weights = w, facecolor='b', alpha=0.75)#, log=True)

	plt.show()
	
		
		
		
main()