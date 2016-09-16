import numpy as np
import scipy.stats as stat
import h5py
import matplotlib.pyplot as plt

def compare(filenameA, filenameB):
    fa = h5py.File(filenameA)
    fb = h5py.File(filenameB)

    ha = fa['halos']
    hb = fb['halos']

    pa = ha['percentile']
    pb = hb['percentile']

    print(pa[:10])
    print(pb[:10])

# compute 2d histogram
    xmin = 0.
    xmax = 1.
    f, xedges, yedges = np.histogram2d(pa, pb, bins=100, range=[[xmin,xmax],[xmin,xmax]],normed=True)
    f = np.log10(f)
    f[f == -np.inf] = -10
    xmid = np.diff(xedges) * 0.5 + xedges[:-1]
    ymid = np.diff(yedges) * 0.5 + yedges[:-1]
    xx, yy = np.meshgrid(xmid,ymid)

    fig = plt.figure()
    ax = fig.gca()
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(xmin,xmax)
    ax.set_xlabel(filenameA)
    ax.set_ylabel(filenameB)
    ax.set_title('correlation between env. density percentiles\nin mass bins of varying width')

    im = ax.pcolormesh(xx,yy,f)
    plt.colorbar(im)

#    plt.scatter(pa,pb)

    plt.savefig('mass_bin_width_comparison.png')

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fileA')
parser.add_argument('fileB')
args = parser.parse_args()

compare(args.fileA, args.fileB)
