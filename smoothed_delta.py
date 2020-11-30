#!/usr/bin/env python

"""
outputs smoothed delta bedgraph, t input == delta bedgraph

to visualize smoothing function operated on delta files
when picking all local maxima

"""

import argparse
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from numpy import *
import numpy
from scipy.signal import argrelextrema
import operator

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def read_coverage(wig_fyle):
    #create list of coverage values, position dependent
    with open(wig_fyle,"r") as inp:
        firstline = inp.readline()
        output_list = [[float(i) for ind,i in enumerate(line.strip().split("\t")) if ind == 3] for line in inp.readlines() if str(line[0]) != '#']
        flat_list = [item for sublist in output_list for item in sublist]
    return flat_list

def smoother(delta_list,window_len):
    return savitzky_golay(numpy.asarray(delta_list),window_len,3,deriv=0, rate=1).tolist()

def writer(wig_fyle,smoothed_list):
    #creates delta wig file, based off values in wig file
    new_name = 'smoothed_'+str(wig_fyle)
    counter = 0
    with open(wig_fyle,'r') as inp:
        firstline = inp.readline()
        with open(new_name,"w") as outp:
            outp.write(firstline)
            while counter < len(smoothed_list)-1:
                for line in inp:
                    if str(line[0]) == '#':
                        outp.write(line)
                    else:
                        line = line.strip().split("\t")
                        outp.write("\t".join(line[0:3]))
                        outp.write("\t")
                        outp.write(str(round(smoothed_list[counter],1)))
                        outp.write("\n")
                        counter += 1


def main():
    parser = argparse.ArgumentParser(description='calculates the delta value at each position [ x upstream - x downstream].')
    parser.add_argument('-window_len',type=int,default=11,help='length of window used for smoothing, note, NOT NUCLEOTIDES, based on x lines of wig file, odd integer, default = 11')
    parser.add_argument('-file', default = None, help='Specific <.wig>s to use', nargs='+',dest='bedgraphs')
    parser.add_argument('-batchdir',action="store_true",default=False,help = 'Use all delta <.bedgraphs> in the directory')
    args = parser.parse_args()

    if args.bedgraphs != None:
        for fyle in args.bedgraphs:
            sub_data = read_coverage(fyle)
            sub_smooothed = smoother(sub_data,args.window_len)
            writer(fyle,sub_smooothed)

    if args.batchdir == True:
        for fyle in sorted(glob.glob('delta_*.bedgraph')):
            sub_data = read_coverage(fyle)
            sub_smooothed = smoother(sub_data,args.window_len)
            writer(fyle,sub_smooothed)


if __name__ == '__main__':
    main()
