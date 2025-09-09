import sys
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
import numpy as np
import scipy.stats as sci_stats
import scipy.optimize as sci_optimize
import scipy.interpolate as sci_int
from astropy.io import fits
import json
import warnings
import pygtc
import re

# Use ordered dictionaries for the observations
from collections import OrderedDict


"""
List of exceptions 
"""
class MissingFileException(Exception):
    pass

class MissingKeywordException(Exception):
    pass

class OutOfRangeException(Exception):
    pass

def difference_utc2tdb(jd):
    ### jd        date                sec        leap diff
    #  2441317.5 1972-01-01T00:00:00 2272060800 10   42.184
    #  2441499.5 1972-07-01T00:00:00 2287785600 11   43.184
    #  2441683.5 1973-01-01T00:00:00 2303683200 12   44.184
    #  2442048.5 1974-01-01T00:00:00 2335219200 13   45.184
    #  2442413.5 1975-01-01T00:00:00 2366755200 14   46.184
    #  2442778.5 1976-01-01T00:00:00 2398291200 15   47.184
    #  2443144.5 1977-01-01T00:00:00 2429913600 16   48.184
    #  2443509.5 1978-01-01T00:00:00 2461449600 17   49.184
    #  2443874.5 1979-01-01T00:00:00 2492985600 18   50.184
    #  2444239.5 1980-01-01T00:00:00 2524521600 19   51.184
    #  2444786.5 1981-07-01T00:00:00 2571782400 20   52.184
    #  2445151.5 1982-07-01T00:00:00 2603318400 21   53.184
    #  2445516.5 1983-07-01T00:00:00 2634854400 22   54.184
    #  2446247.5 1985-07-01T00:00:00 2698012800 23   55.184
    #  2447161.5 1988-01-01T00:00:00 2776982400 24   56.184
    #  2447892.5 1990-01-01T00:00:00 2840140800 25   57.184
    #  2448257.5 1991-01-01T00:00:00 2871676800 26   58.184
    #  2448804.5 1992-07-01T00:00:00 2918937600 27   59.184
    #  2449169.5 1993-07-01T00:00:00 2950473600 28   60.184
    #  2449534.5 1994-07-01T00:00:00 2982009600 29   61.184
    #  2450083.5 1996-01-01T00:00:00 3029443200 30   62.184
    #  2450630.5 1997-07-01T00:00:00 3076704000 31   63.184
    #  2451179.5 1999-01-01T00:00:00 3124137600 32   64.184
    #  2453736.5 2006-01-01T00:00:00 3345062400 33   65.184
    #  2454832.5 2009-01-01T00:00:00 3439756800 34   66.184
    #  2456109.5 2012-07-01T00:00:00 3550089600 35   67.184
    #  2457204.5 2015-07-01T00:00:00 3644697600 36   68.184
    #  2457754.5 2017-01-01T00:00:00 3692217600 37   69.184

    jd_table = np.asarray([2441317.5, 2441499.5, 2441683.5, 2442048.5, 2442413.5, 2442778.5, 2443144.5, 2443509.5, 2443874.5, 2444239.5,
                           2444786.5, 2445151.5, 2445516.5, 2446247.5, 2447161.5, 2447892.5, 2448257.5, 2448804.5, 2449169.5, 2449534.5,
                           2450083.5, 2450630.5, 2451179.5, 2453736.5, 2454832.5, 2456109.5, 2457204.5, 2457754.5])

    df_table = np.asarray([42.184, 43.184, 44.184, 45.184, 46.184, 47.184, 48.184, 49.184, 50.184, 51.184,
                52.184, 53.184, 54.184, 55.184, 56.184, 57.184, 58.184, 59.184, 60.184, 61.184,
                62.184, 63.184, 64.184, 65.184, 66.184, 67.184, 68.184, 69.184])/86400.

    return df_table[(jd_table-jd<0)][-1]
    
    
    
def cleanSpectrumMatrix(ww,ss,plot=False):
    
    from scipy.ndimage import median_filter
    import numpy as np

    ss=ss.astype(np.float64)

    nspikes=3
    
    correctedPixels=list()

    #array where the smoothed spectrum will be stored
    filteredSpectrum=ss.copy()

    #array where the spike-less spectrum will be stored
    cleanedSpectrum=ss.copy()

    for order in range(ss.shape[0]):#working order by order

        median_filter(ss[order,:], size=5, mode='nearest')

        #get a smoothed version of the spectrum where spikes are still there
        filteredSpectrum[order,:] = median_filter(ss[order,:], size=5, mode='nearest')

        #the absolute deviations between original and smoothed order (in principle, they are ~0 except for spikes)
        dev=np.abs(ss[order,:]-filteredSpectrum[order,:])
    
        mad=np.median(np.abs(dev))

        #the pixels sorted by how much they deviate from 0 (last elements are the most deviant ones)
        worstPixels=np.argsort(dev)
        # worstPixels=match[worstPixels]

        cc=list()

        #replacing the last elements with the median of the neighboring pixels
        for b in worstPixels[-nspikes:]:
                        
            if dev[b]>10*mad:
                window=np.where( (np.abs(ww[order,:]-ww[order,b])<.05) & (np.abs(ww[order,:]-ww[order,b])>0) )[0]
         
                cleanedSpectrum[order,b]=np.median(ss[order,window])
                
                cc.append(b)

        correctedPixels.append(cc)
        

    if plot:

        fig,axes=plt.subplots(nrows=4,figsize=(12,12),sharex=True)
        
        for match in range(ss.shape[0]):

            axes[0].plot(ww[match,:],ss[match,:])
            axes[1].plot(ww[match,:],filteredSpectrum[match,:])
            axes[2].plot(ww[match,:],ss[match,:]-filteredSpectrum[match,:])
            axes[3].plot(ww[match,:],cleanedSpectrum[match,:])

            axes[2].scatter(ww[match,correctedPixels[match]],ss[match,correctedPixels[match]]-filteredSpectrum[match,correctedPixels[match]],color='black')
        
        axes[2].legend()
        
        axes[0].set_ylabel('Original sp')
        axes[1].set_ylabel('Smoothed sp')
        axes[2].set_ylabel('Original-Smoothed')
        axes[3].set_ylabel('Cleaned')
        
        # axes[0].set_xlim([5576,5578])
        # axes[0].set_xlim([5535,5536])
    
    return cleanedSpectrum