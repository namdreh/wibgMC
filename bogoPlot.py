#!/usr/bin/python 
#
# bogoRun.py
# Chris Herdman
# 10.31.2014
#
# Plot bogoRun data

import argparse
import math
import numpy as np
import matplotlib.pyplot as plt
import MCstat

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 
    
    # setup the command line parser options 
    parser = argparse.ArgumentParser(description=
                                        'Plot PIMC estimator vs. parameter')
    parser.add_argument('fileNames', help='Scalar estimator files', nargs='+')
    args = parser.parse_args()
                                    
    fileName = args.fileNames[0]
    
    dataType = fileName.split('-')[1]
    
    # Get first line of file
    MCfile = open(fileName,'r')
    MCfile.readline()
    x = np.array(MCfile.readline().split()[1:]).astype(np.float)
    MCfile.close()
    
    # Get MC data
    MCdata = np.loadtxt(fileName)
    mean = np.average(MCdata,axis=0)
    # error = np.std(MCdata,axis=0)/np.sqrt(MCdata.shape[0])
    bins = MCstat.bin(np.array(MCdata))
    error = bins[-1,:]    
    
    if dataType == 'occ':
        N = np.sum(MCdata[0,:])
        y = mean/N
        dy = error/N
        n0 = y[0]
        dn0 = dy[0]
    elif dataType == 'pn0':
        N = 2*x[-1]
        y = mean
        dy = error
        n0 = np.sum( mean*2*x)/N
        
        ai = 1.0-(2.0*x)/N
        dn0 = np.sqrt( np.sum((error[0:-1]*ai[0:-1])**2) )
                
        nzInds = np.nonzero(y[:-1] > 2.0*dy[:-1])[0]        
        ynz = y[nzInds]
        dynz = dy[nzInds]
        
        y = ynz
        dy = dynz
        x = x[nzInds]
        
        Sigma = np.sum(ynz)
        Svn = -np.sum( ynz*np.log(ynz) ) - (1.0-Sigma)*np.log(1.0-Sigma)
        dSdP = -(np.log(ynz) - np.log(1.0-Sigma))
        dSvn = np.sqrt(np.sum( (dSdP*dynz)**2  ))
                
        S2 = -np.log( np.sum( ynz**2 ) + (1.0-Sigma)**2)
        dS2dP = (-1.0)*np.exp(S2)*(2.0*ynz-2.0*(1.0-Sigma))
        dS2 = np.sqrt(np.sum( (dS2dP*dynz)**2  ))
        
        print "\tvon Neumann entropy:\t" + str(Svn) + ' +/- ' + str(dSvn)
        print "\t2nd Renyi entropy:\t" + str(S2) + ' +/- ' + str(dS2) 
         
        
    print "\tCondensate fraction:\t" + str(n0) + ' +/- ' + str(dn0)    
    
    f,ax = plt.subplots()
    ax.errorbar(x,y,dy,linestyle='--',marker='o',markersize=4,capsize=4)
    ax.set_yscale('log')
    plt.xlim([x[0]-0.25,x[-1]+0.25])
    plt.xlabel(r'$m$')
    plt.ylabel(r'$n$')
    plt.show()
    
    
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()