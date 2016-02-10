#!/usr/bin/python 
#
# bogoRun.py
# Chris Herdman
# 09.06.2014
#
# Monte Carlo sample a Bogoliubov ground state

import math
import numpy as np
import argparse
import time
import calendar
import sys


# -----------------------------------------------------------------------------
def getFileParameters(fileName,xval):
    Data = np.loadtxt(fileName)
    pval = np.interp(xval,Data[:,0],Data[:,1])
    return pval


# -----------------------------------------------------------------------------
class Progress:
    ''' Progress Meter '''
    # Constructor for Progress Class
    # --------------------------------------------------------
    def __init__(self,Max,deltaPerc):
        self.Max = float(Max)
        self.deltaPerc = deltaPerc
        self.nextPerc = deltaPerc
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    def update(self,count):
        perc = (count/self.Max)*100
        if perc >= self.nextPerc:
            sys.stdout.write('\t'+str(self.nextPerc)+'%')
            sys.stdout.flush()
            self.nextPerc = self.nextPerc + self.deltaPerc
            if self.nextPerc > 100:
                sys.stdout.write("\n")
    # --------------------------------------------------------
    
    

# -----------------------------------------------------------------------------
class BogoParams:
    ''' BogoMC Parameters Class '''
    # Constructor for BogoParams Class
    # ------------------------------------------------------------------------
    def __init__(self):
        
        # Configuration params
        self.N = 0
        self.M = 0
        self.L = 0
        self.nDim = 1
        
        # Interaction params
        self.lam = 0.0
        self.U0 = 0.0
        self.wfnType = ''
        self.CdFile = ''
        self.n0ad = 0.0
        
        # MC params
        self.seed = 0
        
        # Job tag & id
        self.id = 0
        self.tag = ''
        self.logFileName = ''
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    def parseComLin(self):
        paramDict = {}
        parser = argparse.ArgumentParser(description=
                'Monte Carlo sampler for the Bogoliubov ground state\n\
                ./bogoRun.py -N 4 -L 1 -U 10 -M 8 --Nbins 10 -e 10',
                formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument("-N",type=int, help="number of particles")
        parser.add_argument("-M",type=int, help="number of states")
        parser.add_argument("-L",type=float, help="length",default=1.0)
        parser.add_argument("-l","--lam",type=float, help="h-bar^2/2m",
                                default=1.0)
        parser.add_argument("-U",type=float, help="interaction strenth")
        parser.add_argument("--nad",type=float, help="dimensionless strength")
        parser.add_argument("--type",type=str, help="wavefunction type",
                                default='Ueda')
        parser.add_argument("--Cd_file",type=str, help="C_d parameter file")
        parser.add_argument("--sigma_file",type=str, 
                                                    help="sigma parameter file")                                                                
        parser.add_argument("--Nbins",type=int, help="number of bins")
        parser.add_argument("--binsize",type=int, help="size of bins",
                                default=10)
        parser.add_argument("--equil","-e",type=int,default=10, 
                                help="number of equilibration steps")
        parser.add_argument("--dimensions","-d",type=int,default=1, 
                                  help="number of spatial dimensions")                                
        parser.add_argument("-s","--seed",type=int, help="RNG seed",
                                default=0)
        parser.add_argument("--update",type=str, help="update type",
                                default='modes')                            
        
        args = parser.parse_args()
        
        # set parameters
        self.nDim = args.dimensions        
        self.N = args.N
        self.M = args.M
        self.L = args.L
        self.lam = args.lam
        self.U0 = args.U
        self.nad = args.nad
        self.seed = args.seed
        self.Nbins = args.Nbins
        self.binsize = args.binsize
        self.equil = args.equil
        self.wfnType = args.type
        self.CdFile = args.Cd_file
        self.sigmaFile = args.sigma_file
        self.update = args.update      
        
        # Derived Parameters
        self.n = 2.0*self.N/self.L
        
        # Interaction strength
        if self.nad is not None:
            self.a = (self.nad/self.n)**(1.0/3.0)
            self.U0 = 8.0*math.pi*self.lam*self.a
        elif self.U0 is not None:
            self.a = self.U0/(8.0*math.pi*self.lam)
            self.nad = self.n*(self.n)**(3.0)
        else:
            print "Iteraction strength not defined: set -U or --nad"
            sys.exit()              
        
        if self.wfnType == 'Etto':
            if (self.CdFile is not None) and (self.sigmaFile is not None):
                self.Cd = getFileParameters(self.CdFile,self.nad)
                self.sigma = getFileParameters(self.sigmaFile,self.nad)
            else:
                print "Parameter files not defined: set --Cd_file AND --sigma_File"
                sys.exit()              
            
        
        # set id        
        self.id = str(int(calendar.timegm(time.gmtime())-1416362400)).zfill(9)
                
        # build tag string
        self.tag = "{}-{:8.6f}-{}".format(self.N,self.nad,self.id)
        
        self.createLog()
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    def output(self):
        print "\nStarting bogoMC with the following parameters:"
        print'\t{: <10}{: ^7}{:<12}'.format('bogoID',':',self.id)
        print'\t{: <10}{: ^7}{:<12}'.format('wfn. type',':',self.wfnType)
        print'\t{: <8}{: ^5}{:<12}{: <8}{: ^5}{:<8}'.format( \
                                                 'updates','=',self.update,'seed','=',self.seed)
        print'\t{: <8}{: ^5}{:<12}'.format('dim.','=',self.nDim)        
        print'\t{: <8}{: ^5}{:<12}{: <8}{: ^5}{:<8}'.format( \
                                                 'N','=',self.N,'L','=',self.L)
        print'\t{: <8}{: ^5}{:<12}{: <8}{: ^5}{:<8}'.format( \
                                        'lambda','=',self.lam,'M','=',self.M)                                    
        print'\t{: <8}{: ^5}{:<12.4f}{: <8}{: ^5}{:<8.4f}'.format( \
                                        'U0','=',self.U0,'a','=',self.a)                                    
        print'\t{: <8}{: ^5}{:<12.6f}'.format('na^d','=',self.nad)                                        
        if self.wfnType == 'Etto':
            print'\t{: <8}{: ^5}{:<12.4f}{: <8}{: ^5}{:<12.4}'.format( \
                                        'Cd','=',self.Cd,'sigma','=',self.sigma)
                                    
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    def createLog(self):
        self.logFileName = 'bogo-log-'+self.tag+'.dat'
        file = open(self.logFileName,'w')
        file.write('----------------------------\n')
        file.write('Initial parameters\n')
        file.write('----------------------------\n')
        file.write('{: <12}{: ^7}{:<9}\n'.format('bogoID','=',self.id))
        file.write('{: <12}{: ^7}{: <9}\n'.format('seed','=',self.seed))
        file.write('{: <12}{: ^7}{: <9}\n'.format('dimensions','=',self.nDim))
        file.write('{: <12}{: ^7}{: <9}\n'.format('updates','=',self.update))
        file.write('{: <12}{: ^7}{: <9}\n'.format('wfn. type','=',self.wfnType))
        if self.wfnType == 'Etto':
            file.write('{: <12}{: ^7}{: <9.7f}\n'.format('Cd','=',self.Cd))
            file.write('{: <12}{: ^7}{: <9.7f}\n'.format('sigma','=',self.sigma))
        file.write('{: <12}{: ^7}{: <9}\n'.format('Nb','=',2*self.N))
        file.write('{: <12}{: ^7}{: <9.7f}\n'.format('L','=',self.L))
        file.write('{: <12}{: ^7}{: <9}\n'.format('max M','=',self.M))
        file.write('{: <12}{: ^7}{: <9.7f}\n'.format('lambda','=',self.lam))
        file.write('{: <12}{: ^7}{: <9.7f}\n'.format('U0','=',self.U0))
        file.write('{: <12}{: ^7}{: <9.7f}\n'.format('a','=',self.a))
        file.write('{: <12}{: ^7}{: <9.7f}\n'.format('na^d','=',self.nad))
        file.write('{: <12}{: ^7}{: <9}\n'.format('num. of bins','=',self.Nbins))
        file.write('{: <12}{: ^7}{: <9}\n'.format('bin size','=',self.binsize))
        file.write('{: <12}{: ^7}{: <9}\n'.format('equil. steps','=',self.equil))
        file.write('----------------------------\n')
        file.close()
    
    # --------------------------------------------------------


# -----------------------------------------------------------------------------
class BogoConfig:
    ''' Paired boson Fock state object '''
    
    # Constructor for BogoConfig Class
    # ------------------------------------------------------------------------
    def __init__(self,nDim,N,L,M):
        
        # Define Physical Parameters
        self.nDim = nDim                    # nDim: # of spatial dimensions
        self.N = N                          # N: Number of pairs
        self.L = L                          # L: Length of box
        self.n = 2.0*N/L                    # n: density of bosons
        self.M = M                          # M: Number of states
        
        # Define member data structures
        Mtup = (M,) + tuple((2*M-1) for d in range(1,self.nDim))
        self.k2 = np.zeros(Mtup,dtype=float) # M: Momenta of states
        self.occ = np.zeros(Mtup,dtype=int) # occ: Occupation of states
        
        locStr = ''
        for d in range(self.nDim):
            locStr = locStr + 'i,'
        self.loc = np.zeros(N,dtype=(locStr))
        self.k = np.zeros(Mtup,dtype=(locStr))
        
        # Initialize pair occupations
        self.occ[tuple(0 for d in range(self.nDim))] = N        
        self.Mocc = 0
                
        # Set momenta
        twoPionLsq = (2.0*math.pi/L)**2
        mv = np.zeros(nDim,dtype=int)
        if self.nDim == 1:
            for x in range(self.M):
                self.k[x] = (x,)
                self.k2[x] = twoPionLsq*(x**2)
        elif self.nDim == 2:
            for x in range(self.M):
                for y in range( self.M ):
                    for ysign in ('pos','neg'):
                        if ysign == 'pos':
                            yidx = 2*y
                            my = y
                        else:
                            yidx = 2*y-1
                            my = -y                           
                        self.k[x][yidx] = (x,my)
                        self.k2[x][yidx] = twoPionLsq*(x**2+my**2)
        elif self.nDim == 3:
            for x in range(self.M):
                for y in range( self.M ):
                    for z in range( self.M ):
                        for ysign in ('pos','neg'):
                            if ysign == 'pos':
                                yidx = 2*y
                                my = y
                            else:
                                yidx = 2*y-1
                                my = -y
                            for zsign in ('pos','neg'):
                                if zsign == 'pos':
                                    zidx = 2*z
                                    mz = z
                                else:
                                    zidx = 2*z-1
                                    mz = -z                                
                                self.k[x][yidx][zidx] = (x,my,mz)
                                self.k2[x][yidx][zidx] = twoPionLsq*(x**2+my**2+mz**2)
        else:
            print 'Invalid dimension!!'
            sys.exit()
    
    # --------------------------------------------------------
    
    # Check move function
    # -- returns true if move is valid
    # --------------------------------------------------------
    def checkMove(self,mi,mf):
        allowed = False
        if mi[0] >= 0 and mf[0] >= 0 and mi[0] < self.M and mf[0] < self.M:
            if self.nDim == 1:
                if self.occ[mi] > 0:
                    allowed = True
            else:
                if min(mi) >= 0 and min(mf) >= 0 and max(mi) < (2*self.M-1) and max(mf) < (2*self.M-1):
                    if self.occ[mi] > 0:
                        allowed = True        
        return allowed
    # --------------------------------------------------------
    
    # Move Pair function
    # -- attepmts to move a pair from mi to mf
    # -- returns success boolean      
    # ------------------------------------------------------------------------
    def movePair(self,mi,mf):
        allowed = self.checkMove(mi,mf)
        if allowed:
            self.occ[mi] = self.occ[mi] - 1
            self.occ[mf] = self.occ[mf] + 1                
        return allowed
    # --------------------------------------------------------
    
    # ------------------------------------------------------------------------
    def moveThisPair(self,n,mf):
        mi = tuple(self.loc[n])
        allowed = self.checkMove(mi,mf)
        if allowed:
            self.occ[mi] = self.occ[mi] - 1
            self.occ[mf] = self.occ[mf] + 1
            self.loc[n] = mf
                
        return allowed
    # --------------------------------------------------------
    
    
    # ------------------------------------------------------------------------
    def MoccRatio(self,mi,mf):
        ratio = 1.0
        if mi == self.Mocc:
            if (mf > self.Mocc) or (self.occ[mi] == 1):
                ratio = (mf+1)/(float(self.Mocc+1))
        
        return ratio
    # --------------------------------------------------------
    
    
    # Truncate down to M
    # --------------------------------------------------------
    def truncate(self,M):
        self.M = M
        Mtup = (slice(0,M),) + tuple(slice(0,(2*M-1)) for d in range(1,self.nDim))
        self.k2 = self.k2[Mtup]
        self.k = self.k[Mtup]
        self.occ = self.occ[Mtup]
    # --------------------------------------------------------
    


# -----------------------------------------------------------------------------
class Estimator:
    ''' Estimator class for a BogoConfig'''
    
    # --------------------------------------------------------
    def __init__(self,config,tag,id,M,name):
        self.config = config
        self.M = M
        self.bin = np.zeros(self.M)
        self.fileName = 'bogo-'+name+'-'+tag+'.dat'
        self.file = open(self.fileName,'w')
        self.count = 0
        self.id = id
        
        
        # build header
        self.file.write('# BOGOID: '+str(self.id)+'\n')
    # --------------------------------------------------------
    
    
    # --------------------------------------------------------
    def write(self):
        
        # Compute bin mean
        self.bin = self.bin/self.count
        
        # create data string
        binStr = ''
        for m in range(len(self.bin)):
            binStr = binStr + "{:16.12f}".format(self.bin[m])
        binStr = binStr + "\n"
        
        # Write data to file
        self.file.write(binStr)
        
        # Reset bin data and count
        for m in range(len(self.bin)):
            self.bin[m] = 0.0
        self.count = 0
    # --------------------------------------------------------
    
    
    # --------------------------------------------------------
    def close(self):
        self.file.close()
    # --------------------------------------------------------
    


# -----------------------------------------------------------------------------
class BogoOccEstimator(Estimator):
    ''' Estimator class for a BogoConfig'''
    
    # --------------------------------------------------------
    def __init__(self,config,tag,id):
        Estimator.__init__(self,config,tag,id,config.M,'occ')
        
        # build header
        self.header = '#'+'0'.rjust(15)
        for m in range(1,self.M):
            self.header = self.header + str(m).rjust(16)
        self.header = self.header + '\n'
        
        self.file.write(self.header)
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    def measure(self):
        self.bin = self.bin + self.config.occ
        self.count  = self.count + 1
    # --------------------------------------------------------
    
    
    # --------------------------------------------------------
    # Truncate down to M
    def truncate(self,M):
        self.M = M
        self.bin = np.resize(self.bin,M)
    # --------------------------------------------------------
    

    
# -----------------------------------------------------------------------------
class BogoPN0Estimator(Estimator):
    ''' Estimator class for a BogoConfig'''
    
    # --------------------------------------------------------
    def __init__(self,config,tag,id):
        Estimator.__init__(self,config,tag,id,config.N+1,'pn0')
        
        # build header
        self.header = '#'+'0'.rjust(15)
        for m in range(1,self.M):
            self.header = self.header + str(m).rjust(16)
        self.header = self.header + '\n'
        self.nDim = self.config.nDim
        self.condInd = ()
        for i in range(self.nDim):
            self.condInd = self.condInd+(0,)
        
        self.file.write(self.header)
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    def measure(self):
        self.bin[self.config.occ[self.condInd]] = \
                                self.bin[self.config.occ[self.condInd]] +1
        self.count  = self.count + 1
    # --------------------------------------------------------
    
        



# -----------------------------------------------------------------------------
class BogoMC:
    ''' Monte Carlo object to sample a Bogoliubov ground state '''
    
    # Constructor for BogoM Class
    # ------------------------------------------------------------------------
    def __init__(self,params):
        
        self.id  = params.id      
        self.tag = params.tag
        self.logFileName = params.logFileName
        
        # Instantiate configuration
        self.config = BogoConfig(params.nDim,params.N,params.L,params.M)
        
        # Instantiate estimators
        self.estimator = None
        
        # Define Physical Parameters
        self.nDim = params.nDim
        self.wfnType = params.wfnType
        self.updates = params.update
        self.M = 2
        self.maxMFlag = False
        self.lam = params.lam                      # lam: \habar^2/2M
        self.U0 = params.U0                        # U0: Interaction strength
        self.nad = params.nad                                   # n_0 a^d
        self.a = params.a
        
        self.k0 = np.sqrt(8*np.pi*self.a*self.config.n)
        
        # Define member data structures
        Mtup = (self.M,) + tuple((2*self.M-1) for d in range(1,self.nDim))
        self.ck2 = np.zeros(Mtup,dtype=float)  
                                                # ck2: weight of pair in state         
        
        # Instantiate RNG
        self.seed = params.seed
        self.rng = np.random.RandomState(self.seed)
        
        # MC parameters
        self.MCstep = self.config.N
        self.att = 0
        self.acc = 0
        self.startTime = int(calendar.timegm(time.gmtime()))
                
        # Set wavefunction weights
        if self.wfnType == 'Ueda':
            self.setUedaWeights()
        elif self.wfnType == 'Etto':
            self.Cd = params.Cd
            self.sigma = params.sigma
            self.setEttoWeights()
        elif self.wfnType == 'GaussLegg':
            self.setGaussLeggWeights()
            
    # --------------------------------------------------------
    
    # Set wfn. Ueda weights for bogoMC object
    # ------------------------------------------------------------------------
    def setUedaWeights(self):
        
        # for m in range(0,self.config.M):
        #     ek = self.lam*(self.config.k2[m])
        #     #x = 2.0*np.pi*np.sqrt(self.lam*self.config.L/(self.U0*2.0*self.config.N))*m
        #     x = np.sqrt(ek/(self.U0*2.0*self.config.N/(self.config.L**(3.0))))
        #     #self.ck2[m]  = 1.0+2.0*x**2-2.0*x*np.sqrt(x**2+1.0)
        #     self.ck2[m]  = (1.0+x**2-x*np.sqrt(x**2+2.0))**2.0
        
        k2 = self.config.k2/(self.k0**2)
        self.ck2 = (1.0+k2-np.sqrt(k2)*np.sqrt(k2+2.0))**2
        
        if self.nDim > 1:
            self.ck2[0,1::2] = 0.0
        if self.nDim > 2:
            self.ck2[0,0,1::2] = 0.0
        
    # Set wfn. Leggett weights for bogoMC object
    # ------------------------------------------------------------------------
    def setGaussLeggWeights(self):
        
        k2 = self.config.k2/(self.k0**2)
        exp = np.exp(4*np.pi*self.nad*k2)
        self.ck2 = ( k2*exp + 1.0 - np.sqrt((k2*exp+1)**2-1 ))**2
    
    # Set wfn. Ettouhami weights for bogoMC object
    # ------------------------------------------------------------------------
    def setEttoWeights(self):
        k2 = self.config.k2/(self.k0**2)
        # Qk2 = (1.0/self.Cd)*(k2 + self.sigma)*np.exp(4*np.pi*self.nad*(k2))
        # self.ck2 = (1 + Qk2-np.sqrt((1+Qk2)**2 - 1))**2
        
        exponent = 4*np.pi*self.nad*(k2)
        largeExp = np.nonzero( exponent > 100.0)
        exponent[largeExp] = 0.0
        
        Qk2 = (1.0/self.Cd)*(k2 + self.sigma)*np.exp(exponent)
        self.ck2 = (1 + Qk2-np.sqrt((1+Qk2)**2 - 1))**2
        self.ck2[largeExp] = 0.0
        
        zeroTup = ()
        for d in range(self.nDim):
            zeroTup = zeroTup+(0,)
        self.ck2[zeroTup] = 1 
        
        if self.nDim > 1:
            self.ck2[0,1::2] = 0.0
        if self.nDim > 2:
            self.ck2[0,0,1::2] = 0.0
    
     
    # --------------------------------------------------------
    def initEstimators(self):
        #self.occEst = BogoOccEstimator(self.config,self.tag,self.id)
        self.pn0Est = BogoPN0Estimator(self.config,self.tag,self.id)
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    def measure(self):
        #self.occEst.measure()
        self.pn0Est.measure()
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    def closeEstimators(self):
        #self.occEst.close()
        self.pn0Est.close()
    # --------------------------------------------------------
    
    
    # --------------------------------------------------------
    def resetAcc(self):
        self.att = 0
        self.acc = 0
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    def update2(self):
        
        # increment att counter
        self.att = self.att + 1
        
        # Choose indicies
        miv = ()
        for d in range(self.nDim):
            miv = miv +(self.rng.randint(self.M),)
        
        # choose direction to update
        if self.nDim > 1:
            dir = self.rng.randint(self.nDim)
        else:
            dir = 0
                
        mi = miv[dir]
                
        if mi == 0:
            mf = 1
        elif mi == self.M-1:
            mf = self.M-2
        else:
            mf = mi + 1 - 2*self.rng.randint(2)
        
        mfv = ()
        for d in range(self.nDim):
            if d == dir:
                mfv = mfv + (mf,)
            else:
                mfv = mfv + (miv[d],)
        
        accepted = False
        if self.config.checkMove(miv,mfv):
            weight = self.ck2[mfv]/self.ck2[miv]
            r = self.rng.random_sample()
            
            # Set detailed balance factor to correct for move selection
            if mi == 0 or mi == self.M-1:
                dbFactor = 0.5
            elif mf == 0 or mf == self.M-1:
                dbFactor = 2.0
            else:
                dbFactor = 1.0
            if weight*dbFactor > r:
                accepted = self.config.movePair(miv,mfv)
                self.acc = self.acc + 1
        
        return accepted
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    def update(self):
        
        # increment att counter
        self.att = self.att + 1
        
        # Choose indicies
        # mi = self.rng.randint(self.M)        
        if self.config.Mocc == 0:
            mi = 0
        else:
            mi = self.rng.randint(self.config.Mocc+1)
        
        if mi == 0:
            mf = 1
        elif mi == self.M-1:
            mf = self.M-2
        else:
            mf = mi + 1 - 2*self.rng.randint(2)
        
        accepted = False
        if self.config.checkMove(mi,mf):
            weight = self.ck2[mf]/self.ck2[mi]
            r = self.rng.random_sample()
            
            # Set detailed balance factor to correct for move selection
            if mi == 0 or mi == self.M-1:
                dbFactor = 0.5
            elif mf == 0 or mf == self.M-1:
                dbFactor = 2.0
            else:
                dbFactor = 1.0
                
            #set detail balance factor if Mocc changes     
            MoccFact = self.config.MoccRatio(mi,mf)
                        
            if weight*dbFactor*(1.0/MoccFact) > r:
                accepted = self.config.movePair(mi,mf)
                self.acc = self.acc + 1
        
        return accepted
    
    # --------------------------------------------------------
    def updateThis(self):
        
        # increment att counter
        self.att = self.att + 1
        
        # Choose particle
        n = self.rng.randint(self.config.N)
        
        # Choose indicies
        miv = tuple(self.config.loc[n])+tuple()
        
        # choose direction to update
        if self.nDim > 1:
            dir = self.rng.randint(self.nDim)
        else:
            dir = 0
                
        mi = miv[dir]
                
        if mi == 0:
            mf = 1
        elif (dir == 0) and (mi == self.M-1):
            mf = self.M-2
        elif (dir > 0) and (mi == 2*self.M-2):
            mf = 2*self.M-3
        else:
            mf = mi + 1 - 2*self.rng.randint(2)
            
        mfv = ()
        for d in range(self.nDim):
            if d == dir:
                mfv = mfv + (mf,)
            else:
                mfv = mfv + (miv[d],)
                
        accepted = False
        if self.config.checkMove(miv,mfv):
            weight = self.ck2[mfv]/self.ck2[miv]
            r = self.rng.random_sample()
            
            # Set detailed balance factor to correct for move selection
            dbFactor = 1.0
            if mi == 0 or ((dir == 0) and (mi == self.M-1)) or ((dir > 0) and (mi == 2*self.M-2)):
                dbFactor = 0.5*dbFactor
            if mf == 0 or ((dir == 0) and (mf == self.M-1)) or ((dir > 0) and (mf == 2*self.M-2)):
                dbFactor = 2.0*dbFactor
                             
            dbFactor = dbFactor*(self.config.occ[mfv]+1.0)/float(self.config.occ[miv])
                                    
            if weight*dbFactor > r:
                accepted = self.config.moveThisPair(n,mfv)
                self.acc = self.acc + 1
                
        return accepted
    # --------------------------------------------------------
    
    
    # --------------------------------------------------------
    def MCequil(self):
        for u in range(self.MCstep):
            if self.updates == 'particles':
                accepted = self.updateThis()
            else: 
                acepted = self.update2()
            
            # Adjust M if necessary
            if not self.maxMFlag:
                
                lastOcc = False
                for d in range(self.nDim):
                    tup = ()
                    for d2 in range(self.nDim):
                        if d < 1:
                            if d2 == d:
                                tup = tup+(self.M-1,)
                            else:
                                tup = tup+(slice(0,self.M),)
                        else:
                            if d2 == d:
                                tup = tup+(2*self.M-2,)
                            else:
                                tup = tup+(slice(0,2*self.M-1),)
                    if np.any(self.config.occ[tup] > 0):
                        lastOcc = True
                        break
                if lastOcc:
                    if self.M == self.config.M:
                        print '\n\tLast state occupied! Increase M!'
                        self.maxMFlag = True
                    else:
                        self.M = self.M + 1
                        # self.MCstep = self.M
                
                # if self.config.occ[self.M-1] > 0:
                #     if self.M == self.config.M:
                #         print '\n\tLast state occupied! Increase M!'
                #         self.maxMFlag = True
                #     else:
                #         self.M = self.M + 1
                #         self.MCstep = self.M
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    def MCupdate(self):
        accCount = 0
        for u in range(self.MCstep):
            if self.updates == 'particles':
                accepted = self.updateThis()
            else: 
                accepted = self.update2()
            if accepted:
                accCount = accCount + 1
        self.measure()
        return (self.MCstep,accCount)
    # --------------------------------------------------------
    
    
    # Truncate down to M
    # --------------------------------------------------------
    def truncate(self):
        # self.MCstep = self.M
        Mtup = (slice(0,self.M),) + tuple(slice(0,(2*self.M-1)) for d in range(1,self.nDim))
        self.ck2 = self.ck2[Mtup]        
        self.config.truncate(self.M)
         
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    def writeData(self):
        #self.occEst.write()
        self.pn0Est.write()
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    def writeEquilLog(self):
        file = open(self.logFileName,'a')
        file.write('----------------------------\n')
        file.write('Post equilibration\n')
        file.write('----------------------------\n')
        file.write('{: <12}{: ^7}{: <12}\n'.format('Mequil','=',self.M))
        file.write('----------------------------\n')
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    def writeStatsLog(self):
        file = open(self.logFileName,'a')
        file.write('----------------------------\n')
        file.write('Run statistics\n')
        file.write('----------------------------\n')
        file.write('{: <12}{: ^7}{: <12}\n'.format('updates att.','=',self.att))
        file.write('{: <12}{: ^7}{: <12}\n'.format('updates acc.','=',self.acc))
                                    
        runTime = int(calendar.timegm(time.gmtime())) - self.startTime
        file.write('{: <12}{: ^7}{: <12}\n'.format('runtime (s)','=',runTime))
        file.write('----------------------------\n')
    # --------------------------------------------------------
    
    # --------------------------------------------------------
    def finalize(self):
        self.writeStatsLog()
        self.closeEstimators()
    # --------------------------------------------------------
    

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main():
    
    # -----------------------------------------------
    # Initialize
    params = BogoParams()
    params.parseComLin()
    params.output()
    bogoMC = BogoMC(params)
    
    # -----------------------------------------------
    # Equilibration loop
    progMeter = Progress(params.equil,10)
    print('\nEquilibrating:')
    for b in range(params.equil):
        for MCu in range(params.binsize):
            bogoMC.MCequil()
        progMeter.update(b+1)
    bogoMC.resetAcc()
    bogoMC.truncate()
    print '\tUsing M\t=\t'+str(bogoMC.M)
    bogoMC.initEstimators()
    bogoMC.writeEquilLog()
    
    # -----------------------------------------------
    # Main MC looop
    progMeter = Progress(params.Nbins,10)
    print('\nSampling:')
    for b in range(params.Nbins):
        for MCu in range(params.binsize):
            accData = bogoMC.MCupdate()
        bogoMC.writeData()        
        progMeter.update(b+1)
    
    # -----------------------------------------------
    # Finalize
    bogoMC.finalize()
    print "Statistics:"
    print '\t'+str(bogoMC.att) + ' updates attempted '
    print '\t'+str(bogoMC.acc) + ' updates accepted '
    print '\t'+"Acceptance rate: " + str(bogoMC.acc/(float(bogoMC.att))) + '\n'

# -----------------------------------------------------------------------------
if __name__ == "__main__": 
    main()