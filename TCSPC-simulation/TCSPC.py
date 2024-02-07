import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import convolve
from scipy import stats
from sympy import symbols, solve, Eq, Function
from sympy.matrices import Matrix
import sympy as sp

rng = np.random.default_rng()

class Simulation():
    def __init__(self,amp,tau, run_time=20*60, irfwidth=1e-5,
                 n_bins = 380, window = 20, bg = 10, t0 = 10/19):
        self.amp = amp/np.sum(amp)            #normalized amplitudes array
        self.tau = tau            #lifetimes array (in ns)
        self.run_time = run_time  #data collection time (in s), default =20*60 s
        self.irfwidth = irfwidth  #sigma of Gaussian IRF (in ns), default = 1e-5 ns
        self.n_bins = n_bins      #no. of histogram bins, default = 380
        self.window = window      #decay data time window, ns, default =20
        self.bg = bg              #background count rate 10 per s
        self.t0 = t0              #offset of IRF, 10/19 ns

    def multi_exp_data(self):
        '''Generate TCSPC fluorescence decay data (not Monte Carlo method)
        Inputs: amplitudes - fractional intensities of each lifetime component (1d array)
                lifetimes  - lifetime array (1d array)
                acquisitiontime - in s
                irfwidth   - sigma of Gaussian IRF
                n_bins     - no. of histogram bins, default = 380
                window     - decay data time window, ns, default =20
        Outputs: t (time array), noisydecay (decay data)'''
        amplitudes = self.amp 
        lifetimes  = self.tau
        acquisitiontime = self.run_time
        irfwidth = self.irfwidth
        t0 = self.t0  #IRF offset, ns
        t = np.linspace(0,self.window,self.n_bins)
        if irfwidth == 0:
            irfwidth = 1e-8
 
        # check that each amplitude has a corresponding lifetime
        if len(amplitudes) != len(lifetimes):
            return
        # generate a multiexponential decay starting at 1 at t=0
        # using the supplied amplitudes and lifetimes
        # sum_i^n A_i exp(-t/tau_i)
        puredecay = sum([amplitudes[j] * np.exp(-t / lifetimes[j]) for j in range(len(lifetimes))])
        #IRF
        irf_kernel = stats.norm.pdf(t,loc = t0, scale = irfwidth)
        # convolute the IRF and decay and trim to 381 bins
        Iconvol = convolve(puredecay, irf_kernel, mode='full')[:self.n_bins]

        # we do our measurements at 2500 counts per second
        # calculate how many fluorescence counts per second this corresponds to
        # i.e. subtract background from total counts
        fluorate = 2500 - self.bg
        # calculate total number of fluorescence photons counted in measurement
        totalfluorescence = fluorate * acquisitiontime
        # now scale the multiexponential decay so it contains this many counts
        noiseless = totalfluorescence * Iconvol / np.sum(Iconvol)
        # and add on 'bg' counts per second spread evenly across all bins
        noiseless = noiseless + (self.bg * acquisitiontime /self.n_bins)
        # finally add Poisson noise to each bin
        noisydecay = rng.poisson(noiseless)

        return t,noisydecay

    def MC_exp(self, multi = False):
        '''If multi == False:
            Generate n_tau mono-exponetial decay curves for an array of lifetimes (n_tau) 
            using Monte Carlo method
            Photon count rate is 2500 per s
            Input:  tau       (1d array of lifetimes, size = n_tau)
                    run_time  (in s) data collection time
                    irfwidth  (sigma of Gaussian IRF)
                    n_bins    no. of histogram bins, default = 380
                    window    Decay data time window, ns, default =20
                    multi     generate n_tau mono-exp decays if False, generate 1 multi-exp decay if True
            output: time array, time-domain decay data (381 by n_tau) matrix (380 bins, 20ns window)
         
           If multi == true, generate one multi-exponential decay curves (sum A_i exp(-t/tau_i)'''
        #IRF properties
        t0 = self.t0 # ns, offset
        n_photon = self.run_time*(2500-self.bg) #no. of photon collected, 2500 photons per s
        n_arr = np.ones(n_photon) #array for meshgrid
        N_arr, Tau = np.meshgrid(n_arr,self.tau)
        #Generate time for each photon, sum of normal distribution (IRF) and exponential distribution (decay)
        if multi == False:
            t_tot = rng.normal(t0,self.irfwidth,size = np.shape(Tau)) + rng.exponential(Tau)
            output = np.zeros((len(self.tau),self.n_bins)) #store output data
            for i in range(len(self.tau)):
                output[i],bins = np.histogram(t_tot[i], bins=self.n_bins,range = (0,self.window))
                output[i] += np.full(self.n_bins, int(self.bg*self.run_time/self.n_bins)) # distribute background count uniformly to each bi
        if multi == True:
            # generate an array of n_photon lifetime with weighted probability using amplitude
            tau_arr = rng.choice(self.tau,len(n_arr),p = self.amp)
            t_tot = rng.normal(t0,self.irfwidth,size = np.shape(tau_arr)) + rng.exponential(tau_arr)
            output, bins = np.histogram(t_tot, bins=self.n_bins,range = (0,self.window))
        return bins[:-1],output
    