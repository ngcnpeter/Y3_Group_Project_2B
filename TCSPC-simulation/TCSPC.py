import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import convolve
from scipy import stats
from sympy import symbols, solve, Eq, Function
from sympy.matrices import Matrix
import sympy as sp

rng = np.random.default_rng()

def deconv_fft(signal,kernel):
    '''Deconvolve decay data with IRF kernel using FFT
    Input:  signal - convolved/original signal (1d array)
            kernel - Gaussian kernel (IRF) (same length as signal)
    Output: deconvolved signal (1d array)'''
    deconv_arr =np.fft.ifft(np.fft.fft(signal)/np.fft.fft(kernel))*np.sum(kernel)
    deconv_arr[deconv_arr<1] = 0
    return deconv_arr

class Simulation():
    def __init__(self,amp,tau, run_time=20*60, irfwidth=1e-3,
                 n_bins = 380, window = 20, bg = 10, t0 = 10/19):
        self.amp = amp/np.sum(amp)            #normalized amplitudes array
        self.tau = tau            #lifetimes array (in ns)
        self.run_time = run_time  #data collection time (in s), default =20*60 s
        self.irfwidth = irfwidth  #sigma of Gaussian IRF (in ns), default = 1e-3 ns
        self.n_bins = n_bins      #no. of histogram bins, default = 380
        self.window = window      #decay data time window, ns, default =20
        self.bg = bg              #background count rate 10 per s
        self.t0 = t0              #offset of IRF, 10/19 ns
        self.t = np.linspace(0,window,n_bins) #time array in ns
        self.ker = stats.norm.pdf(self.t,loc = t0,scale = irfwidth) #gaussian kernel

    def multi_exp_data(self,deconv = False):
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
        self.t = np.linspace(0,self.window,self.n_bins)
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
        Iconvol = convolve(puredecay, irf_kernel, mode='full')[:self.n_bins]/np.sum(irf_kernel)

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
        self.noisydecay = rng.poisson(noiseless)

        if deconv == True:
            self.noisydecay = deconv(self.noisydecay,self.ker)

        return self.t,self.noisydecay

    def MC_exp(self, multi = False,deconv = False):
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
                    deconv    deconv with Gaussian IRF if True, default False
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
            self.output = np.zeros((len(self.tau),self.n_bins)) #store output data
            for i in range(len(self.tau)):
                self.output[i],self.bins = np.histogram(t_tot[i], bins=self.n_bins,range = (0,self.window))
                self.output[i] += np.full(self.n_bins, int(self.bg*self.run_time/self.n_bins)) # distribute background count uniformly to each bi
        if multi == True:
            # generate an array of n_photon lifetime with weighted probability using amplitude
            tau_arr = rng.choice(self.tau,len(n_arr),p = self.amp)
            t_tot = rng.normal(t0,self.irfwidth,size = np.shape(tau_arr)) + rng.exponential(tau_arr)
            self.output, self.bins = np.histogram(t_tot, bins=self.n_bins,range = (0,self.window))
        
        if deconv == True:
            self.output = deconv_fft(self.output,self.ker)
        self.bins = self.bins[:-1]
        return self.bins,self.output
    
    def plot(self,ax, MC=True,multi=False,logy = True,deconv = False):
        '''Plot TCSPC decay
           Input: ax      plt axes object
                  MC      default True - >use MC_exp Monte Carlo method, if False -> multi_exp_data
                  logy    default True -> yscale('log'), if False ->yscale('linear') 
           Output: Plot of simulated data'''
        ax.set_title(r'TCSPC Fluorescence Decay ($\tau =$' + str(self.tau) + 'ns)')
        ax.set_xlabel('time/ns')
        ax.set_ylabel('Photon Count')
        if logy == True:
            ax.set_yscale('log')
        #Plot Monte Carlo simulation results
        if MC ==True: 
            self.t,self.y = self.MC_exp(multi = multi,deconv= deconv)
            #arrays of mono exp decay if False, 1 multi-exp decay if True
            if multi == False:
                for i in range(len(self.tau)):
                    ax.plot(self.t,self.y[i],label = str(self.tau[i]) + ' ns')
            else:
                ax.plot(self.t,self.y,label = 'Data')
        #Plot given data generation method            
        else :
            self.t,self.y = self.multi_exp_data(deconv =deconv)
            ax.plot(self.t,self.y,label = 'Data')
        ax.legend()


    
 
