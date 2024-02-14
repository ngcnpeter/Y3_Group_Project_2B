import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import convolve
from scipy import stats
from sympy import *
from sympy.matrices import Matrix
import sympy as sp
from lmfit import Model, Parameters
import lmfit
import inspect
import pandas as pd

rng = np.random.default_rng()
def exp1(t,tau):
    '''return mono-exponential exp(-t/tau)
       t    time array (ns)
       tau  lifetime   (ns)'''
    return np.exp(-t/tau)

def exp2(t,A1,tau1,tau2):
    '''returns bi-exponential A1*exp(-t/tau1) + (1-A1)*exp(-t/tau2)
       t    time array (ns)
       A1   amplitude 1
       tau1 lifetime 1 (ns)
       tau2 lifetime 2 (ns)
    '''
    return A1*np.exp(-t/tau1)+(1-A1)*np.exp(-t/tau2)

def exp_fit(func,tdata,ydata,guess,end = int((15/20*380)),bg = 10, run_time = 20*60,weights = None):
    '''use least-square fit for given exponential function (exp1 or exp2)
       Inputs:
       func      exp function to be fitted 
       tdata     time array (non-trimmed)
       ydata     photon count (non-trimmed)
       guess     guess intial parameters for fitting
       end       trim the end point to avoid low count statistics
       bg        background count per s
       run_time  run_time (s)
       weights   weights for the data points of the fit (1/yerr)
       Outputs:
       result        lmfit result
       params_opt    fitted parameters
       chi2_red      reduced chi2
       fit_report    fit_report from lmfit
       '''
    model = Model(func)
    params = Parameters()
    # Get the parameter names and default values from the input function
    params_name = inspect.signature(func).parameters
    params_name = list(params_name.keys())[1:]  # Exclude 'x' from parameters
    for i,name in enumerate(params_name):
    # Add initial guess value for the parameter
        params.add(name,value=guess[i],min = 0)

    #Trim and scale data for fitting
    ydata = ydata-np.full(len(ydata),int(bg*run_time/len(tdata)))#subtract background from each bin
    max_idx = np.argmax(ydata) #index of data point with maximum photon count N(0)
    tdata = tdata[:end-max_idx] #start from t = 0
    ydata = ydata[max_idx:end]  #start from max.
    ydata = ydata/ydata[0] # scale y data such that the beginning is 1 

    
    result = model.fit(ydata, params, t=tdata,weights = weights) #perform least squares fit
    params_opt = result.params #optimized params
    chi2= result.chisqr #chi squared
    chi2_red = result.chisqr/(len(tdata)-len(params))
    fit_report = result.fit_report()
    return result, params_opt, chi2_red, fit_report

def fit_df(results):
    '''Generates an info_df and par_df to store fitted parameters and chi2
       Input:
       results    list of lmfit ModelResult objects 
       output:
       info_df (information dataframe), par_df (parameter dataframe)
       '''
    
    # Extract the information from the result objects
    par_col = ['_val','init_value','stderr','correl'] #column names for parameter dataframe
    attribute_names = ['chisqr', 'redchi'] + par_col
    info_dict = {attribute_name: [] for attribute_name in attribute_names}
    par_df = pd.DataFrame()

    for result in results:
        #create data frame storing parameters
        par_df_new = pd.concat({k: pd.Series(vars(v)).T for k, v in result.params.items()}, axis=1) #parameter attributes in pd.DataFrame
        par_df_new = par_df_new.loc[par_col] #select value, initial value, error, and correlation 
        par_df_new.loc['correl'] = [{k : f'{float(v):.3g}' for k,v in pair.items() }for pair in  par_df_new.loc['correl'].values] #round correlations
        #append the new df to exisiting df
        par_df = pd.concat([par_df,par_df_new])
        info_dict['chisqr'].append(result.chisqr) #chi2
        info_dict['redchi'].append(result.redchi) #reduced chi2
        for col in par_col[:-1]:
            info_dict[col].append([f'{v:.3g}' for v in par_df_new.T[col].values]) #store as list in this data frame
        info_dict['correl'].append(par_df_new.T['correl'].values) #correlation dictionary
    info_df = pd.DataFrame(info_dict) 
    return info_df, par_df


def deconv_fft(signal,kernel):
    '''Deconvolve decay data with IRF kernel using FFT
    Input:  signal - convolved/original signal (1d array)
            kernel - Gaussian kernel (IRF) (same length as signal)
    Output: deconvolved signal (1d array)'''
    deconv_arr =np.fft.ifft(np.fft.fft(signal)/np.fft.fft(kernel))*np.sum(kernel)
    deconv_arr[deconv_arr<1] = 0
    return deconv_arr

def kernel(t,t0 = 10/19,sigma = 1e-3):
    '''return Gaussian kernel 
       t      time array
       Optional:
       t0     centre
       sigma  standard deviation
       Default values are provided in given script '''
    return stats.norm.pdf(t,loc = t0, scale = sigma)


def phasor_fft(y,ker,dt):
        '''Generate phasor of multi-exponetial decay curves for an array of lifetimes (n_tau) with corresponding amplitudes
        Photon count rate is 2500 per s
        Input:  y   signal to be transformed
                ker IRF kernel
                dt  time interval (inverse of sampling rate)

        output: angular frequency w, phasor (array of complex number, real<->g, -imag <->s coordinate )'''
        if len(np.shape(y)) == 1:
            y_sum = np.sum(y)
        else:
            y_sum = np.sum(y, axis = 1)
        #transpose to allow division for multiple decay curves
        phasor = (np.fft.fft(y).T/y_sum).T/np.fft.fft(ker)*np.sum(ker)  
        phasor.imag*=-1 #times -1 to imaginary component to make it positive
        freq = np.fft.fftfreq(len(ker), d=dt) #frequency
        w = 2*np.pi*freq #angular frequency
        return w, phasor


def phasor_coordinate(w,tau):
    '''returns phasor coordinates in a phasor plot with given
       Inputs: 
       w     angular frequency
       tau   lifetime

       Outputs:
       g (horizontal), and s (vertical) coordinates
       '''
    return 1/(1+(w*tau)**2), w*tau/(1+(w*tau)**2)

def exp_FT(omega,tau):
    '''Analytic solution to Fourier Transform (normalized, i.e. divided by int_0^infty exp(-t/tau)dt) 
    of mono exponential decay with components lifetime tau
    Input:
    omega     angular frequency array
    tau       lifetime'''
    W, Tau = np.meshgrid(omega,tau)
    return 1/(1+(W*Tau)**2) + 1j*W*Tau/(1+(W*Tau)**2)

def multi_exp_FT(omega,A,tau):
    '''Analytic solution to Fourier Transform (normalized, i.e. divided by int_0^infty exp(-t/tau)dt) 
    of multi exponential decay with components lifetime tau
    Input:
    omega     angular frequency array
    A         amplitude array
    tau       lifetime'''
    coeff = A*tau/np.sum(A*tau) #coefficient of the sum of mono_exp_FT
    mono_arr = exp_FT(omega,tau)#array of FT of each lifetime
    return np.dot(coeff,mono_arr)


def phasor_plot(ax,w,phasor):
    '''Create phasor plot for data transformed at a/an array of angular frequency w
       Inputs:
       ax      plk.axes object for plotting
       w       angular frequency (value or array) /GHz chosen to be plotted
       phasor  FFT of a decay curve
       '''
    x = np.linspace(0,1,1000)
    y_circ = np.sqrt(0.5**2-(x-0.5)**2)
    for i in range(len(w)):
        w0 = w[i]
        ax.scatter(np.real(phasor[i]), np.imag(phasor[i]), label = f'f = {w0/2/np.pi:.3f} GHz')
        ax.plot(x,y_circ,'k') #universal circle
        ax.legend()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel('G')
    ax.set_ylabel('S')
    ax.set_title(f'Phasor Plot')
    ax.axis('equal')
    ax.grid()

def phasor_solve(w,phasor,n=2,num = False,guess=None):
    '''Solve for samplitudes and lifetimes from simulated phasor coordinates
       Input: w        angular frequency array
              phasor   output from phasor_fft
              n        number of components (Default 2)
              num      True for numerical solution, False for analytic solution
              guess    guess for numerical solution'''
    # Define the variables and symbols
    A = symbols('A1:%d' % (n+1)) #amplitudes
    t = symbols('t1:%d' % (n+1)) #lifetimes

    equations = [sum([A[j]for j in range(n)])-1]
    At_sum = sum([A[j]*t[j]for j in range(n)])

    # Generate the equations using different angular frequencies
    for i in range(1,2*n):
        equation = sum([A[j]*t[j]/At_sum/ ((w[i] * t[j])**2 + 1) for j in range(n)]) - np.real(phasor)[i] #g coordinate of phasor
        equations.append(equation)

    # Solve the system of equations
    if num == True:
        solution = nsolve(equations,[n for n in A]+[n for n in t],guess, solver='bisect')
        solution = np.concatenate(np.array(solution).astype(float)) #1d array
    else:
        solution = solve(equations)[0]
        solution = {str(k):float(v) for k,v in solution.items()} #convert symbols to string
    return solution

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
        self.t = np.linspace(0,window,n_bins+1)[:-1] #time array in ns
        self.dt = window/n_bins
        self.ker = kernel(self.t,t0,irfwidth) #gaussian kernel
        self.MC_exp()
        t,self.y_arr = self.MC_exp_hist(multi = False) #arrays of mono-exp decays
        t,self.y = self.MC_exp_hist(multi = True)      #array of multi-exponential decays
        t,self.y2 = self.multi_exp_data()              #array of multi-exponential decays
        #store angular frequency and correspnding fourier transform of self.y 
        #in self.w and self.phasor
        self.phasor_fft()  
        #store photon count data  of 100 simulationsin self.sim_data 
        # and corresponding phasor data in self.phasor_data
        self.repeat_sim(100) 

    def multi_exp_data(self,deconv = False):
        '''Generate TCSPC fluorescence decay data (not Monte Carlo method)
        Inputs: amplitudes - fractional intensities of each lifetime component (1d array)
                lifetimes  - lifetime array (1d array)
                acquisitiontime - in s
                irfwidth   - sigma of Gaussian IRF
                n_bins     - no. of histogram bins, default = 380
                window     - decay data time window, ns, default =20
        Outputs: self.t (time array), self.y2 (decay data)'''
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
        puredecay = sum([amplitudes[j] * np.exp(-self.t / lifetimes[j]) for j in range(len(lifetimes))])
        #IRF
        irf_kernel = stats.norm.pdf(self.t,loc = t0, scale = irfwidth)
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
        self.y2 = rng.poisson(noiseless)

        if deconv == True:
            self.y2 = deconv_fft(self.y2,self.ker)

        return self.t,self.y2

    def MC_exp(self, multi = False):
        '''If multi == False:
            Generate n_tau mono-exponetial decay curves for an array of lifetimes (n_tau) 
            using Monte Carlo method
            Photon count rate is 2500 per s
            Input:  
                    n_bins    no. of histogram bins, default = 380 when n_bins is not specified
                    multi     generate n_tau mono-exp decays if False, generate 1 multi-exp decay if True
                    deconv    deconv with Gaussian IRF if True, default False
            output: array of generated total time for each photon
         
           If multi == true, generate one multi-exponential decay curves (sum A_i exp(-t/tau_i)'''
        #IRF properties
        t0 = self.t0 # ns, offset
        n_photon = self.run_time*(2500-self.bg) #no. of photon collected, 2500 photons per s
        n_arr = np.ones(n_photon) #array for meshgrid
        N_arr, Tau = np.meshgrid(n_arr,self.tau)
        #set to default value if not provided
        
        #Generate time for each photon, sum of normal distribution (IRF) and exponential distribution (decay)
        if True:
            self.t_tot_2D = rng.normal(t0,self.irfwidth,size = np.shape(Tau)) + rng.exponential(Tau)
        if True:
            # generate an array of n_photon lifetime with weighted probability using amplitude
            # note that normalized exponential probability density function (pdf) is tau*exp(-t/tau)
            # the weighting of exp pdf of different tau becomes: A_i*tau_i, A_i is the amplitude of decay
            tau_arr = rng.choice(self.tau,len(n_arr),p = self.amp*self.tau/np.sum(self.amp*self.tau))
            self.t_tot = rng.normal(t0,self.irfwidth,size = np.shape(tau_arr)) + rng.exponential(tau_arr)
        if multi == True:
            return self.t_tot
        else:
            return self.t_tot_2D

    def MC_exp_hist(self,n_bins = None, window = None,deconv = False,multi = True):
        '''Histogram the MC_exp generated total time data
        return bins array and photon number (against time) array'''
        n_bins = n_bins or self.n_bins 
        window = window or self.window 
        #multi = True -> multi-exp decay
        if multi == True: 
            self.y, self.bins = np.histogram(self.t_tot, bins=n_bins,range = (0,window))
            self.y += np.full(n_bins, int(self.bg*self.run_time/n_bins))
        else:
        #multi = False -> separate decays
            self.y = np.zeros((len(self.tau),n_bins)) #store output data
            for i in range(len(self.tau)): 
                self.y[i],self.bins = np.histogram(self.t_tot_2D[i], bins=n_bins,range = (0,window))
                self.y[i] += np.full(n_bins, int(self.bg*self.run_time/n_bins)) # distribute background count uniformly to each bin
            self.y_arr = self.y
        self.bins = self.bins[:-1]
        if deconv ==True:
            self.y =deconv_fft(self.y,kernel(self.bins))
        return self.bins,self.y

    def plot(self,ax, y, logy = True,deconv = False):
        '''Plot TCSPC decay
           Input: ax      plt axes object
                  MC      default True - >use MC_exp Monte Carlo method, if False -> multi_exp_data
                  y       photon number array (ydata)
                  logy    default True -> yscale('log'), if False ->yscale('linear') 
           Output: Plot of simulated data'''
        ax.set_title(r'TCSPC Fluorescence Decay ($\tau =$' + str(self.tau) + 'ns)')
        ax.set_xlabel('time/ns')
        ax.set_ylabel('Photon Count')
        if logy == True:
            ax.set_yscale('log')
        if deconv == True:
            y =deconv_fft(y,kernel(self.bins))

        if len(np.shape(y)) == 1:
            #arrays of mono exp decay if shape != 1, 1 multi-exp decay if ==1
            ax.plot(self.t,y,label = 'Data')
        else:
            #plot mono-exp decay curves of different lifetimes
            for i in range(len(self.tau)):
                    ax.plot(self.t,y[i],label = str(self.tau[i]) + ' ns')
        
        ax.legend()
    def fit(self,func,y = None,plot = False,guess = None,end = None, bg = None, run_time = None,ax=None):
        #set default values from object attributes unless specified
        if y is None:
            y = self.y #photon count
        guess = guess or list(self.amp[:-1])+self.tau #initial guess for fit
        end = end or int((self.n_bins*3/4)) #end index
        bg = bg or self.bg
        run_time = run_time or self.run_time
        self.fit_result, self.par, self.chi2_red,self.fit_report = exp_fit(
            func,self.t,y,guess = guess,end = end, bg = bg, run_time = run_time)
        if plot == True:
            #pass an ax object for fitting
            self.fit_result.plot_fit(ax)
            ax.set_yscale('log')
            ax.set_ylabel('Photon Count')
            ax.set_xlabel('time/ns')
            

    def phasor_fft(self,MC=True,multi = True,n_bins = None, window = None):
        '''Generate phasor of multi-exponetial decay curves for an array of lifetimes (n_tau) with corresponding amplitudes
        Photon count rate is 2500 per s
        Input:  amp (1d array of amplitudes of each lifetime component)
                tau (1d array of lifetimes, size = n_tau)
                run_time (in s)
                irfwidth  (sigma of Gaussian IRF)
                n_bins no. of bins. default self.n_bins

        output: angular frequency w, phasor (array of complex number, real<->g, -imag <->s coordinate )'''
        #set default values unless specified
        n_bins = n_bins or self.n_bins
        window = window or self.window 
        self.bins, self.y = self.MC_exp_hist(multi = multi,n_bins=n_bins) #call this function to update y,y_arr
        
        if MC == False:
            y = self.y2
        else:
            y = self.y
        ker = kernel(self.bins)
        self.w, self.phasor = phasor_fft(y,ker,self.window/len(self.bins))
        return self.w, self.phasor

    def repeat_sim(self,n_repeat,MC = False, multi = True):
        '''Store photon count array of n_repeat simulations in sim_data (n_repeat by n_bins) array'''
        #create array to store simulation data
        self.sim_data = np.zeros((n_repeat,self.n_bins))
        #repeat simulation
        for i in range(n_repeat):
            if MC == False:
                bins,y = self.multi_exp_data()
            if MC == True:
                self.MC_exp(multi = multi) #default multi-exponential
                bins,y = self.MC_exp_hist()
            self.sim_data[i] = y             #store simulated data
        self.w, self.phasor_data = phasor_fft(self.sim_data,self.ker,self.dt) #transform stored data to phasor


class Phasor(Simulation):
    def __init__(self, amp, tau, run_time=20 * 60, irfwidth=0.001, n_bins=380, window=20, bg=10, t0=10 / 19):
        super().__init__(amp, tau, run_time, irfwidth, n_bins, window, bg, t0)
        self.A_funcs_list(n = len(self.tau))
        self.generate_df()
    
    def A_funcs_list(self,n=2):    
        '''Return a list of functions to evaluate amplitude  A_i of n-eponential 
        from calculated weighting (f_i) and lifetimes (t_i) in phasor plot
        in the form {'A_i':A_i_sol_func} 
        Input:
        tuple of (f1, ... fn, t1,...,tn)
        Output:
        array of A1,...An
        '''
        A = symbols('A1:%d' % (n+1)) #amplitudes
        f = symbols('f1:%d' % (n+1)) #weighting of each phasor -> the system is more easily solved for f than A
        t = symbols('t1:%d' % (n+1)) #lifetimes
        A_eqs = [A[i]*t[i]/sum([A[j]*t[j] for j in range(n)])-f[i] for i in range(n-1)]
        A_eqs.append(sum([A[i] for i in range(n)])-1 )
        A_soln = solve(A_eqs,[n for n in A]) #solve for A in terms of f and t
        #Dictionary of functions to evaluate A_i {'A_i':A_i_sol_func}
        self.A_funcs = [lambdify([v for v in f]+[v for v in t], expr) for key, expr in A_soln.items()]
        return self.A_funcs

    def A_solve(self,f_t_tuple):
        '''Return an array of amplitude  A_i of n-eponential evaluated at
        calculated weighting (f_i) and lifetimes (t_i) in phasor plot
        using functions in self.A_funcs
        in the form {'A_i':A_i_sol_func} 
        Input:
        tuple of (f1, ... fn, t1,...,tn)
        Output:
        array of A1,...An
        '''
        A_array = [func(*f_t_tuple) for func in self.A_funcs]
        return np.array(A_array)

    def phasor_solve(self,w=None,phasor=None,n=2,num = False,guess=None):
        '''Solve for samplitudes and lifetimes from simulated phasor coordinates
        Input: w        angular frequency array
                phasor   output from phasor_fft
                n        number of components (Default 2)
                num      True for numerical solution, False for analytic solution
                guess    guess for numerical solution'''
        #set default values unless given
        if w is None:
            w = self.w 
        if phasor is None:
            phasor = self.phasor
        # Define the variables and symbols
        A = symbols('A1:%d' % (n+1)) #amplitudes
        f = symbols('f1:%d' % (n+1)) #weighting of each phasor -> the system is more easily solved for f than A
        t = symbols('t1:%d' % (n+1)) #lifetimes
        equations = [sum([f[j]for j in range(n)])-1]

        # Generate the equations using different angular frequencies to solve for f and t
        for i in range(1,2*n):
            #equation = sum([A[j]*t[j]/At_sum/ ((w[i] * t[j])**2 + 1) for j in range(n)]) - np.real(phasor)[i] #g coordinate of phasor
            equation = sum([f[j]/ ((w[i] * t[j])**2 + 1) for j in range(n)]) - np.real(phasor)[i] #g coordinate of phasor
            equations.append(equation)

        # Solve the system of equations
        if num == True:
            self.solution = nsolve(equations,[n for n in f]+[n for n in t],guess, solver='bisect')
            self.solution = np.concatenate(np.array(self.solution).astype(float)) #1d array
            self.solution[:n] = self.A_solve(self.solution)
        else:
            self.solution = solve(equations)[0]
            self.solution = np.array([v for v in self.solution.values()]).astype(float) #turn solution into 1d array
            self.solution[:n] = self.A_solve(self.solution)
        return self.solution

    def generate_df(self,w = None, phasor_data = None,num=True):
        ''' Generate DataFrame of solved A and t values for self.phasor_data (repeated simulations)
        Input: 
        w          angular frequency array
        phasor_data phasor_data array of repeated simulations
        num        True for numerical solution, False for analytic solution
        Output:
        phasor_df DataFrame to store parameters of bi-exponential decays'''
        if w is None:
            w = self.w 
        if phasor_data is None:
            phasor_data = self.phasor_data
        self.df = pd.DataFrame()
        for i in range(len(self.phasor_data)):
            sol = self.phasor_solve(w,phasor_data[i],num = num, guess = list((self.amp*self.tau)/np.sum(self.amp*self.tau))+self.tau) #solution
            sol = {k:v for k,v in zip(['A1','A2','t1','t2'],sol)} #convert solution to dict 
            phasor_dict = {str(round(self.w[n]/np.pi/2,2)):phasor_data[i,n]for n in range(1,5)} #record the phasor positions
            sol.update(phasor_dict)# append dictionary
            self.df = pd.concat([self.df,pd.DataFrame(sol, index=[i])]) #concatenate the results into 1 dataframe
        return self.df



 
