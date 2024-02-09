# Y3_Group_Project_2B
# A Critical Comparison of Time and Frequency Domain Measurements of Fluorescence Lifetimes

## Project Description

#### Time-Correlated Single Photon Counting (TCSPC) Analysis

This Python script provides functions for performing time-correlated single photon counting (TCSPC) analysis. It includes functions for fitting mono-exponential and bi-exponential decay curves, deconvolving decay data with an instrument response function (IRF) kernel using FFT, generating phasor plots, and solving for fractional intensities and lifetimes from phasor coordinates.

## Dependencies

The following Python packages are required to run the script:
- `numpy`
- `matplotlib`
- `scipy`
- `sympy`
- `lmfit`

## Description of folders and files

- `lit_review` : Papers for literature review
- `MC_FLIM_FRET` : MatLab codes developed by WeiYue Chen at University of Cambridge to simulate and analyse fluorescence decay in applications such as Fluorescence Lifetime Imaging Microscopy (FLIM) and F$\"o$rster Resonance Energy Transfer (FRET)$. This folder was downloaded at https://laser.ceb.cam.ac.uk/research/resources.
- `TCSPC-imiage-simulation` : a public repository storing MatLab codes for FLIM images
- `TCSPC-simulation` : Codes and data for UCL Y3 Physics Group Project (Group 2) supervised by Prof. Angus Bain and Dr. Thomas Blacker. Contains pdf documents of tasks and python notebooks and script to aid the comparison of time-domain and frequency domain analysis of muli-exponential fluorescence decay
## Usage

1. Import the required functions from `TCSPC.py` into your Python script:

```python
from TCSPC import *
```

2. Generate or load your time-resolved decay data.

3. Choose the appropriate function based on your analysis needs:

- `exp1(t, tau)`: Returns a mono-exponential decay curve for a given time array `t` and lifetime `tau`.
- `exp_fit(func, tdata, ydata, guess)`: Performs a least-square fit of the provided exponential function (`func`) to the time-resolved data (`tdata` and `ydata`) with an initial guess for the fitting parameters (`guess`).
- `deconv_fft(signal, kernel)`: Deconvolves the decay data (`signal`) with the instrument response function (IRF) kernel using Fast Fourier Transform (FFT).
- `phasor_fft(y, ker, dt)`: Generates phasor coordinates for multi-exponential decay curves given an array of lifetimes (`y`) and corresponding amplitudes, along with the IRF kernel (`ker`) and time interval (`dt`).
- `phasor_plot(ax, w, phasor)`: Creates a phasor plot for data transformed at a/an array of angular frequencies (`w`) using the phasor coordinates (`phasor`).
- `phasor_solve(w, phasor, n=2, num=False, guess=None)`: Solves for fractional intensities and lifetimes from simulated phasor coordinates (`phasor`) using either numerical or analytic solution methods.

4. Perform the desired analysis by calling the appropriate function with the necessary arguments.

5. Visualize the results using matplotlib.

6. `Simulation(tau, amp)`: A class that represents a TCSPC simulation wtih mono or multi-exponential components of lifetimes `tau` and amplitudes `amp`. It allows you to generate synthetic decay data, perform fitting, deconvolution, and phasor analysis, and visualize the results. Demonstration of the use of this class is in `interim_report.ipynb`

## Example

Here's an example of fitting a mono-exponential decay curve using the `exp_fit` function:

```python
import numpy as np
import matplotlib.pyplot as plt
from TCSPC import exp1, exp_fit

# Generate synthetic data
t = np.linspace(0, 10, 1000)  # Time array
tau = 1.5  # Lifetime
data = exp1(t, tau) + np.random.normal(0, 0.05, len(t))  # Simulated data with noise

# Initial guess for fitting parameters
guess = [1.0]  # Guess for tau

# Perform exponential fitting
result, params_opt, chi2_red, fit_report = exp_fit(exp1, t, data, guess)

# Plot the original data and fitted curve
plt.plot(t, data, 'bo', label='Data')
plt.plot(t, result.best_fit, 'r-', label='Fit')
plt.xlabel('Time')
plt.ylabel('Intensity')
plt.legend()
plt.show()

# Print the optimized parameters and fit report
print('Optimized Parameters:')
for name, param in params_opt.items():
    print(f'{name}: {param.value}')

print('\nFit Report:')
print(fit_report)

print('\nReduced Chi-square:', chi2_red)
```

