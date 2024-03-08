# Y3_Group_Project_2B
# A Critical Comparison of Time and Frequency Domain Measurements of Fluorescence Lifetimes

## Project Description

#### Time-Correlated Single Photon Counting (TCSPC) Analysis

This repository contains Jupyter notebooks and a Python module for the analysis of fluorescence decay data using various methods. The notebooks explore different techniques, simulations, and fitting algorithms to analyze the lifetime of fluorescent molecules.

## Dependencies

The following Python packages are required to run the script:
- `numpy`
- `matplotlib`
- `scipy`
- `sympy`
- `lmfit`
- `mtalg`

## Description of folders and files

- `lit_review` : Papers for literature review
- `MC_FLIM_FRET` : MatLab codes developed by WeiYue Chen at University of Cambridge to simulate and analyse fluorescence decay in applications such as Fluorescence Lifetime Imaging Microscopy (FLIM) and F$\"o$rster Resonance Energy Transfer (FRET)$. This folder was downloaded at https://laser.ceb.cam.ac.uk/research/resources.
- `TCSPC-imiage-simulation` : a public repository storing MatLab codes for FLIM images
- `TCSPC-simulation` : Codes and data for UCL Y3 Physics Group Project (Group 2) supervised by Prof. Angus Bain and Dr. Thomas Blacker. Contains pdf documents of tasks and python notebooks and script to aid the comparison of time-domain and frequency domain analysis of muli-exponential fluorescence decay

# Fluorescence Decay Analysis

This repository contains Jupyter notebooks and a Python module for the analysis of fluorescence decay data using various methods. The notebooks explore different techniques, simulations, and fitting algorithms to analyze the lifetime of fluorescent molecules.

## Notebooks

1. [phasor.ipynb](https://github.com/ngcnpeter/Y3_Group_Project_2B/blob/main/TCSPC-simulation/phasor.ipynb): Investigates the effect of noise and blind decomposition of phasor to obtain fluorescence decay parameters.
2. [2Bsimulation.ipynb](https://github.com/ngcnpeter/Y3_Group_Project_2B/blob/main/TCSPC-simulation/2Bsimulation.ipynb): Explores different methods to simulate fluorescence decay, including introducing Poisson noise and Monte Carlo simulation.
3. [interim_report.ipynb](https://github.com/ngcnpeter/Y3_Group_Project_2B/blob/main/TCSPC-simulation/interim_report_ipynb): Provides a temporary comparison of fitting and phasor methods of analysis to fluorescence decay data.
4. [Comparing by critical lifetime ratio (PHASOR part) (1).ipynb](https://github.com/ngcnpeter/Y3_Group_Project_2B/blob/main/TCSPC-simulation/comparing_by_critical_lifetime_ratio_ipynb): Compares separability and ability of the phasor method to resolve close lifetime components in fluorescence decay.
5. [τseparation.ipynb](https://github.com/ngcnpeter/Y3_Group_Project_2B/blob/main/TCSPC-simulation/tau_separation_ipynb): Compares different fitting methods for the separation of lifetimes in fluorescence decay.
6. [numberofphotons.ipynb](https://github.com/ngcnpeter/Y3_Group_Project_2B/blob/main/TCSPC-simulation/_number_of_photons_ipynb): Investigates the effect of total collected photon number on fluorescence decay analysis, especially fitting method.
7. [fitting.ipynb](https://github.com/ngcnpeter/Y3_Group_Project_2B/blob/main/TCSPC-simulation/_fitting_ipynb): Compares different fitting algorithms and methods for analyzing time-domain fluorescence decay.
8. [MLE.ipynb](https://github.com/ngcnpeter/Y3_Group_Project_2B/blob/main/TCSPC-simulation/_MLE_ipynb): Uses Maximum Likelihood Estimator (MLE) to fit fluorescence decay data.

## Python Module

The `TCSPC.py` module contains the `Simulation` and `Phasor` classes utilized across the Jupyter notebooks. It includes functions for generating fluorescence decay data, transforming data to phasor using Fast Fourier Transform (FFT), and producing fit results or phasor decomposition results for repeated simulations. This Python script provides functions for performing time-correlated single photon counting (TCSPC) analysis. It includes functions for fitting mono-exponential and bi-exponential decay curves, deconvolving decay data with an instrument response function (IRF) kernel using FFT, generating phasor plots, and solving for amplitudes and lifetimes from phasor coordinates.

## Usage

To run the notebooks and utilize the Python module, follow these steps:
1. Go to https://github.com/ngcnpeter/Y3_Group_Project_2B.git`
2. Install the required dependencies: `pip install -r requirements.txt`
3. Open the desired Jupyter notebook in your Jupyter environment and execute the cells.
4. Make use of the provided notebook functionalities and the functionalities provided by the `TCSPC` module.


Here are some examples:
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
- `phasor_solve(w, phasor, n=2, num=False, guess=None)`: Solves for amplitudes and lifetimes from simulated phasor coordinates (`phasor`) using either numerical or analytic solution methods.
- `phasor_solve_num()`: Solves for amplitudes and lifetimes from simulated phasor coordinates (`phasor`) using numerical methods and phasor from pure, convolved decay.

4. Perform the desired analysis by calling the appropriate function with the necessary arguments.

5. Visualize the results using matplotlib.

6. `Simulation(tau, amp)`: A class that represents a TCSPC simulation wtih mono or multi-exponential components of lifetimes `tau` and amplitudes `amp`. It allows you to generate synthetic decay data, perform fitting, deconvolution, and phasor analysis, and visualize the results. Demonstration of the use of this class is in `interim_report.ipynb`

