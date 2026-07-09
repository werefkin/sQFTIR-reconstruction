# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 14:28:45 2026

@author: r.gap
"""

import sQFTIR_processing
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

path = r'./data/'

da_sample = np.loadtxt(
    path +
    'fringes_sample_Subt2__0__13-47-47-776.txt',
    skiprows=14)
da_bg = np.loadtxt(path + 'fringes_bg_Subt2__0__13-47-11-375.txt', skiprows=14)
da_ref = np.loadtxt(path + 'ref_Subt2__0__13-48-05-876.txt', skiprows=14)

wavelengths = da_sample[:, 0]
FD_spectrum_sample = da_sample[:, 1]  # FD = Fourier Domain
FD_spectrum_ref = da_ref[:, 1]
FD_spectrum_bg = da_bg[:, 1]

spectrum_sample, wavenumbers = sQFTIR_processing.calculate_ftir_spectrum(
    FD_spectrum_sample, FD_spectrum_ref, wavelengths)
spectrum_bg, wavenumbers = sQFTIR_processing.calculate_ftir_spectrum(
    FD_spectrum_bg, FD_spectrum_ref, wavelengths)

absorbance = -np.log10(spectrum_sample / spectrum_bg)

fig, ax = plt.subplots()
ax.plot(wavenumbers, absorbance)
ax.set_xlabel('Wavenumbers (cm$^{-1}$)')
ax.set_ylabel('Absorbance')
ax.set_xlim(3000, 2400)
ax.set_ylim(0, 1.8)


# calculate spectral resolution
lambda_central = 810 * 1e-7  # in cm
delta_lambda = 70.25 * 1e-7  # in cm -- bandwidth
n = 1044  # number of spectral channels
zmax = 1 / 4 * (lambda_central)**2 / (delta_lambda / 1044)
maxpos = 116
spectral_resolution = 1 / (maxpos / (n / 2) * zmax)
