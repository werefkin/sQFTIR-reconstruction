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


sqftir_data_sa = sQFTIR_processing.calculate_ftir_spectrum(
    FD_spectrum_sample, FD_spectrum_ref, wavelengths)
sqftir_data_bg = sQFTIR_processing.calculate_ftir_spectrum(
    FD_spectrum_bg, FD_spectrum_ref, wavelengths)

absorbance = -np.log10(sqftir_data_sa.spectrum / sqftir_data_bg.spectrum)

fig, ax = plt.subplots()
ax.plot(sqftir_data_sa.corr_wn, absorbance)
ax.set_title('Absorbance spectrum')
ax.set_xlabel('Wavenumbers (cm$^{-1}$)')
ax.set_ylabel('Absorbance')
ax.set_xlim(3000, 2400)
ax.set_ylim(0, 1.8)


# Original spectrum
fig, ax = plt.subplots()
ax.plot(wavelengths,FD_spectrum_sample)
ax.set_title('Original spectrum')
ax.set_xlabel('Wavelengths (nm)')
ax.set_ylabel('Intensity')

# Remapped FD spectrum
fig, ax = plt.subplots()
ax.plot(sqftir_data_sa.wn_corr*1e7,sqftir_data_sa.spectrum_fd)
ax.set_title('FD Spectrum remapped to k-space')
ax.set_xlabel('Wavenumbers (cm$^{-1}$)')
ax.set_ylabel('Intensity')

# Time domain signal
fig, ax = plt.subplots()
ax.plot(sqftir_data_sa.spectrum_td)
ax.set_title('Time domain signal')
ax.set_xlabel('Sample no.')
ax.set_ylabel('Intensity')

# Zero-padded kernel
fig, ax = plt.subplots()
ax.plot(sqftir_data_sa.kernel_pad)
ax.set_title('Zero-centerd & zero-padded kernel')
ax.set_xlabel('Sample no.')
ax.set_ylabel('Intensity')

# Spectrum sample
fig, ax = plt.subplots()
ax.plot(sqftir_data_sa.corr_wn, sqftir_data_sa.spectrum)
ax.set_title('Sample spectrum')
ax.set_xlabel('Wavenumbers (cm$^{-1}$)')
ax.set_ylabel('Intensity')



# calculate spectral resolution
lambda_central = 805 * 1e-7  # in cm
delta_lambda = 70.25 * 1e-7  # in cm -- bandwidth
n = 1044  # number of spectral channels
zmax = 1 / 4 * (lambda_central)**2 / (delta_lambda / 1044)
maxpos = 116
spectral_resolution = 1 / (maxpos / (n / 2) * zmax)
print('spectral resolution (cm^-1):', spectral_resolution)