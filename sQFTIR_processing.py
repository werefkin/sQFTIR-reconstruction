# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 13:01:25 2026

@author: r.gap
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftshift, ifft
from scipy.interpolate import interp1d


def fourier_transform_interferograms(
    interferograms
):
    n_samples = len(interferograms)
    delta = 15  # to avoid detection of the DC component
    maxpos = np.argmax(abs(interferograms[int(n_samples / 2) + delta:]))
    print('maxpos:', maxpos)
    interferogram_half = interferograms[int(n_samples / 2):]

    # optionally hardcode maxpos in order to fix spectral resolution
    # maxpos = 77

    print('maxpos_actual:', maxpos + delta)
    interferogram_re = interferogram_half[:2 * (maxpos + delta)]
    print('len_interferogram_re', len(interferogram_re))
    # plt.figure()

    # plt.plot(interferogram_re)
    # plt.title('interferogram_re')
    # DC removal
    interferograms = interferogram_re - np.mean(interferogram_re, axis=0)

    # Apodization
    window = np.blackman(len(interferogram_re))
    kernel = interferogram_re * window[:]

    # Zero-padding to 1024 samples
    pad_number = int(1024)  # desired length of interferogram
    padding_needed = int(1 / 2 * (pad_number - len(interferogram_re)))
    kernel_pad = np.pad(
        kernel,
        (padding_needed,
         padding_needed),
        mode='constant')

    # plt.figure()
    # plt.plot(kernel)
    # plt.title('interferogram, BMH window and padding')

    # FFT (NO pre-shift!)
    spectrum = fft(kernel_pad, norm="ortho")
    print('shape_spectrum', np.shape(spectrum))

    return abs(spectrum)


def remap_to_k(data, ref_spectrum, xx_lambd):
    data = data - ref_spectrum
    lambda_space_range = xx_lambd[:]
    k_space_range = 1 / lambda_space_range
    wn_corr = k_space_range
    remap_interp_func = interp1d(
        k_space_range,
        data,
        axis=0,
        fill_value="extrapolate")
    spectral_interferograms = remap_interp_func(wn_corr)
    return spectral_interferograms, wn_corr


def calculate_ftir_spectrum(fd_spectrum, fd_reference, xaxis_lambda):

    # remap to k (to be linear in wavenumbers)
    spectrum_fd, wn_corr = remap_to_k(fd_spectrum, fd_reference, xaxis_lambda)
    # Fourier transform to time domain
    spectrum_td = fftshift(ifft(spectrum_fd, norm='ortho'))
    # plt.figure()
    # plt.plot(spectrum_td)
    # plt.title('Spectrum time domain')
    # fringes in k space
    # plt.plot(wn_corr*1e7, spectrum_fd)
    # center interferogram and Fourier transform back
    spectrum = fourier_transform_interferograms(spectrum_td)
    # correct x axis
    xspace = np.arange(np.shape(xaxis_lambda)[0])
    interp_axis_func = interp1d(xspace, xaxis_lambda, axis=0, kind='cubic')
    xax_interp = interp_axis_func(np.linspace(
        xspace[0], xspace[-1], len(spectrum)))
    wn = 1 / xax_interp * 1e7
    pump_wl = 659.75  # pump laser wavelength in nm
    corr_wn = 1 / pump_wl * 1e7 - wn

    return spectrum, corr_wn
