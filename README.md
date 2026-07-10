# sQFTIR Reconstruction Protocol

This repository contains a Python implementation of scanless quantum Fourier-transform infrared spectroscopy. An example provides the post-processing protocol used to reconstruct mid-infrared spectra from near-infrared spectral interferograms acquired at a center wavelength of approximately 805 nm.

The file `sQFTIR_processing.py` contains the processing functions, while `example.py` demonstrates the complete reconstruction workflow using the datasets stored in the `./data` directory.

## Reconstruction workflow

The reconstruction protocol consists of the following steps:

1. **DC correction**

   The reference signal is subtracted from both the background and sample interferograms, `fringes_bg` and `fringes_sample`, respectively.

2. **Remapping to k-space**

   The corrected interferograms are resampled onto a uniformly spaced wavenumber axis to enforce linear sampling in k-space.

3. **Fourier transformation to the conjugate domain**

   An inverse Fourier transform is applied to the remapped interferogram to retrieve time-domain signals. The coherence burst is then identified in one of the two symmetric halves of the transformed signal, isolated, and shifted to the center of the processing window.

4. **Spectral reconstruction**

   Apodization and zero-padding are applied to the isolated coherence burst. A subsequent Fourier transform yields the demodulated transmission spectrum.

## Example

The `example.py` script applies the reconstruction protocol to:

* a background measurement acquired without a sample in the idler beam path; and
* a sample measurement acquired with a polypropylene sample in the idler beam path.

The reconstructed background and sample spectra are then used to calculate the absorbance spectrum.
