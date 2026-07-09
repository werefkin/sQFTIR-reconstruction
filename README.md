# sQFTIR Reconstruction Protocol



Calculates Fourier domain (FD) mid-infrared spectra from FD interferograms in the near-infrared (center wavelength = 805 nm). 



sQFTIR\_processing.py provides the functions and example.py demonstrates the reconstruction protocol on the data (stored in ./data).



The reconstruction protocol works as follows: 

FD interferograms are DC corrected (ref is subtracted from fringes\_bg and fringes\_sample) and then remapped to k-space (enforce linearity in k-space). 

An inverse Fourier transform is performed and the coherence burst of one of the identical halves of the interferogram is located, isolated and re-centered.  

Apodization, zero padding and Fourier Transformation yield the demodulated transmission spectrum.





The example.py file runs the reconstruction protocol on a Background (no sample in the idler path) and a sample measurement (polypropylene sample in the idler path) and computes the absorbance. 



&#x20; 

