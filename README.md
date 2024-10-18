# WD-BASS
White Dwarf Binary And Single Star (WD-BASS) package for spectral fitting. It can do spectroscopic, photometric and hybrid fitting with synthetic spectra for single stars or binaries. It can also do Gaussian/Lorentzian profile fitting. The program excels not only for WDs, but any source where you have a relatively smooth continuum.

[![DOI](https://zenodo.org/badge/786251349.svg)](https://zenodo.org/doi/10.5281/zenodo.11188043)

# How to get started?
There is a manual that has the most detail and describes every input/output, and I have made a flowchart to assist in the flow of the package. I recommend installing the package and taking a look at the example data which can be viewed as tutorials. Input files here are for the tutorial data, but I have made scripts to automate the creation of the input file (read the manual, particularly mentions of Create_yaml.py).

On top of zipping and extracting the package, you should place these two folders in the extracted file location of the folder "WD-BASS" (file size was too large for github)
https://drive.google.com/drive/folders/1w8e7O43tm7-AH5OrbYMe3U9msdy7_xjv?usp=sharing
You should copy the commands from "Bashrc_commands.txt" into your .bashrc file which is found in your home directory.


**Pip install command:**

1) pip install -U astropy pillow natsort matplotlib emcee corner scipy schwimmbad numba dust-extinction astroquery tqdm

2) pip install mpi4py    OR   conda install mpi4py

# Usage of the package and citation
WD-BASS was first introduced in Munday et al 2024 for use with the DBL survey and the outline of the most important technical details of the package is stated in Section 4 of this work. If there are updates, I will include these by default into WD-BASS and include a change log in the manual. 

If you use WD-BASS in your research, you must cite Munday et al 2024 (https://ui.adsabs.harvard.edu/abs/2024MNRAS.532.2534M%2F/abstract) and I ask that you please include mention of where this pipeline is found, either through a link to this github page or through citation of the Zenodo record (https://zenodo.org/records/11188044). If you wish to develop your own code and you put lines of code or subroutines from WD-BASS your own scripts, you must still cite these two references.

As a last request, if you have any suggestions on how to improve the function or efficiency of the code, if you find a bug, or you would like me to include e.g. a different synthetic spectral grid or new cooling sequences into WD-BASS, please contact me via email at james.munday98(at)gmail.com or raise an issue on Github. I have invested many hours into this package and I hope that it assists you in your research! I would love to talk about your research or project ideas with this code, so you are more than welcome to contact me about this at any time. I am happy to collaborate and fit your spectra upon request too. This is especially relevant for if you are fitting binaries (single or double lined) which are more challenging as the parameters are more degenerate.


# Currently supported model grids (and references)
- Obtainable in this package only:

DA WD, 3D-NLTE: Munday et al 2024. With mention of Tremblay et al 2013b, 2015, Kilic et al 2021, Tremblay & Bergeron 2009. See Section 5.1.1 of Munday et al 2024 for details.

- Obtainable at https://warwick.ac.uk/fac/sci/physics/research/astro/people/tremblay/modelgrids/

DA WD, 3D: Tremblay et al 2013b, 2015
DB/DBA/DC (helium-rich atmosphere): Cukanovaite et al 2021

- Gaussian/Lorentzian fitting: Nothing

- (In closed beta testing) Atmospheric models for hot subdwarfs. Details to come


# Change log (YYYYMMDD):
20240706:  Added the ability to fit the photometry only. Command run_single_phot has been added to the bashrc command list to perform this task.
20241009:  Code neatening and speeding improvements, small bug fixes

