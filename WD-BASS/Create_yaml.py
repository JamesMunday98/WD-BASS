import os, sys
sys.path.append(os.environ['WD_BASS_INSTALL_DIR']+"/scripts")
from miscAstro import miscAstro
from natsort import natsorted
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


### Change the below lines to your desired setup. Can be as many or as little as you want
desired_lines = [6562, 4861, 4342, 4000, 2000]
resolution_lines = [6000, 4000, 3000, 2000, 1000]
model_region = [[-90,90], [-70,70], [-60,60], [-50,50]]
norm_region = [[-90,90], [-70,70], [-60,60], [-50,50]]
cut_region = [[-100,100], [-80,80], [-70,70], [-60,60]]
RV_boundaries = [[-200,200], [-200,200], [-200,200], [-200,200]]




if len(sys.argv) > 1:
    single_or_double = sys.argv[1]
    if not (single_or_double=="single" or single_or_double=="double"): raise ValueError("Use is - 'python Create_dbl_yaml.py double'  /  'python Create_dbl_yaml.py single'")
else:
    raise ValueError("Use is - 'python Create_dbl_yaml.py double'  /  'python Create_dbl_yaml.py single'")



all_names, all_refwl, all_mod, all_norm, all_cut, all_hjd, all_resolution, all_RV_boundaries, all_share_rv, all_sigma = [], [], [], [], [], [], [], [], [], []
for fil in natsorted(os.listdir(os.getcwd())):
	if fil.endswith(".dat"):
		print(fil)
		if True:#try:
			wl, flux = np.loadtxt(fil,usecols=[0,1],unpack=True)
			line1=open(fil).readlines()[0]
			if "OrderedDict" in line1:
				HJD = line1.split("HJD")[1].split(")")[0][3:];    RAdeg = line1.split("RA")[1].split(")")[0][3:];    decdeg = line1.split("Dec")[1].split(")")[0][3:]
			elif line1.startswith("#"):
				line1=line1.split(" ")
				vals = [x for x in line1 if x]
				RAdeg=float(vals[1]);    decdeg=float(vals[2]);    HJD=float(vals[3][:-1])
			else:
				print("I can not read the HJD of file:  " + str(fi))
				print("Line 1 is:");    print(line1);    print()
				check=False
		
		
			for refwl, mod, norm, cut, res, rv in zip(desired_lines, model_region, norm_region, cut_region, resolution_lines, RV_boundaries):
				mask = (wl >= refwl+cut[0])  &  (wl <= refwl+cut[1])
				
				if len(wl[mask])>1:
					#plt.plot(wl[mask], flux[mask]);  plt.show()
					
					if fil in all_names:
						all_share_rv.append(np.argwhere(fil==np.asarray(all_names))[0][0])
					else:
						all_share_rv.append(-1)
					
					all_names.append(fil)
					all_refwl.append(refwl)
					all_mod.append(mod)
					all_norm.append(norm)
					all_cut.append(cut)
					all_hjd.append(HJD)
					all_resolution.append(res)
					all_sigma.append(4)
					all_RV_boundaries.append(rv)
					
				
		#except Exception as e: print("not liked, ", e)


all_names, all_refwl, all_mod, all_norm, all_cut, all_hjd, all_resolution, final_all_RV_boundaries, all_share_rv, all_sigma  =  np.asarray(all_names), np.asarray(all_refwl), np.asarray(all_mod), np.asarray(all_norm), np.asarray(all_cut), np.asarray(all_hjd), np.asarray(all_resolution), np.asarray(all_RV_boundaries), np.asarray(all_share_rv), np.asarray(all_sigma)



newargs = np.argsort(all_share_rv)

all_names, all_refwl, all_mod, all_norm, all_cut, all_hjd, all_resolution, all_share_rv = all_names[newargs], all_refwl[newargs], all_mod[newargs], all_norm[newargs], all_cut[newargs], all_hjd[newargs], all_resolution[newargs], all_share_rv[newargs]




mask = np.argsort(all_refwl[all_share_rv!=-1])
all_names[all_share_rv!=-1],       all_refwl[all_share_rv!=-1]     =    all_names[all_share_rv!=-1][mask],        all_refwl[all_share_rv!=-1][mask]
all_mod[all_share_rv!=-1],         all_norm[all_share_rv!=-1]      =    all_mod[all_share_rv!=-1][mask],          all_norm[all_share_rv!=-1][mask]
all_cut[all_share_rv!=-1],         all_hjd[all_share_rv!=-1]       =    all_cut[all_share_rv!=-1][mask],          all_hjd[all_share_rv!=-1][mask]
all_resolution[all_share_rv!=-1],  all_share_rv[all_share_rv!=-1]  =    all_resolution[all_share_rv!=-1][mask],   all_share_rv[all_share_rv!=-1][mask]





final_all_names, final_all_refwl, final_all_hjd, final_all_resolution, final_all_share_rv, final_all_sigma  = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]) #, np.array([]), np.array([]), np.array([])


for re in np.unique(all_refwl)[::-1]:
	mask=all_refwl==re
	final_all_names = np.append(final_all_names, all_names[mask])
	final_all_refwl = np.append(final_all_refwl, all_refwl[mask])
	final_all_hjd = np.append(final_all_hjd, all_hjd[mask])
	final_all_resolution = np.append(final_all_resolution, all_resolution[mask])
	final_all_share_rv = np.append(final_all_share_rv, all_share_rv[mask])
	final_all_sigma = np.append(final_all_sigma, all_sigma[mask])

final_all_mod, final_all_norm, final_all_cut  =  all_mod, all_norm, all_cut


for afile in np.unique(final_all_names):
	minloc=np.amin(np.argwhere(final_all_names==afile))
	final_all_share_rv[(final_all_names==afile) & (final_all_share_rv!=-1)] = minloc

final_all_sharerv = final_all_share_rv.astype(int)








string_to_write="spectraHa: "
string_to_write+="["
for xx in final_all_names:
	string_to_write += str(xx) + ", "
string_to_write = string_to_write[:-2] + "]\n"

string_to_write+="spectra_source_type: 'wl_flux_fluxerr' # need to INSERT your source type here. Options are 'wl_flux_fluxerr', 'Xshooter'. Default is 'wl_flux_fluxerr' \n"
string_to_write+="modelHa: ["
for xx in final_all_mod:
	string_to_write += str(xx) + ", "
string_to_write = string_to_write[:-2] + "]  # in angstroms"


string_to_write+="\nnormaliseHa: ["
for xx in final_all_norm:
	string_to_write += str(xx) + ", "
string_to_write = string_to_write[:-2] + "]  # in angstroms"

string_to_write+="\ncut_Ha_max: ["
for xx in final_all_cut:
	string_to_write += str(xx) + ", "
string_to_write = string_to_write[:-2] + "]  # in angstroms"

string_to_write+="\nresolutions: ["
for xx in final_all_resolution:
	string_to_write += str(xx) + ", "
string_to_write = string_to_write[:-2] + "] # R=lambda/dlambda"


string_to_write+="\nrefwave: [" 
for xx in final_all_refwl:
	string_to_write += str(xx) + ", "
string_to_write = string_to_write[:-2] + "]"

string_to_write+="\nshare_rv: ["
for aval in final_all_share_rv:
	string_to_write += str(int(aval)) + ", "
string_to_write = string_to_write[:-2]  +  "]  # 0 indexed!  needs all shared rvs to be at the end of the line, otherwise won't work properly!!!!!!"

if single_or_double=="double":
	string_to_write+="\nstarType: ["
	for i in range(2):
		string_to_write+='"DA", '
	string_to_write=string_to_write[:-2]+"]  # 'DA' or 'DBA'"
else:
	string_to_write+="\nstarType: ['DA'] # 'DA' or 'DBA"



string_to_write+="\nRV_boundaries1: ["
if single_or_double=="single":
	for aval in final_all_RV_boundaries:
		string_to_write += str(aval) + ", "
else:
	for aval in final_all_RV_boundaries:
		string_to_write += str(aval)  + ", "
string_to_write=string_to_write[:-2] + "] # in kms-1"



if single_or_double=="double":
	string_to_write+="\nRV_boundaries2: [" 
	
	for aval in final_all_RV_boundaries:
		string_to_write += str(aval)  + ", "
	string_to_write = string_to_write[:-2] + "] # in kms-1"



string_to_write+="\nHJD_Values: ["
for aval in final_all_hjd:
	string_to_write += str(aval)  + ", "
string_to_write = string_to_write[:-2] + "]"

RAdeg, decdeg = float(RAdeg), float(decdeg)
RAhr, dechr = miscAstro.ra_dec_deg_to_hr(RAdeg,decdeg)


if single_or_double=="double":
	string_to_write+="\nforced_teff: [0, 0]   # number or 0"
	string_to_write+="\np0teff: [[5000,15000], [5000,15000]]   # number"
	string_to_write+="\nforced_logg: [0, 0]   # number or 0"
	string_to_write+="\np0logg: [[7,9], [7,9]]   # number"
	string_to_write+="\nforced_HoverHe: [0, 0] # only used if DBA is in starType. number or 0"
	string_to_write+="\np0HoverHe: [[0,0], [0,0]]   # number"
else:
	string_to_write+="\nforced_teff: [0]   # number or 0"
	string_to_write+="\np0teff: [[5000,15000]]   # number"
	string_to_write+="\nforced_logg: [0]   # number or 0"
	string_to_write+="\np0logg: [[7,9]]   # number"
	string_to_write+="\nforced_HoverHe: [0] # only used if DBA is in starType. number or 0"
	string_to_write+="\np0HoverHe: [[0,0]]   # number"



string_to_write+="\nforced_scaling: ['WD'] # False, 'WD' or float. If false, will be added to the mcmc parameters"
string_to_write+="\nplot_corner: True # True or False"
string_to_write+="\nforced_ephemeris: [False] # expects [False] or [T0, P0]  (in days)"
string_to_write+="\nforced_K1: [False] # in kms-1. Options [False], [float], ['Fit']"
string_to_write+="\nK1_boundaries: [-300,300] # in kms-1"
string_to_write+="\nforced_Vgamma1: [False] # in kms-1. Options [float], ['Fit']"
string_to_write+="\nVgamma1_boundaries: [-100,100] # in kms-1"
if single_or_double=="double":
	string_to_write+="\nforced_K2: [False] # in kms-1. Options [False], [float], ['Fit']"
	string_to_write+="\nK2_boundaries: [-300,300] # in kms-1"
	string_to_write+="\nforced_Vgamma2: [False] # False, 'Fit' or number"
	string_to_write+="\nVgamma2_boundaries: [-100,100] # in kms-1"
string_to_write+="\nplot_fit: [False] # [False] or e.g. [0, 1, 4] for spectra 1 2 5"
string_to_write+="\nfit_phot_SED: False"
string_to_write+='\nRA: "' + RAhr+'"  # only used when fit_phot_SED is True. Hexadesimal string or a float. Set to None if not necessary'
string_to_write+='\nDec: "' + dechr + '" # only used when fit_phot_SED is True. Hexadesimal string or a float. Set to None if not necessary'
string_to_write+='\nexpected_Gmag: 15 # only used when fit_phot_SED is True. Integer/float/None'
string_to_write+='\nnwalkers: 50'
string_to_write+='\nnburnin: 100'
string_to_write+='\nnsteps: 50'
string_to_write+='\nreddening_Ebv: "lookup" # "lookup" or number'
string_to_write+="\nfile_ignore: []"
string_to_write+="\nsigma_clip: ["
for aval in final_all_sigma:
	string_to_write += str(aval)  + ", "
string_to_write = string_to_write[:-2] + "]"


string_to_write+='\nparallax: "A"'
string_to_write+='\nparallax_uncertainty: "A"'
string_to_write+='\nphot_min_val: 0'
string_to_write+='\nplot_phot_spectrum: False'



string_to_write = string_to_write.replace("  ", ", ")

if single_or_double=="double":
	with open("example_Config_dbl.yaml", "w+") as f:
		f.write(string_to_write)
else:
	with open("example_Config_sgl.yaml", "w+") as f:
		f.write(string_to_write)




print(string_to_write)

