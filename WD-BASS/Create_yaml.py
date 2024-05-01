import os, sys
sys.path.append(os.environ['WD_BASS_INSTALL_DIR']+"/scripts")
from miscAstro import miscAstro
from natsort import natsorted
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


subdwarf=True
### Change the below lines to your desired setup. Can be as many or as little as you want
if not subdwarf:
	desired_lines = [6562.79, 4861.35, 4340.472, 4101.734, 3970.075]
	resolution_lines = [6000, 4000, 3000, 2000, 1000]
	model_region = [[-90,90], [-70,70], [-60,60], [-50,50]]
	norm_region = [[-90,90], [-70,70], [-60,60], [-50,50]]
	cut_region = [[-100,100], [-80,80], [-70,70], [-60,60]]
	RV_boundaries = [[-200,200], [-200,200], [-200,200], [-200,200]]
	p0teff_min, p0teff_max = 5000, 15000
	p0logg_min, p0logg_max = 7, 9
	p0HoverHe_min, p0HoverHe_max = 0, 0
	starType1, starType2 = "DA", "DA"

else:
	desired_lines = [4861.35, 4340.472, 4101.734, 3970.075]
	resolution_lines = [5076, 4532, 4283, 4145]
	model_region = [[-70,70], [-50,50], [-40,40], [-30,30]]
	norm_region = [[-70,70], [-50,50], [-40,40], [-30,30]]
	cut_region = [[-80,80], [-60,60], [-50,50], [-40,40]]
	RV_boundaries = [[-350,350], [-350,350], [-350,350], [-350,350]]
	p0teff_min, p0teff_max = 16000, 54000
	p0logg_min, p0logg_max = 4.6, 7
	p0HoverHe_min, p0HoverHe_max = -5.05, -0.041
	starType1, starType2 = "sd", "sd"




if len(sys.argv) > 1:
    single_or_double = sys.argv[1]
    if not (single_or_double=="single" or single_or_double=="double"): raise ValueError("Use is - 'python Create_dbl_yaml.py double'  /  'python Create_dbl_yaml.py single'")
else:
    raise ValueError("Use is - 'python Create_dbl_yaml.py double'  /  'python Create_dbl_yaml.py single'")


all_names, all_refwl, all_mod, all_norm, all_cut, all_hjd, all_resolution, all_RV_boundaries, all_share_rv, all_sigma = [], [], [], [], [], [], [], [], [], []
for line in desired_lines:
	for fil in natsorted(os.listdir(os.getcwd())):
		if fil.endswith(".dat"):
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
					if line == refwl:
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





def move_minus_ones_left(arr):
    # Find indices of "-1" values
    minus_one_indices = np.where(arr == -1)[0]
    
    # Create a new array with the same shape as the input array
    new_arr = np.zeros_like(arr)
    
    # Move "-1" values to the left
    new_arr[:len(minus_one_indices)] = -1
    
    # Fill the rest of the array with non-"-1" values
    non_minus_one_indices = np.where(arr != -1)[0]
    new_arr[len(minus_one_indices):] = arr[non_minus_one_indices]
    
    # Report the indices of the input array compared to the result
    input_indices = np.arange(len(arr))
    result_indices = np.concatenate((minus_one_indices, non_minus_one_indices))
    
    return new_arr, input_indices, result_indices

# Example usage:
result, input_indices, result_indices = move_minus_ones_left(all_share_rv)


all_names, all_refwl, all_mod, all_norm, all_cut, all_hjd, all_resolution, all_share_rv = all_names[result_indices], all_refwl[result_indices], all_mod[result_indices], all_norm[result_indices], all_cut[result_indices], all_hjd[result_indices], all_resolution[result_indices], all_share_rv[result_indices]




names_share_rv_equals_minus1 = all_names[all_share_rv==-1]
for cn, (aname, an_rv) in enumerate(zip(all_names, all_share_rv)):
	if an_rv!=-1:
		all_share_rv[cn] = np.argwhere(aname==names_share_rv_equals_minus1)[0][0]
		





final_all_names, final_all_refwl, final_all_mod, final_all_norm, final_all_cut, final_all_hjd, final_all_resolution, final_all_share_rv, final_all_sigma = all_names, all_refwl, all_mod, all_norm, all_cut, all_hjd, all_resolution, all_share_rv, all_sigma







string_to_write="spectraHa: "
string_to_write+="["
for xx in final_all_names:
	string_to_write += str(xx) + ", "
string_to_write = string_to_write[:-2] + "]\n"

string_to_write+="spectra_source_type: 'wl_flux_fluxerr' # need to INSERT your source type here. Options are 'wl_flux_fluxerr', 'Xshooter'. Default is 'wl_flux_fluxerr' \n"
string_to_write+="modelHa: ["
for xx in final_all_mod:
	tempxx=str(xx)
	string_to_write += str(tempxx.split("  ")[0]+","+tempxx.split("  ")[1]) + ", "
string_to_write = string_to_write[:-2] + "]  # in angstroms"


string_to_write+="\nnormaliseHa: ["
for xx in final_all_norm:
	tempxx=str(xx)
	string_to_write += str(tempxx.split("  ")[0]+","+tempxx.split("  ")[1]) + ", "
string_to_write = string_to_write[:-2] + "]  # in angstroms"

string_to_write+="\ncut_Ha_max: ["
for xx in final_all_cut:
	tempxx=str(xx)
	string_to_write += str(tempxx.split("  ")[0]+","+tempxx.split("  ")[1]) + ", "
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
	string_to_write+="'" + starType1 + "', '" + starType2 + "'"
	string_to_write=string_to_write+"]  # 'DA' or 'DBA'"
else:
	string_to_write+="\nstarType: ['" + starType1 + "'] # 'DA' or 'DBA"



string_to_write+="\nRV_boundaries1: ["
if single_or_double=="single":
	for aval in final_all_RV_boundaries:
		tempxx=str(aval)
		string_to_write += str(tempxx.split("  ")[0]+","+tempxx.split("  ")[1]) + ", "
else:
	for aval in final_all_RV_boundaries:
		tempxx=str(aval)
		string_to_write += str(tempxx.split("  ")[0]+","+tempxx.split("  ")[1]) + ", "
string_to_write=string_to_write[:-2] + "] # in kms-1"



if single_or_double=="double":
	string_to_write+="\nRV_boundaries2: [" 
	
	for aval in final_all_RV_boundaries:
		tempxx=str(aval)
		string_to_write += str(tempxx.split("  ")[0]+","+tempxx.split("  ")[1]) + ", "
	string_to_write = string_to_write[:-2] + "] # in kms-1"



string_to_write+="\nHJD_Values: ["
for aval in final_all_hjd:
	string_to_write += str(aval)  + ", "
string_to_write = string_to_write[:-2] + "]"

RAdeg, decdeg = float(RAdeg), float(decdeg)
RAhr, dechr = miscAstro.ra_dec_deg_to_hr(RAdeg,decdeg)


if single_or_double=="double":
	string_to_write+="\nforced_teff: [0, 0]   # number or 0"
	string_to_write+="\np0teff: [["+str(p0teff_min)+","+str(p0teff_max)+"], ["+str(p0teff_min)+","+str(p0teff_max)+"]]   # number"
	string_to_write+="\nforced_logg: [0, 0]   # number or 0"
	string_to_write+="\np0logg: [["+str(p0logg_min)+","+str(p0logg_max)+"], ["+str(p0logg_min)+","+str(p0logg_max)+"]]   # number"
	string_to_write+="\nforced_HoverHe: [0, 0] # only used if DBA is in starType. number or 0"
	string_to_write+="\np0HoverHe: [["+str(p0HoverHe_min)+","+str(p0HoverHe_max)+"], ["+str(p0HoverHe_min)+","+str(p0HoverHe_max)+"]]   # number"
else:
	string_to_write+="\nforced_teff: [0]   # number or 0"
	string_to_write+="\np0teff: [["+str(p0teff_min)+","+str(p0teff_max)+"]]   # number"
	string_to_write+="\nforced_logg: [0]   # number or 0"
	string_to_write+="\np0logg: [["+str(p0logg_min)+","+str(p0logg_max)+"]]   # number"
	string_to_write+="\nforced_HoverHe: [0] # only used if DBA is in starType. number or 0"
	string_to_write+="\np0HoverHe: [["+str(p0HoverHe_min)+","+str(p0HoverHe_max)+"]]   # number"



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
if not subdwarf:
	string_to_write+='\nreddening_Ebv: "lookup" # "lookup" or number'
else:
	string_to_write+='\nreddening_Ebv: 0 # "lookup" or number'
string_to_write+="\nfile_ignore: []"
string_to_write+="\nsigma_clip: ["
for aval in final_all_sigma:
	string_to_write += str(aval)  + ", "
string_to_write = string_to_write[:-2] + "]"


string_to_write+='\nparallax: "A"'
string_to_write+='\nparallax_uncertainty: "A"'
string_to_write+='\nphot_min_val: 0'
string_to_write+='\nplot_phot_spectrum: False'



#string_to_write = string_to_write.replace("  ", ", ")

if single_or_double=="double":
	with open("example_Config_dbl.yaml", "w+") as f:
		f.write(string_to_write)
else:
	with open("example_Config_sgl.yaml", "w+") as f:
		f.write(string_to_write)




#print(string_to_write)

