import numpy as np
from numpy import interp, polyfit
import os, sys, yaml, requests, warnings, emcee, corner
sys.path.append(os.environ['WD_BASS_INSTALL_DIR']+"/scripts")
from load_All_Models import load_models
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from schwimmbad import MPIPool
from importable_MTR_function import get_MTR, get_MTR_DB
warnings.filterwarnings('ignore')
from Fit_Photometric_SED import Fit_phot
from checkLocalStars import checkLocalStars
from miscAstro import miscAstro
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.convolution import Gaussian1DKernel,convolve
from numba import njit, jit
from scipy.optimize import curve_fit
from io import StringIO
import pandas as pd


with open('example_Config_sgl.yaml') as file:
    config_info = yaml.load(file, Loader=yaml.FullLoader)

#with open('example_Config_dbl.yaml') as file:
#    config_info = yaml.load(file, Loader=yaml.FullLoader)

speed_of_light = 2.99792458E5

spectra_source_type = config_info["spectra_source_type"]
if spectra_source_type!="wl_flux_fluxerr" and spectra_source_type!="Xshooter": raise ValueError("spectra_source_type is not recognised in the config file")
input_files = np.asarray(config_info["spectraHa"])
reference_wl = np.asarray(config_info["refwave"])
modelHa = np.asarray(config_info["modelHa"])
normaliseHa_all = np.asarray(config_info["normaliseHa"])
cut_Ha_all = np.asarray(config_info["cut_Ha_max"])
resolutions = np.asarray(config_info["resolutions"])
#resolutions*=2
#print("WARNING WARNING resolution doubled")
share_rv = np.asarray(config_info["share_rv"])
starType1 = np.asarray(config_info["starType"])[0]
RV_boundaries1 = np.asarray(config_info["RV_boundaries1"])
HJD_values = np.asarray(config_info["HJD_Values"])
forced_teff1 = np.asarray(config_info["forced_teff"])
forced_logg1 = np.asarray(config_info["forced_logg"])
forced_HoverHe1 = np.asarray(config_info["forced_HoverHe"])
plot_corner=np.asarray(config_info["plot_corner"])
forced_Scaling=np.asarray(config_info["forced_scaling"])[0]
if forced_Scaling!=False and (forced_Scaling!="WD" and forced_Scaling!="sd"): forced_Scaling=float(forced_Scaling)
forced_ephemeris=np.asarray(config_info["forced_ephemeris"]) # expects False or [T0, P0]
if forced_ephemeris[0]==False:  forced_ephemeris=forced_ephemeris[0];  forced_T0=None;  forced_P0=None
else: forced_T0=float(forced_ephemeris[0]);  forced_P0=float(forced_ephemeris[1])
forced_K1=config_info["forced_K1"][0]
if not forced_K1==False and not forced_K1=="Fit":  forced_K1=float(forced_K1)
K1_boundaries=np.asarray(config_info["K1_boundaries"]);             p0K1=np.array([float(K1_boundaries[0]), float(K1_boundaries[1])])
forced_Vgamma1=config_info["forced_Vgamma1"][0]
if not forced_Vgamma1=="Fit":  forced_Vgamma1=float(forced_Vgamma1)
Vgamma1_boundaries=np.asarray(config_info["Vgamma1_boundaries"]);   p0Vgamma1=np.array([float(Vgamma1_boundaries[0]), float(Vgamma1_boundaries[1])])
plot_fit=np.asarray([config_info["plot_fit"]])
fit_phot_SED=np.asarray([config_info["fit_phot_SED"]])
RA=np.asarray([config_info["RA"]])[0]
Dec=np.asarray([config_info["Dec"]])[0]


if len(sys.argv)>2 and "one_by_oneMCMC" in sys.argv[2]:
	with open("out/result.out") as resultfile:
		result_lines = resultfile.readlines()
		for iii, lll in enumerate(result_lines):
			if starType1.startswith("D") or starType1.startswith("sd"):
				if "Teff1:" in lll:       
					T1_med=float(result_lines[iii+1].split("\t")[0])
					if forced_teff1==0:      forced_teff1=[T1_med]
				elif "Logg1:" in lll:
					logg1_med=float(result_lines[iii+1].split("\t")[0])
					if forced_logg1==0:      forced_logg1=[logg1_med]
				elif "H/He1:" in lll:
					HoverHe1_med=float(result_lines[iii+1].split("\t")[0])
					if forced_HoverHe1==0:   forced_HoverHe1=[HoverHe1_med]
				elif "Scaling:" in lll:
					Scaling=float(result_lines[iii+1].split("\t")[0])
					if forced_Scaling==False: forced_Scaling=[Scaling]
			else:
				if "A1_1:" in lll:        A1_1_med=float(result_lines[iii+1].split("\t")[0])
				elif "sig1_1:" in lll:    sig1_1_med=float(result_lines[iii+1].split("\t")[0])
				elif "A1_2:" in lll:      A1_2_med=float(result_lines[iii+1].split("\t")[0])
				elif "sig1_2:" in lll:    sig1_2_med=float(result_lines[iii+1].split("\t")[0])
	
	fit_phot_SED=False
	




## gaussian/lorentzian things
if "GG" in starType1 or "LL" in starType1 or "GL" in starType1:
	try:  
		A1_1_boundaries = np.asarray([config_info["A1_1_boundaries"]])[0]
		A1_2_boundaries = np.asarray([config_info["A1_2_boundaries"]])[0]
		sigma1_1_boundaries = np.asarray([config_info["sigma1_1_boundaries"]])[0]
		sigma1_2_boundaries = np.asarray([config_info["sigma1_2_boundaries"]])[0]
	except: raise ValueError("Missing boundaries for star1: A1_1_boundaries/sigma1_1_boundaries/A1_2_boundaries/sigma1_2_boundaries")


if type(RA)==np.str_ and type(Dec)==np.str_:  RAdeg, Decdeg= miscAstro.ra_dec_hr_to_deg(RA,Dec)
elif type(RA)==float and type(Dec)==float:  RAdeg, Decdeg= RA, Dec;  RA, Dec = miscAstro.ra_dec_deg_to_hr(RAdeg,Decdeg)
else: print(type(RA), type(Dec));  raise ValueError

expected_Gmag=np.asarray([config_info["expected_Gmag"]])[0]
nwalkers=np.asarray([config_info["nwalkers"]])[0]
burnin=np.asarray([config_info["nburnin"]])[0]
nsteps=np.asarray([config_info["nsteps"]])[0]
try:  reddening_Ebv=float(np.asarray([config_info["reddening_Ebv"]])[0])
except: reddening_Ebv=np.asarray([config_info["reddening_Ebv"]])[0]

if reddening_Ebv=="lookup":
	try:
		star_plx=float(np.asarray(config_info["parallax"])) # in mas
		star_e_plx=float(np.asarray(config_info["parallax_uncertainty"])) # in mas
	except:
		RAdeg, Decdeg, star_mag, star_plx, star_e_plx, warning = checkLocalStars.find_star_in_gaia_edr3(RAdeg, Decdeg, 10, predicted_Gmag=expected_Gmag)
		#raise ValueError(star_mag)
		#print(star_mag.value[0], star_plx.value[0], star_e_plx.value[0], warning)
	c = SkyCoord(ra=RAdeg*u.degree, dec=Decdeg*u.degree, frame='icrs')
	gal_coords= c.galactic
	gal_l=gal_coords.l.value; gal_b=gal_coords.b.value
	
	#complete_with_stilism("Stilism.csv")
	url = "http://stilism.obspm.fr/reddening?frame=galactic&vlong={}&ulong=deg&vlat={}&ulat=deg&distance={}"
	
	res = requests.get(url.format(gal_l, gal_b, 1000/star_plx), allow_redirects=True)
	files = StringIO(res.content.decode("utf-8"))
	dfstilism = pd.read_csv(files)
	dist_stilism = dfstilism["distance[pc]"][0]
	red_stilism = dfstilism["reddening[mag]"][0]
	dist_uncstilism = dfstilism["distance_uncertainty[pc]"][0]
	red_unc_min = dfstilism["reddening_uncertainty_min[mag]"][0]
	red_unc_max = dfstilism["reddening_uncertainty_max[mag]"][0]
	
	reddening_Ebv=red_stilism
else:
	try: float(reddening_Ebv)
	except Exception as e: print(e); reddening_Ebv=0
try:  sigma_clip=np.asarray(config_info["sigma_clip"])
except: sigma_clip=np.full((len(modelHa),),-1)
try: plot_burnin=np.asarray(config_info["plot_burnin"])
except: plot_burnin=False
try: want_gaiadr3=np.asarray(config_info["want_gaiadr3"])
except: want_gaiadr3=True
if fit_phot_SED:
	try:
		plax=float(np.asarray(config_info["parallax"])) # in mas
		plax_unc=float(np.asarray(config_info["parallax_uncertainty"])) # in mas
		found_parallax=True
	except:
		RAdeg, Decdeg, star_mag, star_plx, star_e_plx, warning = checkLocalStars.find_star_in_gaia_edr3(RAdeg, Decdeg, 10, predicted_Gmag=expected_Gmag)
		try: star_mag, plax, plax_unc = star_mag.value, star_plx.value, star_e_plx.value
		except: star_mag, plax, plax_unc = star_mag, star_plx, star_e_plx
		
		try: 
			float(plax); found_parallax=True
		except: found_parallax=False
	try:  phot_min_val=float(np.asarray(config_info["phot_min_val"]))
	except: phot_min_val=0
	try:  want_wise1=np.asarray(config_info["want_wise1"])
	except: want_wise1=False
	try:  want_wise2=np.asarray(config_info["want_wise2"])
	except: want_wise2=False
	try:  want_2mass=np.asarray(config_info["want_2mass"])
	except: want_2mass=False
	try: ignore_filt = np.asarray(config_info["ignore_filt"])
	except: ignore_filt = []
	try: want_galex = np.asarray(config_info["want_galex"])
	except: want_galex=False
		
		
else: found_parallax=False


sys_args = sys.argv

if sys.argv[1]=="ATM":

	file_ignore=np.asarray([config_info["file_ignore"]])[0]

	if len(file_ignore) > 0:  dowhile=True
	else: dowhile=False
	
	
	if len(sys_args)==4 and "one_by_oneMCMC" in sys.argv[2]:
		list_indexes=[]
		for i in range(len(share_rv)):
			list_indexes.append(i)
		list_indexes=np.asarray(list_indexes)
		

	while dowhile:
		stopwhile=True
		original_input_files = input_files
		original_share_rv = share_rv
		for aaaa in file_ignore:
			try:
				original_file_index = np.argwhere(original_input_files==aaaa)[0][0]
			except: raise ValueError(aaaa)
			for zcnt, aaa_sharerv in enumerate(original_share_rv):
				if aaa_sharerv>=original_file_index:
					share_rv[zcnt]-=1
				
				if aaa_sharerv>=original_file_index and len(sys_args)==4 and "one_by_oneMCMC" in sys.argv[2]:    list_indexes[zcnt]-=1
			
			
			#print(len(normaliseHa_all),len(cut_Ha_all),len(resolutions),len(reference_wl),len(RV_boundaries1),len(RV_boundaries2), len(HJD_values), len(input_files))
			mask_ignore_files = input_files!=aaaa
			modelHa = modelHa[mask_ignore_files]
			normaliseHa_all=normaliseHa_all[mask_ignore_files]
			cut_Ha_all=cut_Ha_all[mask_ignore_files]
			resolutions=resolutions[mask_ignore_files]
			reference_wl=reference_wl[mask_ignore_files]
			share_rv=share_rv[mask_ignore_files]
			RV_boundaries1=RV_boundaries1[mask_ignore_files]
			HJD_values=HJD_values[mask_ignore_files]
			input_files=input_files[mask_ignore_files]
			sigma_clip=sigma_clip[mask_ignore_files]
			if len(sys_args)==4 and "one_by_oneMCMC" in sys.argv[2]:
				list_indexes=list_indexes[mask_ignore_files]
				removed_minus1_pos = np.argwhere(mask_ignore_files==False)[0][0]
				list_indexes[removed_minus1_pos:]-=1
			
			file_ignore_mask=file_ignore!=aaaa
			file_ignore=file_ignore[file_ignore_mask]
			
			
			stopwhile=False
			break
		if stopwhile: dowhile=False
		
	
	if len(sys_args)==4 and "one_by_oneMCMC" in sys.argv[2]:
		wanted_index = (share_rv==int(sys.argv[3])) | (list_indexes==int(sys.argv[3]))
		modelHa = modelHa[wanted_index]
		normaliseHa_all=normaliseHa_all[wanted_index]
		cut_Ha_all=cut_Ha_all[wanted_index]
		resolutions=resolutions[wanted_index]
		reference_wl=reference_wl[wanted_index]
		share_rv=share_rv[wanted_index]
		RV_boundaries1=RV_boundaries1[wanted_index]
		HJD_values=HJD_values[wanted_index]
		input_files=input_files[wanted_index]
		sigma_clip=sigma_clip[wanted_index]


if len(sys_args)==4 and "one_by_oneMCMC" in sys.argv[2]:
	while True:
		if np.amin(share_rv[share_rv!=-1])!=0:   share_rv[share_rv!=-1]-=1
		else:   break


if len(share_rv)==0:
	print("No RVs passed. Exiting script")
	sys.exit()
	
	
	


print("Processing :   ", len(input_files), "   spectra")
for afi, aha, anorm, acut in zip(input_files, modelHa, normaliseHa_all, cut_Ha_all):
	print(afi, aha, anorm, acut)
try: os.mkdir("out")
except: None

if fit_phot_SED:
	if type(RA)==np.str_ and type(Dec)==np.str_:  RAdeg, Decdeg= miscAstro.ra_dec_hr_to_deg(RA,Dec)
	elif type(RA)==float and type(Dec)==float:  RAdeg, Decdeg= RA, Dec;  RA, Dec = miscAstro.ra_dec_deg_to_hr(RAdeg,Decdeg)
	else: print(type(RA), type(Dec));  raise ValueError
	#checkLocalStars.find_star_in_gaia_edr3(RAdeg, Decdeg, 5, predicted_Gmag=expected_Gmag)


if (starType1 != "DA" and starType1 != "DBA" and starType1!="DB") and forced_Scaling == "WD" and fit_phot_SED:   raise ValueError("WD scaling only allowed for DA stars")
if (forced_K1=="Fit") or ((forced_K1==False and not (forced_Vgamma1=="Fit" or isinstance(forced_Vgamma1, float)))):  raise ValueError("Error - knowing the ephemeris and fitting K1 + RV2s is weird. Not allowed")
if (forced_K1==False and forced_Vgamma1=="Fit"):  raise ValueError("Can't fit Vgamma with K free")



if not len(input_files)==len(reference_wl)==len(normaliseHa_all)==len(cut_Ha_all)==len(resolutions)==len(share_rv)==len(RV_boundaries1)==len(sigma_clip)==len(HJD_values):
	print("ALERT");  print("ALERT"); print("ALERT"); print("ALERT"); print("ALERT")
	print("Mismatch of length of required arrays")
	print("input_files", len(input_files));  print("reference_wl", len(reference_wl))
	print("normaliseHa", len(normaliseHa_all));  print("cut_Ha_all", len(cut_Ha_all))
	print("resolutions", len(resolutions));  print("share_rv", len(share_rv))
	print("RV_boundaries1", len(RV_boundaries1))
	print("HJD_values", len(HJD_values))
	print("sigma_clip", len(sigma_clip))
	raise ValueError


theminww_loadgrid, themaxww_loadgrid=1,1500000
if fit_phot_SED:
	filterdir=os.environ['WD_BASS_INSTALL_DIR']+"/Filters"
	filter_dict = {"GAIA/GAIA3:Grp":np.array([np.loadtxt(filterdir+"/GAIA_GAIA3.Grp.dat",unpack=True), 7619.96, "PHOTON"], dtype=object),
		"GAIA/GAIA3:Gbp":np.array([np.loadtxt(filterdir+"/GAIA_GAIA3.Gbp.dat",unpack=True), 5035.75, "PHOTON"], dtype=object), 
		"GAIA/GAIA3:G": np.array([np.loadtxt(filterdir+"/GAIA_GAIA3.G.dat",unpack=True), 5822.39, "PHOTON"], dtype=object),
		"SDSS:z": np.array([np.loadtxt(filterdir+"/SLOAN_SDSS.z.dat",unpack=True), 8922.78, "PHOTON"], dtype=object),
		"SDSS:i": np.array([np.loadtxt(filterdir+"/SLOAN_SDSS.i.dat",unpack=True), 7457.89, "PHOTON"], dtype=object),
		"SDSS:r": np.array([np.loadtxt(filterdir+"/SLOAN_SDSS.r.dat",unpack=True), 6141.12, "PHOTON"], dtype=object),
		"SDSS:g": np.array([np.loadtxt(filterdir+"/SLOAN_SDSS.g.dat",unpack=True), 4671.78, "PHOTON"], dtype=object),
		"SDSS:u": np.array([np.loadtxt(filterdir+"/SLOAN_SDSS.u.dat",unpack=True), 3608.04, "PHOTON"], dtype=object),
		"PAN-STARRS/PS1:y": np.array([np.loadtxt(filterdir+"/PAN-STARRS_PS1.y.dat",unpack=True),9613.5, "PHOTON"], dtype=object),
		"PAN-STARRS/PS1:z": np.array([np.loadtxt(filterdir+"/PAN-STARRS_PS1.z.dat",unpack=True),8668.5, "PHOTON"], dtype=object),
		"PAN-STARRS/PS1:i": np.array([np.loadtxt(filterdir+"/PAN-STARRS_PS1.i.dat",unpack=True), 7503.03, "PHOTON"], dtype=object),
		"PAN-STARRS/PS1:r": np.array([np.loadtxt(filterdir+"/PAN-STARRS_PS1.r.dat",unpack=True), 6156.4, "PHOTON"], dtype=object),
		"PAN-STARRS/PS1:g": np.array([np.loadtxt(filterdir+"/PAN-STARRS_PS1.g.dat",unpack=True),4810.9, "PHOTON"], dtype=object),
		"Johnson:V": np.array([np.loadtxt(filterdir+"/Generic_Johnson.V.dat",unpack=True),5466.1, "ENERGY"], dtype=object),
		#"BarrosU":np.array([np.loadtxt(filterdir+"/WHT_ULTRACAM.u.dat",unpack=True), 3481.95, "PHOTON"], dtype=object),
		#"BarrosG":np.array([np.loadtxt(filterdir+"/WHT_ULTRACAM.g.dat",unpack=True), 4762.3, "PHOTON"], dtype=object),
		#"BarrosR":np.array([np.loadtxt(filterdir+"/WHT_ULTRACAM.r.dat",unpack=True), 6256.2, "PHOTON"], dtype=object),
		#"BarrosI":np.array([np.loadtxt(filterdir+"/WHT_ULTRACAM.i.dat",unpack=True), 7585.9, "PHOTON"], dtype=object),
		"HST_F775W":np.array([np.loadtxt(filterdir+"/HST_WFC3_UVIS1.F775W.dat",unpack=True), 7612.8, "PHOTON"], dtype=object),
		"HST_F390W":np.array([np.loadtxt(filterdir+"/HST_WFC3_UVIS1.F390W.dat",unpack=True), 4022.2, "PHOTON"], dtype=object),
		"HST_F225W":np.array([np.loadtxt(filterdir+"/HST_WFC3_UVIS1.F225W.dat",unpack=True), 2372.8, "PHOTON"], dtype=object),
		"XMM-OT:UVW1":np.array([np.loadtxt(filterdir+"/XMM_OM.UVW1_filter.dat",unpack=True), 2971.0, "PHOTON"], dtype=object),
		"SWIFT:U":np.array([np.loadtxt(filterdir+"/Swift_UVOT.U_trn.dat",unpack=True), 3520.95, "PHOTON"], dtype=object),
		"GALEX:NUV":np.array([np.loadtxt(filterdir+"/GALEX_GALEX.NUV.dat",unpack=True), 2303.37, "PHOTON"], dtype=object),
		"GALEX:FUV":np.array([np.loadtxt(filterdir+"/GALEX_GALEX.FUV.dat",unpack=True), 1548.85, "PHOTON"], dtype=object),
		"GAIA/GAIA2:Grp":np.array([np.loadtxt(filterdir+"/GAIA_GAIA2.Grp.dat",unpack=True), 7592.04, "ENERGY"], dtype=object),
		"GAIA/GAIA2:Gbp":np.array([np.loadtxt(filterdir+"/GAIA_GAIA2.Gbp.dat",unpack=True), 5014.05, "ENERGY"], dtype=object), 
		"GAIA/GAIA2:G": np.array([np.loadtxt(filterdir+"/GAIA_GAIA2.G.dat",unpack=True), 5838.81, "ENERGY"], dtype=object), 
		"XMM-OT:UVM2": np.array([np.loadtxt(filterdir+"/XMM_OM.UVM2_filter.dat",unpack=True), 2329.54, "PHOTON"], dtype=object),
		"XMM-OT:UVW2": np.array([np.loadtxt(filterdir+"/XMM_OM.UVW2_filter.dat",unpack=True), 2143.84, "PHOTON"], dtype=object),
		"UVOT:UVW2": np.array([np.loadtxt(filterdir+"/XMM_OM.UVW2_filter.dat",unpack=True), 2143.84, "PHOTON"], dtype=object),
		"UVOT:UVW1": np.array([np.loadtxt(filterdir+"/XMM_OM.UVW1_filter.dat",unpack=True), 2971.72, "PHOTON"], dtype=object),
		"XMM-OT:U": np.array([np.loadtxt(filterdir+"/XMM_OM.U_filter.dat",unpack=True), 3515.73, "PHOTON"], dtype=object),
		"XMM-OT:V": np.array([np.loadtxt(filterdir+"/XMM_OM.V_filter.dat",unpack=True), 5437.81, "PHOTON"], dtype=object),
		"Johnson:B": np.array([np.loadtxt(filterdir+"/Generic_Johnson.B.dat",unpack=True), 4369.53, "ENERGY"], dtype=object),
		"WISE:W1": np.array([np.loadtxt(filterdir+"/WISE_WISE.W1.dat",unpack=True), 33154.27, "ENERGY"], dtype=object),
		"WISE:W2": np.array([np.loadtxt(filterdir+"/WISE_WISE.W2.dat",unpack=True), 45644.77, "ENERGY"], dtype=object),
		"WISE:W3": np.array([np.loadtxt(filterdir+"/WISE_WISE.W3.dat",unpack=True), 107866.13, "ENERGY"], dtype=object),
		"2MASS:J": np.array([np.loadtxt(filterdir+"/2MASS_2MASS.J.dat",unpack=True), 12285.64, "ENERGY"], dtype=object),
		"2MASS:H": np.array([np.loadtxt(filterdir+"/2MASS_2MASS.H.dat",unpack=True), 16385.40, "ENERGY"], dtype=object),
		"2MASS:Ks": np.array([np.loadtxt(filterdir+"/2MASS_2MASS.Ks.dat",unpack=True), 21521.61, "ENERGY"], dtype=object)
		}





	url = Fit_phot.get_url_CDS(RA + " " + Dec,5)
	sedfilter, sed_wl, sedflux, sedfluxe = np.array([]), np.array([]), np.array([]), np.array([])
	try:
		sedfilter, sed_wl, sedflux, sedfluxe, catalogue = Fit_phot.get_data(url, False)
		if len(sedfilter)<=3:    url2 = Fit_phot.get_url_CDS(RA + " " + Dec,10);    sedfilter2, sed_wl2, sedflux2, sedfluxe2, catalogue2 = Fit_phot.get_data(url, False)
		try: 
			if len(sedfilter2)>len(sedfilter):    sedfilter, sed_wl, sedflux, sedfluxe, catalogue = sedfilter2, sed_wl2, sedflux2, sedfluxe2, catalogue2
		except: None
		
		
		if len(sedfilter)>=1:
			list_gaia=[]
			if "II/328/allwise" in catalogue: allwise_found=True
			else: allwise_found=False
			for cnt, (ff, cat) in enumerate(zip(sedfilter, catalogue)):
				if not ff in ignore_filt:
					if "SDSS" in ff:
						if str(cat)=="J/MNRAS/508/3877/maincat": list_gaia.append(True)
						else: list_gaia.append(False)
					elif "PAN" in ff:
						if cat=="II/349/ps1": list_gaia.append(True)
						else: list_gaia.append(False)
					elif "GAIA3" in ff:
						if cat=="I/350/gaiaedr3" and want_gaiadr3:  list_gaia.append(True)
						#if cat=="I/350/gaiaedr3" and ("GAIA3:Gbp" in ff or "GAIA3:Grp" in ff):  list_gaia.append(True)
						#elif cat=="I/350/gaiaedr3" and want_gaiadr3==True:  list_gaia.append(True)
						else: list_gaia.append(False)
					elif "GALEX:NUV" in ff and want_galex:
						if cat=="II/312/ais": list_gaia.append(True)
						else: list_gaia.append(False)
					elif want_wise1 and "WISE" in ff and "allwise" in cat and "W1" in ff: list_gaia.append(True)
					elif want_wise1 and allwise_found==False and "WISE" in ff and str(cat)=="II/311/wise" and "W1" in ff: list_gaia.append(True)
					elif want_wise2 and "WISE" in ff and "allwise" in cat and "W2" in ff: list_gaia.append(True)
					elif want_wise2 and allwise_found==False and "WISE" in ff and str(cat)=="II/311/wise" and "W2" in ff: list_gaia.append(True)
					elif want_2mass and "2MASS" in ff and not str(cat)=="J/MNRAS/472/4173/table1": list_gaia.append(True)
					#elif "GAIA" in ff:   list_gaia.append(False)
					#elif "XMM" in ff or "UVOT" in ff or "GALEX" in ff: list_gaia.append(False)
					else:   list_gaia.append(False)
				else: list_gaia.append(False)
			
			
			list_gaia=np.asarray(list_gaia)
			
			
			if len(sedfilter[list_gaia])==0:
				for cnt, (ff, cat) in enumerate(zip(sedfilter, catalogue)):
					if not ff in ignore_filt:
						if "SDSS" in ff and str(cat)=="IV/38/tic":
							list_gaia[cnt]=True
			
			sedfilter, sed_wl, sedflux, sedfluxe = sedfilter[list_gaia], sed_wl[list_gaia], sedflux[list_gaia], sedfluxe[list_gaia]
			
			
			
			
			
			masksome=sedflux>phot_min_val
			sedfilter, sed_wl, sedflux, sedfluxe = sedfilter[masksome], sed_wl[masksome], sedflux[masksome], sedfluxe[masksome]
			
			
			if want_2mass:
				mask_nan_fluxes = ( ~((sedfilter=="2MASS:H") | (sedfilter=="2MASS:J") | (sedfilter=="2MASS:Ks"))  ) | (((sedfilter=="2MASS:H") | (sedfilter=="2MASS:J") | (sedfilter=="2MASS:Ks"))  &    (~np.isnan(sedfluxe.astype(float))))
				sedfilter, sed_wl, sedflux, sedfluxe = sedfilter[mask_nan_fluxes], sed_wl[mask_nan_fluxes], sedflux[mask_nan_fluxes], sedfluxe[mask_nan_fluxes]
				
				
				for filtername in ["2MASS:H", "2MASS:J", "2MASS:Ks"]:
					if not filtername in ignore_filt:
						mask1=sedfilter==filtername
						m_wl, m_f, m_fe = np.mean(sed_wl[mask1]), np.mean(sedflux[mask1]), np.mean(sedfluxe[mask1])
						sedfilter, sed_wl, sedflux, sedfluxe = sedfilter[~mask1], sed_wl[~mask1], sedflux[~mask1], sedfluxe[~mask1]
						if m_f!=0 and not np.isnan(m_fe):
							sedfilter, sed_wl, sedflux, sedfluxe = np.append(sedfilter, filtername), np.append(sed_wl, m_wl), np.append(sedflux, m_f), np.append(sedfluxe, m_fe)
				
			
			sedfluxe[sedfluxe==0] = np.amin(sedfluxe[sedfluxe!=0])

			
			#sedflux = sedflux / 1E23;   sedfluxe = sedfluxe / 1E23
			sed_wl*=10
		
		else:
			theminww_loadgrid, themaxww_loadgrid=1000,15000
			fit_phot_SED=False
			print("ALERT:")
			print("No available photometry of the source in this RA and Dec in a 5'' radius")
		
		try:
			External_AB_mag = np.asarray(config_info["External_AB_mag"]).astype(float)
			External_AB_mag_err = np.asarray(config_info["External_AB_mag_err"]).astype(float)
			External_Flux_Jy=10**(23-(External_AB_mag+48.594)/2.5)  #  erg/cm2/s/Hz
			External_Flux_Jy_err=np.sqrt(  ((-(0.9210340371976184*10**((External_AB_mag+243/5)*-0.4+23)))*External_AB_mag_err)**2   )
			External_Filter = np.asarray(config_info["External_Filter"])
			External_wl = []
			for an_external_filt in External_Filter:
				try:
					External_wl.append(filter_dict[str(an_external_filt)][1])
				except Exception as e: 
					print("External filter with name", an_external_filt, "not recognised. I am removing this filter from the list")
					mask = External_Filter!=an_external_filt
					External_Filter, External_wl, External_Flux_Jy, External_Flux_Jy_err = External_Filter[mask], External_wl[mask], External_Flux_Jy[mask], External_Flux_Jy_err[mask]
					
			sedfilter, sed_wl, sedflux, sedfluxe  =  np.concatenate((sedfilter, External_Filter)), np.concatenate((sed_wl, External_wl)), np.concatenate((sedflux, External_Flux_Jy)), np.concatenate((sedfluxe, External_Flux_Jy_err))
		except: None
		
	
	except:
		External_AB_mag = np.asarray(config_info["External_AB_mag"]).astype(float)
		External_AB_mag_err = np.asarray(config_info["External_AB_mag_err"]).astype(float)
		External_Flux_Jy=10**(23-(External_AB_mag+48.594)/2.5)  #  erg/cm2/s/Hz
		External_Flux_Jy_err=np.sqrt(  ((-(0.9210340371976184*10**((External_AB_mag+243/5)*-0.4+23)))*External_AB_mag_err)**2   )
		External_Filter = np.asarray(config_info["External_Filter"])
		External_wl = []
		for an_external_filt in External_Filter:
			try:
				External_wl.append(filter_dict[str(an_external_filt)][1])
			except Exception as e: 
				print("External filter with name", an_external_filt, "not recognised. I am removing this filter from the list")
				mask = External_Filter!=an_external_filt
				External_Filter, External_wl, External_Flux_Jy, External_Flux_Jy_err = External_Filter[mask], External_wl[mask], External_Flux_Jy[mask], External_Flux_Jy_err[mask]
				
		sedfilter, sed_wl, sedflux, sedfluxe  =  np.concatenate((sedfilter, External_Filter)), np.concatenate((sed_wl, External_wl)), np.concatenate((sedflux, External_Flux_Jy)), np.concatenate((sedfluxe, External_Flux_Jy_err))
	
	
		
	try: 
		if np.asarray(config_info["plot_phot_spectrum"]):
			plt.errorbar(sed_wl, sedflux, sedfluxe,fmt='.k');   plt.show();   plt.clf()
	except: None
		
		
		
		
	theminww, themaxww = 999999, -999999
	for afilt in sedfilter:
		filter_dict[afilt]
		ww, flfl  =  filter_dict[afilt][0]
		mask=flfl>=0
		ww, flfl = ww[mask], flfl[mask]
		minww, maxww = np.amin(ww), np.amax(ww)
		
		if minww<theminww: theminww=minww
		if maxww>themaxww: themaxww=maxww
			
		print(afilt)
		
	theminww_loadgrid = theminww-160
	themaxww_loadgrid = themaxww+160
	
	
	if theminww_loadgrid>3646:  theminww_loadgrid=3646
	if themaxww_loadgrid<6800:  themaxww_loadgrid=6800
	
	
	theminww-=60
	themaxww+=60
	print(theminww, themaxww)
	if theminww== 999999 and themaxww==-999999:  raise ValueError("The filters are likely not recognised in the script. Add them to 'filter_dict'")
	
	













#Uniform priors over which the variables will be allowed to vary
p0T1=np.asarray(config_info["p0teff"])[0].astype(float)  # max teff is 14000 for DA
p0logg1=np.asarray(config_info["p0logg"])[0].astype(float)
p0HoverHe1=np.asarray(config_info["p0HoverHe"])[0].astype(float)
try:  p0scaling=np.asarray(config_info["p0scaling"])[0].astype(float)
except: p0scaling = None

if found_parallax:
	p0parallax=[plax-3*plax_unc, plax+3*plax_unc]  # only used for plotting. gaussian prior
if starType1.startswith("sd") and fit_phot_SED:
	want_R=True
	Dguess = 3.086E16 * 1000/plax
	try: p0R = np.asarray(config_info["p0R"])[0].astype(float)
	except: p0R = np.array([0.1,0.3])
else:	want_R=False


ndim, num_DBA = 0, 0
p0range = np.array([]);   p0labels = np.array([])
pos_min = np.array([]);   pos_max = np.array([])
used_RV_boundaries, star_DBA = [], []


if starType1.startswith("D") or starType1.startswith("sd"):
	if sys.argv[1] == "ATM":

		for cn, (fteff, flogg, fHHe, starType) in enumerate(zip([forced_teff1], [forced_logg1], [forced_HoverHe1], [starType1])):
			if fteff == 0 and flogg == 0:
				ndim+=2 # teff + logg
				p0range = np.array([p0T1, p0logg1]);         p0labels = np.array(["T1", "logg1"])
				pos_min = np.array([p0T1[0], p0logg1[0]]);   pos_max = np.array([p0T1[1], p0logg1[1]])
			elif fteff == 0:
				ndim+=1 # teff
				p0range = np.array([p0T1]);         p0labels = np.array(["T1"])
				pos_min = np.array([p0T1[0]]);   pos_max = np.array([p0T1[1]])
			elif flogg == 0:
				ndim+=1 # logg
				p0range = np.array([p0logg1]);         p0labels = np.array(["logg"])
				pos_min = np.array([p0logg1[0]]);   pos_max = np.array([p0logg1[1]])
			
			if starType=="DBA" or starType.startswith("sd") and fHHe==0:
				num_DBA+=1
				star_DBA.append(str(cn+1))

		for xx in range(num_DBA):
			ndim+=1 # teff + logg + H/He
			if star_DBA[xx]=="1":
				p0range = np.concatenate((p0range, np.array([p0HoverHe1])));         p0labels = np.concatenate((p0labels, np.array(["H/He1"])))
				pos_min = np.concatenate((pos_min, np.array([p0HoverHe1[0]])));   pos_max = np.concatenate((pos_max, np.array([p0HoverHe1[1]])))



	if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
		if forced_K1==False:
			#print(len(input_files), len(share_rv), len(RV_boundaries1))
			for i, (fi, arv, arvBounds1) in enumerate(zip(input_files, share_rv, RV_boundaries1)):
				if arv==-1:
					# this is to configure the unique RVs in the mcmc
					pos_min = np.append(pos_min, arvBounds1[0]);   pos_max = np.append(pos_max, arvBounds1[1])
					ndim+=1 #Number of individual RVs to vary in the mcmc run
					
					# this is for corner plotting
					try:      p0range = np.concatenate((p0range, np.array([arvBounds1])))
					except:   p0range = np.array([arvBounds1])
					p0labels = np.append(p0labels,"RV1_"+str(i))
					
					
					used_RV_boundaries.append(arvBounds1)
				else:   print("ignoring entry of index (zero-indexed) ", i, " as a separate set of RVs, it is paired with the RV of index ", arv)

		
	elif isinstance(forced_P0, float) and isinstance(forced_T0, float):
		if forced_K1=="Fit":
			pos_min = np.append(pos_min, p0K1[0]);   pos_max = np.append(pos_max, p0K1[1])
			ndim+=1
			
			# this is for corner plotting
			try:     p0range = np.concatenate((p0range, np.array([p0K1])))
			except:  p0range = np.array([p0K1])
			p0labels = np.append(p0labels,"K1")
		
		if forced_Vgamma1=="Fit":
			pos_min = np.append(pos_min, p0Vgamma1[0]);   pos_max = np.append(pos_max, p0Vgamma1[1])
			ndim+=1
			
			# this is for corner plotting
			try:     p0range = np.concatenate((p0range, np.array([p0Vgamma1])))
			except:  p0range = np.array([p0Vgamma1])
			p0labels = np.append(p0labels,"Vg1")

	else: raise ValueError



	if sys.argv[1] == "ATM" and fit_phot_SED:
		if found_parallax:
			ndim+=1
			pos_min = np.append(pos_min, p0parallax[0]);   pos_max = np.append(pos_max, p0parallax[1])
			p0range = np.concatenate((p0range, np.array([p0parallax])))
			p0labels = np.append(p0labels,"Parallax")

		if forced_Scaling==False:
			ndim+=1
			pos_min = np.append(pos_min, p0scaling[0]);   pos_max = np.append(pos_max, p0scaling[1])
			p0range = np.concatenate((p0range, np.array([p0scaling])))
			p0labels = np.append(p0labels,"Scaling")
		
		if want_R:
			ndim+=1
			pos_min = np.append(pos_min, p0R[0]);   pos_max = np.append(pos_max, p0R[1])
			p0range = np.concatenate((p0range, np.array([p0R])))
			p0labels = np.append(p0labels,"R")
			

elif (starType1.startswith("LL") or starType1.startswith("GG") or starType1.startswith("GL")):
	#Uniform priors over which the variables will be allowed to vary

	ndim = 0
	p0range = np.array([]);   p0labels = np.array([])
	pos_min = np.array([]);   pos_max = np.array([])
	used_RV_boundaries = []


	if sys.argv[1] == "ATM":
		ndim+=4
		p0range = np.array([A1_1_boundaries, sigma1_1_boundaries, A1_2_boundaries, sigma1_2_boundaries]);   p0labels = np.array(["A1_1", "sig1_1", "A1_2", "sig1_2"])
		pos_min = np.array([A1_1_boundaries[0], sigma1_1_boundaries[0], A1_2_boundaries[0], sigma1_2_boundaries[0]])
		pos_max = np.array([A1_1_boundaries[1], sigma1_1_boundaries[1], A1_2_boundaries[1], sigma1_2_boundaries[1]])
			

		

	if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
		if forced_K1==False:
			for i, (fi, arv, arvBounds1) in enumerate(zip(input_files, share_rv, RV_boundaries1)):
				if arv==-1:
					# this is to configure the unique RVs in the mcmc
					pos_min = np.append(pos_min, arvBounds1[0]);   pos_max = np.append(pos_max, arvBounds1[1])
					ndim+=1 #Number of individual RVs to vary in the mcmc run
					
					# this is for corner plotting
					try:    p0range = np.concatenate((p0range, np.array([arvBounds1])))
					except: p0range = np.array([arvBounds1])
					p0labels = np.append(p0labels,"RV1_"+str(i))
					
					
					used_RV_boundaries.append(arvBounds1)
				else:   print("ignoring entry of index (zero-indexed) ", i, " as a separate set of RVs, it is paired with the RV of index ", arv)

		
	elif isinstance(forced_P0, float) and isinstance(forced_T0, float):
		if forced_K1=="Fit":
			pos_min = np.append(pos_min, p0K1[0]);   pos_max = np.append(pos_max, p0K1[1])
			ndim+=1
			
			# this is for corner plotting
			try:     p0range = np.concatenate((p0range, np.array([p0K1])))
			except:  p0range = np.array([p0K1])
			p0labels = np.append(p0labels,"K1")
		
		if forced_Vgamma1=="Fit":
			pos_min = np.append(pos_min, p0Vgamma1[0]);   pos_max = np.append(pos_max, p0Vgamma1[1])
			ndim+=1
			
			# this is for corner plotting
			try:     p0range = np.concatenate((p0range, np.array([p0Vgamma1])))
			except:  p0range = np.array([p0Vgamma1])
			p0labels = np.append(p0labels,"Vg1")

	else: raise ValueError


else: raise ValueError("WD-BASS Error: Star types must be two degenerate stars (Dx) or both be of: GG, LL")





psize = pos_max - pos_min
p0 = [pos_min + psize*np.random.rand(ndim) for i in range(nwalkers)]

pier_or_antoine="mixed"#"pier3Dphot_antoine1Dspec"

if (starType1=="DA") and pier_or_antoine=="pier3Dphot_antoine1Dspec":
	# load Pier
	wl_all, flux_all, Teff_all, Grav_all = load_models.load_models_DA(minwl=theminww_loadgrid,maxwl=themaxww_loadgrid)
	logg_all=np.round(np.log10(Grav_all),1)
	Teff_all_logg_lt_7 = Teff_all[logg_all<7];     wl_all_logg_lt_7 = wl_all[logg_all<7];     flux_all_logg_lt_7 = flux_all[logg_all<7];      logg_all_logg_lt_7 = logg_all[logg_all<7]
	Teff_all_logg_gteq_7 = Teff_all[logg_all>=7];  wl_all_logg_gteq_7 = wl_all[logg_all>=7];  flux_all_logg_gteq_7 = flux_all[logg_all>=7];   logg_all_logg_gteq_7 = logg_all[logg_all>=7]
	
	Teffs_allowed_mask = ((Teff_all == 5000) | (Teff_all == 5500) | (Teff_all == 6000) | (Teff_all == 6500) | (Teff_all == 7000) | (Teff_all == 7500) | (Teff_all == 8000) | (Teff_all == 8500) | (Teff_all == 9000) | (Teff_all == 9500) | (Teff_all == 10000) | (Teff_all == 10500) | (Teff_all == 11000) | (Teff_all == 11500) | (Teff_all == 12000) | (Teff_all == 12500) | (Teff_all == 13000) | (Teff_all == 13500) | (Teff_all == 14000))
	a6p5_7p0_mask = Teffs_allowed_mask & ((logg_all==6.5) | (logg_all==7))
	Teff_all_logg_6p5_7p0 = Teff_all[a6p5_7p0_mask];
	wl_all_logg_6p5_7p0 = wl_all[a6p5_7p0_mask]
	flux_all_logg_6p5_7p0 = flux_all[a6p5_7p0_mask]
	logg_all_logg_6p5_7p0 = logg_all[a6p5_7p0_mask]
	
	
	
	unique_Teffs_logg_gteq_7, unique_Teffs_logg_lt_7, unique_Teffs_logg_6p5_7p0  =  np.unique(Teff_all_logg_gteq_7), np.unique(Teff_all_logg_lt_7), np.unique(Teff_all_logg_6p5_7p0)
	
	# Load in the 1D DA models to cover the out-of-bounds region
	wl_all_1D, flux_all_1D, Teff_all_1D, logg_all_1D = load_models.load_model_1D_DA(minwl=theminww_loadgrid,maxwl=themaxww_loadgrid)
	unique_Teffs_1D = np.unique(Teff_all_1D)
	
	
	
	# LOAD ANTOINE SPECTROSCOPY
	
	wl_all_ant1D, flux_all_ant1D, Teff_all_ant1D, Grav_all_ant1D = load_models.load_models_DA_AntoineNLTE(minwl=np.amin(reference_wl)-150,maxwl=np.amax(reference_wl)+250);   mask_wl_5100_6200 = ((wl_all_ant1D<=6200) & (wl_all_ant1D>=5100));   wl_all_ant1D, flux_all_ant1D, Teff_all_ant1D, Grav_all_ant1D = wl_all_ant1D[~mask_wl_5100_6200], flux_all_ant1D[~mask_wl_5100_6200], Teff_all_ant1D[~mask_wl_5100_6200], Grav_all_ant1D[~mask_wl_5100_6200]
	
	logg_all_ant1D=np.round(np.log10(Grav_all_ant1D),2)

	unique_Teffs_ant1D  =  np.unique(Teff_all_ant1D)
	
	
	

elif (starType1=="DA") and pier_or_antoine=="pier":
	# Load in the 3D DA models
	if fit_phot_SED:         wl_all, flux_all, Teff_all, Grav_all = load_models.load_models_DA(minwl=theminww_loadgrid,maxwl=themaxww_loadgrid)
	else:                    wl_all, flux_all, Teff_all, Grav_all = load_models.load_models_DA(minwl=np.amin(reference_wl)-150,maxwl=np.amax(reference_wl)+250);   mask_wl_5100_6200 = ((wl_all<=6200) & (wl_all>=5100));   wl_all, flux_all, Teff_all, Grav_all = wl_all[~mask_wl_5100_6200], flux_all[~mask_wl_5100_6200], Teff_all[~mask_wl_5100_6200], Grav_all[~mask_wl_5100_6200]
	
	logg_all=np.round(np.log10(Grav_all),1)
	Teff_all_logg_lt_7 = Teff_all[logg_all<7];     wl_all_logg_lt_7 = wl_all[logg_all<7];     flux_all_logg_lt_7 = flux_all[logg_all<7];      logg_all_logg_lt_7 = logg_all[logg_all<7]
	Teff_all_logg_gteq_7 = Teff_all[logg_all>=7];  wl_all_logg_gteq_7 = wl_all[logg_all>=7];  flux_all_logg_gteq_7 = flux_all[logg_all>=7];   logg_all_logg_gteq_7 = logg_all[logg_all>=7]



	Teffs_allowed_mask = ((Teff_all == 5000) | (Teff_all == 5500) | (Teff_all == 6000) | (Teff_all == 6500) | (Teff_all == 7000) | (Teff_all == 7500) | (Teff_all == 8000) | (Teff_all == 8500) | (Teff_all == 9000) | (Teff_all == 9500) | (Teff_all == 10000) | (Teff_all == 10500) | (Teff_all == 11000) | (Teff_all == 11500) | (Teff_all == 12000) | (Teff_all == 12500) | (Teff_all == 13000) | (Teff_all == 13500) | (Teff_all == 14000))
	a6p5_7p0_mask = Teffs_allowed_mask & ((logg_all==6.5) | (logg_all==7))
	Teff_all_logg_6p5_7p0 = Teff_all[a6p5_7p0_mask];
	wl_all_logg_6p5_7p0 = wl_all[a6p5_7p0_mask]
	flux_all_logg_6p5_7p0 = flux_all[a6p5_7p0_mask]
	logg_all_logg_6p5_7p0 = logg_all[a6p5_7p0_mask]


	unique_Teffs_logg_gteq_7, unique_Teffs_logg_lt_7, unique_Teffs_logg_6p5_7p0  =  np.unique(Teff_all_logg_gteq_7), np.unique(Teff_all_logg_lt_7), np.unique(Teff_all_logg_6p5_7p0)
	
	
	
	
	# Load in the 1D DA models to cover the out-of-bounds region
	if fit_phot_SED:        wl_all_1D, flux_all_1D, Teff_all_1D, logg_all_1D = load_models.load_model_1D_DA(minwl=theminww_loadgrid,maxwl=themaxww_loadgrid)
	else:                   wl_all_1D, flux_all_1D, Teff_all_1D, logg_all_1D = load_models.load_model_1D_DA(minwl=np.amin(reference_wl)-150,maxwl=np.amax(reference_wl)+250);   mask_wl_5100_6200 = ((wl_all_1D<=6200) & (wl_all_1D>=5100));   wl_all_1D, flux_all_1D, Teff_all_1D, logg_all_1D = wl_all_1D[~mask_wl_5100_6200], flux_all_1D[~mask_wl_5100_6200], Teff_all_1D[~mask_wl_5100_6200], logg_all_1D[~mask_wl_5100_6200]
	unique_Teffs_1D = np.unique(Teff_all_1D)
	
	
	
	

elif (starType1=="DA") and pier_or_antoine=="antoine":
	if fit_phot_SED:         wl_all_ant1D, flux_all_ant1D, Teff_all_ant1D, Grav_all_ant1D = load_models.load_models_DA_AntoineNLTE(minwl=theminww_loadgrid,maxwl=themaxww_loadgrid)
	else:                    wl_all_ant1D, flux_all_ant1D, Teff_all_ant1D, Grav_all_ant1D = load_models.load_models_DA_AntoineNLTE(minwl=np.amin(reference_wl)-150,maxwl=np.amax(reference_wl)+250);   mask_wl_5100_6200 = ((wl_all_ant1D<=6200) & (wl_all_ant1D>=5100));   wl_all_ant1D, flux_all_ant1D, Teff_all_ant1D, Grav_all_ant1D = wl_all_ant1D[~mask_wl_5100_6200], flux_all_ant1D[~mask_wl_5100_6200], Teff_all_ant1D[~mask_wl_5100_6200], Grav_all_ant1D[~mask_wl_5100_6200]
	
	logg_all_ant1D=np.round(np.log10(Grav_all_ant1D),2)


	unique_Teffs_ant1D  =  np.unique(Teff_all_ant1D)
	
elif starType1=="DA" and pier_or_antoine=="mixed":
	if fit_phot_SED:         wl_all_synth, flux_all_synth, Teff_all_synth, Grav_all_synth = load_models.load_models_DA_3D_NLTE(minwl=theminww_loadgrid,maxwl=themaxww_loadgrid)
	else:                    wl_all_synth, flux_all_synth, Teff_all_synth, Grav_all_synth = load_models.load_models_DA_3D_NLTE(minwl=np.amin(reference_wl)-150,maxwl=np.amax(reference_wl)+250);   mask_wl_5100_6200 = ((wl_all_synth<=6200) & (wl_all_synth>=5100));   wl_all_synth, flux_all_synth, Teff_all_synth, Grav_all_synth = wl_all_synth[~mask_wl_5100_6200], flux_all_synth[~mask_wl_5100_6200], Teff_all_synth[~mask_wl_5100_6200], Grav_all_synth[~mask_wl_5100_6200]
	
	if False:
		## Trim the grid to only include the values within p0Teff and p0logg
		min_p0T = np.amin(np.array([p0T1[0], p0T1[1]]));    max_p0T = np.amax(np.array([p0T1[0], p0T1[1]]))
		
		mask_TEFF = (np.unique(Teff_all_synth) >= min_p0T)  &  (np.unique(Teff_all_synth) <= max_p0T)
		try:  mask_TEFF[np.amin(np.argwhere(mask_TEFF==True)) - 1] = True
		except: None
		try: mask_TEFF[np.amax(np.argwhere(mask_TEFF==True)) + 1] = True
		except: None
		
		min_p0logg = np.amin(np.array([p0logg1[0], p0logg1[1]]));    max_p0logg = np.amax(np.array([p0logg1[0], p0logg1[1]]))
		
		mask_LOGG = (np.unique(Grav_all_synth) >= min_p0logg)  &  (np.unique(Grav_all_synth) <= max_p0logg)
		try: mask_LOGG[np.amin(np.argwhere(mask_LOGG==True)) - 1] = True
		except: None
		try: mask_LOGG[np.amax(np.argwhere(mask_LOGG==True)) + 1] = True
		except: None
		
		
		newMinT, newMaxT = np.amin(np.unique(Teff_all_synth)[mask_TEFF]), np.amax(np.unique(Teff_all_synth)[mask_TEFF])
		newMinLogg, newMagLogg = np.amin(np.unique(Grav_all_synth)[mask_LOGG]), np.amax(np.unique(Grav_all_synth)[mask_LOGG])
		
		
		mask_atm_grid = (Teff_all_synth <= newMaxT) & (Teff_all_synth >= newMinT)  &  (Grav_all_synth <= newMagLogg) & (Grav_all_synth >= newMinLogg)
		
		wl_all_synth, flux_all_synth, Teff_all_synth, Grav_all_synth = wl_all_synth[mask_atm_grid], flux_all_synth[mask_atm_grid], Teff_all_synth[mask_atm_grid], Grav_all_synth[mask_atm_grid]
	
	
	logg_all_synth=Grav_all_synth
	unique_Teffs_synth  =  np.unique(Teff_all_synth)
	
	
	
	
if starType1=="DBA" or starType1=="DB":
	# Load in the 3D DB/DBA/DC LTE spectra
	if fit_phot_SED:         wl_all, flux_all, Teff_all, Grav_all, H_over_He_all = load_models.load_models_DBA(minwl=theminww_loadgrid,maxwl=themaxww_loadgrid)
	else:                    wl_all, flux_all, Teff_all, Grav_all, H_over_He_all = load_models.load_models_DBA(minwl=np.amin(reference_wl)-150,maxwl=np.amax(reference_wl)+250)#;   mask_wl_5100_6200 = ((wl_all<=6200) & (wl_all>=5100));   wl_all, flux_all, Teff_all, Grav_all, H_over_He_all = wl_all[~mask_wl_5100_6200], flux_all[~mask_wl_5100_6200], Teff_all[~mask_wl_5100_6200], Grav_all[~mask_wl_5100_6200], H_over_He_all[~mask_wl_5100_6200]
	#H_over_He_all=np.log10(H_over_He_all)
	#logg_all=np.round(np.log10(Grav_all),1)
	logg_all=Grav_all

	#Teff_all=np.round(Teff_all,-3) # this is naughty, but decreases run time with a small loss of accuracy to normalised profiles because teff grids aren't exactly the same
	Teff_all_logg = Teff_all;   wl_all_logg = wl_all;     flux_all_logg = flux_all;      logg_all_logg = logg_all;     HoverHe_all_logg_DBA = H_over_He_all
	unique_HoverHe_DBA=np.unique(HoverHe_all_logg_DBA)
	
	if starType1=="DB":
		maskDBonly=HoverHe_all_logg_DBA==30
		Teff_all_logg = Teff_all[maskDBonly];   wl_all_logg = wl_all[maskDBonly];     flux_all_logg = flux_all[maskDBonly];      logg_all_logg = logg_all[maskDBonly]
		del HoverHe_all_logg_DBA, H_over_He_all
	
	
	#print(np.amin(H_over_He_all), np.amax(H_over_He_all))
	
if starType1=="DC":
	if fit_phot_SED:         wl_all, flux_all, Teff_all, Grav_all = load_models.load_models_DBA_1eMinus5(minwl=theminww_loadgrid,maxwl=themaxww_loadgrid)
	else:                    wl_all, flux_all, Teff_all, Grav_all = load_models.load_models_DBA_1eMinus5(minwl=np.amin(reference_wl)-150,maxwl=np.amax(reference_wl)+250)
	
	
	## Trim the grid to only include the values within p0Teff and p0logg
	min_p0T = np.amin(np.array([p0T1[0], p0T1[1]]));    max_p0T = np.amax(np.array([p0T1[0], p0T1[1]]))
	
	mask_TEFF = (np.unique(Teff_all) >= min_p0T)  &  (np.unique(Teff_all) <= max_p0T)
	mask_TEFF[np.amin(np.argwhere(mask_TEFF==True)) - 1] = True;      mask_TEFF[np.amax(np.argwhere(mask_TEFF==True)) + 1] = True
	
	min_p0logg = np.amin(np.array([p0logg1[0], p0logg1[1]]));    max_p0logg = np.amax(np.array([p0logg1[0], p0logg1[1]]))
	
	mask_LOGG = (np.unique(Grav_all) >= min_p0logg)  &  (np.unique(Grav_all) <= max_p0logg)
	mask_LOGG[np.amin(np.argwhere(mask_LOGG==True)) - 1] = True;    mask_LOGG[np.amax(np.argwhere(mask_LOGG==True)) + 1] = True
	
	
	newMinT, newMaxT = np.amin(np.unique(Teff_all)[mask_TEFF]), np.amax(np.unique(Teff_all)[mask_TEFF])
	newMinLogg, newMagLogg = np.amin(np.unique(Grav_all)[mask_LOGG]), np.amax(np.unique(Grav_all)[mask_LOGG])
	
	
	mask_atm_grid = (Teff_all <= newMaxT) & (Teff_all >= newMinT)  &  (Grav_all <= newMagLogg) & (Grav_all >= newMinLogg)
	
	wl_all, flux_all, Teff_all, Grav_all = wl_all[mask_atm_grid], flux_all[mask_atm_grid], Teff_all[mask_atm_grid], Grav_all[mask_atm_grid]
	
	
	
	logg_all=np.round(Grav_all,1)
	Teff_all_logg = Teff_all;   wl_all_logg = wl_all;     flux_all_logg = flux_all;      logg_all_logg = logg_all


elif starType1=="sd":
	minP0H, maxP0H = np.amin(p0HoverHe1), np.amax(p0HoverHe1)
	if maxP0H>-0.041 or minP0H<-5.05:
		raise ValueError("p0HoverHe goes outside of grid bounds that are [-5.05, -0.041]")

	if fit_phot_SED:         wl_all, flux_all, Teff_all, Grav_all, H_over_He_all = load_models.load_models_Subdwarfs(minwl=theminww_loadgrid,maxwl=themaxww_loadgrid)
	else:                    wl_all, flux_all, Teff_all, Grav_all, H_over_He_all = load_models.load_models_Subdwarfs(minwl=np.amin(reference_wl)-150,maxwl=np.amax(reference_wl)+250)
	
	#raise ValueError(np.amin(Teff_all), np.amax(Teff_all), np.unique(Teff_all))
	
	if forced_teff1==0:
		## Trim the grid to only include the values within p0Teff and p0logg
		min_p0T = np.amin(np.array([p0T1[0], p0T1[1]]));    max_p0T = np.amax(np.array([p0T1[0], p0T1[1]]))
		mask_TEFF = (np.unique(Teff_all) >= min_p0T)  &  (np.unique(Teff_all) <= max_p0T)
		try: mask_TEFF[np.amin(np.argwhere(mask_TEFF==True)) - 1] = True
		except: None
		try: mask_TEFF[np.amax(np.argwhere(mask_TEFF==True)) + 1] = True
		except: None
		
		min_p0logg = np.amin(np.array([p0logg1[0], p0logg1[1]]));    max_p0logg = np.amax(np.array([p0logg1[0], p0logg1[1]]))
		mask_LOGG = (np.unique(Grav_all) >= min_p0logg)  &  (np.unique(Grav_all) <= max_p0logg)
		try: mask_LOGG[np.amin(np.argwhere(mask_LOGG==True)) - 1] = True
		except: None
		try: mask_LOGG[np.amax(np.argwhere(mask_LOGG==True)) + 1] = True
		except: None
		
		min_p0HoverHe = np.amin(np.array([p0HoverHe1[0], p0HoverHe1[1]]));    max_p0HoverHe = np.amax(np.array([p0HoverHe1[0], p0HoverHe1[1]]))
		mask_HOVERHE = (np.unique(H_over_He_all) >= min_p0HoverHe)  &  (np.unique(H_over_He_all) <= max_p0HoverHe)
		try: mask_HOVERHE[np.amin(np.argwhere(mask_HOVERHE==True)) - 1] = True
		except: None
		try: mask_HOVERHE[np.amax(np.argwhere(mask_HOVERHE==True)) + 1] = True
		except: None
		
		
		newMinT, newMaxT = np.amin(np.unique(Teff_all)[mask_TEFF]), np.amax(np.unique(Teff_all)[mask_TEFF])
		newMinLogg, newMaxLogg = np.amin(np.unique(Grav_all)[mask_LOGG]), np.amax(np.unique(Grav_all)[mask_LOGG])
		newMinHoverHe, newMaxHoverHe = np.amin(np.unique(H_over_He_all)[mask_HOVERHE]), np.amax(np.unique(H_over_He_all)[mask_HOVERHE])
		
		
		mask_atm_grid = (Teff_all <= newMaxT) & (Teff_all >= newMinT)  &  (Grav_all <= newMaxLogg) & (Grav_all >= newMinLogg)  &  (H_over_He_all <= newMaxHoverHe) & (H_over_He_all >= newMinHoverHe)
		
		wl_all, flux_all, Teff_all, Grav_all, H_over_He_all = wl_all[mask_atm_grid], flux_all[mask_atm_grid], Teff_all[mask_atm_grid], Grav_all[mask_atm_grid], H_over_He_all[mask_atm_grid]
		
		

	Teff_all_logg = Teff_all;   wl_all_logg = wl_all;     flux_all_logg = flux_all;      logg_all_logg = Grav_all;     HoverHe_all_logg = H_over_He_all
	unique_HoverHe=np.unique(HoverHe_all_logg)
	
	
	



#print(np.amin(Teff_all), np.amax(Teff_all))
#print(np.amin(np.log10(Grav_all)), np.amax(np.log10(Grav_all)))




list_norm_wl_grids, list_normalised_flux, list_normalised_err = [], [], []
mask_out_Ha_min_all, mask_out_Ha_max_all, cut_limits_min_all, cut_limits_max_all = [], [], [], []

min_wl_all_files, max_wl_all_files = 9E9, 0
# Normalise the input data
for files, normaliseHa, cut_Ha, ref_wl in zip(input_files, normaliseHa_all, cut_Ha_all, reference_wl):
	if spectra_source_type=="wl_flux_fluxerr":
		try:  wl_data, flux_data, flux_e_data = np.loadtxt(os.getcwd()+"/"+files, skiprows=1, unpack=True)
		except:  
			##### this part of the code is for when you do not have flux errors. I do my best to guess the SNR of the data given an expected continuum SNR under the assumption that your spectral line follows a gaussian profile. In the "gauss" function below, mmm and ccc are included to fit the general trend of the continuum as a linear fit. 
			
			##### I used this method to fit to FIES observations as errors are not output in its reduction pipeline
			wl_data, flux_data = np.loadtxt(os.getcwd()+"/"+files, skiprows=2, unpack=True)
			
			cut_limits_min, cut_limits_max = ref_wl+cut_Ha[0],  ref_wl+cut_Ha[1]
			norm_limits_min, norm_limits_max = ref_wl+normaliseHa[0],  ref_wl+normaliseHa[1]
			
			
			
			mask_cut = (wl_data>=cut_limits_min) & (wl_data<=cut_limits_max)
			mask_norm = ((wl_data>=cut_limits_min) & (wl_data<=norm_limits_min)) | ((wl_data<=cut_limits_max) & (wl_data>=norm_limits_max))
			
			
			from scipy.optimize import curve_fit
			popt, pcov = curve_fit(gauss_if_no_flux_error, wl_data[mask_cut], flux_data[mask_cut], p0 = [-0.8, ref_wl, 10, 1, 1, 1], bounds=[[-np.inf,ref_wl-5,0,0,-np.inf,-np.inf], [0,ref_wl+5,80,np.inf,np.inf,np.inf]])
			
			# plt.plot(wl_data[mask_cut], flux_data[mask_cut], c='k');  plt.plot(wl_data[mask_cut], gauss_if_no_flux_error(wl_data[mask_cut], *popt));  plt.show()
			
			
			med_flux = np.median(flux_data[mask_norm])   # this gets the median flux of the normalised region
			
			
			predicted_SNR_of_normalise_region_per_pixel = 10  # this is the predicted flux of the normalised part of the spectrum.  I use this to predict what the SNR is for all other pixels below
			
			flux_e_data = med_flux /predicted_SNR_of_normalise_region_per_pixel   *  np.sqrt(med_flux/gauss_if_no_flux_error(wl_data, *popt))  # propogate SNR to all parts of the spectrum. I assume that there is no detector readout noise and so your SNR only depends on the root(#photons)
			
			#plt.plot(wl_data[mask_cut], flux_e_data[mask_cut]);   plt.plot(wl_data[mask_cut], flux_data[mask_cut],c='k');   plt.show()
			
			
			
			flux_e_data[flux_e_data<0] = np.amax(flux_e_data[mask_cut][flux_e_data[mask_cut]>0]) #  if negative flux, the error on these measurements is taken as the maximum positive error. Bit arbitrary... better to have actual flux errors!! I do this in case there is a big negative outlier (e.g. poor sky subtraction)
			
		args = np.argsort(wl_data)
		wl_data=wl_data[args]
		flux_data=flux_data[args]
		flux_e_data=flux_e_data[args]
		#print(wl_data)
	elif spectra_source_type=="Xshooter":
		hdul=fits.open(files)
		header=hdul[0].header
		#plt.title(header["SPEC_RES"])
		arm=header["HIERARCH ESO SEQ ARM"]
		data=hdul[1].data
		
		if False:
			MJDstart=header["MJD-OBS"]
			MJDend=header["MJD-END"]
			MJDmid = (MJDstart+MJDend)/2
			
			RA, dec = header["RA"], header["DEC"]
			
			c=SkyCoord(str(RA) + " " + str(dec), frame="icrs", unit=(u.deg,u.deg))
			ra=c.ra.deg; dec=c.dec.deg
			print(ra, dec)
			jd=Time(MJDmid,format='mjd',scale='utc',location=VLT)
			
			
			corr=jd.light_travel_time(c,kind="barycentric")
			MJDmid=(jd+corr).value
			BJDtime=MJDmid
		
		
		qual=data["qual"][0]
		snr=data["snr"][0]
		
		mask=((qual==0) & (snr>2))
		
		wl_data=data["wave"][0][mask]*10
		flux_data=data["flux"][0][mask]   # erg cm-2 s-1 AA-1
		flux_e_data=data["err"][0][mask]  # erg cm-2 s-1 AA-1
		
		flux_data = 3.33564095E+04 * flux_data * np.square(wl_data) / 1E23   # convert to  erg/cm^2/s/Hz
		flux_e_data = 3.33564095E+04 * flux_e_data * np.square(wl_data) / 1E23 # convert to erg/cm^2/s/Hz
        
	plt.clf()
	# I deredden the observed spectra (spectra, not flux calibrated photometry - this is modelled with a reddened synthetic spectrum) so that the synthetic spectra do not need to be redenned every iteration. This may become an issue if you have very low resolution (R<100) data, but otherwise it's not something to worry about
	from dust_extinction.parameter_averages import G23
	ext = G23(Rv=3.1)
	flux_data /= ext.extinguish(wl_data*u.AA, Ebv=reddening_Ebv)
	flux_e_data /= ext.extinguish(wl_data*u.AA, Ebv=reddening_Ebv)
	
	
	#plt.plot(wl_data, (flux_data-orig_f)/flux_data,c='g');  plt.show()
        
	#plt.plot(wl_data,flux_data);  plt.title(files);  plt.show();  plt.clf()
	
	minwlhere=np.amin(wl_data);        maxwlhere=np.amax(wl_data)
	
	if not (ref_wl < maxwlhere and ref_wl > minwlhere) and minwlhere<ref_wl-100: raise ValueError("File named: " + str(files) + "  does not have a valid reference wl (" + str(ref_wl) + ", while the minimum wl in the file was "+str(minwlhere)+")")
	if minwlhere<min_wl_all_files:     min_wl_all_files = minwlhere -50
	if maxwlhere>max_wl_all_files:     max_wl_all_files = maxwlhere +50
	
	
	
	cut_limits_min, cut_limits_max = ref_wl+cut_Ha[0],  ref_wl+cut_Ha[1]
	cut_limits_min_all.append(cut_limits_min);    cut_limits_max_all.append(cut_limits_max)
	
	
	mask_out_Ha_min, mask_out_Ha_max = ref_wl+normaliseHa[0],  ref_wl+normaliseHa[1]
	mask_out_Ha_min_all.append(mask_out_Ha_min);    mask_out_Ha_max_all.append(mask_out_Ha_max)
	
	mask_out = (wl_data>cut_limits_min) & (wl_data<cut_limits_max)
	wl_data=wl_data[mask_out];   flux_data = flux_data[mask_out];   flux_e_data = flux_e_data[mask_out]
	
	mask_norm = (wl_data>cut_limits_min) & (wl_data<cut_limits_max) & ((wl_data<mask_out_Ha_min) | (wl_data>mask_out_Ha_max))
	
	
	wl_data_normmask = wl_data[mask_norm]; flux_data_normmask = flux_data[mask_norm];   fluxe_data_normmask = flux_e_data[mask_norm]
	
	
	
	
	good=flux_data_normmask==flux_data_normmask
	while True:
		current_length = len(flux_data_normmask[good])
		m,c = polyfit(wl_data_normmask[good], flux_data_normmask[good], w=1/fluxe_data_normmask[good], deg=1)
		
		flux_norm_sigma_clip = flux_data_normmask - (m*wl_data_normmask+c)
		
		std = np.std(flux_data_normmask[good])
		good = good & ((np.abs(flux_norm_sigma_clip)  < 1.5*np.abs(std)))
		
		#plt.scatter(wl_data_normmask[~good], flux_data_normmask[~good],c='r')
		
		if len(flux_data_normmask[good]) == current_length:
			break
	#plt.plot(wl_data[mask_norm], flux_data[mask_norm]);   plt.show();   plt.clf()
	normalised_wl_grid = wl_data;    normalised_flux = flux_data/(m*wl_data+c);    normalised_err = flux_e_data/(m*wl_data+c)
	list_norm_wl_grids.append(normalised_wl_grid);     list_normalised_flux.append(normalised_flux);     list_normalised_err.append(normalised_err)
	
	
try:
	if fit_phot_SED==False:
		mask = (wl_all_logg_6p5_7p0 > min_wl_all_files) & (wl_all_logg_6p5_7p0 < max_wl_all_files)
		Teff_all_logg_6p5_7p0 = Teff_all_logg_6p5_7p0[mask];    wl_all_logg_6p5_7p0 = wl_all_logg_6p5_7p0[mask];    flux_all_logg_6p5_7p0 = flux_all_logg_6p5_7p0[mask];    logg_all_logg_6p5_7p0 = logg_all_logg_6p5_7p0[mask]
except: None

try:
	if fit_phot_SED==False:
		mask = (wl_all_logg_lt_7 > min_wl_all_files) & (wl_all_logg_lt_7 < max_wl_all_files)
		Teff_all_logg_lt_7 = Teff_all_logg_lt_7[mask];     wl_all_logg_lt_7 = wl_all_logg_lt_7[mask];     flux_all_logg_lt_7 = flux_all_logg_lt_7[mask];      logg_all_logg_lt_7 = logg_all_logg_lt_7[mask]
except: None

try: 
	if fit_phot_SED==False:
		mask = (wl_all_logg_gteq_7 > min_wl_all_files) & (wl_all_logg_gteq_7 < max_wl_all_files)
		Teff_all_logg_gteq_7 = Teff_all_logg_gteq_7[mask];  wl_all_logg_gteq_7 = wl_all_logg_gteq_7[mask]; flux_all_logg_gteq_7 = flux_all_logg_gteq_7[mask];  logg_all_logg_gteq_7 = logg_all_logg_gteq_7[mask]
except: None

try:
	if fit_phot_SED==False:
		mask = (wl_all_1D > min_wl_all_files) & (wl_all_1D < max_wl_all_files)
		wl_all_1D = wl_all_1D[mask]; flux_all_1D = flux_all_1D[mask]; Teff_all_1D = Teff_all_1D[mask]; logg_all_1D = logg_all_1D[mask]
except: None

try:
	if fit_phot_SED==False:
		mask = (wl_all_logg > min_wl_all_files) & (wl_all_logg < max_wl_all_files)
		Teff_all_logg = Teff_all_logg[mask];   wl_all_logg = wl_all_logg[mask];     flux_all_logg = flux_all_logg[mask];      logg_all_logg = logg_all_logg[mask];     HoverHe_all_logg_DBA = HoverHe_all_logg_DBA[mask]
except: None



if len(p0labels) == 1:
	p0labels = np.flip(np.append(p0labels, "Dummy"))
	dummybound = [-1, 1]
	pos_min = np.flip(np.append(pos_min, dummybound[0]))
	pos_max = np.flip(np.append(pos_max, dummybound[-1]))
	p0range = np.flip(np.concatenate((p0range, np.array([dummybound]))),axis=0)
	ndim+=1
	psize = pos_max - pos_min
	p0 = [pos_min + psize*np.random.rand(ndim) for i in range(nwalkers)]


@njit
def return_DAgrids(temperature_star, logg_star):
	pier_or_antoine="mixed"#"pier3Dphot_antoine1Dspec"
	minval=10000;   maxval=-10000
	#if temperature_star>25000: raise ValueError
	
	if pier_or_antoine=="pier" or pier_or_antoine=="pier3Dphot_antoine1Dspec":
		if not (logg_star<7 and (temperature_star>14000 or temperature_star<5000)):   # 3D DA grid (unique structure for each of the below categories). Limits Logg = 5-9 . Teff 5000-14000
			# start by finding the nearest 2 Teffs, then take the two for interpolation limits
			if logg_star>6.5 and logg_star<=7.0:
				Teff_all = Teff_all_logg_6p5_7p0;    wl_all = wl_all_logg_6p5_7p0;    flux_all = flux_all_logg_6p5_7p0;    logg_all = logg_all_logg_6p5_7p0
				
				
				temdiff = temperature_star - unique_Teffs_logg_6p5_7p0
				Teff_min=unique_Teffs_logg_6p5_7p0[np.argwhere(temdiff==np.amin(temdiff[temdiff>0]))[0][0]]
				Teff_max=unique_Teffs_logg_6p5_7p0[np.argwhere(temdiff==np.amax(temdiff[temdiff<0]))[0][0]]
				
				
			elif logg_star<7:
				Teff_all = Teff_all_logg_lt_7;    logg_all = logg_all_logg_lt_7;    flux_all = flux_all_logg_lt_7;    wl_all = wl_all_logg_lt_7
				
				temdiff = temperature_star - unique_Teffs_logg_lt_7
				Teff_min=unique_Teffs_logg_lt_7[np.argwhere(temdiff==np.amin(temdiff[temdiff>0]))[0][0]]
				Teff_max=unique_Teffs_logg_lt_7[np.argwhere(temdiff==np.amax(temdiff[temdiff<0]))[0][0]]
				
				
			else:
				Teff_all = Teff_all_logg_gteq_7;    logg_all = logg_all_logg_gteq_7;    flux_all = flux_all_logg_gteq_7;    wl_all = wl_all_logg_gteq_7
				
				temdiff = temperature_star - unique_Teffs_logg_gteq_7
				Teff_min=unique_Teffs_logg_gteq_7[np.argwhere(temdiff==np.amin(temdiff[temdiff>0]))[0][0]]
				Teff_max=unique_Teffs_logg_gteq_7[np.argwhere(temdiff==np.amax(temdiff[temdiff<0]))[0][0]]
				
				
			#if Teff_min==Teff_max==5000:
			#	if logg_star>=7:  Teff_max=5250
			#	else: Teff_max=5500
			#elif Teff_min==Teff_max==14000:     Teff_min=13500
			
			
			
		    	
		    	# then find the nearest 2 loggs
			if logg_star<7:   list_search = [5,5.5,6,6.5]
			else:         list_search = [7,7.5,8,8.5,9]
			
			if not (logg_star>6.5 and logg_star<7.0):
				for logg_opts in list_search:
					if logg_opts-logg_star<np.abs(logg_opts-maxval):    maxval = logg_opts
					if logg_star-logg_opts>0:   minval = logg_opts
			else:
				minval=6.5;   maxval=7.0
			
			mask_logg = (logg_all<=maxval) & (logg_all>=minval) & (Teff_all<=Teff_max) & (Teff_all>=Teff_min)
			
			
			
			Grav1_N = logg_all[mask_logg];    wl_all1_N=wl_all[mask_logg];    flux1_N=flux_all[mask_logg];    Teff1_N=Teff_all[mask_logg]
			
			if len(np.unique(Teff1_N))==3:
				un_teffs = np.unique(Teff1_N)
				if temperature_star<un_teffs[1]:   newmask = (Teff1_N!=un_teffs[2])
				else:   newmask = (Teff1_N!=un_teffs[0])
				
				Grav1_N, wl_all1_N, flux1_N, Teff1_N = Grav1_N[newmask], wl_all1_N[newmask], flux1_N[newmask], Teff1_N[newmask]
			
			
			if pier_or_antoine=="pier":
				return Grav1_N, wl_all1_N, flux1_N, Teff1_N
			elif pier_or_antoine=="pier3Dphot_antoine1Dspec":
				Teff_all = Teff_all_ant1D;    wl_all = wl_all_ant1D;    flux_all = flux_all_ant1D;    logg_all = logg_all_ant1D
				
				
				temdiff = temperature_star - unique_Teffs_ant1D
				Teff_min=unique_Teffs_ant1D[np.argwhere(temdiff==np.amin(temdiff[temdiff>0]))[0][0]]
				Teff_max=unique_Teffs_ant1D[np.argwhere(temdiff==np.amax(temdiff[temdiff<0]))[0][0]]
				
				


				minval=10000;   maxval=-10000
				# then find the nearest 2 loggs
				list_search = [6.5,6.75,7,7.25,7.5,7.75,8,8.25,8.5,8.75,9]

				for logg_opts in list_search:
					if logg_opts-logg_star<np.abs(logg_opts-maxval):    maxval = logg_opts
					if logg_star-logg_opts>0:   minval = logg_opts


				mask_logg = (logg_all_ant1D<=maxval) & (logg_all_ant1D>=minval) & (Teff_all_ant1D<=Teff_max) & (Teff_all_ant1D>=Teff_min)
				

				antGrav1_N = logg_all_ant1D[mask_logg];    antwl_all1_N=wl_all_ant1D[mask_logg];    antflux1_N=flux_all_ant1D[mask_logg];    antTeff1_N=Teff_all_ant1D[mask_logg]


				if len(np.unique(antTeff1_N))==3:
					un_teffs = np.unique(antTeff1_N)
					if temperature_star<un_teffs[1]:   newmask = (antTeff1_N!=un_teffs[2])
					else:   newmask = (antTeff1_N!=un_teffs[0])
					
					
					antGrav1_N, antwl_all1_N, antflux1_N, antTeff1_N = antGrav1_N[newmask], antwl_all1_N[newmask], antflux1_N[newmask], antTeff1_N[newmask]

				return antGrav1_N, antwl_all1_N, antflux1_N, antTeff1_N, Grav1_N, wl_all1_N, flux1_N, Teff1_N
				


		else:
			# find the nearest 2 Teffs and loggs for the 1D DA grid also for when it goes out of bounds, then take the two for interpolation limits
			Teff_all = Teff_all_1D;    wl_all = wl_all_1D;    flux_all = flux_all_1D;    logg_all = logg_all_1D
			
			temdiff = temperature_star - unique_Teffs_1D
			Teff_min=unique_Teffs_1D[np.argwhere(temdiff==np.amin(temdiff[temdiff>0]))[0][0]]
			Teff_max=unique_Teffs_1D[np.argwhere(temdiff==np.amax(temdiff[temdiff<0]))[0][0]]
			
			
				
			list_search = [6.5,7,7.5,8,8.5,9,9.5]
			
			for logg_opts in list_search:
				if logg_opts-logg_star<np.abs(logg_opts-maxval):    maxval = logg_opts
				if logg_star-logg_opts>0:   minval = logg_opts
			
			mask_logg = (logg_all<=maxval) & (logg_all>=minval) & (Teff_all<=Teff_max) & (Teff_all>=Teff_min)
			
			
			if len(logg_all[mask_logg])==0:
				raise ValueError("ERROR - Grid of 1D models is logg 6.5 - 9.5... you are probably searching a logg that is not allowed")
			
			Grav1_N = logg_all[mask_logg];    wl_all1_N=wl_all[mask_logg];    flux1_N=flux_all[mask_logg];    Teff1_N=Teff_all[mask_logg]
			
			#print(temperature_star, logg_star, np.amin(Grav1_N), np.amax(Grav1_N), np.amin(Teff1_N), np.amax(Teff1_N), np.amin(wl_all1_N), np.amax(wl_all1_N))
			#wl_grid = np.arange(6562-100, 6562+100,0.07)
			#model_spectrum=griddata(np.array([Teff1_N, wl_all1_N, Grav1_N]).T, flux1_N, np.array([np.full((len(wl_grid),),temperature_star), wl_grid, np.full((len(wl_grid),),logg_star)]).T, method="linear")
			#plt.plot(wl_grid, model_spectrum)
			#plt.title(str(temperature_star) + "  " + str(logg_star))
			#plt.show();   plt.clf()
			
			if len(np.unique(Teff1_N))==3:
				un_teffs = np.unique(Teff1_N)
				if temperature_star<un_teffs[1]:   newmask = (Teff1_N!=un_teffs[2])
				else:   newmask = (Teff1_N!=un_teffs[0])
				
				Grav1_N, wl_all1_N, flux1_N, Teff1_N = Grav1_N[newmask], wl_all1_N[newmask], flux1_N[newmask], Teff1_N[newmask]
			
			
			
			
			if pier_or_antoine=="pier":
				return Grav1_N, wl_all1_N, flux1_N, Teff1_N
			
			
			elif pier_or_antoine=="pier3Dphot_antoine1Dspec":
				Teff_all = Teff_all_ant1D;    wl_all = wl_all_ant1D;    flux_all = flux_all_ant1D;    logg_all = logg_all_ant1D
				
				
				temdiff = temperature_star - unique_Teffs_ant1D
				Teff_min=unique_Teffs_ant1D[np.argwhere(temdiff==np.amin(temdiff[temdiff>0]))[0][0]]
				Teff_max=unique_Teffs_ant1D[np.argwhere(temdiff==np.amax(temdiff[temdiff<0]))[0][0]]
				
				
				minval=10000;   maxval=-10000
				# then find the nearest 2 loggs
				list_search = [6.5,6.75,7,7.25,7.5,7.75,8,8.25,8.5,8.75,9]

				for logg_opts in list_search:
					if logg_opts-logg_star<np.abs(logg_opts-maxval):    maxval = logg_opts
					if logg_star-logg_opts>0:   minval = logg_opts


				mask_logg = (logg_all_ant1D<=maxval) & (logg_all_ant1D>=minval) & (Teff_all_ant1D<=Teff_max) & (Teff_all_ant1D>=Teff_min)


				antGrav1_N = logg_all_ant1D[mask_logg];    antwl_all1_N=wl_all_ant1D[mask_logg];    antflux1_N=flux_all_ant1D[mask_logg];    antTeff1_N=Teff_all_ant1D[mask_logg]
				


				if len(np.unique(antTeff1_N))==3:
					un_teffs = np.unique(antTeff1_N)
					if temperature_star<un_teffs[1]:   newmask = (antTeff1_N!=un_teffs[2])
					else:   newmask = (antTeff1_N!=un_teffs[0])
					
					antGrav1_N, antwl_all1_N, antflux1_N, antTeff1_N = antGrav1_N[newmask], antwl_all1_N[newmask], antflux1_N[newmask], antTeff1_N[newmask]

				return antGrav1_N, antwl_all1_N, antflux1_N, antTeff1_N, Grav1_N, wl_all1_N, flux1_N, Teff1_N
			
			
			
			
	elif pier_or_antoine=="antoine":
		if logg_star>=6.5:
			Teff_all = Teff_all_ant1D;    wl_all = wl_all_ant1D;    flux_all = flux_all_ant1D;    logg_all = logg_all_ant1D
			
			
			temdiff = temperature_star - unique_Teffs_ant1D
			Teff_min=unique_Teffs_ant1D[np.argwhere(temdiff==np.amin(temdiff[temdiff>0]))[0][0]]
			Teff_max=unique_Teffs_ant1D[np.argwhere(temdiff==np.amax(temdiff[temdiff<0]))[0][0]]
			
			
		else: raise ValueError
			
		
		
	    	
	    	# then find the nearest 2 loggs
		list_search = [6.5,6.75,7,7.25,7.5,7.75,8,8.25,8.5,8.75,9]
		
		for logg_opts in list_search:
			if logg_opts-logg_star<np.abs(logg_opts-maxval):    maxval = logg_opts
			if logg_star-logg_opts>0:   minval = logg_opts
		
		
		mask_logg = (logg_all_ant1D<=maxval) & (logg_all_ant1D>=minval) & (Teff_all_ant1D<=Teff_max) & (Teff_all_ant1D>=Teff_min)
		
		
		Grav1_N = logg_all_ant1D[mask_logg];    wl_all1_N=wl_all_ant1D[mask_logg];    flux1_N=flux_all_ant1D[mask_logg];    Teff1_N=Teff_all_ant1D[mask_logg]
		
		
		if len(np.unique(Teff1_N))==3:
			un_teffs = np.unique(Teff1_N)
			if temperature_star<un_teffs[1]:   newmask = (Teff1_N!=un_teffs[2])
			else:   newmask = (Teff1_N!=un_teffs[0])
			
			Grav1_N, wl_all1_N, flux1_N, Teff1_N = Grav1_N[newmask], wl_all1_N[newmask], flux1_N[newmask], Teff1_N[newmask]
		
		return Grav1_N, wl_all1_N, flux1_N, Teff1_N
		
	elif pier_or_antoine == "mixed":
		if logg_star>=6.5:
			Teff_all = Teff_all_synth;    wl_all = wl_all_synth;    flux_all = flux_all_synth;    logg_all = logg_all_synth
			
			
			temdiff = temperature_star - unique_Teffs_synth
			Teff_min=unique_Teffs_synth[np.argwhere(temdiff==np.amin(temdiff[temdiff>0]))[0][0]]
			Teff_max=unique_Teffs_synth[np.argwhere(temdiff==np.amax(temdiff[temdiff<0]))[0][0]]
			
		else: raise ValueError
			
		
		
	    	
	    	# then find the nearest 2 loggs
		list_search = [6.5,7,7.5,8,8.5,9]
		
		for logg_opts in list_search:
			if logg_opts-logg_star<np.abs(logg_opts-maxval):    maxval = logg_opts
			if logg_star-logg_opts>0:   minval = logg_opts
		
		
		mask_logg = (logg_all<=maxval) & (logg_all>=minval) & (Teff_all<=Teff_max) & (Teff_all>=Teff_min)
		
		
		Grav1_N = logg_all[mask_logg];    wl_all1_N=wl_all[mask_logg];    flux1_N=flux_all[mask_logg];    Teff1_N=Teff_all[mask_logg]
		
		
		if len(np.unique(Teff1_N))==3:
			un_teffs = np.unique(Teff1_N)
			if temperature_star<un_teffs[1]:   newmask = (Teff1_N!=un_teffs[2])
			else:   newmask = (Teff1_N!=un_teffs[0])
			
			Grav1_N, wl_all1_N, flux1_N, Teff1_N = Grav1_N[newmask], wl_all1_N[newmask], flux1_N[newmask], Teff1_N[newmask]
		
		return Grav1_N, wl_all1_N, flux1_N, Teff1_N



@njit
def return_DBgrids(temperature_star, logg_star):
	# find nearest Teff in models in the DA/DBA/DC grid, then take the two for interpolation limits
	minval=10000;   maxval=-10000
	if logg_star>=7.5 and logg_star<8: minval=7.5; maxval=8
	elif logg_star>=8 and logg_star<8.5: minval=8; maxval=8.5
	elif logg_star>=8.5 and logg_star<=9: minval=8.5; maxval=9
	elif logg_star<7.5: raise ValueError("Minimum logg for grid is 7.5. You attempted with " + str(logg_star))
	else: print(logg_star); print(logg_star); print(logg_star);  raise ValueError
	
	
	
	Teff_all = Teff_all_logg;  wl_all = wl_all_logg;     flux_all = flux_all_logg;      logg_all = logg_all_logg
	
	

	mask_logg = (logg_all<=maxval) & (logg_all>=minval)
	Grav_N = logg_all[mask_logg];    wl_all_N=wl_all[mask_logg];    flux_N=flux_all[mask_logg];    Teff_N=Teff_all[mask_logg]
	
	
	unique_Teffs_DC=np.unique(Teff_N)
	
	temdiff = temperature_star - unique_Teffs_DC
	Teff_min=unique_Teffs_DC[np.argwhere(temdiff==np.amin(temdiff[temdiff>0]))[0][0]]
	Teff_max=unique_Teffs_DC[np.argwhere(temdiff==np.amax(temdiff[temdiff<0]))[0][0]]
	
	
	mask_logg_then_T = (Teff_N<=Teff_max) & (Teff_N>=Teff_min)
	
	
	#Grav_N = Grav_N[mask_logg_then_T]
	
	if len(Grav_N[mask_logg_then_T])==0:
		raise ValueError("ERROR - Grid of DBA/DB/DC")
	
	Grav1_N = Grav_N[mask_logg_then_T];    wl_all1_N=wl_all_N[mask_logg_then_T];    flux1_N=flux_N[mask_logg_then_T];    Teff1_N=Teff_N[mask_logg_then_T]
	return Grav1_N, wl_all1_N, flux1_N, Teff1_N



@njit
def return_DBAgrids(temperature_star, logg_star, HoverHestar):
	# find nearest Teff in models in the DA/DBA/DC grid, then take the two for interpolation limits
	if starType1=="DBA":
		minval=10000;   maxval=-10000
		if logg_star>=7.5 and logg_star<8: minval=7.5; maxval=8
		elif logg_star>=8 and logg_star<8.5: minval=8; maxval=8.5
		elif logg_star>=8.5 and logg_star<=9: minval=8.5; maxval=9
		elif logg_star<7.5: raise ValueError("Minimum logg for grid is 7.5. You attempted with " + str(logg_star))
		else: print(logg_star); print(logg_star); print(logg_star);  raise ValueError
		
		Teff_all = Teff_all_logg;  wl_all = wl_all_logg;     flux_all = flux_all_logg;      logg_all = logg_all_logg;   HoverHe_all = HoverHe_all_logg_DBA
		
		#if HoverHestar>=2 and HoverHestar<3: HoverHe_min=2; HoverHe_max=3
		#elif HoverHestar>=3 and HoverHestar<4: HoverHe_min=3; HoverHe_max=4
		#elif HoverHestar>=4 and HoverHestar<5: HoverHe_min=4; HoverHe_max=5
		#elif HoverHestar>=5 and HoverHestar<7: HoverHe_min=5; HoverHe_max=7
		#elif HoverHestar>=7 and HoverHestar<=30: HoverHe_min=7; HoverHe_max=30
		#else: print(HoverHestar);  print(HoverHestar);  raise ValueError
		
		if HoverHestar>=2 and HoverHestar<=5: HoverHe_min=2; HoverHe_max=5
		elif HoverHestar>=5 and HoverHestar<=8: HoverHe_min=5; HoverHe_max=8
		elif HoverHestar>=8 and HoverHestar<=30: HoverHe_min=8; HoverHe_max=30
		else: print(HoverHestar);  print(HoverHestar);  raise ValueError
		
		mask_logg = (logg_all<=maxval) & (logg_all>=minval) & (HoverHe_all>=HoverHe_min) & (HoverHe_all<=HoverHe_max)
		Grav_N = logg_all[mask_logg];    wl_all_N=wl_all[mask_logg];    flux_N=flux_all[mask_logg];    Teff_N=Teff_all[mask_logg];     HoverHe_N = HoverHe_all[mask_logg]
		
		
		unique_Teffs_DBA=np.unique(Teff_N)
	
		temdiff = temperature_star - unique_Teffs_DBA
		if np.amax(temdiff)>0:    Teff_min=unique_Teffs_DBA[np.argwhere(temdiff==np.amin(temdiff[temdiff>0]))[0][0]]
		else: raise ValueError("DB grid temperature grid out of bounds")
		if np.amin(temdiff)<0:    Teff_max=unique_Teffs_DBA[np.argwhere(temdiff==np.amax(temdiff[temdiff<0]))[0][0]]
		else: raise ValueError("DB grid temperature grid out of bounds")
		
		
		#mask_logg = (Teff_N<=Teff_max+1000) & (Teff_N>=Teff_min-1000)
		mask_logg_then_T = (Teff_N<=Teff_max) & (Teff_N>=Teff_min)
		
		
		#Grav_N = Grav_N[mask_logg_then_T]
		
		if len(Grav_N[mask_logg_then_T])==0:   raise ValueError("ERROR - Grid of DBA/DB/DC")
		
		return Grav_N[mask_logg_then_T], wl_all_N[mask_logg_then_T], flux_N[mask_logg_then_T], Teff_N[mask_logg_then_T], HoverHe_N[mask_logg_then_T]
		
		
	else:
		Teff_all = Teff_all_logg;  wl_all = wl_all_logg;     flux_all = flux_all_logg;      logg_all = logg_all_logg;   HoverHe_all = HoverHe_all_logg
		
		
		Teff_min = np.unique(Teff_all)[np.unique(Teff_all)<=temperature_star][-1]
		Teff_max = np.unique(Teff_all)[np.unique(Teff_all)>=temperature_star][0]
		logg_min = np.unique(logg_all)[np.unique(logg_all)<=logg_star][-1]
		logg_max = np.unique(logg_all)[np.unique(logg_all)>=logg_star][0]
		HoverHe_min = np.unique(HoverHe_all)[np.unique(HoverHe_all)<=HoverHestar][-1]
		HoverHe_max = np.unique(HoverHe_all)[np.unique(HoverHe_all)>=HoverHestar][0]
		
		mask = (Teff_all<=Teff_max) & (Teff_all>=Teff_min)  &  (logg_all<=logg_max) & (logg_all>=logg_min)  &  (HoverHe_all<=HoverHe_max) & (HoverHe_all>=HoverHe_min)
		
		
		if len(logg_all[mask])==0:   raise ValueError("ERROR - Grid of subdwarfs")
		
		return logg_all[mask], wl_all[mask], flux_all[mask], Teff_all[mask], HoverHe_all[mask]
	
	
	
	
	
@njit
def return_DCgrids(temperature_star, logg_star):
	# find nearest Teff in models in the DA/DBA/DC grid, then take the two for interpolation limits
	minval=10000;   maxval=-10000
	if logg_star>=7.5 and logg_star<8: minval=7.5; maxval=8
	elif logg_star>=8 and logg_star<8.5: minval=8; maxval=8.5
	elif logg_star>=8.5 and logg_star<=9: minval=8.5; maxval=9
	elif logg_star<7.5: raise ValueError("Minimum logg for grid is 7.5. You attempted with " + str(logg_star))
	else: print(logg_star); print(logg_star); print(logg_star);  raise ValueError
	
	
	
	Teff_all = Teff_all_logg;  wl_all = wl_all_logg;     flux_all = flux_all_logg;      logg_all = logg_all_logg
	
	

	mask_logg = (logg_all<=maxval) & (logg_all>=minval)
	Grav_N = logg_all[mask_logg];    wl_all_N=wl_all[mask_logg];    flux_N=flux_all[mask_logg];    Teff_N=Teff_all[mask_logg]
	
	
	unique_Teffs_DC=np.unique(Teff_N)
	
	temdiff = temperature_star - unique_Teffs_DC
	Teff_min=unique_Teffs_DC[np.argwhere(temdiff==np.amin(temdiff[temdiff>0]))[0][0]]
	Teff_max=unique_Teffs_DC[np.argwhere(temdiff==np.amax(temdiff[temdiff<0]))[0][0]]
	
	
	mask_logg_then_T = (Teff_N<=Teff_max) & (Teff_N>=Teff_min)
	
	
	#Grav_N = Grav_N[mask_logg_then_T]
	
	if len(Grav_N[mask_logg_then_T])==0:
		raise ValueError("ERROR - Grid of DBA/DB/DC")
	
	Grav1_N = Grav_N[mask_logg_then_T];    wl_all1_N=wl_all_N[mask_logg_then_T];    flux1_N=flux_N[mask_logg_then_T];    Teff1_N=Teff_N[mask_logg_then_T]
	return Grav1_N, wl_all1_N, flux1_N, Teff1_N




@njit
def initiate_interp(ww, wl_all1_N, Grav1_N, Teff1_N,flux1_N, unique_Ts, unique_Gs): # return the fluxes that correspond to the wavelength desired for each T/logg grid
	ffffsss=[]
	maska = (wl_all1_N==ww)
	for tt in unique_Ts:
		maskb= maska & (tt==Teff1_N)
		for gg in unique_Gs:
			ffffsss.append(flux1_N[maskb & (Grav1_N==gg)][0])
	
	return ffffsss



@jit
def bilinear_interpolation(x, y, points):
    # https://stackoverflow.com/questions/8661537/how-to-perform-bilinear-interpolation-in-python
    # gives near enough the same result as scipy griddata, but is faster.
    '''Interpolate (x,y) from values associated with four points.

    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.

        >>> bilinear_interpolation(12, 5.5,
        ...                        [(10, 4, 100),
        ...                         (20, 4, 200),
        ...                         (10, 6, 150),
        ...                         (20, 6, 300)])
        165.0

    '''
    # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

    #  points = sorted(points)               # order points by x, then by y    # JM commented in my script because it's already passed this way
    
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points
    
    
    # JM this is a test: not necessary, but bugged fixed for me!
    #if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
    #    raise ValueError('points do not form a rectangle', x1, x2, y1, y2, _x1, _x2, _y1, _y2, x, y)
    #if not x1 <= x <= x2 or not y1 <= y <= y2:
    #    raise ValueError('(x, y) not within the rectangle', x1, x2, y1, y2, _x1, _x2, _y1, _y2, x, y)

    return (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)
           ) / ((x2 - x1) * (y2 - y1) + 0.0)



@njit
def return_model_spectrum_DA(wl_all1_N, ref_wl, cut_limits_min, cut_limits_max, Grav1_N, flux1_N, Teff1_N, temperature_star, logg_star):
	if ref_wl>6000:    excess_slack=10 # angstroms to allow variation with RV of star. Here, Halpha goes +-450kms-1. Increase this if RV diff larger
	else:  excess_slack=7.5 # angstroms to allow variation with RV of star. Here, Hbeta goes +-500kms-1. Increase this if RV diff larger
	
	mask_logg_wl = ((wl_all1_N > ref_wl+cut_limits_min-excess_slack) & (wl_all1_N < ref_wl+cut_limits_max+excess_slack))
	Grav_N_N = Grav1_N[mask_logg_wl];   wl_all_N_N=wl_all1_N[mask_logg_wl];    flux_N_N=flux1_N[mask_logg_wl];    Teff_N_N=Teff1_N[mask_logg_wl]
	
	# interpolate for a model at the reference wavelength with this mcmc interation	
	wl_grid, unique_Ts, unique_Gs = np.unique(wl_all_N_N), np.unique(Teff_N_N), np.unique(Grav_N_N)
	
	#ffffsss1 = flux_N_N[(Teff_N_N==unique_Ts[0]) & (Grav_N_N==unique_Gs[0])]
	#ffffsss2 = flux_N_N[(Teff_N_N==unique_Ts[0]) & (Grav_N_N==unique_Gs[1])]
	#ffffsss3 = flux_N_N[(Teff_N_N==unique_Ts[1]) & (Grav_N_N==unique_Gs[0])]
	#ffffsss4 = flux_N_N[(Teff_N_N==unique_Ts[1]) & (Grav_N_N==unique_Gs[1])]
	
	
	maskT0 = Teff_N_N==unique_Ts[0]
	maskT1 = Teff_N_N==unique_Ts[1]
	maskG0 = Grav_N_N==unique_Gs[0]
	maskG1 = Grav_N_N==unique_Gs[1]
	ffffsss1 = flux_N_N[maskT0 & maskG0]
	ffffsss2 = flux_N_N[maskT0 & maskG1]
	ffffsss3 = flux_N_N[maskT1 & maskG0]
	ffffsss4 = flux_N_N[maskT1 & maskG1]
	
	
	try:
		model_spectrum = (ffffsss1 * (unique_Ts[1] - temperature_star) * (unique_Gs[1] - logg_star) +           ffffsss3 * (temperature_star - unique_Ts[0]) * (unique_Gs[1] - logg_star) +            ffffsss2 * (unique_Ts[1] - temperature_star) * (logg_star - unique_Gs[0]) +            ffffsss4 * (temperature_star - unique_Ts[0]) * (logg_star - unique_Gs[0])           ) / ((unique_Ts[1] - unique_Ts[0]) * (unique_Gs[1] - unique_Gs[0]))
	except:	raise ValueError(unique_Ts, unique_Gs, unique_Ts[0], unique_Ts[1], unique_Gs[0], unique_Gs[1], len(ffffsss1), len(ffffsss2), len(ffffsss3), len(ffffsss4), temperature_star, logg_star)

	#plt.plot(wl_grid, model_spectrum,c='r', ls='--');  plt.show()
	
	
	#plt.plot(wl_grid, model_spectrum,c='r'); plt.plot(wl_grid, model_spectrum_starx,c='k'); plt.show()
	
	if ref_wl>6500:     fine_grid=np.linspace(np.amin(wl_grid),np.amax(wl_grid),int((np.amax(wl_grid) - np.amin(wl_grid))*50)) # 0.02AA spacing
	elif ref_wl>4500:   fine_grid=np.linspace(np.amin(wl_grid),np.amax(wl_grid),int((np.amax(wl_grid) - np.amin(wl_grid))*33)) # 0.03AA spacing
	elif ref_wl>4200:   fine_grid=np.linspace(np.amin(wl_grid),np.amax(wl_grid),int((np.amax(wl_grid) - np.amin(wl_grid))*20)) # 0.05AA spacing
	else:               fine_grid=np.linspace(np.amin(wl_grid),np.amax(wl_grid),int((np.amax(wl_grid) - np.amin(wl_grid))*10)) # 0.1AA spacing
	
	
	#fine_grid=np.linspace(np.amin(wl_grid),np.amax(wl_grid),int((np.amax(wl_grid) - np.amin(wl_grid))*400)) # 0.0025AA spacing
	model_spectrum=interp(fine_grid, wl_grid, model_spectrum)
	
	return fine_grid, model_spectrum
	


def trilinear_interpolation(wavelengths_model, temps_model, loggs_model, HoverHes_model, fluxes_model, t_star, logg_star, HoverHe_star):
	# interpolate for a model at the reference wavelength with this mcmc interation	
	unique_Ts, unique_Gs, unique_HoverHes = np.unique(temps_model), np.unique(loggs_model), np.unique(HoverHes_model)
	
	Temp_d = (t_star-unique_Ts[0]) / (unique_Ts[1]-unique_Ts[0])
	logg_d = (logg_star-unique_Gs[0]) / (unique_Gs[1]-unique_Gs[0])
	HoverHe_d = (HoverHe_star-unique_HoverHes[0]) / (unique_HoverHes[1]-unique_HoverHes[0])
	
	
	# https://en.wikipedia.org/wiki/Trilinear_interpolation
	#### worked example:
	# t1, t2 = 27000, 28000
	# logg1, logg2 = 5, 5.2
	# HoverHe1, HoverHe2 = -4.55, -4.05
	
	# c000 ->   t1 = 27000, logg1 = 5, HoverHe1 = -4.55         ->   t1 = unique_Ts[0], logg1 = unique_Gs[0], HoverHe1=unique_HoverHes[0]
	# c001 ->   t1 = 27000, logg1 = 5, HoverHe1 = -4.05         ->   t1 = unique_Ts[0], logg1 = unique_Gs[0], HoverHe1=unique_HoverHes[1]
	# c010 ->   t1 = 27000, logg1 = 5.2, HoverHe1 = -4.55       ->   t1 = unique_Ts[0], logg1 = unique_Gs[1], HoverHe1=unique_HoverHes[0]
	# c011 ->   t1 = 27000, logg1 = 5.2, HoverHe1 = -4.05       ->   t1 = unique_Ts[0], logg1 = unique_Gs[1], HoverHe1=unique_HoverHes[1]
	
	
	# c100 ->   t1 = 28000, logg1 = 5, HoverHe1 = -4.55         ->   t1 = unique_Ts[1], logg1 = unique_Gs[0], HoverHe1=unique_HoverHes[0]
	# c101 ->   t1 = 28000, logg1 = 5, HoverHe1 = -4.05         ->   t1 = unique_Ts[1], logg1 = unique_Gs[0], HoverHe1=unique_HoverHes[1]
	# c110 ->   t1 = 28000, logg1 = 5.2, HoverHe1 = -4.55       ->   t1 = unique_Ts[1], logg1 = unique_Gs[1], HoverHe1=unique_HoverHes[0]
	# c111 ->   t1 = 28000, logg1 = 5.2, HoverHe1 = -4.05       ->   t1 = unique_Ts[1], logg1 = unique_Gs[1], HoverHe1=unique_HoverHes[1]
	
	
	flux000 = fluxes_model[(temps_model==unique_Ts[0]) & (loggs_model==unique_Gs[0]) & (HoverHes_model==unique_HoverHes[0])]
	flux001 = fluxes_model[(temps_model==unique_Ts[0]) & (loggs_model==unique_Gs[0]) & (HoverHes_model==unique_HoverHes[1])]
	flux010 = fluxes_model[(temps_model==unique_Ts[0]) & (loggs_model==unique_Gs[1]) & (HoverHes_model==unique_HoverHes[0])]
	flux011 = fluxes_model[(temps_model==unique_Ts[0]) & (loggs_model==unique_Gs[1]) & (HoverHes_model==unique_HoverHes[1])]
	
	flux100 = fluxes_model[(temps_model==unique_Ts[1]) & (loggs_model==unique_Gs[0]) & (HoverHes_model==unique_HoverHes[0])]
	flux101 = fluxes_model[(temps_model==unique_Ts[1]) & (loggs_model==unique_Gs[0]) & (HoverHes_model==unique_HoverHes[1])]
	flux110 = fluxes_model[(temps_model==unique_Ts[1]) & (loggs_model==unique_Gs[1]) & (HoverHes_model==unique_HoverHes[0])]
	flux111 = fluxes_model[(temps_model==unique_Ts[1]) & (loggs_model==unique_Gs[1]) & (HoverHes_model==unique_HoverHes[1])]
	
	
	
	c00 = flux000 * (1-Temp_d)  +  flux100*Temp_d
	c01 = flux001 * (1-Temp_d)  +  flux101*Temp_d
	c10 = flux010 * (1-Temp_d)  +  flux110*Temp_d
	c11 = flux011 * (1-Temp_d)  +  flux111*Temp_d
	
	
	c0 = c00*(1-logg_d) + (c10*logg_d)
	c1 = c01*(1-logg_d) + (c11*logg_d)
	
	c = c0*(1-HoverHe_d) + c1*HoverHe_d
	
	
	return c



#@njit
# temperature grid for DBs is irregular. Don't try to speed it up and accept the griddata way
def return_model_spectrum_DBA(wl_all1_N, ref_wl, cut_limits_min, cut_limits_max, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, temperature_star, logg_star, HoverHe_star):
	if ref_wl>6000:    excess_slack=10 # angstroms to allow variation with RV of star. Here, Halpha goes +-450kms-1. Increase this if RV diff larger
	else:  excess_slack=7.5 # angstroms to allow variation with RV of star. Here, Hbeta goes +-500kms-1. Increase this if RV diff larger
	
	mask_logg_wl = ((wl_all1_N > ref_wl+cut_limits_min-excess_slack) & (wl_all1_N < ref_wl+cut_limits_max+excess_slack))
	Grav_N_N = Grav1_N[mask_logg_wl];   wl_all_N_N=wl_all1_N[mask_logg_wl];    flux_N_N=flux1_N[mask_logg_wl];    Teff_N_N=Teff1_N[mask_logg_wl];    HoverHe_N_N = HoverHe1_N[mask_logg_wl]
	
	
	
	
	# interpolate for a model at the reference wavelength with this mcmc interation	
	wl_grid, unique_Ts, unique_Gs, unique_HoverHes = np.unique(wl_all_N_N), np.unique(Teff_N_N), np.unique(Grav_N_N), np.unique(HoverHe_N_N)
	model_spectrum=griddata(np.array([Teff_N_N, wl_all_N_N, Grav_N_N, HoverHe_N_N]).T, flux_N_N, np.array([np.full((len(wl_grid),),temperature_star), wl_grid, np.full((len(wl_grid),),logg_star), np.full((len(wl_grid),),HoverHe_star)]).T, method="linear")
	
	if ref_wl>6500:     fine_grid=np.linspace(np.amin(wl_grid),np.amax(wl_grid),int((np.amax(wl_grid) - np.amin(wl_grid))*10)) # 0.1AA spacing
	elif ref_wl>4500:   fine_grid=np.linspace(np.amin(wl_grid),np.amax(wl_grid),int((np.amax(wl_grid) - np.amin(wl_grid))*10)) # 0.1AA spacing
	elif ref_wl>4200:   fine_grid=np.linspace(np.amin(wl_grid),np.amax(wl_grid),int((np.amax(wl_grid) - np.amin(wl_grid))*10)) # 0.1AA spacing
	else:               fine_grid=np.linspace(np.amin(wl_grid),np.amax(wl_grid),int((np.amax(wl_grid) - np.amin(wl_grid))*10)) # 0.1AA spacing
	
	model_spectrum=np.interp(fine_grid, wl_grid, model_spectrum)
	
	return fine_grid, model_spectrum
	
	
@njit
def return_model_spectrum_subdwarf(wl_all1_N, ref_wl, cut_limits_min, cut_limits_max, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, temperature_star, logg_star, HoverHe_star):
	if ref_wl>6000:    excess_slack=10 # angstroms to allow variation with RV of star. Here, Halpha goes +-450kms-1. Increase this if RV diff larger
	else:  excess_slack=7.5 # angstroms to allow variation with RV of star. Here, Hbeta goes +-500kms-1. Increase this if RV diff larger
	
	mask_logg_wl = ((wl_all1_N > ref_wl+cut_limits_min-excess_slack) & (wl_all1_N < ref_wl+cut_limits_max+excess_slack))
	Grav_N_N = Grav1_N[mask_logg_wl];   wl_all_N_N=wl_all1_N[mask_logg_wl];    flux_N_N=flux1_N[mask_logg_wl];    Teff_N_N=Teff1_N[mask_logg_wl];    HoverHe_N_N = HoverHe1_N[mask_logg_wl]
	
	
	# interpolate for a model at the reference wavelength with this mcmc interation	
	wl_grid, unique_Ts, unique_Gs, unique_HoverHes = np.unique(wl_all_N_N), np.unique(Teff_N_N), np.unique(Grav_N_N), np.unique(HoverHe_N_N)
	
	Temp_d = (temperature_star-unique_Ts[0]) / (unique_Ts[1]-unique_Ts[0])
	logg_d = (logg_star-unique_Gs[0]) / (unique_Gs[1]-unique_Gs[0])
	HoverHe_d = (HoverHe_star-unique_HoverHes[0]) / (unique_HoverHes[1]-unique_HoverHes[0])
	
	
	# https://en.wikipedia.org/wiki/Trilinear_interpolation
	#### worked example:
	# t1, t2 = 27000, 28000
	# logg1, logg2 = 5, 5.2
	# HoverHe1, HoverHe2 = -4.55, -4.05
	
	# c000 ->   t1 = 27000, logg1 = 5, HoverHe1 = -4.55         ->   t1 = unique_Ts[0], logg1 = unique_Gs[0], HoverHe1=unique_HoverHes[0]
	# c001 ->   t1 = 27000, logg1 = 5, HoverHe1 = -4.05         ->   t1 = unique_Ts[0], logg1 = unique_Gs[0], HoverHe1=unique_HoverHes[1]
	# c010 ->   t1 = 27000, logg1 = 5.2, HoverHe1 = -4.55       ->   t1 = unique_Ts[0], logg1 = unique_Gs[1], HoverHe1=unique_HoverHes[0]
	# c011 ->   t1 = 27000, logg1 = 5.2, HoverHe1 = -4.05       ->   t1 = unique_Ts[0], logg1 = unique_Gs[1], HoverHe1=unique_HoverHes[1]
	
	
	# c100 ->   t1 = 28000, logg1 = 5, HoverHe1 = -4.55         ->   t1 = unique_Ts[1], logg1 = unique_Gs[0], HoverHe1=unique_HoverHes[0]
	# c101 ->   t1 = 28000, logg1 = 5, HoverHe1 = -4.05         ->   t1 = unique_Ts[1], logg1 = unique_Gs[0], HoverHe1=unique_HoverHes[1]
	# c110 ->   t1 = 28000, logg1 = 5.2, HoverHe1 = -4.55       ->   t1 = unique_Ts[1], logg1 = unique_Gs[1], HoverHe1=unique_HoverHes[0]
	# c111 ->   t1 = 28000, logg1 = 5.2, HoverHe1 = -4.05       ->   t1 = unique_Ts[1], logg1 = unique_Gs[1], HoverHe1=unique_HoverHes[1]
	
	
	
	mask_Ts_equal_uniqueT0 = Teff_N_N==unique_Ts[0]
	mask_Ts_equal_uniqueT1 = Teff_N_N==unique_Ts[1]
	mask_Gs_equal_uniqueG0 = Grav_N_N==unique_Gs[0]
	mask_Gs_equal_uniqueG1 = Grav_N_N==unique_Gs[1]
	mask_Hs_equal_uniqueH0 = HoverHe_N_N==unique_HoverHes[0]
	mask_Hs_equal_uniqueH1 = HoverHe_N_N==unique_HoverHes[1]
	
	
	
	c00 = flux_N_N[mask_Ts_equal_uniqueT0 & mask_Gs_equal_uniqueG0 & mask_Hs_equal_uniqueH0]   * (1-Temp_d)  +  flux_N_N[mask_Ts_equal_uniqueT1 & mask_Gs_equal_uniqueG0 & mask_Hs_equal_uniqueH0] *Temp_d
	c01 = flux_N_N[mask_Ts_equal_uniqueT0 & mask_Gs_equal_uniqueG0 & mask_Hs_equal_uniqueH1]   * (1-Temp_d)  +  flux_N_N[mask_Ts_equal_uniqueT1 & mask_Gs_equal_uniqueG0 & mask_Hs_equal_uniqueH1] *Temp_d
	c10 = flux_N_N[mask_Ts_equal_uniqueT0 & mask_Gs_equal_uniqueG1 & mask_Hs_equal_uniqueH0]   * (1-Temp_d)  +  flux_N_N[mask_Ts_equal_uniqueT1 & mask_Gs_equal_uniqueG1 & mask_Hs_equal_uniqueH0] *Temp_d
	c11 = flux_N_N[mask_Ts_equal_uniqueT0 & mask_Gs_equal_uniqueG1 & mask_Hs_equal_uniqueH1]   * (1-Temp_d)  +  flux_N_N[mask_Ts_equal_uniqueT1 & mask_Gs_equal_uniqueG1 & mask_Hs_equal_uniqueH1] *Temp_d
	
	
	c0 = c00*(1-logg_d) + (c10*logg_d)
	c1 = c01*(1-logg_d) + (c11*logg_d)
	
	model_spectrum = c0*(1-HoverHe_d) + c1*HoverHe_d
	
	
	fine_grid=np.linspace(np.amin(wl_grid),np.amax(wl_grid),int((np.amax(wl_grid) - np.amin(wl_grid))*10)) # 0.1AA spacing
	model_spectrum=np.interp(fine_grid, wl_grid, model_spectrum)
	
	
	return fine_grid, model_spectrum



if forced_Scaling=="WD":
	install_path = os.environ['WD_BASS_INSTALL_DIR']
	loaded_Althaus = np.load(install_path + "/saved_MTR/Althaus_2013_full_nomasses.npy")
	loaded_Istrate = np.load(install_path + "/saved_MTR/Istrate_Z0p02_diffusion_nomasses.npy")
	loaded_CO = np.load(install_path + "/saved_MTR/table_valuesCO.npy")


checkT1, checklogg1, checkHoverHe1, checkscaling, checkR, checkDummy = False, False, False, False, False, False
if forced_teff1==0: checkT1=True
if forced_logg1==0: checklogg1=True
if forced_HoverHe1==0 and (starType1=="DBA" or starType1.startswith("sd")): checkHoverHe1=True
if forced_Scaling==False: checkscaling=True
if fit_phot_SED and starType1.startswith("sd"): checkR = True
if "Dummy" in p0labels: checkDummy=True





#Test that next step falls within boundaries allowed by the priors
def lnprior(theta, arguments):
	input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, normalised_wavelength, normalised_flux, normalised_err,  inp_resolution, used_RV_boundaries, HJD_values, sigma_clip, the_unique_wavelengths, spec1wl, spec1flux = arguments
	if "T1" in p0labels:            args = np.argwhere(p0labels=="T1")[0][0];          T1 = theta[args]
	if "logg1" in p0labels:         args = np.argwhere(p0labels=="logg1")[0][0];       logg1 = theta[args]
	if "H/He1" in p0labels:         args = np.argwhere(p0labels=="H/He1")[0][0];       HoverHe1 = theta[args]
	if "K1" in p0labels:            args = np.argwhere(p0labels=="K1")[0][0];          mcmc_K1 = theta[args]
	if "Vg1" in p0labels:           args = np.argwhere(p0labels=="Vg1")[0][0];         mcmc_Vgamma1 = theta[args]
	if "Scaling" in p0labels:       args = np.argwhere(p0labels=="Scaling")[0][0];     Scaling = theta[args]
	if "RV1_0" in p0labels:         num_start_RVs = np.argwhere(p0labels=="RV1_0")[0][0]
	elif "RV1_1" in p0labels:       num_start_RVs = np.argwhere(p0labels=="RV1_1")[0][0]
	elif "RV1_2" in p0labels:       num_start_RVs = np.argwhere(p0labels=="RV1_2")[0][0]
	elif "RV1_3" in p0labels:       num_start_RVs = np.argwhere(p0labels=="RV1_3")[0][0]
	elif "RV1_4" in p0labels:       num_start_RVs = np.argwhere(p0labels=="RV1_4")[0][0]
	if "Parallax" in p0labels:      args = np.argwhere(p0labels=="Parallax")[0][0];     mcmc_parallax = theta[args]
	if "R" in p0labels:       args = np.argwhere(p0labels=="R")[0][0];     mcmc_R = theta[args]
	if "Dummy" in p0labels:   args = np.argwhere(p0labels=="Dummy")[0][0];     mcmc_Dummy = theta[args]
	
	
	if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
		if forced_Scaling==False: RV=theta[num_start_RVs:-1];  Scaling=theta[-1]
		else:                     RV=theta[num_start_RVs:]
		
		
		
	passed = True
	if checkT1:  
		if not p0T1[0]<T1<p0T1[1]: passed = False
	if checklogg1:  
		if not p0logg1[0]<logg1<p0logg1[1]: passed = False
	if checkHoverHe1:  
		if not p0HoverHe1[0]<HoverHe1<p0HoverHe1[1]: passed = False
	if checkscaling:
		if not p0scaling[0]<Scaling<p0scaling[1]: passed = False
	if checkR:
		if not p0R[0]<mcmc_R<p0R[1]: passed = False
	if checkDummy:
		if not dummybound[0]<mcmc_Dummy<dummybound[1]: passed=False
	if passed==False:  return -np.inf
	
	
	if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
		for anRV, anRVbound in zip(RV,used_RV_boundaries):  # go through all individual RVs and make sure they are in boundaries
			if not anRVbound[0] < anRV < anRVbound[1]:        return -np.inf
	else:
		if forced_K1=="Fit":
			if forced_Vgamma1=="Fit":
				if not (p0K1[0] < mcmc_K1 < p0K1[1] and p0Vgamma1[0] < mcmc_Vgamma1 < p0Vgamma1[1]):    return -np.inf
			else:
				if not p0K1[0] < mcmc_K1 < p0K1[1]:       return -np.inf
		
		if forced_K1!="Fit" and forced_Vgamma1=="Fit":
			if not p0Vgamma1[0] < mcmc_Vgamma1 < p0Vgamma1[1]:    return -np.inf
		
	
	if "Parallax" in p0labels: 
		parallax_prior = np.log(1.0/(np.sqrt(2*np.pi)*plax_unc))-0.5*(mcmc_parallax-plax)**2/plax_unc**2 + np.log(1.0/(np.sqrt(2*np.pi)*plax_unc))-0.5*(mcmc_parallax-plax)**2/plax_unc**2
		return parallax_prior
	return 0.0
	



def lnlike(theta, arguments):
	input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, normalised_wavelength, normalised_flux, normalised_err,  inp_resolution, used_RV_boundaries, HJD_values, sigma_clip, the_unique_wavelengths, spec1wl, spec1flux = arguments
	if type(input_files)==np.str_:
		input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, normalised_wavelength, normalised_flux, normalised_err,  inp_resolution, used_RV_boundaries, HJD_values, sigma_clip = [input_files], [share_rv], [reference_wl], [cut_Ha_all], [normaliseHa_all], [normalised_wavelength], [normalised_flux], [normalised_err],  [inp_resolution], [used_RV_boundaries], [HJD_values], [sigma_clip]
	
	if "T1" in p0labels:            args = np.argwhere(p0labels=="T1")[0][0];          T1 = theta[args]
	else:                           T1 = forced_teff1
	
	if "logg1" in p0labels:         args = np.argwhere(p0labels=="logg1")[0][0];       logg1 = theta[args]
	else:                           logg1 = forced_logg1

	if "H/He1" in p0labels:         args = np.argwhere(p0labels=="H/He1")[0][0];       HoverHe1 = theta[args]
	else:  
		try: HoverHe1 = forced_HoverHe1
		except: None
	
	if "K1" in p0labels:            args = np.argwhere(p0labels=="K1")[0][0];          mcmc_K1 = theta[args]
	else:                           mcmc_K1 = forced_K1
	
	if "Vg1" in p0labels:           args = np.argwhere(p0labels=="Vg1")[0][0];         mcmc_Vgamma1 = theta[args]
	else:                           mcmc_Vgamma1 = forced_Vgamma1
	
	if "Parallax" in p0labels:      args = np.argwhere(p0labels=="Parallax")[0][0];    mcmc_parallax = theta[args]
	
	if "Scaling" in p0labels:       args = np.argwhere(p0labels=="Scaling")[0][0];     Scaling = theta[args]
	else:                           Scaling = forced_Scaling

	if "R" in p0labels:       args = np.argwhere(p0labels=="R")[0][0];     mcmc_R = theta[args]
	
	if "RV1_0" in p0labels:         num_start_RVs = np.argwhere(p0labels=="RV1_0")[0][0]
	
	
	if not type(spec1flux).__module__ == np.__name__:
		if starType1=="DA":
			if not pier_or_antoine=="pier3Dphot_antoine1Dspec":
				Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DAgrids(T1, logg1)
			else:
				Grav1_N, wl_all1_N, flux1_N, Teff1_N, photGrav1_N, photwl_all1_N, photflux1_N, photTeff1_N = return_DAgrids(T1, logg1)
				if len(wl_all1_N)==0:   raise ValueError(len(wl_all1_N), len(Grav1_N), len(flux1_N), len(Teff1_N), T1, logg1)
			
		elif starType1=="DBA" or starType1.startswith("sd"):
			Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N = return_DBAgrids(T1, logg1, HoverHe1)
		elif starType1=="DB":
			Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DBgrids(T1, logg1)
		elif starType=="DC":
			Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DCgrids(T1, logg1)
		
		
	
	if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
		if forced_Scaling==False: RV=theta[num_start_RVs:-1];  Scaling=theta[-1]
		else:                     RV=theta[num_start_RVs:]
	else:
		if forced_Scaling==False: Scaling=theta[-1]
	
	
	if fit_phot_SED:
		if forced_Scaling=="WD" and starType1=="DA":
			R1 = get_MTR(T1, logg=logg1, return_R_from_T_logg=True, loaded_Istrate=loaded_Istrate, loaded_CO=loaded_CO, loaded_Althaus=loaded_Althaus)
			
			if np.isnan(R1): raise ValueError(T1, logg1)
		if forced_Scaling=="WD" and starType1=="DBA" or starType1=="DB" or starType1=="DC":
			R1 = get_MTR_DB(T1, logg=logg1)
			
			if np.isnan(R1): raise ValueError(T1, logg1)
		elif isinstance(forced_Scaling, float): Scaling=forced_Scaling
		
		if starType1=="DA" or starType1=="DB" or starType1=="DC":
			smeared_wl, smeared_flux = Fit_phot.fit_phot_SED_single(Grav1_N, wl_all1_N, flux1_N, Teff1_N, None, T1, logg1, None, min_wl=theminww, max_wl=themaxww, starType1=starType1, R1=R1, parallax=mcmc_parallax, red=reddening_Ebv)
		elif starType1=="DBA":
			smeared_wl, smeared_flux = Fit_phot.fit_phot_SED_single(Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N, T1, logg1, HoverHe1, min_wl=theminww, max_wl=themaxww, starType1=starType1, R1=R1, parallax=mcmc_parallax, red=reddening_Ebv)
		elif starType1.startswith("sd"):
			smeared_wl, smeared_flux = Fit_phot.fit_phot_SED_single(Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N, T1, logg1, HoverHe1, min_wl=theminww, max_wl=themaxww, starType1=starType1, R1=mcmc_R, parallax=mcmc_parallax, red=reddening_Ebv)
		
		
		rchisq_phot, chisq_phot = Fit_phot.process_photometry_in_each_pb(smeared_wl, smeared_flux, sedfilter, sed_wl, sedflux, sedfluxe, filter_dict=filter_dict, theminww_plot=theminww+50, themaxww_plot=themaxww+50, single_or_double="single", return_points_for_phot_model=False)
		
		if np.round(chisq_phot,7)==0:
			plt.plot(smeared_wl, smeared_flux)
			plt.title(str(T1) + "  "+ str(logg1))
			plt.show()
			raise ValueError
	
	
	
	
	chisq_spec=0
	
	index = np.arange(0,len(input_files),1)
	
	if fit_phot_SED:
		no_need_to_recompute = np.amin(smeared_wl) < np.amin(the_unique_wavelengths)+np.amin(cut_Ha_all.T) - 10 and np.amax(smeared_wl) > np.amax(the_unique_wavelengths)+np.amax(cut_Ha_all.T) + 10
	elif type(spec1flux).__module__ == np.__name__: no_need_to_recompute=True
	else:   no_need_to_recompute = False
	
	for ref_wl in the_unique_wavelengths:
		mask_ref_wl = reference_wl==ref_wl
		ref_wl_cut_lims_min, ref_wl_cut_lims_max = cut_Ha_all[mask_ref_wl].T
		cut_wl_min, cut_wl_max = np.amin(ref_wl_cut_lims_min), np.amax(ref_wl_cut_lims_max)
		
		
		
		if no_need_to_recompute:
			if fit_phot_SED:
				mask =  (smeared_wl>ref_wl+cut_wl_min-10) & (smeared_wl<ref_wl+cut_wl_max+10)
				model_wl1, model_spectrum_star1 = smeared_wl[mask], smeared_flux[mask]
			elif type(spec1flux).__module__ == np.__name__:
				mask1 =  (spec1wl>ref_wl+cut_wl_min-10) & (spec1wl<ref_wl+cut_wl_max+10)
				model_wl1, model_spectrum_star1 = spec1wl[mask1], spec1flux[mask1]

		
		else:
			if starType1=="DA" or starType1=="DC" or starType1=="DB":  model_wl1, model_spectrum_star1 = return_model_spectrum_DA(wl_all1_N, ref_wl, cut_wl_min, cut_wl_max, Grav1_N, flux1_N, Teff1_N, T1, logg1)
			elif starType1=="DBA": model_wl1, model_spectrum_star1 = return_model_spectrum_DBA(wl_all1_N, ref_wl, cut_wl_min, cut_wl_max, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, T1, logg1, HoverHe1)
			elif starType1.startswith("sd"): model_wl1, model_spectrum_star1 = return_model_spectrum_subdwarf(wl_all1_N, ref_wl, cut_wl_min, cut_wl_max, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, T1, logg1, HoverHe1)
			
		
		
		
		# go through all files one by one and do the same as above, but recalculating the interpolated spectrum, which makes things slower
		for (ii, in_fi, sh_rv, modha, cut_lim, norm_lim, aHJD, sigclipspec) in zip(index[mask_ref_wl], input_files[mask_ref_wl], share_rv[mask_ref_wl], modelHa[mask_ref_wl], cut_Ha_all[mask_ref_wl], normaliseHa_all[mask_ref_wl], HJD_values[mask_ref_wl], sigma_clip[mask_ref_wl]):  # the mcmc RVs are gone through 1 by one, uses the ii counter.
			modelHa_min, modelHa_max = modha
			cut_limits_min, cut_limits_max = cut_lim
			norm_limits_min, norm_limits_max = norm_lim
			
			
			
			
			if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
				if sh_rv == -1:  mcmc_rv1 = RV[ii]
				else:            mcmc_rv1 = RV[sh_rv]
			else:
				phi = ((aHJD - forced_T0)%forced_P0)/forced_P0
				if forced_K1=="Fit":     
					if forced_Vgamma1=="Fit":    mcmc_rv1 = mcmc_Vgamma1 + mcmc_K1*np.sin(2*np.pi*phi)
					else:    mcmc_rv1 = forced_Vgamma1 + mcmc_K1*np.sin(2*np.pi*phi)
				else:                   
					if forced_Vgamma1=="Fit":    mcmc_rv1 = mcmc_Vgamma1 + forced_K1*np.sin(2*np.pi*phi)
					else: mcmc_rv1 = forced_Vgamma1 + forced_K1*np.sin(2*np.pi*phi)
				if forced_K1!="Fit" and forced_Vgamma1=="Fit":     mcmc_rv1 = mcmc_Vgamma1 + forced_K1*np.sin(2*np.pi*phi)
					
			
			
			
			normalised_wavelength = list_norm_wl_grids[ii];    normalised_flux = list_normalised_flux[ii];    normalised_err = list_normalised_err[ii];    inp_resolution = resolutions[ii]
			
			
			dlam1 = ref_wl*mcmc_rv1/speed_of_light


			
				
			# Smear to the resolution desired
			interparmodel = convolve(model_spectrum_star1, Gaussian1DKernel(stddev=0.5*(ref_wl/inp_resolution)/(model_wl1[10]-model_wl1[9])), boundary = 'extend')
			interparr = interp(normalised_wavelength, model_wl1+dlam1, interparmodel)
			
			
			
			mask = ((normalised_wavelength>ref_wl+cut_limits_min) & (normalised_wavelength<ref_wl+norm_limits_min)) | ((normalised_wavelength>ref_wl+norm_limits_max)   & (normalised_wavelength<ref_wl+cut_limits_max))
			
			
					
			try:
				m,c =  polyfit(normalised_wavelength[mask], interparr[mask], deg=1)
				interparr /= (m*normalised_wavelength + c)
			except: return -np.inf
			
			
			
			desired_range = (normalised_wavelength>ref_wl+modelHa_min)  &  (normalised_wavelength<ref_wl+modelHa_max)
			if sigclipspec!=-1:
				resid = normalised_flux-(interparr)
				resid_in_sig = resid/np.std(np.abs(resid[desired_range]))
				#plt.scatter(normalised_wavelength[desired_range], resid_in_sig[desired_range]);   plt.show()
				
				clip_mask = (np.abs(resid_in_sig[desired_range])>sigclipspec) & ((normalised_wavelength[desired_range] > ref_wl+5) | (normalised_wavelength[desired_range] < ref_wl-5))
				#plt.axhline(0, c='k');   plt.axhline(-1*sigclipspec, c='grey', ls='--');  plt.axhline(sigclipspec, c='grey', ls='--');  plt.scatter(normalised_wavelength[desired_range][clip_mask], resid_in_sig[desired_range][clip_mask],c='r');  plt.show();  plt.clf()
			else: clip_mask=normalised_flux[desired_range]!=normalised_flux[desired_range]
			
			
			off =  polyfit(normalised_wavelength[desired_range][~clip_mask], normalised_flux[desired_range][~clip_mask] - interparr[desired_range][~clip_mask], w=1/normalised_err[desired_range][~clip_mask], deg=0)[0]
			
			
			if False:
				plt.plot(normalised_wavelength[desired_range], normalised_flux[desired_range], c='r')
				plt.plot(normalised_wavelength[desired_range], interparr[desired_range], c='orange')
				plt.plot(normalised_wavelength[desired_range], interparr[desired_range] + off, c='g')
				plt.show()
			
			
			
			
			#plt.plot(normalised_wavelength[desired_range] - ref_wl, normalised_flux[desired_range], c='k');   plt.plot(normalised_wavelength[desired_range] - ref_wl, interparr[desired_range], c='orange')
			#plt.title(str(len(desired_range[desired_range==True])) + "   " + str(ref_wl));   plt.show();  plt.clf()
			if False:
				plt.plot(normalised_wavelength[desired_range], normalised_flux[desired_range],c='b')
				plt.plot(normalised_wavelength[desired_range][~clip_mask], normalised_flux[desired_range][~clip_mask],c='k')
				plt.scatter(normalised_wavelength[desired_range][clip_mask], normalised_flux[desired_range][clip_mask],c='r')
				plt.show()

			## Now on to calculating chisq
			chisq_indiv = -0.5*np.sum((np.square(normalised_flux[desired_range][~clip_mask]-(off+interparr[desired_range][~clip_mask])))/np.square(normalised_err[desired_range][~clip_mask]))
			chisq_spec += chisq_indiv
	try: rchisq_spec=chisq_spec/(len(share_rv[share_rv==-1]) - 1)
	except: rchisq_spec=chisq_spec  #  when 1 spectrum is used, enters here
	
	
	if fit_phot_SED:
		#if np.isnan(chisq_spec): plt.plot(normalised_wavelength, interparr1); plt.plot(normalised_wavelength, interparr2); plt.title(str(T1) + "  "+ str(logg1)+"  "+ str(T2)+"  "+ str(logg2));  plt.show()
		if np.isnan(chisq_spec) == True or np.isnan(chisq_phot) == True:    return -np.inf
		else:  print(chisq_spec,chisq_phot);  return chisq_spec + chisq_phot # (1/np.abs(rchisq_spec)) * chisq_spec + (1/np.abs(rchisq_phot)) * chisq_phot   # this is not right, but it's something. i really need to caclulate what chisq would be corresponds to the results when rchisq=1
	else:
		if np.isnan(chisq_spec) == True :    return -np.inf
		else:   return chisq_spec
		
		


def lnprob(theta, arguments):
	lp = lnprior(theta, arguments)
	if not np.isfinite(lp):
		return -np.inf
	return lp + lnlike(theta, arguments)






#Test that next step falls within boundaries allowed by the priors
def lnprior_gauss_lorentz(theta, arguments):
	input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, normalised_wavelength, normalised_flux, normalised_err,  inp_resolution, used_RV_boundaries, HJD_values, sigma_clip, the_unique_wavelengths, spec1wl, spec1flux, inputScaling = arguments
	if "K1" in p0labels:            args = np.argwhere(p0labels=="K1")[0][0];          mcmc_K1 = theta[args]
	if "Vg1" in p0labels:           args = np.argwhere(p0labels=="Vg1")[0][0];         mcmc_Vgamma1 = theta[args]
	if "RV1_0" in p0labels:         num_start_RVs = np.argwhere(p0labels=="RV1_0")[0][0]
	elif "RV1_1" in p0labels:       num_start_RVs = np.argwhere(p0labels=="RV1_1")[0][0]
	elif "RV1_2" in p0labels:       num_start_RVs = np.argwhere(p0labels=="RV1_2")[0][0]
	elif "RV1_3" in p0labels:       num_start_RVs = np.argwhere(p0labels=="RV1_3")[0][0]
	elif "RV1_4" in p0labels:       num_start_RVs = np.argwhere(p0labels=="RV1_4")[0][0]
	
	args = np.argwhere(p0labels=="A1_1")[0][0];   gaussLorentz1_A1 = theta[args]
	args = np.argwhere(p0labels=="A1_2")[0][0];   gaussLorentz1_A2 = theta[args]
	args = np.argwhere(p0labels=="sig1_1")[0][0];   gaussLorentz1_sigma1 = theta[args]
	args = np.argwhere(p0labels=="sig1_2")[0][0];   gaussLorentz1_sigma2 = theta[args]
	
	
	if not (A1_1_boundaries[0] < gaussLorentz1_A1 < A1_1_boundaries[1] or A1_2_boundaries[0] < gaussLorentz1_A2 < A1_2_boundaries[1]):
		return -np.inf
	
	if not (sigma1_1_boundaries[0] < gaussLorentz1_sigma1 < sigma1_1_boundaries[1] or sigma1_2_boundaries[0] < gaussLorentz1_sigma2 < sigma1_2_boundaries[1]):
		return -np.inf
	
	
	
	if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
		RV=theta[num_start_RVs:]
		
		
		
	passed = True
	
	if passed==False:  return -np.inf
	
	
	if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
		for anRV, anRVbound in zip(RV,used_RV_boundaries):  # go through all individual RVs and make sure they are in boundaries
			if not anRVbound[0] < anRV < anRVbound[1]:        return -np.inf
	else:
		if forced_K1=="Fit":
			if forced_Vgamma1=="Fit":
				if not (p0K1[0] < mcmc_K1 < p0K1[1] and p0Vgamma1[0] < mcmc_Vgamma1 < p0Vgamma1[1]):    return -np.inf
			else:
				if not p0K1[0] < mcmc_K1 < p0K1[1]:       return -np.inf
		
		if forced_K1!="Fit" and forced_Vgamma1=="Fit":
			if not p0Vgamma1[0] < mcmc_Vgamma1 < p0Vgamma1[1]:    return -np.inf
	
	return 0.0
	
	

@njit
def gauss_2x_mcmc(x, a1, sigma1, a2, sigma2):#  no offset 'b' because I fit offset the observations to best match the profile anyway. Ignoring b1 b2 also means that b1=b2=0.  No need to include a mean here
	return a1*np.exp(-0.5 * np.square(x/sigma1))   +   a2*np.exp(-0.5 * np.square(x/sigma2))

@njit
def lorentz_2x_mcmc(x, a1, sigma1, a2, sigma2):#  no offset 'b' because I fit offset the observations to best match the profile anyway. Ignoring b1 b2 also means that b1=b2=0.  No need to include a mean here
	return (a1/np.pi)  *  (sigma1/(x**2 + sigma1**2))  +  (a2/np.pi)  *  (sigma2/(x**2 + sigma2**2))

@njit
def gauss_lorentz_mcmc(x, a1, sigma1, a2, sigma2):#  no offset 'b' because I fit offset the observations to best match the profile anyway. Ignoring b1 b2 also means that b1=b2=0.  No need to include a mean here
	return a1*np.exp(-0.5 * np.square(x/sigma1))  +  (a2/np.pi)  *  (sigma2/(x**2 + sigma2**2))
	




def lnlike_gauss_lorentz(theta, arguments):
	input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, normalised_wavelength, normalised_flux, normalised_err,  inp_resolution, used_RV_boundaries, HJD_values, sigma_clip, the_unique_wavelengths, spec1wl, spec1flux, inputScaling = arguments
	if type(input_files)==np.str_:
		input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, normalised_wavelength, normalised_flux, normalised_err,  inp_resolution, used_RV_boundaries, HJD_values, sigma_clip = [input_files], [share_rv], [reference_wl], [cut_Ha_all], [normaliseHa_all], [normalised_wavelength], [normalised_flux], [normalised_err],  [inp_resolution], [used_RV_boundaries], [HJD_values], [sigma_clip]
	
	args = np.argwhere(p0labels=="A1_1")[0][0];   gaussLorentz1_A1 = theta[args]
	args = np.argwhere(p0labels=="A1_2")[0][0];   gaussLorentz1_A2 = theta[args]
	args = np.argwhere(p0labels=="sig1_1")[0][0];   gaussLorentz1_sigma1 = theta[args]
	args = np.argwhere(p0labels=="sig1_2")[0][0];   gaussLorentz1_sigma2 = theta[args]
	

	
	if "K1" in p0labels:            args = np.argwhere(p0labels=="K1")[0][0];          mcmc_K1 = theta[args]
	else:                           mcmc_K1 = forced_K1
	
	if "Vg1" in p0labels:           args = np.argwhere(p0labels=="Vg1")[0][0];         mcmc_Vgamma1 = theta[args]
	else:                           mcmc_Vgamma1 = forced_Vgamma1
	
	if "RV1_0" in p0labels:         num_start_RVs = np.argwhere(p0labels=="RV1_0")[0][0]
	
	

	
	if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
		RV=theta[num_start_RVs:]
	
	
	
	
	chisq_spec=0
	
	index = np.arange(0,len(input_files),1)
	
	the_unique_wavelengths=[np.amax(the_unique_wavelengths)]
	for ref_wl in the_unique_wavelengths:
		mask_ref_wl = reference_wl==ref_wl
		ref_wl_cut_lims_min, ref_wl_cut_lims_max = cut_Ha_all[mask_ref_wl].T
		cut_wl_min, cut_wl_max = np.amin(ref_wl_cut_lims_min), np.amax(ref_wl_cut_lims_max)
		model_wl1 = np.linspace(cut_wl_min, cut_wl_max, int(cut_wl_max-cut_wl_min)*10)
		if starType1=="GG":     model_spectrum_star1 = gauss_2x_mcmc(model_wl1, gaussLorentz1_A1, gaussLorentz1_sigma1, gaussLorentz1_A2, gaussLorentz1_sigma2)
		elif starType1=="LL":   model_spectrum_star1 = lorentz_2x_mcmc(model_wl1, gaussLorentz1_A1, gaussLorentz1_sigma1, gaussLorentz1_A2, gaussLorentz1_sigma2)
		elif starType1=="GL":   model_spectrum_star1 = gauss_lorentz_mcmc(model_wl1, gaussLorentz1_A1, gaussLorentz1_sigma1, gaussLorentz1_A2, gaussLorentz1_sigma2)
		
		model_wl1 = model_wl1 + ref_wl
		
		#plt.plot(model_wl1, model_spectrum_star1,c='k');   plt.plot(model_wl2, model_spectrum_star2,c='r');   plt.show()
		
		
		
		
		# go through all files one by one and do the same as above, but recalculating the interpolated spectrum, which makes things slower
		for (ii, in_fi, sh_rv, modha, cut_lim, norm_lim, aHJD, sigclipspec) in zip(index[mask_ref_wl], input_files[mask_ref_wl], share_rv[mask_ref_wl], modelHa[mask_ref_wl], cut_Ha_all[mask_ref_wl], normaliseHa_all[mask_ref_wl], HJD_values[mask_ref_wl], sigma_clip[mask_ref_wl]):  # the mcmc RVs are gone through 1 by one, uses the ii counter.
			modelHa_min, modelHa_max = modha
			cut_limits_min, cut_limits_max = cut_lim
			norm_limits_min, norm_limits_max = norm_lim
			
			
			
			
			if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
				if sh_rv == -1:  mcmc_rv1 = RV[ii]
				else:            mcmc_rv1 = RV[sh_rv]
			else:
				phi = ((aHJD - forced_T0)%forced_P0)/forced_P0
				if forced_K1=="Fit":     
					if forced_Vgamma1=="Fit":    mcmc_rv1 = mcmc_Vgamma1 + mcmc_K1*np.sin(2*np.pi*phi)
					else:    mcmc_rv1 = forced_Vgamma1 + mcmc_K1*np.sin(2*np.pi*phi)
				else:                   
					if forced_Vgamma1=="Fit":    mcmc_rv1 = mcmc_Vgamma1 + forced_K1*np.sin(2*np.pi*phi)
					else: mcmc_rv1 = forced_Vgamma1 + forced_K1*np.sin(2*np.pi*phi)
				if forced_K1!="Fit" and forced_Vgamma1=="Fit":     mcmc_rv1 = mcmc_Vgamma1 + forced_K1*np.sin(2*np.pi*phi)
					
			
			
			
			normalised_wavelength = list_norm_wl_grids[ii];    normalised_flux = list_normalised_flux[ii];    normalised_err = list_normalised_err[ii];    inp_resolution = resolutions[ii]
			
			
			dlam1 = ref_wl*mcmc_rv1/speed_of_light
			
			
			# interpolate the model onto the observation's wavelength grid
			interparr = interp(model_wl1, model_wl1+dlam1, model_spectrum_star1)   +   1
			
			
			# Smear to the resolution desired
			
			resstd=0.5*(ref_wl/inp_resolution)/(model_wl1[10]-model_wl1[9])
			
			interparr = convolve(interparr, Gaussian1DKernel(stddev=resstd), boundary = 'extend')
			
			
			
			interparr = interp(normalised_wavelength, model_wl1, interparr)
				
			#plt.plot(normalised_wavelength, normalised_flux);   plt.plot(normalised_wavelength, interparr);   plt.show()
			
			
			mask = ((normalised_wavelength>ref_wl+cut_limits_min) & (normalised_wavelength<ref_wl+norm_limits_min)) | ((normalised_wavelength>ref_wl+norm_limits_max)   & (normalised_wavelength<ref_wl+cut_limits_max))
			
			
			try:
				m,c =  polyfit(normalised_wavelength[mask], interparr[mask], deg=1)
				interparr /= (m*normalised_wavelength + c)
			except: return -np.inf
			
			
			if False:
				mask_rv=True
				lam_rv_remove = ref_wl*50/speed_of_light
				if mask_rv:   aaamask = (normalised_wavelength > ref_wl + dlam2 + lam_rv_remove) | (normalised_wavelength < ref_wl + dlam2 - lam_rv_remove)
				else:         aaamask = normalised_wavelength==normalised_wavelength
			else:   aaamask = normalised_wavelength==normalised_wavelength
			
			
			
			desired_range = (normalised_wavelength>ref_wl+modelHa_min)  &  (normalised_wavelength<ref_wl+modelHa_max) & aaamask
			if sigclipspec!=-1:
				resid = normalised_flux-(interparr)
				resid_in_sig = resid/np.std(np.abs(resid[desired_range]))
				
				clip_mask = (np.abs(resid_in_sig[desired_range])>sigclipspec) & ((normalised_wavelength[desired_range] > ref_wl+5) | (normalised_wavelength[desired_range] < ref_wl-5))
			else: clip_mask=normalised_flux[desired_range]!=normalised_flux[desired_range]
			
			
			
			off = polyfit(normalised_wavelength[desired_range][~clip_mask], normalised_flux[desired_range][~clip_mask] - interparr[desired_range][~clip_mask], w=1/normalised_err[desired_range][~clip_mask], deg=0)[0]
			

			## Now on to calculating chisq
			chisq_indiv = -0.5*np.sum((np.square(normalised_flux[desired_range & aaamask][~clip_mask]-(off+interparr[desired_range & aaamask][~clip_mask])))/np.square(normalised_err[desired_range & aaamask][~clip_mask]))
			chisq_spec += chisq_indiv
			
			
	
	
	if np.isnan(chisq_spec) == True :    return -np.inf
	else:   return chisq_spec
		
		


def lnprob_gauss_lorentz(theta, arguments):
	lp = lnprior_gauss_lorentz(theta, arguments)
	if not np.isfinite(lp):
		return -np.inf
	return lp + lnlike_gauss_lorentz(theta, arguments)















if sys.argv[1]=="RV" or sys.argv[1]=="RV_gauss":
	desired_wl=6562.81
	# load in atmospheric result
	found_out_scaling=False
	with open("out/result.out") as resultfile:
		result_lines = resultfile.readlines()
		for iii, lll in enumerate(result_lines):
			if starType1.startswith("D") or starType1.startswith("sd"):
				if "Teff1:" in lll:
					T1_med=float(result_lines[iii+1].split("\t")[0])
				elif "Logg1:" in lll:
					logg1_med=float(result_lines[iii+1].split("\t")[0])
				elif "H/He1:" in lll:
					HoverHe1_med=float(result_lines[iii+1].split("\t")[0])
				elif "Scaling:" in lll:
					Scaling=float(result_lines[iii+1].split("\t")[0])
					found_out_scaling=True
			else:
				if "A1_1:" in lll:
					A1_1_med=float(result_lines[iii+1].split("\t")[0])
				elif "sig1_1:" in lll:
					sig1_1_med=float(result_lines[iii+1].split("\t")[0])
				elif "A1_2:" in lll:
					A1_2_med=float(result_lines[iii+1].split("\t")[0])
				elif "sig1_2:" in lll:
					sig1_2_med=float(result_lines[iii+1].split("\t")[0])
					
		resultfile.close()
	
	
	if starType1.startswith("D") or starType1.startswith("sd"):  #  missing starType1.startswith("sd") bit
		try: T1_med
		except: T1_med = forced_teff1
		try: logg1_med
		except: logg1_med = forced_logg1
		
		if found_out_scaling==False:
			if forced_Scaling==False:	Scaling_med = Scaling_med
			elif forced_Scaling=="WD" and starType1=="DA":
				R1 = get_MTR(T1_med, logg=logg1_med, return_R_from_T_logg=True, loaded_Istrate=loaded_Istrate, loaded_CO=loaded_CO, loaded_Althaus=loaded_Althaus)
			elif forced_Scaling=="WD" and starType1=="DBA"  or starType1=="DB" or starType1=="DC":
				R1 = get_MTR_DB(T1_med, logg=logg1_med)
			else: Scaling_med=forced_Scaling
	
	
		pier_or_antoine="mixed"#"pier3Dphot_antoine1Dspec"
		for tempcnt, (temperature_star, logg_star, starType) in enumerate(zip([T1_med],[logg1_med], [starType1])):
			minval=10000;   maxval=-10000
			if starType=="DA":
				if not pier_or_antoine=="pier3Dphot_antoine1Dspec":
					if tempcnt==0:     Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DAgrids(temperature_star, logg_star)
				else:
					if tempcnt==0:     Grav1_N, wl_all1_N, flux1_N, Teff1_N, photGrav1_N, photwl_all1_N, photflux1_N, photTeff1_N = return_DAgrids(temperature_star, logg_star)
				

			elif starType=="DBA" or starType.startswith("sd"):
				if tempcnt==0:    HoverHestar=HoverHe1_med;    Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N = return_DBAgrids(temperature_star, logg_star, HoverHestar)
			
			
			elif starType=="DB":
				if tempcnt==0:    Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DBgrids(temperature_star, logg_star)
			
			elif starType=="DC":
				if tempcnt==0:    Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DCgrids(temperature_star, logg_star)
		
	else:   Scaling_med=1

	
	
	#raise ValueError(np.amin(wl_all1_N), np.amax(wl_all1_N), 6562.81)
	cut_limits_min, cut_limits_max = np.amin(cut_Ha_all.T), np.amax(cut_Ha_all.T)
	
	if starType1.startswith("D")  or starType1.startswith("sd"):
		if starType1=="DA" or starType1=="DB" or starType1=="DC": model_wl1, model_spectrum_star1 = return_model_spectrum_DA(wl_all1_N, desired_wl, cut_limits_min-20, cut_limits_max+20, Grav1_N, flux1_N, Teff1_N, T1_med, logg1_med)
		elif starType1=="DBA": model_wl1, model_spectrum_star1 = return_model_spectrum_DBA(wl_all1_N, desired_wl, cut_limits_min-20, cut_limits_max+20, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, T1_med, logg1_med, HoverHe1_med)
		elif starType1.startswith("sd"):  model_wl1, model_spectrum_star1 = return_model_spectrum_subdwarf(wl_all1_N, desired_wl, cut_limits_min-20, cut_limits_max+20, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, T1_med, logg1_med, HoverHe1_med)
	
	else:
		model_wl1 = np.linspace(cut_limits_min, cut_limits_max, int(cut_limits_max-cut_limits_min)*10) #+ ref_wl
		if starType1=="GG":     model_spectrum_star1 = gauss_2x_mcmc(model_wl1, A1_1_med, sig1_1_med, A1_2_med, sig1_2_med) + 1
		elif starType1=="LL":   model_spectrum_star1 = lorentz_2x_mcmc(model_wl1, A1_1_med, sig1_1_med, A1_2_med, sig1_2_med) + 1
		elif starType1=="GL":   model_spectrum_star1 = gauss_lorentz_mcmc(model_wl1, A1_1_med, sig1_1_med, A1_2_med, sig1_2_med) + 1
		model_wl1 = model_wl1 + np.amax(reference_wl)
			
	
	
	from lmfit.models import GaussianModel
	from lmfit import Model
	from scipy.optimize import least_squares
	
	
	if sys.argv[1]=="RV_gauss":
		desired_refwl = np.amax(reference_wl)
		mask_ref_wl = reference_wl==desired_refwl
		
		
		#########################################################################################################################################################################
		##### ADD EXTRA GAUSSIAN IF DESIRED
		#########################################################################################################################################################################
		try: file_ignore = np.asarray([config_info["file_ignore"]])[0]
		except: file_ignore = []
		
		
		
		
		
		ndim=2
		p0A1 = [-0.1, 0]
		p0std_dev1 = [0, 1]
		
		pos_min = np.array([p0A1[0], p0std_dev1[0]])
		pos_max = np.array([p0A1[1], p0std_dev1[1]])
		p0range = np.array([p0A1, p0std_dev1])
		p0labels = np.array(["A1", "std1"])
		for ii, in_fi in enumerate(input_files[mask_ref_wl]):
			if not in_fi in file_ignore:
				rvbound1 = RV_boundaries1[ii]
				
				ndim+=1
				pos_min = np.append(pos_min, rvbound1[0]);   pos_max = np.append(pos_max, rvbound1[1])
				p0range = np.concatenate((p0range, np.array([rvbound1])))
				p0labels = np.append(p0labels,"RV1_"+str(ii))
				
		
		nwalkers=len(p0labels)+(len(p0labels)-2)*2 + 2
		burnin=50
		nsteps=50
		
		def lnprior_gauss(theta, arguments):
			A1, std_dev1 = theta[0], theta[1]
			RVs = theta[2:]
			RVs_star1 = RVs
			
			
			if p0A1[0]<A1<p0A1[1] and p0std_dev1[0]<std_dev1<p0std_dev1[1]:
				return 0.0
			else:   return np.inf
		
		def agauss(x, A, cen, std_dev):
			return A*np.exp(-np.square(x-cen)/(2*std_dev**2))

		def lnlike_gauss(theta, arguments):
			A1, std_dev1 = theta[0], theta[1]
			RVs = theta[2:]
			RVs_star1 = RVs
			
			input_files, reference_wl, cut_Ha_all, normaliseHa_all, list_norm_wl_grids, list_normalised_flux, list_normalised_err,  resolutions, used_RV_boundaries, wl1, spec1, desired_refwl = arguments
			
			
			gauss1 = agauss(wl1-desired_refwl, A1, 0, std_dev1)
			
			mask_ref_wl = reference_wl==desired_refwl
			chisq=0
			for ii, (in_fi, modha, cut_lim, norm_lim, RV1) in enumerate(zip(input_files[mask_ref_wl], modelHa[mask_ref_wl], cut_Ha_all[mask_ref_wl], normaliseHa_all[mask_ref_wl], RVs_star1)):
				if not in_fi in file_ignore:
					wl_RV1 = RV1/speed_of_light * desired_refwl
					normalised_wavelength, normalised_flux, normalised_err, inp_resolution = list_norm_wl_grids[ii], list_normalised_flux[ii], list_normalised_err[ii], resolutions[ii]
					rvbound1 = RV_boundaries1[ii]
					modelHa_min, modelHa_max = modha
					cut_limits_min, cut_limits_max = cut_lim
					norm_limits_min, norm_limits_max = norm_lim
					
					
					resAA = desired_refwl/inp_resolution
					
					
					smear_model_spectrum_star1 = convolve(spec1, Gaussian1DKernel(stddev=0.5*resAA/(wl1[10]-wl1[9])), boundary = 'extend')
					
					smear_model_spectrum_star1 = interp(wl1, wl1+wl_RV1, smear_model_spectrum_star1)
					
					
					interparr = smear_model_spectrum_star1
					
					
					mask = ((wl1>desired_refwl+cut_limits_min) & (wl1<desired_refwl+norm_limits_min)) | ((wl1>desired_refwl+norm_limits_max)   & (wl1<desired_refwl+cut_limits_max))
					
					m,c =  polyfit(wl1[mask], interparr[mask], deg=1)
					interparr /= (m*wl1 + c)
					interparr = interp(normalised_wavelength, wl1, interparr)
					
					
					smear_gauss1 = convolve(gauss1, Gaussian1DKernel(stddev=0.5*resAA/(wl1[10]-wl1[9])), boundary = 'extend')
					
					smear_gauss1 = interp(normalised_wavelength, wl1+wl_RV1, smear_gauss1)
					
					
					#plt.plot(normalised_wavelength, interparr+smear_gauss1+smear_gauss2,c='g')
					#plt.plot(normalised_wavelength, interparr,c='r')
					#plt.plot(normalised_wavelength, normalised_flux, c='k')
					##plt.plot(normalised_wavelength, smear_gauss1+1, c='orange');   plt.plot(normalised_wavelength, smear_gauss2+1, c='b')
					##plt.plot(normalised_wavelength, smear_gauss1+smear_gauss2+1, c='purple')
					#plt.title(str(A1) + "   " + str(A2) + "   " + str(std_dev1) + "   " + str(std_dev2));   plt.show()
					
					
					desired_range = (normalised_wavelength>desired_wl-20)  &  (normalised_wavelength<desired_wl+20)
					chisq += -0.5*np.sum((np.square(normalised_flux[desired_range]-interparr[desired_range]))/np.square(normalised_err[desired_range]))
					
					
			# raise ValueError(chisq, in_fi, "@@@@", file_ignore, input_files, ii)
			return chisq
			
		
		
		def lnprob_gauss(theta, arguments):
			lp = lnprior_gauss(theta, arguments)
			if not np.isfinite(lp):
				return -np.inf
			return lp + lnlike_gauss(theta, arguments)
		
					
		
				
		#raise ValueError(nwalkers, p0labels)
		pool = MPIPool()
		if not pool.is_master():
			pool.wait()
			sys.exit(0)
		print("Starting MCMC")
		


		psize = pos_max - pos_min
		p0 = [pos_min + psize*np.random.rand(ndim) for i in range(nwalkers)]
		
		
		
		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_gauss, pool=pool, args=[[input_files, reference_wl, cut_Ha_all, normaliseHa_all, list_norm_wl_grids, list_normalised_flux, list_normalised_err,  resolutions, used_RV_boundaries, model_wl1, model_spectrum_star1, desired_refwl]])
		pos, prob, state = sampler.run_mcmc(p0,burnin,progress=True)
		samples = sampler.flatchain
		np.savetxt("RVfits/Gauss_MCMC_samples_burnin.dat",samples)
		np.savetxt("RVfits/Gauss_MCMC_burninpos.dat",pos)
		print("Finished burn-in")



		#Make a corner plot which shows values
		samples=np.loadtxt("RVfits/Gauss_MCMC_samples_burnin.dat")

		sampler.reset()
		pos=np.loadtxt("RVfits/Gauss_MCMC_burninpos.dat")
		

		pos, prob, state = sampler.run_mcmc(pos,nsteps,progress=True)
		print("MCMC completed")
		pool.close()

		samples = sampler.flatchain
		np.savetxt("RVfits/Gauss_MCMC_samples.dat",samples)

		#Make a corner plot which shows values
		samples=np.loadtxt("RVfits/Gauss_MCMC_samples.dat")
		#raise ValueError(len(p0range), len(samples))
		print(samples)
		
		A1_result = samples[::,0];    std_dev1_result = samples[::,1]
		A1_med = np.percentile(A1_result, 50, axis=0);  A1_min = np.percentile(A1_result, 16, axis=0);   A1_max = np.percentile(A1_result, 84, axis=0)
		std_dev1_med = np.percentile(std_dev1_result, 50, axis=0);  std_dev1_min = np.percentile(std_dev1_result, 16, axis=0);   std_dev1_max = np.percentile(std_dev1_result, 84, axis=0)
		
		
		#raise ValueError(A1_med, std_dev1_med, A2_med, std_dev2_med)
		#if plot_corner==True:
		#	fig2 = corner.corner(samples, labels=p0labels, range=p0range, bins=30,smooth=True, quantiles=[0.16, 0.5, 0.84], show_titles=True, labels_args={"fontsize": 40}, title_fmt = '.3f')
		#	fig2.savefig("RVfits/Gauss_corner.pdf")
		#	del fig2
		#plt.close();  plt.clf()
		
		
		
		
		
		
		
		
		
		
	#########################################################################################################################################################################################
	##### RUN RVS FOR A SOLUTION
	#########################################################################################################################################################################################

	#@jit
	def WDmodel(x, RV1, wl1, spec1, gauss1=np.array([1,1])):
		wl_RV1 = RV1/speed_of_light * desired_refwl
		
		
		# interpolate the model onto the observation's wavelength grid
		interparr = interp(x, wl1+wl_RV1, spec1)
		
		
		mask = ((x>desired_refwl+cut_limits_min) & (x<desired_refwl+norm_limits_min)) | ((x>desired_refwl+norm_limits_max)   & (x<desired_refwl+cut_limits_max))
		
		m,c =  polyfit(x[mask], interparr[mask], deg=1)
		interparr /= (m*x + c)
		
		if gauss1[0]!=1:
			interparr += interp(x,wl1+wl_RV1,gauss1)
		
		return interparr

	
	def calc_model(params, x, wl1, spec1,gauss1):
		'''Uses a cosine to model the pulsations'''
		RV1, extra_norm, coff = params
		#plt.scatter(x, (WDmodel(x, RV1, RV2, wl1, spec1, wl2, spec2,gauss1,gauss2) + offset)/(x*extra_norm + coff))
		#plt.title(str(extra_norm) + "  " + str(coff))	
		#plt.show()
		return (WDmodel(x, RV1, wl1, spec1,gauss1))/(x*extra_norm + coff)


	def residual(params, xxx, yyy, yyyerr,wl1,spec1,gauss1):
		'''Calculates residuals for a given model'''
		model = calc_model(params, xxx, wl1, spec1,gauss1)

		#mask_narrow = (xxx > desired_refwl-20) & (xxx < desired_refwl+20)
		
		#return (yyy[mask_narrow]-model[mask_narrow]) / yyyerr[mask_narrow]
		return (yyy-model) / yyyerr


	def fit_bootstrap(params, xxx, yyy, yyyerr, bounds, num_its, wl1, spec1, gauss1=np.array([1,1])):
		'''Samples the data a thousand times, fits each sample using least-squares,
		then reports median and standard deviation of the result'''
		RV1_guess, extra_norm_guess, coff_guess = params
		RV1, offset, extra_norm, coff = [], [], [], []
		data = np.array((xxx, yyy, yyyerr)).T
		for i in range(num_its):
			# Select 100% of our sample
			sample = data[np.random.choice(len(data),size=int(1.0*len(data)), replace=True)]
			# get the vectors back
			x = sample.T[0]
			y = sample.T[1]
			err = sample.T[2]
			lsm_fit = least_squares(residual, [RV1_guess, extra_norm_guess, coff_guess], method='trf',
		                        args=(x,y,err,wl1,spec1,gauss1), ftol=1E-10, bounds=((bounds[0], bounds[2], bounds[4]), (bounds[1], bounds[3], bounds[5])))
			
			
			
			RV1.append(lsm_fit.x[0])
			extra_norm.append(lsm_fit.x[1])
			coff.append(lsm_fit.x[2])
		
		
		return np.median(RV1), np.std(RV1), np.median(extra_norm), np.std(extra_norm), np.median(coff), np.std(coff)
	
	
	
	def calc_model2(params, x, wl1, spec1,gauss1):
		'''Uses a cosine to model the pulsations'''
		RV1, offset = params
		#plt.scatter(x, WDmodel(x, RV1, RV2, wl1, spec1, wl2, spec2,gauss1,gauss2) + offset)
		#plt.show()
		return WDmodel(x, RV1, wl1, spec1, gauss1) + offset


	def residual2(params, xxx, yyy, yyyerr,wl1,spec1,gauss1):
		'''Calculates residuals for a given model'''
		model = calc_model2(params, xxx, wl1, spec1, gauss1)

		mask_narrow = (xxx > desired_refwl-20) & (xxx < desired_refwl+20)
		
		return (yyy[mask_narrow]-model[mask_narrow]) / yyyerr[mask_narrow]


	def fit_bootstrap2(params, xxx, yyy, yyyerr, bounds, num_its, wl1, spec1, gauss1=np.array([1,1])):
		'''Samples the data a thousand times, fits each sample using least-squares,
		then reports median and standard deviation of the result'''
		RV1_guess, off_guess = params
		RV1, offset = [], []
		data = np.array((xxx, yyy, yyyerr)).T
		for i in range(num_its):
			sample = data[np.random.choice(len(data),size=int(1.0*len(data)), replace=True)]
			x = sample.T[0];    y = sample.T[1];    err = sample.T[2]
			lsm_fit = least_squares(residual2, [RV1_guess, off_guess], method='trf',
		                        args=(x,y,err,wl1,spec1,gauss1), ftol=1E-10, bounds=((bounds[0], bounds[2]), (bounds[1], bounds[3])))
			
			
			RV1.append(lsm_fit.x[0]);    offset.append(lsm_fit.x[1])
		
		
		return np.median(RV1), np.std(RV1), np.median(offset), np.std(offset)
	
	
	




	try: os.mkdir("RVfits")
	except: None
	
		
	rvfilename=[]
	desired_refwl = np.amax(reference_wl)
	mask_ref_wl = reference_wl==desired_refwl
	for ii, (in_fi, modha, cut_lim, norm_lim, aHJD, sigclipspec) in enumerate(zip(input_files[mask_ref_wl], modelHa[mask_ref_wl], cut_Ha_all[mask_ref_wl], normaliseHa_all[mask_ref_wl], HJD_values[mask_ref_wl], sigma_clip[mask_ref_wl])):
		try:
			normalised_wavelength, normalised_flux, normalised_err, inp_resolution = list_norm_wl_grids[ii], list_normalised_flux[ii], list_normalised_err[ii], resolutions[ii]
			rvbound1 = RV_boundaries1[ii]
			modelHa_min, modelHa_max = modha
			cut_limits_min, cut_limits_max = cut_lim
			norm_limits_min, norm_limits_max = norm_lim
			
			
			resAA = desired_refwl/inp_resolution
			
			
			smear_model_spectrum_star1 = convolve(model_spectrum_star1, Gaussian1DKernel(stddev=0.5*resAA/(model_wl1[10]-model_wl1[9])), boundary = 'extend')
			
			
			if sys.argv[1]=="RV_gauss":
				gauss1 = agauss(model_wl1, A1_med, desired_wl, std_dev1_med)
				
				smear_gauss1 = convolve(gauss1, Gaussian1DKernel(stddev=0.5*resAA/(model_wl1[10]-model_wl1[9])), boundary = 'extend')
			
			
			
			
			Low_SNR=True  # note that this improvement to normalisation only makes a noticeable difference when the SNR is low (because the normalisation of the data is not good). For speed, it can absolutely be ignored, but if you're analysing a small sample or on a case by case basis I recommend to leave it on. I noticed an accuracy of a couple 100ms-1 improved in SNR=25 in the wings for Halpha fits and a big improvement in precision for SNR=10-15 data.
			if Low_SNR:
				# First try to improve the accuracy of the normalisation using the model and the full wavelength range supplied
				RV1 = np.mean(rvbound1)   # this is the initial guess for the bootstrap
				if sys.argv[1]=="RV_gauss":
					RV1, RV1err, extra_norm, extra_norm_err, coff, coff_err = fit_bootstrap([RV1, 0, 1],  normalised_wavelength,  normalised_flux,  normalised_err, bounds=[rvbound1[0], rvbound1[1], -0.02,0.02, -20, 20], num_its=25, wl1=model_wl1, spec1=smear_model_spectrum_star1, gauss1=smear_gauss1)
				else:
					RV1, RV1err, extra_norm, extra_norm_err, coff, coff_err = fit_bootstrap([RV1, 0, 1],  normalised_wavelength,  normalised_flux,  normalised_err, bounds=[rvbound1[0], rvbound1[1], -0.02,0.02, -20, 20], num_its=25, wl1=model_wl1, spec1=smear_model_spectrum_star1)
				
				# clip out the bad points
				if sigclipspec!=-1:
					resid = normalised_flux - (WDmodel(normalised_wavelength, RV1, wl1=model_wl1, spec1=smear_model_spectrum_star1))/(normalised_wavelength*extra_norm + coff)
					resid_in_sig = resid/np.std(np.abs(resid))
					clip_mask = (np.abs(resid_in_sig)<sigclipspec) | ((normalised_wavelength < desired_refwl+5) & (normalised_wavelength > desired_refwl-5))
				else: clip_mask=normalised_flux==normalised_flux
				
				
				
				# Try to improve the accuracy of the normalisation again with the full wavelength range, this time with the bad points sigma clipped
				if sys.argv[1]=="RV_gauss":
					RV1, RV1err, extra_norm, extra_norm_err, coff, coff_err = fit_bootstrap([RV1, extra_norm, coff],  normalised_wavelength[clip_mask],  normalised_flux[clip_mask],  normalised_err[clip_mask], bounds=[rvbound1[0], rvbound1[1], -0.02,0.02, -20, 20], num_its=100, wl1=model_wl1, spec1=smear_model_spectrum_star1, gauss1=smear_gauss1)
				else:
					RV1, RV1err, extra_norm, extra_norm_err, coff, coff_err = fit_bootstrap([RV1, extra_norm, coff],  normalised_wavelength[clip_mask],  normalised_flux[clip_mask],  normalised_err[clip_mask], bounds=[rvbound1[0], rvbound1[1], -0.02,0.02, -20, 20], num_its=100, wl1=model_wl1, spec1=smear_model_spectrum_star1) 
				
				
				
				# Apply the tweaked normalisation to the data
				#plt.plot(normalised_wavelength, normalised_flux, c='k')
				normalised_flux*=(normalised_wavelength*extra_norm + coff)
				mask = ((normalised_wavelength>desired_refwl+cut_limits_min) & (normalised_wavelength<desired_refwl+norm_limits_min)) | ((normalised_wavelength>desired_refwl+norm_limits_max)   & (normalised_wavelength<desired_refwl+cut_limits_max))
				# and lastly make the spectrum normalised to 1 again using the desired region with 1 sigma clipping iteration to the normalisation region
				m_, c_ = polyfit(normalised_wavelength[mask]-desired_refwl, normalised_flux[mask], deg=1)
				resid=normalised_flux[mask] - (normalised_wavelength[mask]*m_ + c)
				try:
					mask2=np.abs(resid)<2.5
					m_, c_ = polyfit(normalised_wavelength[mask][mask2]-desired_refwl, normalised_flux[mask][mask2], deg=1)
					normalised_flux/=c_
				except: 
					normalised_flux/=c_
				
				
			
			if sigclipspec!=-1:
				resid = normalised_flux - (WDmodel(normalised_wavelength, RV1, wl1=model_wl1, spec1=smear_model_spectrum_star1))
				resid_in_sig = resid/np.std(np.abs(resid))
				clip_mask = (np.abs(resid_in_sig)<sigclipspec) | ((normalised_wavelength < desired_refwl+5) & (normalised_wavelength > desired_refwl-5))
			else: clip_mask=normalised_flux==normalised_flux
			
			
			
			mask_narrow = ((normalised_wavelength > desired_refwl-25) & (normalised_wavelength < desired_refwl+25)) | ((normalised_wavelength<desired_refwl+norm_limits_min+5) | (normalised_wavelength>desired_refwl+norm_limits_max-5))
			if sigclipspec!=-1:
				# fit for RV1 RV2 with a narrowed range around the line cores
				RV1, RV1err, off, off_err = fit_bootstrap2([RV1, 0],  normalised_wavelength[clip_mask],  normalised_flux[clip_mask],  normalised_err[clip_mask], bounds=[rvbound1[0], rvbound1[1], -0.2, 0.2], num_its=50, wl1=model_wl1, spec1=smear_model_spectrum_star1)
				
				
				# clip the narrowed area
				mask_narrow = ((normalised_wavelength > desired_refwl-25) & (normalised_wavelength < desired_refwl+25)) | ((normalised_wavelength<desired_refwl+norm_limits_min+5) | (normalised_wavelength>desired_refwl+norm_limits_max-5))
				resid = normalised_flux[mask_narrow] - (WDmodel(normalised_wavelength, RV1, wl1=model_wl1, spec1=smear_model_spectrum_star1) + off)[mask_narrow]
				resid_in_sig = resid/np.std(np.abs(resid))
				clip_mask2 = ((np.abs(resid_in_sig)<sigclipspec) | ((normalised_wavelength[mask_narrow] < desired_refwl+5) & (normalised_wavelength[mask_narrow] > desired_refwl-5))) | ((normalised_wavelength[mask_narrow]>desired_refwl+norm_limits_max-5) | (normalised_wavelength[mask_narrow]<desired_refwl+norm_limits_min+5))
			else: clip_mask2=normalised_flux[mask_narrow]==normalised_flux[mask_narrow]; off=0
			
			
			
			# fit for RV1 RV2 with a narrowed range around the line cores for a final result
			RV1, RV1err, off, off_err = fit_bootstrap2([RV1, off],  normalised_wavelength[mask_narrow][clip_mask2],  normalised_flux[mask_narrow][clip_mask2],  normalised_err[mask_narrow][clip_mask2], bounds=[rvbound1[0], rvbound1[1], -0.2, 0.2], num_its=200, wl1=model_wl1, spec1=smear_model_spectrum_star1)
			
			
			if True:
				print(in_fi)
				fig, (ax, ax2) = plt.subplots(1,2,gridspec_kw={'width_ratios': [2, 2]})#,sharey=True)
				plt.subplots_adjust(wspace=.0)
				#fine_grid = np.linspace(np.amin(normalised_wavelength[clip_mask]), np.amax(normalised_wavelength[clip_mask]), 100000)
				fine_grid = normalised_wavelength
				if sys.argv[1]=="RV_gauss":
					ax.plot(normalised_wavelength,  normalised_flux,c='k');    ax.plot(fine_grid, WDmodel(fine_grid, RV1, wl1=model_wl1, spec1=smear_model_spectrum_star1, gauss1=smear_gauss1) + off, c='orange'); ax.scatter(normalised_wavelength[~clip_mask],  normalised_flux[~clip_mask],c='pink'); ax.scatter(normalised_wavelength[mask_narrow][~clip_mask2],  normalised_flux[mask_narrow][~clip_mask2],c='r')
				else:
					ax.plot(normalised_wavelength,  normalised_flux,c='k');    ax.plot(fine_grid, WDmodel(fine_grid, RV1, wl1=model_wl1, spec1=smear_model_spectrum_star1) + off, c='orange'); ax.scatter(normalised_wavelength[~clip_mask],  normalised_flux[~clip_mask],c='pink'); ax.scatter(normalised_wavelength[mask_narrow][~clip_mask2],  normalised_flux[mask_narrow][~clip_mask2],c='r')
				mask_small_region = (normalised_wavelength>desired_refwl-17) & (normalised_wavelength<desired_refwl+17)
				if sys.argv[1]=="RV_gauss":
					ax2.plot(speed_of_light*(normalised_wavelength[mask_small_region]-desired_refwl)/desired_refwl,  normalised_flux[mask_small_region],c='k');    ax2.plot(speed_of_light*(fine_grid[mask_small_region]-desired_refwl)/desired_refwl, WDmodel(fine_grid, RV1, wl1=model_wl1, spec1=smear_model_spectrum_star1, gauss1=smear_gauss1)[mask_small_region] + off, c='orange')
				else:
					ax2.plot(speed_of_light*(normalised_wavelength[mask_small_region]-desired_refwl)/desired_refwl,  normalised_flux[mask_small_region],c='k');    ax2.plot(speed_of_light*(fine_grid[mask_small_region]-desired_refwl)/desired_refwl, WDmodel(fine_grid, RV1, wl1=model_wl1, spec1=smear_model_spectrum_star1)[mask_small_region]  + off, c='orange')
				ax2.set_yticks([])
				ax2.axvline(-100,c='grey', ls='dotted')
				#ax2.axvline(0,c='grey')
				ax2.axvline(100,c='grey', ls='dotted')
				ax.set_title("RV1=" + str(np.round(RV1,1)) + " " + str(np.round(RV1err,1)));  ax2.set_title("R=" + str(inp_resolution) + "  " + str(rvbound1)); plt.savefig("RVfits/"+str(in_fi)+".png", dpi=300) # plt.show()
				plt.clf()
				rvfilename.append([in_fi, aHJD, RV1, RV1err])
		except Exception as e: print(e); rvfilename.append([in_fi, aHJD, "ERROR", "ERROR"])
	np.savetxt("RVfits/RVfit_results.dat", rvfilename, fmt="%s", delimiter="\t",header="Observed velocities. For relativistic correction to source velocities (~0.14kms-1 for 200kms-1),  vsource  =  vobs/np.sqrt((1+vobs/c) / (1-vobs/c) )")
	
	
	

elif sys.argv[1]=="ATM":
	pool = MPIPool()
	if not pool.is_master():
		pool.wait()
		sys.exit(0)
	print("Starting MCMC")
	if isinstance(forced_teff1, list):   forced_teff1 = forced_teff1[0]
	if isinstance(forced_logg1, list):   forced_logg1 = forced_logg1[0]
	try:
		if isinstance(forced_HoverHe1, list):   forced_HoverHe1 = forced_HoverHe1[0]
	except: None
	
		
	if starType1.startswith("D")  or starType1.startswith("sd"):
		if forced_teff1!=0 and forced_logg1!=0:
			if starType1=="DA":
				Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DAgrids(forced_teff1, forced_logg1)
				model_wl1, model_spectrum_star1 = return_model_spectrum_DA(wl_all1_N, 5000, -5000, 5000, Grav1_N, flux1_N, Teff1_N, forced_teff1, forced_logg1)
			elif starType1.startswith("sd"):
				Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N = return_DBAgrids(forced_teff1[0], forced_logg1[0], forced_HoverHe1[0])
				model_wl1, model_spectrum_star1 = return_model_spectrum_subdwarf(wl_all1_N, 5000, -5000, 5000, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, T1_med, logg1_med, HoverHe1_med)
			
			
			sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool, args=[[input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, list_norm_wl_grids, list_normalised_flux, list_normalised_err,  resolutions, used_RV_boundaries, HJD_values, sigma_clip, np.unique(reference_wl), model_wl1, model_spectrum_star1]])
			
			
		else:
			sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool, args=[[input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, list_norm_wl_grids, list_normalised_flux, list_normalised_err,  resolutions, used_RV_boundaries, HJD_values, sigma_clip, np.unique(reference_wl), None, None]])
	
	else:
		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_gauss_lorentz, pool=pool, args=[[input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, list_norm_wl_grids, list_normalised_flux, list_normalised_err,  resolutions, used_RV_boundaries, HJD_values, sigma_clip, np.unique(reference_wl), None, None, None]])
	
	
	
	
	
	pos, prob, state = sampler.run_mcmc(p0,burnin,progress=True)
	samples = sampler.flatchain
	np.savetxt("out/MCMC_samples_burnin.dat",samples)
	np.savetxt("out/MCMC_burninpos.dat",pos)
	print("Finished burn-in")



	#Make a corner plot which shows values
	samples=np.loadtxt("out/MCMC_samples_burnin.dat")

	if plot_burnin:
		fig = corner.corner(samples, labels=p0labels, range=p0range, bins=30,smooth=True, quantiles=[0.16, 0.5, 0.84], show_titles=True, labels_args={"fontsize": 40}, title_fmt = '.3f')
		fig.savefig("out/corner_burnin.pdf")
		del fig
	plt.close();  plt.clf()

	sampler.reset()
	pos=np.loadtxt("out/MCMC_burninpos.dat")
	

	pos, prob, state = sampler.run_mcmc(pos,nsteps,progress=True)
	print("MCMC completed")
	pool.close()

	samples = sampler.flatchain
	np.savetxt("out/MCMC_samples.dat",samples)

	#Make a corner plot which shows values
	samples=np.loadtxt("out/MCMC_samples.dat")
	if plot_corner:
		fig2 = corner.corner(samples, labels=p0labels, range=p0range, bins=30,smooth=True, quantiles=[0.16, 0.5, 0.84], show_titles=True, labels_args={"fontsize": 40}, title_fmt = '.3f')
		fig2.savefig("out/corner.pdf")
		del fig2
	plt.close();  plt.clf()
	
	
	lines_to_write = []
	print(p0labels)
	if starType1.startswith("D") or starType1.startswith("sd"):
		if "T1" in p0labels:
			args = np.argwhere(p0labels=="T1")[0][0];          Teff1_result = samples[::,args]
			T1_med = np.percentile(Teff1_result, 50, axis=0);  T1_min = np.percentile(Teff1_result, 16, axis=0);   T1_max = np.percentile(Teff1_result, 84, axis=0)
			lines_to_write.append("Teff1:\n")
			lines_to_write.append(str(T1_med) + "\t" + str(np.percentile(Teff1_result, 16, axis=0) - T1_med) + "\t" + str(np.percentile(Teff1_result, 84, axis=0) - T1_med) + "\n")
		else:   T1_med = forced_teff1;   T1_min = 0;   T1_max=0
			

		if "logg1" in p0labels:
			args = np.argwhere(p0labels=="logg1")[0][0];       logg1_result = samples[::,args]
			logg1_med = np.percentile(logg1_result, 50, axis=0);   logg1_min = np.percentile(logg1_result, 16, axis=0);   logg1_max = np.percentile(logg1_result, 84, axis=0)
			lines_to_write.append("Logg1:\n")
			lines_to_write.append(str(logg1_med) + "\t" + str(np.percentile(logg1_result, 16, axis=0) - logg1_med) + "\t" + str(np.percentile(logg1_result, 84, axis=0) - logg1_med) + "\n")
		else:   logg1_med = forced_logg1;   logg1_min = 0;   logg1_max=0


		if "H/He1" in p0labels:
			args = np.argwhere(p0labels=="H/He1")[0][0];       HoverHe1_result = samples[::,args]
			HoverHe1_med = np.percentile(HoverHe1_result, 50, axis=0);   HoverHe1_min = np.percentile(HoverHe1_result, 16, axis=0);   HoverHe1_max = np.percentile(HoverHe1_result, 84, axis=0)
			lines_to_write.append("H/He1:\n")
			lines_to_write.append(str(HoverHe1_med) + "\t" + str(np.percentile(HoverHe1_result, 16, axis=0) - HoverHe1_med) + "\t" + str(np.percentile(HoverHe1_result, 84, axis=0) - HoverHe1_med) + "\n")
		else:   HoverHe1_med = forced_HoverHe1;   HoverHe1_min = 0;   HoverHe1_max=0
	
	else:
		if "A1_1" in p0labels:
			args = np.argwhere(p0labels=="A1_1")[0][0];          A1_1_result = samples[::,args]
			A1_1_med = np.percentile(A1_1_result, 50, axis=0);  A1_1_min = np.percentile(A1_1_result, 16, axis=0);   A1_1_max = np.percentile(A1_1_result, 84, axis=0)
			lines_to_write.append("A1_1:\n")
			lines_to_write.append(str(A1_1_med) + "\t" + str(np.percentile(A1_1_result, 16, axis=0) - A1_1_med) + "\t" + str(np.percentile(A1_1_result, 84, axis=0) - A1_1_med) + "\n")
		
		if "sig1_1" in p0labels:
			args = np.argwhere(p0labels=="sig1_1")[0][0];          sig1_1_result = samples[::,args]
			sig1_1_med = np.percentile(sig1_1_result, 50, axis=0);  sig1_1_min = np.percentile(sig1_1_result, 16, axis=0);   sig1_1_max = np.percentile(sig1_1_result, 84, axis=0)
			lines_to_write.append("sig1_1:\n")
			lines_to_write.append(str(sig1_1_med) + "\t" + str(np.percentile(sig1_1_result, 16, axis=0) - sig1_1_med) + "\t" + str(np.percentile(sig1_1_result, 84, axis=0) - sig1_1_med) + "\n")
		if "A1_2" in p0labels:
			args = np.argwhere(p0labels=="A1_2")[0][0];          A1_2_result = samples[::,args]
			A1_2_med = np.percentile(A1_2_result, 50, axis=0);  A1_2_min = np.percentile(A1_2_result, 16, axis=0);   A1_2_max = np.percentile(A1_2_result, 84, axis=0)
			lines_to_write.append("A1_2:\n")
			lines_to_write.append(str(A1_2_med) + "\t" + str(np.percentile(A1_2_result, 16, axis=0) - A1_2_med) + "\t" + str(np.percentile(A1_2_result, 84, axis=0) - A1_2_med) + "\n")
		if "sig1_2" in p0labels:
			args = np.argwhere(p0labels=="sig1_2")[0][0];          sig1_2_result = samples[::,args]
			sig1_2_med = np.percentile(sig1_2_result, 50, axis=0);  sig1_2_min = np.percentile(sig1_2_result, 16, axis=0);   sig1_2_max = np.percentile(sig1_2_result, 84, axis=0)
			lines_to_write.append("sig1_2:\n")
			lines_to_write.append(str(sig1_2_med) + "\t" + str(np.percentile(sig1_2_result, 16, axis=0) - sig1_2_med) + "\t" + str(np.percentile(sig1_2_result, 84, axis=0) - sig1_2_med) + "\n")



	if "K1" in p0labels:
		args = np.argwhere(p0labels=="K1")[0][0];          K1_result = samples[::,args]
		K1_med = np.percentile(K1_result, 50, axis=0);     K1_min = np.percentile(K1_result, 16, axis=0);   K1_max = np.percentile(K1_result, 84, axis=0)
		lines_to_write.append("K1:\n")
		lines_to_write.append(str(K1_med) + "\t" + str(np.percentile(K1_result, 16, axis=0) - K1_med) + "\t" + str(np.percentile(K1_result, 84, axis=0) - K1_med) + "\n")
	else:   K1_med = forced_K1;   K1_min = 0;   K1_max=0


	if "Vg1" in p0labels:
		args = np.argwhere(p0labels=="Vg1")[0][0];         Vgamma1_result = samples[::,args]
		Vgamma1_med = np.percentile(Vgamma1_result, 50, axis=0);   Vgamma1_min = np.percentile(Vgamma1_result, 16, axis=0);   Vgamma1_max = np.percentile(Vgamma1_result, 84, axis=0)
		lines_to_write.append("Vg1:\n")
		lines_to_write.append(str(Vgamma1_med) + "\t" + str(np.percentile(Vgamma1_result, 16, axis=0) - Vgamma1_med) + "\t" + str(np.percentile(Vgamma1_result, 84, axis=0) - Vgamma1_med) + "\n")
	else:   Vgamma1_med = forced_Vgamma1;   Vgamma1_min = 0;   Vgamma1_max=0


	if "Scaling" in p0labels:
		args = np.argwhere(p0labels=="Scaling")[0][0];     Scaling_result = samples[::,args]
		Scaling_med = np.percentile(Scaling_result, 50, axis=0);   Scaling_min = np.percentile(Scaling_result, 16, axis=0);   Scaling_max = np.percentile(Scaling_result, 84, axis=0)
		lines_to_write.append("Scaling:\n")
		lines_to_write.append(str(Scaling_med) + "\t" + str(np.percentile(Scaling_result, 16, axis=0) - Scaling_med) + "\t" + str(np.percentile(Scaling_result, 84, axis=0) - Scaling_med) + "\n")
	else:   Scaling_med = forced_Scaling;   Scaling_min = 0;   Scaling_max=0
	
	if "Parallax" in p0labels:
		args = np.argwhere(p0labels=="Parallax")[0][0];     parallax_result = samples[::,args]
		parallax_med = np.percentile(parallax_result, 50, axis=0);   parallax_min = np.percentile(parallax_result, 16, axis=0);   parallax_max = np.percentile(parallax_result, 84, axis=0)
		lines_to_write.append("Parallax:\n")
		lines_to_write.append(str(parallax_med) + "\t" + str(np.percentile(parallax_result, 16, axis=0) - parallax_med) + "\t" + str(np.percentile(parallax_result, 84, axis=0) - parallax_med) + "\n")



	if "R" in p0labels:
		args = np.argwhere(p0labels=="R")[0][0];     R_result = samples[::,args]
		R_med = np.percentile(R_result, 50, axis=0);   R_min = np.percentile(R_result, 16, axis=0);   R_max = np.percentile(R_result, 84, axis=0)
		lines_to_write.append("R:\n")
		lines_to_write.append(str(R_med) + "\t" + str(np.percentile(R_result, 16, axis=0) - R_med) + "\t" + str(np.percentile(R_result, 84, axis=0) - R_med) + "\n")





	# save just the median for individual epoch pngs with fit and observations
	if fit_phot_SED: # no starType1.startswith("sd") yet
		if forced_Scaling==False:	Scaling_med = Scaling_med
		elif forced_Scaling=="WD" and starType1=="DA": 
			R1 = get_MTR(T1_med, logg=logg1_med, return_R_from_T_logg=True, loaded_Istrate=loaded_Istrate, loaded_CO=loaded_CO, loaded_Althaus=loaded_Althaus)
		elif forced_Scaling=="WD" and starType1=="DBA" or starType1=="DB" or starType1=="DC":
			R1 = get_MTR_DB(T1_med, logg=logg1_med)
		else: Scaling_med=forced_Scaling




	allRV1s = []
	allRV1s_minerr, allRV1s_maxerr = [], []
	if "RV1_0" in p0labels:
		num_start_RVs = np.argwhere(p0labels=="RV1_0")[0][0]
		if "Scaling" in p0labels or "Parallax" in p0labels:   RVvals_result = samples[::,num_start_RVs:-1]
		else: RVvals_result = samples[::,num_start_RVs:]
		for cnt, item in enumerate(RVvals_result.T):
			allRV1s.append(np.percentile(item, 50, axis=0));   allRV1s_minerr.append(np.percentile(item, 16, axis=0));   allRV1s_maxerr.append(np.percentile(item, 84, axis=0))
			
	else:
		for aHJD in HJD_values:
			phi = ((aHJD - forced_T0)%forced_P0)/forced_P0
			allRV1s.append(Vgamma1_med + K1_med*np.sin(2*np.pi*phi));   allRV1s_minerr.append(0);   allRV1s_maxerr.append(0)
			


	lines_to_write.append("RV1:\n")
	for arvresult1, arvresult1_min, arvresult1_max, ahjd in zip(allRV1s, allRV1s_minerr, allRV1s_maxerr, HJD_values):
		if forced_K1==False:
			lines_to_write.append(str(ahjd) + "\t" + str(arvresult1) + "\t" + str(arvresult1_min) + "\t" + str(arvresult1_max) + "\n")
		else:
			lines_to_write.append(str(ahjd) + "\t" + str(arvresult1) + "\n")

	

	try: os.mkdir("out")
	except: None

	if len(sys_args)==4 and "one_by_oneMCMC" in sys.argv[2]:  spectrum_number = str(sys.argv[-1])
	else:  spectrum_number = ""
	
	with open("out/result"+spectrum_number+".out", 'w') as output_file:
		for aline in lines_to_write:
			output_file.write(aline)


		

if sys.argv[1]=="ATM" or sys.argv[1]=="plotOnly":
	result = open(os.getcwd()+"/out/result.out").readlines()
	allRV1s, allRV2s = [], []
	allRV1sminus, allRV2sminus =[], []
	allRV1smax, allRV2smax =[], []
	HoverHe1_med, HoverHe2_med = 0, 0
	for iii, lll in enumerate(result):
		linenew=lll.split("\n")[0]
		if "T0 " in linenew:
			T0=float(result[iii+1].split("\t")[0])
			T0_minus=float(result[iii+1].split("\t")[1])
			T0_plus=float(result[iii+1].split("\t")[2])
		elif "P " in linenew:
			Porb=float(result[iii+1].split("\t")[0])  * 60 * 60
			Porb_minus=float(result[iii+1].split("\t")[1])  * 60 * 60
			Porb_plus=float(result[iii+1].split("\t")[2])  * 60 * 60
		elif "K1" in linenew:
			K1=float(result[iii+1].split("\t")[0]) * 1000
			K1_minus=float(result[iii+1].split("\t")[1]) * 1000
			K1_plus=float(result[iii+1].split("\t")[2]) * 1000
		elif "K2" in linenew:
			K2=float(result[iii+1].split("\t")[0]) * 1000
			K2_minus=float(result[iii+1].split("\t")[1]) * 1000
			K2_plus=float(result[iii+1].split("\t")[2]) * 1000
		elif "Vg1" in linenew:
			vgamm1=float(result[iii+1].split("\t")[0]) * 1000
			vgamm1_minus=float(result[iii+1].split("\t")[1]) * 1000
			vgamm1_plus=float(result[iii+1].split("\t")[2]) * 1000
		elif "Vg2" in linenew:
			vgamm2=float(result[iii+1].split("\t")[0]) * 1000
			vgamm2_minus=float(result[iii+1].split("\t")[1]) * 1000
			vgamm2_plus=float(result[iii+1].split("\t")[2]) * 1000
		elif "RV1:" in lll:
			for cccc, ffff in enumerate(result[iii+1:]):
				if "RV2:" in ffff: break
				rrrv1_med=float(result[iii+1+cccc].split("\t")[1])
				allRV1s.append(rrrv1_med)
				allRV1sminus.append(float(result[iii+1+cccc].split("\t")[2])-rrrv1_med)
				allRV1smax.append(float(result[iii+1+cccc].split("\t")[3])-rrrv1_med)
		elif "RV2:" in lll:
			for cccc, ffff in enumerate(result[iii+1:]):
				rrrv2_med=float(result[iii+1+cccc].split("\t")[1])
				allRV2s.append(rrrv2_med)
				allRV2sminus.append(float(result[iii+1+cccc].split("\t")[2])-rrrv2_med)
				allRV2smax.append(float(result[iii+1+cccc].split("\t")[3])-rrrv2_med)
		elif "Scaling:" in lll:
			Scaling_med=float(result[iii+1].split("\t")[0])
			Scaling_minus=float(result[iii+1].split("\t")[1])
			Scaling_plus=float(result[iii+1].split("\t")[2])
			found_out_scaling=True



		if starType1.startswith("D") or starType1.startswith("sd") :
			if "Teff1:" in lll:
				T1_med=float(result[iii+1].split("\t")[0])
				T1_minus=float(result[iii+1].split("\t")[1])
				T1_plus=float(result[iii+1].split("\t")[2])
			elif "Logg1:" in lll:
				logg1_med=float(result[iii+1].split("\t")[0])
				logg1_minus=float(result[iii+1].split("\t")[1])
				logg1_plus=float(result[iii+1].split("\t")[2])
			elif "H/He1:" in lll:
				HoverHe1_med=float(result[iii+1].split("\t")[0])
				HoverHe1_minus=float(result[iii+1].split("\t")[1])
				HoverHe1_plus=float(result[iii+1].split("\t")[2])
			elif "Parallax:" in lll:
				parallax_med=float(result[iii+1].split("\t")[0])
				parallax_minus=float(result[iii+1].split("\t")[1])
				parallax_plus=float(result[iii+1].split("\t")[2])
				found_out_parallax=True
		else:
			if "A1_1:" in lll:
				A1_1_med=float(result[iii+1].split("\t")[0])
				A1_1_minus=float(result[iii+1].split("\t")[1])
				A1_1_plus=float(result[iii+1].split("\t")[2])
			elif "sig1_1:" in lll:
				sig1_1_med=float(result[iii+1].split("\t")[0])
				sig1_1_minus=float(result[iii+1].split("\t")[1])
				sig1_1_plus=float(result[iii+1].split("\t")[2])
			elif "A1_2:" in lll:
				A1_2_med=float(result[iii+1].split("\t")[0])
				A1_2_minus=float(result[iii+1].split("\t")[1])
				A1_2_plus=float(result[iii+1].split("\t")[2])
			elif "sig1_2:" in lll:
				sig1_2_med=float(result[iii+1].split("\t")[0])
				sig1_2_minus=float(result[iii+1].split("\t")[1])
				sig1_2_plus=float(result[iii+1].split("\t")[2])
			
	
	
	if starType1.startswith("D") or starType1.startswith("sd"):
		if forced_Scaling==False:	Scaling_med = Scaling_med
		elif forced_Scaling=="WD": 
			if starType1=="DA":      R1 = get_MTR(T1_med, logg=logg1_med, return_R_from_T_logg=True, loaded_Istrate=loaded_Istrate, loaded_CO=loaded_CO, loaded_Althaus=loaded_Althaus)
			elif starType1=="DBA" or starType1=="DC"  or starType1=="DB":   R1 = get_MTR_DB(T1_med, logg=logg1_med)
		else: Scaling_med=forced_Scaling
	
	else:
		Scaling_med = 1



	if starType1.startswith("D")  or starType1.startswith("sd"): # no starType1.startswith("sd") yet
		for tempcnt, (temperature_star, logg_star, starType) in enumerate(zip([T1_med],[logg1_med], [starType1])):
			minval=10000;   maxval=-10000
			if starType=="DA":
				if not pier_or_antoine=="pier3Dphot_antoine1Dspec":
					if tempcnt==0:     Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DAgrids(temperature_star, logg_star)
				else:
					if tempcnt==0:     Grav1_N, wl_all1_N, flux1_N, Teff1_N, photGrav1_N, photwl_all1_N, photflux1_N, photTeff1_N = return_DAgrids(temperature_star, logg_star)
					

			elif starType=="DBA" or starType1.startswith("sd"):
				if tempcnt==0:    HoverHestar=HoverHe1_med;    Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N = return_DBAgrids(temperature_star, logg_star, HoverHestar)
			
			elif starType=="DB":
				if tempcnt==0:     Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DBgrids(temperature_star, logg_star)
				
			elif starType=="DC":
				if tempcnt==0:     Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DCgrids(temperature_star, logg_star)
			




	# if relevant, plot the photometric solution
	if fit_phot_SED:  # no starType1.startswith("sd") yet
		if starType1=="DA" or starType1=="DB" or starType1=="DC":
			smeared_wl, smeared_flux = Fit_phot.fit_phot_SED_single(Grav1_N, wl_all1_N, flux1_N, Teff1_N, None, T1_med, logg1_med, None, min_wl=theminww, max_wl=themaxww, starType1=starType1, R1=R1, parallax=parallax_med, red=reddening_Ebv)
		elif starType1=="DBA":
			smeared_wl, smeared_flux = Fit_phot.fit_phot_SED_single(Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N, T1_med, logg1_med, HoverHe1_med, min_wl=theminww, max_wl=themaxww,starType1=starType1, R1=R1, parallax=parallax_med, red=reddening_Ebv)
		elif starType1.startswith("sd"):
			smeared_wl, smeared_flux = Fit_phot.fit_phot_SED_single(Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N, T1_med, logg1_med, HoverHe1_med, min_wl=theminww, max_wl=themaxww,starType1=starType1, R1=R_med, parallax=parallax_med, red=reddening_Ebv)
		
		
		
		rchisq_phot, chisq_phot = Fit_phot.process_photometry_in_each_pb(smeared_wl, smeared_flux, sedfilter, sed_wl, sedflux, sedfluxe, filter_dict=filter_dict, plot_solution=True, theminww_plot=theminww+50, themaxww_plot=themaxww+50, single_or_double="single", return_points_for_phot_model=False)










	# go and plot the spectra with the fits, don't plot anything with a shared wavelength, unless it's the gaussian/lorentzian fits
	for ii, (in_fi, sh_rv, ref_wl, modha, cut_lim, norm_lim, mask_out_Ha_min, mask_out_Ha_max, aHJD, inp_resolution, sigclipspec) in enumerate(zip(input_files, share_rv, reference_wl, modelHa, cut_Ha_all, normaliseHa_all, mask_out_Ha_min_all, mask_out_Ha_max_all, HJD_values, resolutions, sigma_clip)):
		if ((starType1.startswith("GG") or starType1.startswith("LL") or starType1.startswith("GL")) and sh_rv==-1) or (sh_rv==-1 and not ii in share_rv):
			norm_limits_min, norm_limits_max = norm_lim
			cut_limits_min, cut_limits_max = cut_lim
			modelHa_min, modelHa_max = modha
			
			fig, (ax, ax2) = plt.subplots(2)
			
			if starType1.startswith("D")  or starType1.startswith("sd"):
				if starType1=="DA" or starType1=="DB" or starType1=="DC": model_wl1, model_spectrum_star1 = return_model_spectrum_DA(wl_all1_N, ref_wl, cut_limits_min, cut_limits_max, Grav1_N, flux1_N, Teff1_N, T1_med, logg1_med)
				elif starType1=="DBA": model_wl1, model_spectrum_star1 = return_model_spectrum_DBA(wl_all1_N, ref_wl, cut_limits_min, cut_limits_max, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, T1_med, logg1_med, HoverHe1_med)
				elif starType1.startswith("sd"): model_wl1, model_spectrum_star1 = return_model_spectrum_subdwarf(wl_all1_N, ref_wl, cut_limits_min, cut_limits_max, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, T1_med, logg1_med, HoverHe1_med)
				
			
			else:
				model_wl1 = np.linspace(cut_limits_min, cut_limits_max, int(cut_limits_max-cut_limits_min)*10) #+ ref_wl
				if starType1=="GG":     model_spectrum_star1 = gauss_2x_mcmc(model_wl1, A1_1_med, sig1_1_med, A1_2_med, sig1_2_med)
				elif starType1=="LL":   model_spectrum_star1 = lorenz_2x_mcmc(model_wl1, A1_1_med, sig1_1_med, A1_2_med, sig1_2_med)
				elif starType1=="GL":   model_spectrum_star1 = gauss_lorentz_mcmc(model_wl1, A1_1_med, sig1_1_med, A1_2_med, sig1_2_med)
				model_wl1 = model_wl1 + ref_wl
			
			
			sh_wl = share_rv[ii]
			
			
			
			if sh_wl == -1:  rv_med1 = allRV1s[ii]
			else:            rv_med1 = allRV1s[sh_wl]
			
			
			
			filename = input_files[ii]
			if "/" in filename:    filename=filename.split("/")[-1]
			normalised_wavelength = list_norm_wl_grids[ii];    normalised_flux = list_normalised_flux[ii];    normalised_err = list_normalised_err[ii]

			
			
			
			dlam1 = ref_wl*rv_med1/speed_of_light
			desired_range = (normalised_wavelength>ref_wl+modelHa_min)  &  (normalised_wavelength<ref_wl+modelHa_max)

			
			#####old: DO NOT DELETE BECAUSE IT IS A TINY TINY TINY IMPROVEMENT THAT MAY BE WANTED FOR LINE CORE STUFF
			## interpolate the model onto the observation's wavelength grid
			#interparr = interp(model_wl1, model_wl1+dlam1, model_spectrum_star1)
			## Smear to the resolution desired
			#interparr = convolve(interparr, Gaussian1DKernel(stddev=0.5*(ref_wl/inp_resolution)/(model_wl1[10]-model_wl1[9])), boundary = 'extend')
			#interparr = interp(normalised_wavelength, model_wl1, interparr)
				
			# Smear to the resolution desired
			interparmodel = convolve(model_spectrum_star1, Gaussian1DKernel(stddev=0.5*(ref_wl/inp_resolution)/(model_wl1[10]-model_wl1[9])), boundary = 'extend')
			interparr = interp(normalised_wavelength, model_wl1+dlam1, interparmodel)
			
			if starType1=="GG" or starType1=="LL" or starType1=="GL":   interparr += 1
			
			
			
			
			mask = ((normalised_wavelength>ref_wl+cut_limits_min) & (normalised_wavelength<ref_wl+norm_limits_min)) | ((normalised_wavelength>ref_wl+norm_limits_max)   & (normalised_wavelength<ref_wl+cut_limits_max))
				
			try:   m,c =  polyfit(normalised_wavelength[mask], interparr[mask], deg=1)
			except Exception as eee:
				plt.clf()
				plt.plot(normalised_wavelength, normalised_flux);    plt.plot(normalised_wavelength, interparr)
				plt.axvline(ref_wl+cut_limits_min, c='k');     plt.axvline(ref_wl+cut_limits_max, c='k')
				plt.axvline(ref_wl+norm_limits_min, c='r');    plt.axvline(ref_wl+norm_limits_max, c='r')
				plt.title(str() + in_fi + " " + str(ref_wl) + "    "+str(cut_limits_min) + "  " + str(norm_limits_min) + "  " + str(norm_limits_max) + "  " + str(cut_limits_max))
				plt.show();    plt.clf()
				m,c =  polyfit(normalised_wavelength[mask], interparr[mask], deg=1)
			interparr /= (m*normalised_wavelength + c)
			#interparr/=np.amax(interparr)
			
			
			mask1 = (normalised_wavelength>ref_wl+cut_limits_min) & (normalised_wavelength<ref_wl+norm_limits_min)
			mask2 = (normalised_wavelength>ref_wl+norm_limits_max) & (normalised_wavelength<ref_wl+cut_limits_max)
			
			if False:
				plt.clf()
				plt.plot(normalised_wavelength[mask1], interparr[mask1], c='g')
				plt.plot(normalised_wavelength[mask2], interparr[mask2], c='r')
				plt.show()
			
			
			
			normwl_to_rv = speed_of_light * (normalised_wavelength-ref_wl)/ref_wl 

			mask_norm = ((normalised_wavelength>ref_wl+cut_limits_min) & (normalised_wavelength<ref_wl+norm_limits_min) ) | (normalised_wavelength<ref_wl+cut_limits_max) & ((normalised_wavelength>ref_wl+norm_limits_max))
			

			ax.plot(normalised_wavelength - ref_wl,normalised_flux, c='k')
			ax2.plot(speed_of_light * (normalised_wavelength - ref_wl)/ref_wl,normalised_flux, c='k')
			x,y =normalised_wavelength[~desired_range],  interparr[~desired_range]
			ax.plot(x[x<ref_wl] - ref_wl, y[x<ref_wl], c='b', label="_blue = Cont. model fit");    ax.plot(x[x>ref_wl] - ref_wl, y[x>ref_wl], c='b')
			ax2.plot(speed_of_light * (x[x<ref_wl] - ref_wl)/ref_wl, y[x<ref_wl], c='b', label="_blue = Cont. model fit");    ax2.plot(speed_of_light * (x[x>ref_wl] - ref_wl)/ref_wl, y[x>ref_wl], c='b')
			x, y = normalised_wavelength[mask_norm], normalised_flux[mask_norm]
			ax.plot(x[x<ref_wl] - ref_wl, y[x<ref_wl], c='r', label="_red = For normalising");    ax.plot(x[x>ref_wl] - ref_wl, y[x>ref_wl], c='r')
			ax2.plot(speed_of_light * (x[x<ref_wl] - ref_wl)/ref_wl, y[x<ref_wl], c='r', label="_red = For normalising");    ax2.plot(speed_of_light* (x[x>ref_wl] - ref_wl)/ref_wl, y[x>ref_wl], c='r',zorder=100)
			ax.plot(normalised_wavelength[desired_range] - ref_wl, interparr[desired_range], c='orange');    ax.set_title(filename)
			ax2.plot(speed_of_light * (normalised_wavelength[desired_range] - ref_wl)/ref_wl, interparr[desired_range], c='orange')
			
			if sigclipspec!=-1:
				resid = normalised_flux-interparr
				resid_in_sig = resid/np.std(np.abs(resid[desired_range]))
				clip_mask = (np.abs(resid_in_sig)>sigclipspec) & desired_range & ((normalised_wavelength > ref_wl+5) | (normalised_wavelength < ref_wl-5))
				ax.scatter(normalised_wavelength[clip_mask] - ref_wl, normalised_flux[clip_mask], c='r', s=4)
				ax2.scatter(speed_of_light * (normalised_wavelength[clip_mask] - ref_wl)/ref_wl, normalised_flux[clip_mask], c='r', s=4)

			if starType1=="DBA" or starType1.startswith("sd"):    ax.scatter(0,1,alpha=0, label="H/He1 = " + str(HoverHe1_med)); ax2.scatter(0,1,alpha=0, label="H/He1 = " + str(HoverHe1_med))
			
			
			ax.scatter(0,1,alpha=0, label="RV1 = " + str(np.round(rv_med1,3)));        ax2.scatter(0,1,alpha=0, label="RV1 = " + str(np.round(rv_med1,3)))
			ax.scatter(0,1,alpha=0, label="refwave = " + str(ref_wl));                 ax2.scatter(0,1,alpha=0, label="refwave = " + str(ref_wl))
			if starType1.startswith("D")  or starType1.startswith("sd"):
				ax.scatter(0,1,alpha=0, label="T1 = " + str(np.round(T1_med,0)));          ax2.scatter(0,1,alpha=0, label="T1 = " + str(np.round(T1_med,0)))
				ax.scatter(0,1,alpha=0, label="Logg1 = " + str(np.round(logg1_med,3)));    ax2.scatter(0,1,alpha=0, label="Logg1 = " + str(np.round(logg1_med,3)))
			ax2.set_xlabel("RV [kms-1]");
			ax.legend(loc="lower right", handlelength=0)
			#ax2.legend(loc="lower right", handlelength=0)
			trim_plot2_xvals_min, trim_plot2_xvals_max = -2500, 2500
			ax2.set_xlim(trim_plot2_xvals_min, trim_plot2_xvals_max)
			ax.set_ylim(0.3,1.2);    ax2.set_ylim(0.35,1.1)
			plt.savefig("out/"+filename.split(".dat")[0]+"_"+str(ref_wl)+".png", dpi=300)
			if ii in plot_fit and not False in plot_fit:   plt.show()
			plt.clf();    plt.close()







	#### Now plot just the shared rv combinations
	where_share_rv = (share_rv != -1)
	plt.clf()
	if starType1.startswith("D")  or starType1.startswith("sd"):
		for unique_shared_rv in np.unique(share_rv[where_share_rv]):
			norm_wl_combinations, normalised_flux_combinations, normalised_err_combinations = [], [], []
			fig = plt.figure()
			gs = fig.add_gridspec(2, 2)
			ax = fig.add_subplot(gs[0, :])
			ax2 = fig.add_subplot(gs[1, 0])
			ax3 = fig.add_subplot(gs[1, 1])
			#fig, (ax, ax2) = plt.subplots(2)
			select_shared_mask = (share_rv==unique_shared_rv) | (range(len(share_rv))==unique_shared_rv)
			
			
			mask_out_Ha_min_all_combinations, mask_out_Ha_max_all_combinations  =  np.asarray(mask_out_Ha_min_all)[select_shared_mask], np.asarray(mask_out_Ha_max_all)[select_shared_mask]
			
			reference_wl_combinations = reference_wl[select_shared_mask]
			cut_Ha_combinations = cut_Ha_all[select_shared_mask]
			normaliseHa_combinations = normaliseHa_all[select_shared_mask]
			resolutions_combinations = resolutions[select_shared_mask]
			input_files_combinations = input_files[select_shared_mask]
			modelHa_combinations = modelHa[select_shared_mask]
			sigma_clip_combinations = sigma_clip[select_shared_mask]
			
			for cn_t_or_f, t_or_f in enumerate(select_shared_mask):	
				if t_or_f==True:
					norm_wl_combinations.append(list_norm_wl_grids[cn_t_or_f]); normalised_flux_combinations.append(list_normalised_flux[cn_t_or_f]), normalised_err_combinations.append(list_normalised_err[cn_t_or_f])


			for cn_shared, (in_fi, ref_wl, modha, cut_lim, norm_lim, normalised_wavelength, normalised_flux, normalised_err, inp_resolution, mask_out_Ha_min, mask_out_Ha_max, sigclipspec) in enumerate(zip(input_files_combinations, reference_wl_combinations, modelHa_combinations, cut_Ha_combinations, normaliseHa_combinations, norm_wl_combinations, normalised_flux_combinations, normalised_err_combinations, resolutions_combinations, mask_out_Ha_min_all_combinations, mask_out_Ha_max_all_combinations, sigma_clip_combinations)):
				cut_limits_min, cut_limits_max = cut_lim
				norm_limits_min, norm_limits_max = norm_lim
				modelHa_min, modelHa_max = modha
				
				
				if starType1=="DA" or starType1=="DB" or starType1=="DC": model_wl1, model_spectrum_star1 = return_model_spectrum_DA(wl_all1_N, ref_wl, cut_limits_min, cut_limits_max, Grav1_N, flux1_N, Teff1_N, T1_med, logg1_med)
				elif starType1=="DBA": model_wl1, model_spectrum_star1 = return_model_spectrum_DBA(wl_all1_N, ref_wl, cut_limits_min, cut_limits_max, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, T1_med, logg1_med, HoverHe1_med)
				elif starType1.startswith("sd"):  model_wl1, model_spectrum_star1 = return_model_spectrum_subdwarf(wl_all1_N, ref_wl, cut_limits_min, cut_limits_max, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, T1_med, logg1_med, HoverHe1_med)
				
					
				


				offset=cn_shared*-0.3 * 1.05**cn_shared
				
				
				rv_med1 = allRV1s[unique_shared_rv]
				
				
				
				filename = in_fi
				if "/" in filename:    filename=filename.split("/")[-1]




				dlam1 = ref_wl*rv_med1/speed_of_light
				desired_range = (normalised_wavelength>ref_wl+modelHa_min)  &  (normalised_wavelength<ref_wl+modelHa_max)

				
				#####old: DO NOT DELETE BECAUSE IT IS A TINY TINY TINY IMPROVEMENT THAT MAY BE WANTED FOR LINE CORE STUFF
				## interpolate the model onto the observation's wavelength grid
				#interparr = interp(model_wl1, model_wl1+dlam1, model_spectrum_star1)
				## Smear to the resolution desired
				#interparr = convolve(interparr, Gaussian1DKernel(stddev=0.5*(ref_wl/inp_resolution)/(model_wl1[10]-model_wl1[9])), boundary = 'extend')
				#interparr = interp(normalised_wavelength, model_wl1, interparr)
					
				# Smear to the resolution desired
				interparmodel = convolve(model_spectrum_star1, Gaussian1DKernel(stddev=0.5*(ref_wl/inp_resolution)/(model_wl1[10]-model_wl1[9])), boundary = 'extend')
				interparr = interp(normalised_wavelength, model_wl1+dlam1, interparmodel)
				
				
				
				
				mask = ((normalised_wavelength>ref_wl+cut_limits_min) & (normalised_wavelength<ref_wl+norm_limits_min)) | ((normalised_wavelength>ref_wl+norm_limits_max)   & (normalised_wavelength<ref_wl+cut_limits_max))
				#mask = mask & (interparr[mask]/np.median(interparr[mask]) > 0.9)
				
				#mask1 = ((normalised_wavelength>ref_wl+norm_limits_max)   & (normalised_wavelength<ref_wl+cut_limits_max))
				#mask2 = ((normalised_wavelength>ref_wl+cut_limits_min) & (normalised_wavelength<ref_wl+norm_limits_min))
				#plt.plot(normalised_wavelength, interparr, c='g', ls='--')
				#plt.plot(normalised_wavelength[mask], interparr[mask1], c='r');   plt.plot(normalised_wavelength[mask], interparr[mask2], c='r')
				#plt.title(str(ref_wl) + "  " + str(ref_wl+cut_limits_min) + "  " + str(ref_wl+norm_limits_min) + "  " + str(ref_wl+norm_limits_max) + "  " + str(ref_wl+cut_limits_max))
				#plt.show();   plt.clf()
				
				
				try:
					m,c =  polyfit(normalised_wavelength[mask], interparr[mask], deg=1)
				except Exception as eee:
					print(eee)
					plt.clf()
					plt.plot(normalised_wavelength, normalised_flux)
					plt.plot(normalised_wavelength, interparr)
					plt.axvline(ref_wl+cut_limits_min, c='k')
					plt.axvline(ref_wl+cut_limits_max, c='k')
					plt.axvline(ref_wl+norm_limits_min, c='r')
					plt.axvline(ref_wl+norm_limits_max, c='r')
					plt.title(str() + in_fi + " " + str(ref_wl) + "    "+str(cut_limits_min) + "  " + str(norm_limits_min) + "  " + str(norm_limits_max) + "  " + str(cut_limits_max))
					plt.show()
					plt.clf()
					
					m,c =  polyfit(normalised_wavelength[mask], interparr[mask], deg=1)
				
				interparr /= (m*normalised_wavelength + c)
				

				mask_norm = ((normalised_wavelength>ref_wl+cut_limits_min) & (normalised_wavelength<ref_wl+norm_limits_min) ) | (normalised_wavelength<ref_wl+cut_limits_max) & ((normalised_wavelength>ref_wl+norm_limits_max))
				
				

				ax.plot(normalised_wavelength - ref_wl,normalised_flux + offset, c='k')
				ax2.plot(speed_of_light * (normalised_wavelength - ref_wl)/ref_wl,normalised_flux + offset, c='k')
				if ref_wl==np.amax(reference_wl_combinations):
					an_x = speed_of_light * (normalised_wavelength - ref_wl)/ref_wl
					amask = (an_x>-1750) & (an_x<1750)
					ax3.plot(an_x[amask], normalised_flux[amask] + offset, c='k')
				x, y = normalised_wavelength[mask_norm], normalised_flux[mask_norm]
				x,y =normalised_wavelength[~desired_range],  interparr[~desired_range]
				ax.plot(x[x<ref_wl] - ref_wl, y[x<ref_wl] + offset, c='b', label="_blue = Continued model fit")
				ax2.plot(speed_of_light * (x[x<ref_wl] - ref_wl)/ref_wl, y[x<ref_wl] + offset, c='b', label="_blue = Continued model fit")
				if ref_wl==np.amax(reference_wl_combinations):
					an_x = speed_of_light * (x[x<ref_wl] - ref_wl)/ref_wl
					a_y = y[x<ref_wl] + offset
					amask = (an_x>-1750) & (an_x<1750)
					ax3.plot(an_x[amask], a_y[amask], c='b', label="_blue = Continued model fit")
				ax.plot(x[x>ref_wl] - ref_wl, y[x>ref_wl] + offset, c='b')
				ax2.plot(speed_of_light * (x[x>ref_wl] - ref_wl)/ref_wl, y[x>ref_wl] + offset, c='b')
				if ref_wl==np.amax(reference_wl_combinations):
					an_x = speed_of_light * (x[x>ref_wl] - ref_wl)/ref_wl
					amask = (an_x>-1750) & (an_x<1750)
					a_y = y[x>ref_wl] + offset
					ax3.plot(an_x[amask], a_y[amask], c='b')
				ax.plot(x[x<ref_wl] - ref_wl, y[x<ref_wl] + offset, c='r', label="_red = For normalising")
				ax2.plot(speed_of_light * (x[x<ref_wl] - ref_wl)/ref_wl, y[x<ref_wl] + offset, c='r', label="_red = For normalising")
				if ref_wl==np.amax(reference_wl_combinations):
					an_x = speed_of_light * (x[x<ref_wl] - ref_wl)/ref_wl
					amask = (an_x>-1750) & (an_x<1750)
					a_y = y[x<ref_wl] + offset
					ax2.plot(an_x[amask], a_y[amask], c='r', label="_red = For normalising")
				ax.plot(x[x>ref_wl] - ref_wl, y[x>ref_wl] + offset, c='r')
				ax2.plot(speed_of_light * (x[x>ref_wl] - ref_wl)/ref_wl, y[x>ref_wl] + offset, c='r')
				if ref_wl==np.amax(reference_wl_combinations):
					an_x = speed_of_light * (x[x>ref_wl] - ref_wl)/ref_wl
					amask = (an_x>-1750) & (an_x<1750)
					a_y = y[x>ref_wl] + offset
					ax3.plot(an_x[amask], a_y[amask], c='r')
				ax.plot(normalised_wavelength[desired_range] - ref_wl, interparr[desired_range] + offset, c='orange')
				ax2.plot(speed_of_light * (normalised_wavelength[desired_range] - ref_wl)/ref_wl, interparr[desired_range] + offset, c='orange')
				if ref_wl==np.amax(reference_wl_combinations):
					an_x = speed_of_light * (normalised_wavelength[desired_range] - ref_wl)/ref_wl
					amask = (an_x>-1750) & (an_x<1750)
					a_y = interparr[desired_range] + offset
					ax3.plot(an_x[amask], a_y[amask], c='orange');  ax3.text(-400,1, str(ref_wl))
					
				
				
				
				if sigclipspec!=-1:
					resid = normalised_flux-interparr
					resid_in_sig = resid/np.std(np.abs(resid[desired_range]))
					clip_mask = (np.abs(resid_in_sig)>sigclipspec) & desired_range & ((normalised_wavelength > ref_wl+5) | (normalised_wavelength < ref_wl-5))
					ax.scatter(normalised_wavelength[clip_mask] - ref_wl, normalised_flux[clip_mask] + offset, c='r', s=4)
					ax2.scatter(speed_of_light * (normalised_wavelength[clip_mask] - ref_wl)/ref_wl, normalised_flux[clip_mask] + offset, c='r', s=4)
					if ref_wl==np.amax(reference_wl_combinations):
						xx=speed_of_light * (normalised_wavelength[clip_mask] - ref_wl)/ref_wl
						amask = (xx>-1750) & (xx<1750)
						yy=normalised_flux[clip_mask] + offset
						ax3.scatter(xx[amask], yy[amask], c='r', s=8)
				
				
			
			
			if starType1=="DBA" or starType1.startswith("sd"):    plt.scatter(0,1,alpha=0, label="H/He1 = " + str(HoverHe1_med))
			
			
			ax.scatter(0,1,alpha=0, label="RV1 = " + str(np.round(rv_med1,3)));       ax2.scatter(0,1,alpha=0, label="RV1 = " + str(np.round(rv_med1,3)))
			ax.scatter(0,1,alpha=0, label="refwave = " + str(ref_wl));                ax2.scatter(0,1,alpha=0, label="refwave = " + str(ref_wl))
			if starType1.startswith("D")  or starType1.startswith("sd"):
				ax.scatter(0,1,alpha=0, label="T1 = " + str(np.round(T1_med,0)));         ax2.scatter(0,1,alpha=0, label="T1 = " + str(np.round(T1_med,0)))
				ax.scatter(0,1,alpha=0, label="Logg1 = " + str(np.round(logg1_med,3)));   ax2.scatter(0,1,alpha=0, label="Logg1 = " + str(np.round(logg1_med,3)))
			ax.set_title(filename)
			ax3.set_xlim(-1750,1750)
			ax2.set_xlabel("RV [kms-1]");   ax3.set_xlabel("RV [kms-1]")
			ax.legend(loc='lower right', handlelength=0, fontsize=10)
			plt.savefig("out/"+filename.split(".dat")[0]+"_all_shared_rvs"+".png", dpi=300)
			if unique_shared_rv in plot_fit and not False in plot_fit:   plt.show()
			plt.clf();   plt.close()





