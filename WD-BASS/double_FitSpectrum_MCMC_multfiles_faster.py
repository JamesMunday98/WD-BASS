import numpy as np
from numpy import interp, polyfit
import os, sys, yaml, emcee, corner, warnings
sys.path.append(os.environ['WD_BASS_INSTALL_DIR']+"/scripts")
from load_All_Models import load_models
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from schwimmbad import MPIPool
from importable_MTR_function import get_MTR, get_MTR_DB
warnings.filterwarnings('ignore')
from Fit_Photometric_SED import Fit_phot
from checkLocalStars import checkLocalStars
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.convolution import Gaussian1DKernel,convolve
from numba import njit, jit
from scipy.optimize import curve_fit
from numpy import amin as npamin, amax as npamax, unique as npunique, argwhere as npargwhere, loadtxt as nploadtxt, linspace as nplinspace, pi as np_pi, inf as np_inf, sin as np_sin, square as np_square, sqrt as npsqrt, mean as npmean, power as np_power, sum as npsum
from dust_extinction.parameter_averages import G23
from mpi4py import MPI
from functools import partial


if sys.argv[1]=="ATM" or sys.argv[1]=="photometry_only":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
else: rank=1

sys_args = sys.argv

def print_exception_details(e):
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, "\t", exc_tb.tb_lineno, "\t", e)


with open('example_Config_dbl.yaml') as file:
    config_info = yaml.load(file, Loader=yaml.FullLoader)

speed_of_light = 2.99792458E5

spectra_source_type = config_info["spectra_source_type"]
if spectra_source_type!="wl_flux_fluxerr" and spectra_source_type!="Xshooter": raise ValueError("spectra_source_type is not recognised in the config file")
input_files = np.asarray(config_info["spectraHa"])
reference_wl = np.asarray(config_info["refwave"])
modelHa = np.asarray(config_info["modelHa"])
normaliseHa_all = np.asarray(config_info["normaliseHa"])
cut_Ha_all = np.asarray(config_info["cut_Ha_max"])
resolutions = np.asarray(config_info["resolutions"])
#resolutions*=2;   print("WARNING WARNING resolution doubled")
share_rv = np.asarray(config_info["share_rv"])
starType1, starType2 = np.asarray(config_info["starType"])
RV_boundaries1 = np.asarray(config_info["RV_boundaries1"])
RV_boundaries2 = np.asarray(config_info["RV_boundaries2"])
HJD_values = np.asarray(config_info["HJD_Values"])
forced_teff1, forced_teff2 = np.asarray(config_info["forced_teff"])
forced_logg1, forced_logg2 = np.asarray(config_info["forced_logg"])
forced_HoverHe1, forced_HoverHe2 = np.asarray(config_info["forced_HoverHe"])
plot_corner=np.asarray(config_info["plot_corner"])
forced_Scaling=np.asarray(config_info["forced_scaling"])
if forced_Scaling!=False and forced_Scaling!="WD": forced_Scaling=float(forced_Scaling)
forced_ephemeris=np.asarray(config_info["forced_ephemeris"]) # expects False or [T0, P0]
if forced_ephemeris[0]==False:  forced_ephemeris=forced_ephemeris[0];  forced_T0=None;  forced_P0=None
else: forced_T0=float(forced_ephemeris[0]);  forced_P0=float(forced_ephemeris[1])
forced_K1=config_info["forced_K1"][0]
if not forced_K1==False and not forced_K1=="Fit":  forced_K1=float(forced_K1)
forced_K2=config_info["forced_K2"][0]
if not forced_K2==False and not forced_K2=="Fit":  forced_K2=float(forced_K2)
K1_boundaries=np.asarray(config_info["K1_boundaries"]);             p0K1=np.array([float(K1_boundaries[0]), float(K1_boundaries[1])])
K2_boundaries=np.asarray(config_info["K2_boundaries"]);             p0K2=np.array([float(K2_boundaries[0]), float(K2_boundaries[1])])
forced_Vgamma1=config_info["forced_Vgamma1"][0]
if not forced_Vgamma1=="Fit":  forced_Vgamma1=float(forced_Vgamma1)
Vgamma1_boundaries=np.asarray(config_info["Vgamma1_boundaries"]);   p0Vgamma1=np.array([float(Vgamma1_boundaries[0]), float(Vgamma1_boundaries[1])])
forced_Vgamma2=config_info["forced_Vgamma2"][0]
if not forced_Vgamma2=="Fit":  forced_Vgamma2=float(forced_Vgamma2)
Vgamma2_boundaries=np.asarray(config_info["Vgamma2_boundaries"]);   p0Vgamma2=np.array([float(Vgamma2_boundaries[0]), float(Vgamma2_boundaries[1])])
plot_fit=np.asarray([config_info["plot_fit"]])[0]
fit_phot_SED=np.asarray([config_info["fit_phot_SED"]])
RA=np.asarray([config_info["RA"]])[0]
Dec=np.asarray([config_info["Dec"]])[0]
try: high_RV_amp=bool(np.asarray([config_info["high_RV_amp"]])[0])
except: high_RV_amp=False


if "run_double_one_by_oneMCMC" in sys.argv:
    file_ignore=np.array([])


## gaussian/lorentzian things
if "GG" in starType1 or "LL" in starType1 or "GL" in starType1:
    try:  
        A1_1_boundaries = np.asarray([config_info["A1_1_boundaries"]])[0]
        A1_2_boundaries = np.asarray([config_info["A1_2_boundaries"]])[0]
        sigma1_1_boundaries = np.asarray([config_info["sigma1_1_boundaries"]])[0]
        sigma1_2_boundaries = np.asarray([config_info["sigma1_2_boundaries"]])[0]
    except: raise ValueError("Missing boundaries for star1: A1_1_boundaries/sigma1_1_boundaries/A1_2_boundaries/sigma1_2_boundaries")
if "GG" in starType2 or "LL" in starType2 or "GL" in starType2:
    try:  
        A2_1_boundaries = np.asarray([config_info["A2_1_boundaries"]])[0]
        A2_2_boundaries = np.asarray([config_info["A2_2_boundaries"]])[0]
        sigma2_1_boundaries = np.asarray([config_info["sigma2_1_boundaries"]])[0]
        sigma2_2_boundaries = np.asarray([config_info["sigma2_2_boundaries"]])[0]
    except: raise ValueError("Missing boundaries for star2: A2_1_boundaries/sigma2_1_boundaries/A2_2_boundaries/sigma2_2_boundaries")
if starType1=="quadLorentz" and starType2=="quadLorentz":
    try:  
        A1_1_boundaries = np.asarray([config_info["A1_1_boundaries"]])[0]
        A2_1_boundaries = np.asarray([config_info["A2_1_boundaries"]])[0]
        sigma1_1_boundaries = np.asarray([config_info["sigma1_1_boundaries"]])[0]
        sigma2_1_boundaries = np.asarray([config_info["sigma2_1_boundaries"]])[0]
    except: raise ValueError("Missing boundaries for star2: A2_1_boundaries/sigma2_1_boundaries/A2_2_boundaries/sigma2_2_boundaries")
    

try:  reddening_Ebv=float(np.asarray([config_info["reddening_Ebv"]])[0])
except: reddening_Ebv=0

if reddening_Ebv=="lookup" or fit_phot_SED:
	from miscAstro import miscAstro
	if type(RA)==np.str_ and type(Dec)==np.str_:  RAdeg, Decdeg= miscAstro.ra_dec_hr_to_deg(RA,Dec)
	elif type(RA)==float and type(Dec)==float:  RAdeg, Decdeg= RA, Dec;  RA, Dec = miscAstro.ra_dec_deg_to_hr(RAdeg,Decdeg)
	else: print(type(RA), type(Dec));  raise ValueError

expected_Gmag=np.asarray([config_info["expected_Gmag"]])[0]
nwalkers=np.asarray([config_info["nwalkers"]])[0]
burnin=np.asarray([config_info["nburnin"]])[0]
nsteps=np.asarray([config_info["nsteps"]])[0]


if reddening_Ebv=="lookup":
    try:
        star_plx=float(np.asarray(config_info["parallax"])) # in mas
        star_e_plx=float(np.asarray(config_info["parallax_uncertainty"])) # in mas
    except:
        RAdeg, Decdeg, star_mag, star_plx, star_e_plx, warning = checkLocalStars.find_star_in_gaia_edr3(RAdeg, Decdeg, 10, predicted_Gmag=expected_Gmag)
        #raise ValueError(star_mag)
        #print(star_mag.value[0], star_plx.value[0], star_e_plx.value[0], warning)
    gal_coords = SkyCoord(ra=RAdeg*u.degree, dec=Decdeg*u.degree, frame='icrs').galactic
    gal_l=gal_coords.l.value; gal_b=gal_coords.b.value
    
    #complete_with_stilism("Stilism.csv")
    url = "http://stilism.obspm.fr/reddening?frame=galactic&vlong={}&ulong=deg&vlat={}&ulat=deg&distance={}"
    
    import requests
    from io import StringIO
    import pandas as pd
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
    try: want_galex_FUV = np.asarray(config_info["want_galex_FUV"])
    except: want_galex_FUV=False
else: found_parallax=False




if sys.argv[1]=="ATM" or sys.argv[1]=="plotOnly":

    
    if "run_double_one_by_oneMCMC" in sys.argv:   file_ignore=np.array([])
    else:   file_ignore=np.asarray([config_info["file_ignore"]])[0]
    

    if len(file_ignore) > 0:  dowhile=True
    else: dowhile=False
    
    
    if "run_double_one_by_oneMCMC" in sys.argv:
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
                original_file_index = npargwhere(original_input_files==aaaa)[0][0]
            except: raise ValueError(aaaa)
            for zcnt, aaa_sharerv in enumerate(original_share_rv):
                if aaa_sharerv>=original_file_index:
                    share_rv[zcnt]-=1
                
                if aaa_sharerv>=original_file_index and "run_double_one_by_oneMCMC" in sys.argv:    list_indexes[zcnt]-=1
            
            
            mask_ignore_files = input_files!=aaaa
            modelHa = modelHa[mask_ignore_files]
            normaliseHa_all=normaliseHa_all[mask_ignore_files]
            cut_Ha_all=cut_Ha_all[mask_ignore_files]
            resolutions=resolutions[mask_ignore_files]
            reference_wl=reference_wl[mask_ignore_files]
            share_rv=share_rv[mask_ignore_files]
            RV_boundaries1=RV_boundaries1[mask_ignore_files]
            RV_boundaries2=RV_boundaries2[mask_ignore_files]
            HJD_values=HJD_values[mask_ignore_files]
            input_files=input_files[mask_ignore_files]
            
            
            sigma_clip=sigma_clip[mask_ignore_files]
            if "run_double_one_by_oneMCMC" in sys.argv:
                list_indexes=list_indexes[mask_ignore_files]
                removed_minus1_pos = npargwhere(mask_ignore_files==False)[0][0]
                list_indexes[removed_minus1_pos:]-=1
            
            file_ignore_mask=file_ignore!=aaaa
            file_ignore=file_ignore[file_ignore_mask]
            
            
            stopwhile=False
            break
        if stopwhile: dowhile=False

    
    if "run_double_one_by_oneMCMC" in sys.argv:
        wanted_index = (share_rv==int(sys.argv[3])) | (list_indexes==int(sys.argv[3]))
        modelHa = modelHa[wanted_index]
        normaliseHa_all=normaliseHa_all[wanted_index]
        cut_Ha_all=cut_Ha_all[wanted_index]
        resolutions=resolutions[wanted_index]
        reference_wl=reference_wl[wanted_index]
        share_rv=share_rv[wanted_index]
        RV_boundaries1=RV_boundaries1[wanted_index]
        RV_boundaries2=RV_boundaries2[wanted_index]
        HJD_values=HJD_values[wanted_index]
        input_files=input_files[wanted_index]
        sigma_clip=sigma_clip[wanted_index]
        

if "run_double_one_by_oneMCMC" in sys.argv:
    while True:
        try: # need a try here for the times where there is only 1 spectrum for one by one which of course has share_rv=-1
            if npamin(share_rv[share_rv!=-1])!=0:   share_rv[share_rv!=-1]-=1
            else:   break
        except: break


if len(share_rv)==0:
    print("No RVs passed. Exiting script")
    sys.exit()




if rank==1:
	print("Processing :   ", len(input_files), "   spectra")
	for afi, aha, anorm, acut, ar1, ar2 in zip(input_files, modelHa, normaliseHa_all, cut_Ha_all, RV_boundaries1, RV_boundaries2):
	    print(afi, aha, anorm, acut, "\t", "RVbounds:", ar1, ar2)
try: os.mkdir("out")
except: None

if fit_phot_SED:
    if type(RA)==np.str_ and type(Dec)==np.str_:  RAdeg, Decdeg= miscAstro.ra_dec_hr_to_deg(RA,Dec)
    elif type(RA)==float and type(Dec)==float:  RAdeg, Decdeg= RA, Dec;  RA, Dec = miscAstro.ra_dec_deg_to_hr(RAdeg,Decdeg)
    else: print(type(RA), type(Dec));  raise ValueError
    #checkLocalStars.find_star_in_gaia_edr3(RAdeg, Decdeg, 5, predicted_Gmag=expected_Gmag)


if ((starType1 != "DA" and starType1 != "DBA" and starType1 != "DC") or (starType2 != "DA" and starType2 != "DBA" and starType2 != "DC")) and forced_Scaling == "WD":   raise ValueError("WD scaling only allowed for DA/DB/DC stars")
if (forced_K1=="Fit" and (forced_K2==False and not (forced_Vgamma2=="Fit" or isinstance(forced_Vgamma2, float)))) or ((forced_K1==False and not (forced_Vgamma1=="Fit" or isinstance(forced_Vgamma1, float))) and forced_K2=="Fit"):  raise ValueError("Error - knowing the ephemeris and fitting K1 + RV2s is weird. Not allowed")
if (forced_K1==False and forced_Vgamma1=="Fit") or (forced_K2==False and forced_Vgamma2=="Fit"):  raise ValueError("Can't fit Vgamma with K free")



if not len(input_files)==len(modelHa)==len(reference_wl)==len(normaliseHa_all)==len(cut_Ha_all)==len(resolutions)==len(share_rv)==len(RV_boundaries1)==len(RV_boundaries2)==len(sigma_clip)==len(HJD_values):
    print("ALERT");  print("ALERT"); print("ALERT"); print("ALERT"); print("ALERT")
    print("Mismatch of length of required arrays")
    print("input_files", len(input_files));  print("reference_wl", len(reference_wl))
    print("modelHa", len(modelHa))
    print("normaliseHa", len(normaliseHa_all));  print("cut_Ha_all", len(cut_Ha_all))
    print("resolutions", len(resolutions));  print("share_rv", len(share_rv))
    print("RV_boundaries1", len(RV_boundaries1));  print("RV_boundaries2", len(RV_boundaries2))
    print("HJD_values", len(HJD_values))
    print("sigma_clip", len(sigma_clip))
    raise ValueError


theminww_loadgrid, themaxww_loadgrid=1,1500000
if fit_phot_SED:

    url = Fit_phot.get_url_CDS(RA + " " + Dec,5)
    sedfilter, sed_wl, sedflux, sedfluxe, catalogue = Fit_phot.get_data(url, False)  # wl [nm], flux and error [Jy]
    if len(sedfilter)<=3:    url2 = Fit_phot.get_url_CDS(RA + " " + Dec,10);    sedfilter2, sed_wl2, sedflux2, sedfluxe2, catalogue2 = Fit_phot.get_data(url)
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
                elif "GALEX:FUV" in ff and want_galex_FUV:
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
                    m_wl, m_f, m_fe = npmean(sed_wl[mask1]), npmean(sedflux[mask1]), npmean(sedfluxe[mask1])
                    sedfilter, sed_wl, sedflux, sedfluxe = sedfilter[~mask1], sed_wl[~mask1], sedflux[~mask1], sedfluxe[~mask1]
                    if m_f!=0 and not np.isnan(m_fe):
                        sedfilter, sed_wl, sedflux, sedfluxe = np.append(sedfilter, filtername), np.append(sed_wl, m_wl), np.append(sedflux, m_f), np.append(sedfluxe, m_fe)
            
        
        sedfluxe[sedfluxe==0] = npamin(sedfluxe[sedfluxe!=0])

        try: 
            if np.asarray(config_info["plot_phot_spectrum"]):
                plt.errorbar(sed_wl, sedflux, sedfluxe,fmt='.k'); plt.xlabel("Wavelength [AA"); plt.ylabel("Flux [Jy]");  plt.show();   plt.clf()
        except: None

        sed_wl*=10
        
        
        filterdir=os.environ['WD_BASS_INSTALL_DIR']+"/Filters"
        filter_dict = {"GAIA/GAIA3:Grp":np.array([nploadtxt(filterdir+"/GAIA_GAIA3.Grp.dat",unpack=True), 7619.96, "PHOTON"], dtype=object),
            "GAIA/GAIA3:Gbp":np.array([nploadtxt(filterdir+"/GAIA_GAIA3.Gbp.dat",unpack=True), 5035.75, "PHOTON"], dtype=object), 
            "GAIA/GAIA3:G": np.array([nploadtxt(filterdir+"/GAIA_GAIA3.G.dat",unpack=True), 5822.39, "PHOTON"], dtype=object),
            "SDSS:z": np.array([nploadtxt(filterdir+"/SLOAN_SDSS.z.dat",unpack=True), 8922.78, "PHOTON"], dtype=object),
            "SDSS:i": np.array([nploadtxt(filterdir+"/SLOAN_SDSS.i.dat",unpack=True), 7457.89, "PHOTON"], dtype=object),
            "SDSS:r": np.array([nploadtxt(filterdir+"/SLOAN_SDSS.r.dat",unpack=True), 6141.12, "PHOTON"], dtype=object),
            "SDSS:g": np.array([nploadtxt(filterdir+"/SLOAN_SDSS.g.dat",unpack=True), 4671.78, "PHOTON"], dtype=object),
            "SDSS:u": np.array([nploadtxt(filterdir+"/SLOAN_SDSS.u.dat",unpack=True), 3608.04, "PHOTON"], dtype=object),
            "PAN-STARRS/PS1:y": np.array([nploadtxt(filterdir+"/PAN-STARRS_PS1.y.dat",unpack=True),9613.5, "PHOTON"], dtype=object),
            "PAN-STARRS/PS1:z": np.array([nploadtxt(filterdir+"/PAN-STARRS_PS1.z.dat",unpack=True),8668.5, "PHOTON"], dtype=object),
            "PAN-STARRS/PS1:i": np.array([nploadtxt(filterdir+"/PAN-STARRS_PS1.i.dat",unpack=True), 7503.03, "PHOTON"], dtype=object),
            "PAN-STARRS/PS1:r": np.array([nploadtxt(filterdir+"/PAN-STARRS_PS1.r.dat",unpack=True), 6156.4, "PHOTON"], dtype=object),
            "PAN-STARRS/PS1:g": np.array([nploadtxt(filterdir+"/PAN-STARRS_PS1.g.dat",unpack=True),4810.9, "PHOTON"], dtype=object),
            "Johnson:V": np.array([nploadtxt(filterdir+"/Generic_Johnson.V.dat",unpack=True),5466.1, "ENERGY"], dtype=object),
            #"BarrosU":np.array([nploadtxt(filterdir+"/WHT_ULTRACAM.u.dat",unpack=True), 3481.95, "PHOTON"], dtype=object),
            #"BarrosG":np.array([nploadtxt(filterdir+"/WHT_ULTRACAM.g.dat",unpack=True), 4762.3, "PHOTON"], dtype=object),
            #"BarrosR":np.array([nploadtxt(filterdir+"/WHT_ULTRACAM.r.dat",unpack=True), 6256.2, "PHOTON"], dtype=object),
            #"BarrosI":np.array([nploadtxt(filterdir+"/WHT_ULTRACAM.i.dat",unpack=True), 7585.9, "PHOTON"], dtype=object),
            "HST_F775W":np.array([nploadtxt(filterdir+"/HST_WFC3_UVIS1.F775W.dat",unpack=True), 7612.8, "PHOTON"], dtype=object),
            "HST_F390W":np.array([nploadtxt(filterdir+"/HST_WFC3_UVIS1.F390W.dat",unpack=True), 4022.2, "PHOTON"], dtype=object),
            "HST_F225W":np.array([nploadtxt(filterdir+"/HST_WFC3_UVIS1.F225W.dat",unpack=True), 2372.8, "PHOTON"], dtype=object),
            "XMM-OT:UVW1":np.array([nploadtxt(filterdir+"/XMM_OM.UVW1_filter.dat",unpack=True), 2971.0, "PHOTON"], dtype=object),
            "SWIFT:U":np.array([nploadtxt(filterdir+"/Swift_UVOT.U_trn.dat",unpack=True), 3520.95, "PHOTON"], dtype=object),
            "GALEX:NUV":np.array([nploadtxt(filterdir+"/GALEX_GALEX.NUV.dat",unpack=True), 2303.37, "PHOTON"], dtype=object),
            "GALEX:FUV":np.array([nploadtxt(filterdir+"/GALEX_GALEX.FUV.dat",unpack=True), 1548.85, "PHOTON"], dtype=object),
            "GAIA/GAIA2:Grp":np.array([nploadtxt(filterdir+"/GAIA_GAIA2.Grp.dat",unpack=True), 7592.04, "ENERGY"], dtype=object),
            "GAIA/GAIA2:Gbp":np.array([nploadtxt(filterdir+"/GAIA_GAIA2.Gbp.dat",unpack=True), 5014.05, "ENERGY"], dtype=object), 
            "GAIA/GAIA2:G": np.array([nploadtxt(filterdir+"/GAIA_GAIA2.G.dat",unpack=True), 5838.81, "ENERGY"], dtype=object), 
            "XMM-OT:UVM2": np.array([nploadtxt(filterdir+"/XMM_OM.UVM2_filter.dat",unpack=True), 2329.54, "PHOTON"], dtype=object),
            "XMM-OT:UVW2": np.array([nploadtxt(filterdir+"/XMM_OM.UVW2_filter.dat",unpack=True), 2143.84, "PHOTON"], dtype=object),
            "UVOT:UVW2": np.array([nploadtxt(filterdir+"/XMM_OM.UVW2_filter.dat",unpack=True), 2143.84, "PHOTON"], dtype=object),
            "UVOT:UVW1": np.array([nploadtxt(filterdir+"/XMM_OM.UVW1_filter.dat",unpack=True), 2971.72, "PHOTON"], dtype=object),
            "XMM-OT:U": np.array([nploadtxt(filterdir+"/XMM_OM.U_filter.dat",unpack=True), 3515.73, "PHOTON"], dtype=object),
            "XMM-OT:V": np.array([nploadtxt(filterdir+"/XMM_OM.V_filter.dat",unpack=True), 5437.81, "PHOTON"], dtype=object),
            "Johnson:B": np.array([nploadtxt(filterdir+"/Generic_Johnson.B.dat",unpack=True), 4369.53, "ENERGY"], dtype=object),
            "WISE:W1": np.array([nploadtxt(filterdir+"/WISE_WISE.W1.dat",unpack=True), 33154.27, "ENERGY"], dtype=object),
            "WISE:W2": np.array([nploadtxt(filterdir+"/WISE_WISE.W2.dat",unpack=True), 45644.77, "ENERGY"], dtype=object),
            "WISE:W3": np.array([nploadtxt(filterdir+"/WISE_WISE.W3.dat",unpack=True), 107866.13, "ENERGY"], dtype=object),
            "2MASS:J": np.array([nploadtxt(filterdir+"/2MASS_2MASS.J.dat",unpack=True), 12285.64, "ENERGY"], dtype=object),
            "2MASS:H": np.array([nploadtxt(filterdir+"/2MASS_2MASS.H.dat",unpack=True), 16385.40, "ENERGY"], dtype=object),
            "2MASS:Ks": np.array([nploadtxt(filterdir+"/2MASS_2MASS.Ks.dat",unpack=True), 21521.61, "ENERGY"], dtype=object),
            "OAJ_JPLUS.gSDSS": np.array([nploadtxt(filterdir+"/OAJ_JPLUS.gSDSS.dat",unpack=True), 4748.47, "PHOTON"], dtype=object),
            "OAJ_JPLUS.J0410": np.array([nploadtxt(filterdir+"/OAJ_JPLUS.J0410.dat",unpack=True), 4107.98, "PHOTON"], dtype=object),
            "OAJ_JPLUS.J0861": np.array([nploadtxt(filterdir+"/OAJ_JPLUS.J0861.dat",unpack=True), 8610.16, "PHOTON"], dtype=object),
            "OAJ_JPLUS.iSDSS": np.array([nploadtxt(filterdir+"/OAJ_JPLUS.iSDSS.dat",unpack=True), 7613.86, "PHOTON"], dtype=object),
            "OAJ_JPLUS.J0430": np.array([nploadtxt(filterdir+"/OAJ_JPLUS.J0430.dat",unpack=True), 4298.36, "PHOTON"], dtype=object),
            "OAJ_JPLUS.rSDSS": np.array([nploadtxt(filterdir+"/OAJ_JPLUS.rSDSS.dat",unpack=True), 6206.11, "PHOTON"], dtype=object),
            "OAJ_JPLUS.J0378": np.array([nploadtxt(filterdir+"/OAJ_JPLUS.J0378.dat",unpack=True), 3793.38, "PHOTON"], dtype=object),
            "OAJ_JPLUS.J0515": np.array([nploadtxt(filterdir+"/OAJ_JPLUS.J0515.dat",unpack=True), 5139.67, "PHOTON"], dtype=object),
            "OAJ_JPLUS.uJAVA": np.array([nploadtxt(filterdir+"/OAJ_JPLUS.uJAVA.dat",unpack=True), 3542.20, "PHOTON"], dtype=object),
            "OAJ_JPLUS.J0395": np.array([nploadtxt(filterdir+"/OAJ_JPLUS.J0395.dat",unpack=True), 3938.55, "PHOTON"], dtype=object),
            "OAJ_JPLUS.J0660": np.array([nploadtxt(filterdir+"/OAJ_JPLUS.J0660.dat",unpack=True), 6606.67, "PHOTON"], dtype=object),
            "OAJ_JPLUS.zSDSS": np.array([nploadtxt(filterdir+"/OAJ_JPLUS.zSDSS.dat",unpack=True), 8940.28, "PHOTON"], dtype=object)

            }
        
        
        try:
            
            External_AB_mag = np.asarray(config_info["External_AB_mag"]).astype(float)
            External_AB_mag_err = np.asarray(config_info["External_AB_mag_err"]).astype(float)
            External_Flux_Jy=10**(23-(External_AB_mag+48.594)/2.5)  #  erg/cm2/s/Hz
            External_Flux_Jy_err=npsqrt(  ((-(0.9210340371976184*10**((External_AB_mag+243/5)*-0.4+23)))*External_AB_mag_err)**2   )
            External_Filter = np.asarray(config_info["External_Filter"])
            External_wl = []
            for an_external_filt in External_Filter:
                try:
                    External_wl.append(filter_dict[str(an_external_filt)][1])
                except Exception as e: 
                    print_exception_details(e)
                    print("External filter with name", an_external_filt, "not recognised. I am removing this filter from the list")
                    mask = External_Filter!=an_external_filt
                    External_Filter, External_wl, External_Flux_Jy, External_Flux_Jy_err = External_Filter[mask], External_wl[mask], External_Flux_Jy[mask], External_Flux_Jy_err[mask]
            sedfilter, sed_wl, sedflux, sedfluxe  =  np.concatenate((sedfilter, External_Filter)), np.concatenate((sed_wl, External_wl)), np.concatenate((sedflux, External_Flux_Jy)), np.concatenate((sedfluxe, External_Flux_Jy_err))
        except: None
        
        
        theminww, themaxww = 999999, -999999
        for afilt in sedfilter:
            ww, flfl  =  filter_dict[afilt][0]
            mask=flfl>=0
            ww, flfl = ww[mask], flfl[mask]
            minww, maxww = npamin(ww), npamax(ww)
            
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
    else:
        theminww_loadgrid, themaxww_loadgrid=1000,15000
        fit_phot_SED=False
        print("ALERT:")
        print("No available photometry of the source in this RA and Dec in a 5'' radius")
    


#theminww, themaxww = 3300, 11400
#theminww_loadgrid, themaxww_loadgrid = 3000, 11400







## Now loading the atmospheric grids if necessary
if starType1.startswith("D") or starType2.startswith("D"):

    pier_or_antoine="mixed"#"pier3Dphot_antoine1Dspec"

    if (starType1=="DA" or starType2=="DA") and pier_or_antoine=="mixed":
        if fit_phot_SED:         wl_all_synth, flux_all_synth, Teff_all_synth, Grav_all_synth = load_models.load_models_DA_3D_NLTE(minwl=theminww_loadgrid,maxwl=themaxww_loadgrid)
        else:                    wl_all_synth, flux_all_synth, Teff_all_synth, Grav_all_synth = load_models.load_models_DA_3D_NLTE(minwl=npamin(reference_wl)-150,maxwl=npamax(reference_wl)+250);   mask_wl_5100_6200 = ((wl_all_synth<=6200) & (wl_all_synth>=5100));   wl_all_synth, flux_all_synth, Teff_all_synth, Grav_all_synth = wl_all_synth[~mask_wl_5100_6200], flux_all_synth[~mask_wl_5100_6200], Teff_all_synth[~mask_wl_5100_6200], Grav_all_synth[~mask_wl_5100_6200]
        
        
        if False:
            if forced_teff1==0 and forced_teff2==0 and forced_logg1==0 and forced_logg2==0:
                # Trim the grid to only include the values within p0Teff and p0logg
                min_p0T = npamin(np.array([p0T1[0], p0T1[1], p0T2[0], p0T2[1]]));    max_p0T = npamax(np.array([p0T1[0], p0T1[1], p0T2[0], p0T2[1]]))
                
                mask_TEFF = (npunique(Teff_all_synth) >= min_p0T)  &  (npunique(Teff_all_synth) <= max_p0T)
                mask_TEFF[npamin(npargwhere(mask_TEFF==True)) - 1] = True;      mask_TEFF[npamax(npargwhere(mask_TEFF==True)) + 1] = True
                
                min_p0logg = npamin(np.array([p0logg1[0], p0logg1[1], p0logg2[0], p0logg2[1]]));    max_p0logg = npamax(np.array([p0logg1[0], p0logg1[1], p0logg2[0], p0logg2[1]]))
                
                mask_LOGG = (npunique(Grav_all_synth) >= min_p0logg)  &  (npunique(Grav_all_synth) <= max_p0logg)
                mask_LOGG[npamin(npargwhere(mask_LOGG==True)) - 1] = True;    mask_LOGG[npamax(npargwhere(mask_LOGG==True)) + 1] = True
                
                
                newMinT, newMaxT = npamin(npunique(Teff_all_synth)[mask_TEFF]), npamax(npunique(Teff_all_synth)[mask_TEFF])
                newMinLogg, newMaxLogg = npamin(npunique(Grav_all_synth)[mask_LOGG]), npamax(npunique(Grav_all_synth)[mask_LOGG])
                
                
                mask_atm_grid = (Teff_all_synth <= newMaxT) & (Teff_all_synth >= newMinT)  &  (Grav_all_synth <= newMaxLogg) & (Grav_all_synth >= newMinLogg)
                
                wl_all_synth, flux_all_synth, Teff_all_synth, Grav_all_synth = wl_all_synth[mask_atm_grid], flux_all_synth[mask_atm_grid], Teff_all_synth[mask_atm_grid], Grav_all_synth[mask_atm_grid]
        
        logg_all_synth=Grav_all_synth
        unique_Teffs_synth  =  npunique(Teff_all_synth)
        
        
        
        
    if starType1=="DBA" or starType2=="DBA" or starType1=="DB" or starType2=="DB":
        # Load in the 3D DB/DBA/DC LTE spectra
        if fit_phot_SED:         wl_all, flux_all, Teff_all, Grav_all, H_over_He_all = load_models.load_models_DBA(minwl=theminww_loadgrid,maxwl=themaxww_loadgrid)
        else:                    wl_all, flux_all, Teff_all, Grav_all, H_over_He_all = load_models.load_models_DBA(minwl=npamin(reference_wl)-150,maxwl=npamax(reference_wl)+250)
        
        if False:
            if forced_teff1==0 and forced_teff2==0 and forced_logg1==0 and forced_logg2==0 and forced_HoverHe1==0 and forced_HoverHe2==0:
                # Trim the grid to only include the values within p0Teff and p0logg
                min_p0T = npamin(np.array([p0T1[0], p0T1[1], p0T2[0], p0T2[1]]));    max_p0T = npamax(np.array([p0T1[0], p0T1[1], p0T2[0], p0T2[1]]))
                
                mask_TEFF = (npunique(Teff_all) >= min_p0T)  &  (npunique(Teff_all) <= max_p0T)
                try: mask_TEFF[npamin(npargwhere(mask_TEFF==True)) - 1] = True
                except: None
                try:  mask_TEFF[npamax(npargwhere(mask_TEFF==True)) + 1] = True
                except: None
                
                min_p0logg = npamin(np.array([p0logg1[0], p0logg1[1], p0logg2[0], p0logg2[1]]));    max_p0logg = npamax(np.array([p0logg1[0], p0logg1[1], p0logg2[0], p0logg2[1]]))
                
                mask_LOGG = (npunique(Grav_all) >= min_p0logg)  &  (npunique(Grav_all) <= max_p0logg)
                try: mask_LOGG[npamin(npargwhere(mask_LOGG==True)) - 1] = True
                except: None
                try: mask_LOGG[npamax(npargwhere(mask_LOGG==True)) + 1] = True
                except: None
                
                
                newMinT, newMaxT = npamin(npunique(Teff_all)[mask_TEFF]), npamax(npunique(Teff_all)[mask_TEFF])
                newMinLogg, newMaxLogg = npamin(npunique(Grav_all)[mask_LOGG]), npamax(npunique(Grav_all)[mask_LOGG])
                
                
                mask_atm_grid = (Teff_all <= newMaxT) & (Teff_all >= newMinT)  &  (Grav_all <= newMaxLogg) & (Grav_all >= newMinLogg)
                
                wl_all, flux_all, Teff_all, Grav_all, H_over_He_all = wl_all[mask_atm_grid], flux_all[mask_atm_grid], Teff_all[mask_atm_grid], Grav_all[mask_atm_grid], H_over_He_all[mask_atm_grid]
        
        
        
        logg_all=Grav_all

        #Teff_all=np.round(Teff_all,-3) # this is naughty, but decreases run time with a small loss of accuracy to normalised profiles because teff grids aren't exactly the same
        Teff_all_logg = Teff_all;   wl_all_logg = wl_all;     flux_all_logg = flux_all;      logg_all_logg = logg_all;     HoverHe_all_logg_DBA = H_over_He_all
        unique_HoverHe_DBA=npunique(HoverHe_all_logg_DBA)
        
        
        if starType1=="DB":
            maskDBonly=HoverHe_all_logg_DBA==30
            Teff_all_logg = Teff_all[maskDBonly];   wl_all_logg = wl_all[maskDBonly];     flux_all_logg = flux_all[maskDBonly];      logg_all_logg = logg_all[maskDBonly]
            del HoverHe_all_logg_DBA, H_over_He_all





    
    if starType1=="DC" or starType2=="DC":
        if fit_phot_SED:         wl_all, flux_all, Teff_all, Grav_all = load_models.load_models_DBA_1eMinus5(minwl=theminww_loadgrid,maxwl=themaxww_loadgrid)
        else:                    wl_all, flux_all, Teff_all, Grav_all = load_models.load_models_DBA_1eMinus5(minwl=npamin(reference_wl)-150,maxwl=npamax(reference_wl)+250)
        
        
        if False:
            if forced_teff1==0 and forced_teff2==0 and forced_logg1==0 and forced_logg2==0:
                # Trim the grid to only include the values within p0Teff and p0logg
                min_p0T = npamin(np.array([p0T1[0], p0T1[1], p0T2[0], p0T2[1]]));    max_p0T = npamax(np.array([p0T1[0], p0T1[1], p0T2[0], p0T2[1]]))
                
                mask_TEFF = (npunique(Teff_all) >= min_p0T)  &  (npunique(Teff_all) <= max_p0T)
                mask_TEFF[npamin(npargwhere(mask_TEFF==True)) - 1] = True;      mask_TEFF[npamax(npargwhere(mask_TEFF==True)) + 1] = True
                
                min_p0logg = npamin(np.array([p0logg1[0], p0logg1[1], p0logg2[0], p0logg2[1]]));    max_p0logg = npamax(np.array([p0logg1[0], p0logg1[1], p0logg2[0], p0logg2[1]]))
                
                mask_LOGG = (npunique(Grav_all) >= min_p0logg)  &  (npunique(Grav_all) <= max_p0logg)
                mask_LOGG[npamin(npargwhere(mask_LOGG==True)) - 1] = True;    mask_LOGG[npamax(npargwhere(mask_LOGG==True)) + 1] = True
                
                
                newMinT, newMaxT = npamin(npunique(Teff_all)[mask_TEFF]), npamax(npunique(Teff_all)[mask_TEFF])
                newMinLogg, newMaxLogg = npamin(npunique(Grav_all)[mask_LOGG]), npamax(npunique(Grav_all)[mask_LOGG])
                
                
                mask_atm_grid = (Teff_all <= newMaxT) & (Teff_all >= newMinT)  &  (Grav_all <= newMaxLogg) & (Grav_all >= newMinLogg)
                
                wl_all, flux_all, Teff_all, Grav_all = wl_all[mask_atm_grid], flux_all[mask_atm_grid], Teff_all[mask_atm_grid], Grav_all[mask_atm_grid]
            
        
        
        
        logg_all=np.round(Grav_all,1)
        Teff_all_logg = Teff_all;   wl_all_logg = wl_all;     flux_all_logg = flux_all;      logg_all_logg = logg_all
    





## Now we input all the observed data

list_norm_wl_grids, list_normalised_flux, list_normalised_err = [], [], []
mask_out_Ha_min_all, mask_out_Ha_max_all, cut_limits_min_all, cut_limits_max_all = [], [], [], []

min_wl_all_files, max_wl_all_files = 9E9, 0
spec_wl, spec_fl, spec_fle = [], [], []

# Normalise the input data
@njit
def gauss_if_no_flux_error(x, a, mu, sigma, b):
    return (a * np.exp(-np_square(x-mu) / sigma) + b )#   /  (x*mmm+ccc)
    #return ((a/sigma/npsqrt(2*np_pi)) * np.exp(-0.5 * np_square((x-mu)/sigma)) + b)  /  (x*mmm+ccc)

try:    stack_spectra=np.asarray(config_info["stack_spectra"]).astype(bool)
except: stack_spectra=False

for files, normaliseHa, cut_Ha, ref_wl in zip(input_files, normaliseHa_all, cut_Ha_all, reference_wl):
    if spectra_source_type=="wl_flux_fluxerr":
        try:  wl_data, flux_data, flux_e_data = nploadtxt(os.getcwd()+"/"+files, skiprows=1, unpack=True)  #  if you have flux errors, files are loaded in normally. Great.
        except:
            ##### this part of the code is for when you do not have flux errors. I do my best to guess the SNR of the data given an expected continuum SNR under the assumption that your spectral line follows a gaussian profile. In the "gauss" function below, mmm and ccc are included to fit the general trend of the continuum as a linear fit. 
            
            ##### I used this method to fit to FIES observations as errors are not output in its reduction pipeline
            wl_data, flux_data = nploadtxt(os.getcwd()+"/"+files, skiprows=2, unpack=True)
            if True:
                
                cut_limits_min, cut_limits_max = ref_wl+cut_Ha[0],  ref_wl+cut_Ha[1]
                norm_limits_min, norm_limits_max = ref_wl+normaliseHa[0],  ref_wl+normaliseHa[1]
                
                
                
                mask_cut = (wl_data>=cut_limits_min) & (wl_data<=cut_limits_max)
                mask_norm = ((wl_data>=cut_limits_min) & (wl_data<=norm_limits_min)) | ((wl_data<=cut_limits_max) & (wl_data>=norm_limits_max))
                mask_linecore = ((wl_data<ref_wl-10) & (wl_data>ref_wl-30)) | ((wl_data>ref_wl+10) & (wl_data<ref_wl+30))
                
                
                m6, m5, m4, m3, m2, m, c = polyfit(wl_data[mask_cut] - ref_wl, flux_data[mask_cut], deg=6)
                
                xs = wl_data - ref_wl
                #plt.plot(wl_data[mask_cut], flux_data[mask_cut]);  plt.plot(wl_data[mask_cut], m6*xs**6 + m5*xs**5 + m4*xs**4 + m3*xs**3 + m2*xs**2 + m*xs + c);  plt.show()
                
                
                med_flux = np.median(flux_data[mask_norm])   # this gets the median flux of the normalised region
                
                
                predicted_SNR_of_normalise_region_per_pixel = 25  # this is the predicted flux of the normalised part of the spectrum.  I use this to predict what the SNR is for all other pixels below
                
                flux_e_data = med_flux /predicted_SNR_of_normalise_region_per_pixel   *  npsqrt(med_flux/(m6*xs**6 + m5*xs**5 + m4*xs**4 + m3*xs**3 + m2*xs**2 + m*xs + c))  # propogate SNR to all parts of the spectrum. I assume that there is no detector readout noise and so your SNR only depends on the root(#photons)
                
                #plt.plot(wl_data[mask_cut], flux_e_data[mask_cut]);   plt.plot(wl_data[mask_cut], flux_data[mask_cut],c='k');   plt.show()
                
                
                
                flux_e_data[flux_e_data<0] = npamax(flux_e_data[mask_cut][flux_e_data[mask_cut]>0]) #  if negative flux, the error on these measurements is taken as the maximum positive error. Bit arbitrary... better to have actual flux errors!! I do this in case there is a big negative outlier (e.g. poor sky subtraction)
            else:
                flux_e_data = flux_data/10
                
                #maybe instead of a gaussian it's better to just fit an e.g. 5order polynomial to everything?
            
            
            
        args = np.argsort(wl_data)
        wl_data=wl_data[args]
        flux_data=flux_data[args]
        flux_e_data=flux_e_data[args]
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
        flux_data=data["flux"][0][mask]
        flux_e_data=data["err"][0][mask]
        
        flux_data=data["flux"][0][mask]   # erg cm-2 s-1 AA-1
        flux_e_data=data["err"][0][mask]  # erg cm-2 s-1 AA-1
        
        flux_data = 3.33564095E+04 * flux_data * np_square(wl_data) / 1E23   # convert to  erg/cm^2/s/Hz
        flux_e_data = 3.33564095E+04 * flux_e_data * np_square(wl_data) / 1E23 # convert to erg/cm^2/s/Hz
        
    plt.clf()
    
    
    # I deredden the observed spectra (spectra, not flux calibrated photometry - this is modelled with a reddened synthetic spectrum) so that the synthetic spectra do not need to be redenned every iteration. This may become an issue if you have very low resolution (R<100) data, but otherwise it's not something to worry about
    
    ext = G23(Rv=3.1)
    flux_data /= ext.extinguish(wl_data*u.AA, Ebv=reddening_Ebv)
    flux_e_data /= ext.extinguish(wl_data*u.AA, Ebv=reddening_Ebv)
    
        
    #plt.plot(wl_data,flux_data);  plt.title(files);  plt.show();  plt.clf()
    
    minwlhere=npamin(wl_data);        maxwlhere=npamax(wl_data)
    
    if not (ref_wl < maxwlhere and ref_wl > minwlhere) and minwlhere<ref_wl-100: raise ValueError("File named: " + str(files) + "  does not have a valid reference wl (" + str(ref_wl) + ", while the minimum wl in the file was "+str(minwlhere)+")")
    if minwlhere<min_wl_all_files:     min_wl_all_files = minwlhere -50
    if maxwlhere>max_wl_all_files:     max_wl_all_files = maxwlhere +50
    
    
    
    if stack_spectra==False:
        cut_limits_min, cut_limits_max = ref_wl+cut_Ha[0],  ref_wl+cut_Ha[1]
        cut_limits_min_all.append(cut_limits_min);    cut_limits_max_all.append(cut_limits_max)
        
        
        mask_out_Ha_min, mask_out_Ha_max = ref_wl+normaliseHa[0],  ref_wl+normaliseHa[1]
        mask_out_Ha_min_all.append(mask_out_Ha_min);    mask_out_Ha_max_all.append(mask_out_Ha_max)
        
        if high_RV_amp:   mask_out = (wl_data>cut_limits_min-20) & (wl_data<cut_limits_max+20) #  if I want to allow high RV amps where the cuts are dynamic, extend the range that is snipped to allow for it. DO NOT INCLUDE THIS AS A NEW LIMIT THOUGH! (line above) I am not changing the yaml file limits globally. 20A is arbitrary
        else:             mask_out = (wl_data>cut_limits_min) & (wl_data<cut_limits_max)
        wl_data=wl_data[mask_out];   flux_data = flux_data[mask_out];   flux_e_data = flux_e_data[mask_out]
        
        mask_norm = (wl_data>cut_limits_min) & (wl_data<cut_limits_max) & ((wl_data<mask_out_Ha_min) | (wl_data>mask_out_Ha_max))
        
        
        wl_data_normmask = wl_data[mask_norm]; flux_data_normmask = flux_data[mask_norm];   fluxe_data_normmask = flux_e_data[mask_norm]
        
        
        
        
        good=flux_data_normmask==flux_data_normmask
        num_its = 0
        try:
            while num_its<100:
                current_length = len(flux_data_normmask[good])
                m,c = polyfit(wl_data_normmask[good], flux_data_normmask[good], w=1/fluxe_data_normmask[good], deg=1)
                
                flux_norm_sigma_clip = flux_data_normmask - (m*wl_data_normmask+c)
                
                std = np.std(flux_data_normmask[good])
                good = good & ((np.abs(flux_norm_sigma_clip)  < 1.5*np.abs(std)))
                
                if len(flux_data_normmask[good]) == current_length:
                    break
        except:
            print(files)
            print("The input wavelengths that I am selecting to normaliise the spectra are:")
            print(mask_out_Ha_min, mask_out_Ha_max, cut_limits_min, cut_limits_max)
            m,c = polyfit(wl_data_normmask, flux_data_normmask, w=1/fluxe_data_normmask, deg=1)
            print("DIDN'T LIKE THE NORMALISTATION OF THE OBSERVED SPECTRA, NO CLIPPING ON NORMALISATION APPLIED")
            
        #plt.plot(wl_data[mask_norm], flux_data[mask_norm]);   plt.show();   plt.clf()
        normalised_wl_grid = wl_data;    normalised_flux = flux_data/(m*wl_data+c);    normalised_err = flux_e_data/(m*wl_data+c)
        list_norm_wl_grids.append(normalised_wl_grid);     list_normalised_flux.append(normalised_flux);     list_normalised_err.append(normalised_err)
    else:
        spec_wl.append(wl_data);  spec_fl.append(flux_data);  spec_fle.append(flux_e_data)  #  store the values pre normalisation, allow for stacking



if "run_double_one_by_oneMCMC" in sys.argv and (starType1=="quadLorentz" and starType2=="quadLorentz"):
    argmax = np.argwhere(reference_wl==np.amax(reference_wl))[0][0]
    input_files, share_rv, modelHa, cut_Ha_all, normaliseHa_all, HJD_values, sigma_clip, list_norm_wl_grids, list_normalised_flux, list_normalised_err,  resolutions, RV_boundaries1, RV_boundaries2, reference_wl = input_files[argmax], share_rv[argmax], modelHa[argmax], cut_Ha_all[argmax], normaliseHa_all[argmax], HJD_values[argmax], sigma_clip[argmax], list_norm_wl_grids[argmax], list_normalised_flux[argmax], list_normalised_err[argmax],  resolutions[argmax], RV_boundaries1[argmax], RV_boundaries2[argmax], reference_wl[argmax]
    
    input_files, share_rv, modelHa, cut_Ha_all, normaliseHa_all, HJD_values, sigma_clip, list_norm_wl_grids, list_normalised_flux, list_normalised_err,  resolutions, RV_boundaries1, RV_boundaries2, reference_wl = np.array([input_files]), np.array([share_rv]), np.array([modelHa]), np.array([cut_Ha_all]), np.array([normaliseHa_all]), np.array([HJD_values]), np.array([sigma_clip]), np.array([list_norm_wl_grids]), np.array([list_normalised_flux]), np.array([list_normalised_err]), np.array([resolutions]), np.array([RV_boundaries1]), np.array([RV_boundaries2]), np.array([reference_wl])
    




try:    stack_spectra=np.asarray(config_info["stack_spectra"]).astype(bool)
except: stack_spectra=False

if stack_spectra:
    ## here I am going to stack spectra in X phase bins if desired. It will be useful for WDJ022558 and maybe others
    #want_stacking_bins, want_stacking_period_T0, want_stacking_period, want_stacking_pdot = np.asarray(config_info["want_galex"])  # period in seconds, pdot in s\,s-1
    
    stack_spectra_bins=np.asarray(config_info["stack_spectra_bins"]).astype(int)
    try:  stack_spectra_T0=np.asarray(config_info["stack_spectra_T0"]).astype(float)
    except: print("No stack_spectra_T0 included in the input file. I am assuming 0");  stack_spectra_T0=0
    try: stack_spectra_P0=np.asarray(config_info["stack_spectra_P0"]).astype(float)
    except: raise ValueError("If stacking, you need to say a period to fold on - stack_spectra_P0. The HJD_values time unit needs to be days. You may also supply a pdot - stack_spectra_Pdot. If you're doing this, I assume that seconds per second (s s-1) are the units of Pdot even though P0 is in days.")
    try:     stack_spectra_Pdot=np.asarray(config_info["stack_spectra_Pdot"]).astype(float);  stack_spectra_P0seconds=stack_spectra_P0*86400
    except:  stack_spectra_Pdot=0
    

    arr = np.linspace(0,1,stack_spectra_bins+1);    diff=np.diff(arr);    midpoints=diff/2
    new_share_rv, list_cut, list_norm, list_mod, list_res, list_refwl, list_input_files, list_HJD_values, list_sig, list_RVbounds1, list_RVbounds2 = [], [], [], [], [], [], [], [], [], [], []
    for cn, val in enumerate(share_rv):
        if val==-1: new_share_rv.append(cn)
        else: new_share_rv.append(val)
    new_share_rv = np.asarray(new_share_rv)


    if stack_spectra_Pdot!=0:
        E=(HJD_values-stack_spectra_T0)/stack_spectra_P0
        OminusC = 0.5*stack_spectra_P0seconds*stack_spectra_Pdot*np_square(E)
        HJD_values-=OminusC/86400
    
    
    #for xx,yy,zz, aa, bb in zip(input_files, HJD_values, ((HJD_values-stack_spectra_T0)%stack_spectra_P0)/stack_spectra_P0, RV_boundaries1, RV_boundaries2):
    #	print(xx,yy,zz,aa,bb)
    #raise ValueError
    
    binned_details=[]
    for midpoint in np.round(arr[:-1]+midpoints,15):  #  i round it to 15dp here because without it the floating point precision goes funny
        process_shareRV_val=[]
        for s, hjd, fname in zip(new_share_rv, HJD_values-stack_spectra_T0, input_files):
            if (hjd%stack_spectra_P0)/stack_spectra_P0 >= midpoint-midpoints[0] and (hjd%stack_spectra_P0)/stack_spectra_P0 < midpoint+midpoints[0]:
                process_shareRV_val.append(s)
                binned_details.append([midpoint, hjd+stack_spectra_T0, fname])
        
        bin_shareRV_val=npunique(process_shareRV_val)
        
        bin_wls, bin_flux, bin_fluxe = [], [], []
        if len(bin_shareRV_val)>1:
            mask=new_share_rv!=new_share_rv
            for vv in bin_shareRV_val:
                mask = mask | (new_share_rv==vv)
            
            for index in npargwhere(mask==True):
                index=index[0]
                bin_wls.append(spec_wl[index]);  bin_flux.append(spec_fl[index]);  bin_fluxe.append(spec_fle[index])
            
            
            process_modelHa, process_normaliseHa_all, process_cut_Ha_all, process_resolutions, process_reference_wl, process_RV_boundaries1, process_RV_boundaries2, process_sigma_clip = modelHa[mask], normaliseHa_all[mask], cut_Ha_all[mask], resolutions[mask], reference_wl[mask], RV_boundaries1[mask], RV_boundaries2[mask], sigma_clip[mask]
            
            for un_refwl in npunique(process_reference_wl):
                un_resolution = npunique(process_resolutions[process_reference_wl==un_refwl])
                if len(un_resolution) > 1: raise ValueError("Trying to stack spectra that have different spectral resolutions for equal refwl")
                else: un_resolution = un_resolution[0]
                all_flux, all_fluxe = [], []
                min_cut, max_cut, min_norm, max_norm, min_model, max_model, min_asig, min_RVbound1, max_RVbound1, min_RVbound2, max_RVbound2 = 9999, -9999, -9999, 9999, -9999, 9999, 9999, -9999, 9999, -9999, 9999
                first_time=True
                for cn, (a_refwl, a_mod, a_norm, a_cut, a_sig, a_RVbound1, a_RVbound2) in enumerate(zip(process_reference_wl, process_modelHa, process_normaliseHa_all, process_cut_Ha_all, process_sigma_clip, process_RV_boundaries1, process_RV_boundaries2)):
                    if un_refwl==a_refwl:
                        wl, fl, fle = bin_wls[cn], bin_flux[cn], bin_fluxe[cn]
                        
                        mask_wave = ((wl-un_refwl)>a_cut[0]) & ((wl-un_refwl)<a_cut[1])
                        
                        wl, fl, fle = wl[mask_wave], fl[mask_wave], fle[mask_wave]
                        
                        if len(all_flux)==0: wl_axis = wl
                        else:                     fl = interp(wl_axis, wl, fl,left=9999,right=9999);  fle = interp(wl_axis, wl, fle,left=9999,right=9999)
                        
                        
                        def apply_vertical_shift(x,a):
                            return stored_flux[mask]/a
                        
                        if first_time:   # store first spectrum, use it to scale all other spectra
                            stored_flux, stored_fluxe = fl, fle;     first_time=False
                        else:
                            mask = (stored_flux != 9999) & (fl != 9999)
                            popt, _ = curve_fit(apply_vertical_shift, wl_axis[mask], fl[mask],sigma=npsqrt(np_square(fle[mask]) + np_square(stored_fluxe[mask])))
                            fl*=popt;  fle*=popt


                        all_flux.append(fl);  all_fluxe.append(fle)
                        
                        if a_norm[0] > min_norm:  min_norm = a_norm[0]
                        if a_norm[1] < max_norm:  max_norm = a_norm[1]
                        if a_cut[0]  < min_cut:   min_cut  = a_cut[0]
                        if a_cut[1]  > max_cut:   max_cut  = a_cut[1]
                        if a_mod[0]  > min_model: min_model= a_mod[0]
                        if a_mod[1]  < max_model: max_model= a_mod[1]
                        if a_sig < min_asig:            min_asig = a_sig
                        if a_RVbound1[0] > min_RVbound1:  min_RVbound1 = a_RVbound1[0]
                        if a_RVbound1[1] < max_RVbound1:  max_RVbound1 = a_RVbound1[1]
                        if a_RVbound2[0] > min_RVbound2:  min_RVbound2 = a_RVbound2[0]
                        if a_RVbound2[1] < max_RVbound2:  max_RVbound2 = a_RVbound2[1]
                        
                        
                        
                all_flux, all_fluxe = np.asarray(all_flux), np.asarray(all_fluxe)
                
                stacked_spec_flux, stacked_spec_err = [], []
                for a,b in zip(all_flux.T, all_fluxe.T):
                    mask=(a!=9999) & (b!=9999) & (b>0)
                    a, b = a[mask], b[mask]
                    if len(a)>1:  stacked_spec_flux.append(np.average(a, weights=1/b));     stacked_spec_err.append(npsqrt(1/npsum(1/np_square(b))))
                    else:         stacked_spec_flux.append(a[0]);     stacked_spec_err.append(b[0])
                    
                
                
                stacked_spec_flux, stacked_spec_err = np.asarray(stacked_spec_flux), np.asarray(stacked_spec_err)
                
                
                resid=np.abs(stacked_spec_flux-np.median(stacked_spec_flux))
                mask_again = resid/np.std(resid) < 10 # not sure why, but sometimes the np.average(a, weights=b) comes out to be a stupidly high number. I'll clip it out in case of a spurious point
                wl_axis, stacked_spec_flux, stacked_spec_err = wl_axis[mask_again], stacked_spec_flux[mask_again], stacked_spec_err[mask_again]
                
                
                mask_norm = (wl_axis>un_refwl+min_cut) & (wl_axis<un_refwl+max_cut) & ((wl_axis<un_refwl+min_norm) | (wl_axis>un_refwl+max_norm))
                wl_data_normmask = wl_axis[mask_norm]; flux_data_normmask = stacked_spec_flux[mask_norm];   fluxe_data_normmask = stacked_spec_err[mask_norm]
                
                
                group_left_of_centre = wl_data_normmask<un_refwl;                           group_right_of_centre = wl_data_normmask>un_refwl
                std_left, std_right = np.std(flux_data_normmask[group_left_of_centre]),     np.std(flux_data_normmask[group_right_of_centre])
                med_left, med_right = np.median(flux_data_normmask[group_left_of_centre]),  np.median(flux_data_normmask[group_right_of_centre])
                
                
                
                mask = ((flux_data_normmask > med_left-2.5*std_left) & (flux_data_normmask < med_left+2.5*std_left)) | ((flux_data_normmask > med_right-2.5*std_right) & (flux_data_normmask < med_right+2.5*std_right))
                wl_data_normmask = wl_data_normmask[mask]; flux_data_normmask = flux_data_normmask[mask];   fluxe_data_normmask = fluxe_data_normmask[mask]
                
            
                good=flux_data_normmask==flux_data_normmask
                while True:
                    current_length = len(flux_data_normmask[good])
                    m,c = polyfit(wl_data_normmask[good], flux_data_normmask[good], w=1/fluxe_data_normmask[good], deg=1)
                    
                    flux_norm_sigma_clip = flux_data_normmask - (m*wl_data_normmask+c)
                    
                    std = np.std(flux_data_normmask[good])
                    good = good & ((np.abs(flux_norm_sigma_clip)  < 2*np.abs(std)))
                    
                    if len(flux_data_normmask[good]) == current_length:    break
                
                normalised_wl_grid = wl_axis;    normalised_flux = stacked_spec_flux/(m*wl_axis+c);    normalised_err = stacked_spec_err/(m*wl_axis+c)
                
                
                mask = normalised_flux< 10*np.median(normalised_flux)
                normalised_wl_grid, normalised_flux, normalised_err = normalised_wl_grid[mask], normalised_flux[mask], normalised_err[mask]
                
                    
                list_norm_wl_grids.append(normalised_wl_grid);     list_normalised_flux.append(normalised_flux);     list_normalised_err.append(normalised_err);   list_HJD_values.append(midpoint)
                list_cut.append([min_cut, max_cut]);  list_norm.append([min_norm, max_norm]);   list_mod.append([min_model, max_model]);  list_res.append(un_resolution);   list_refwl.append(un_refwl);  list_input_files.append("Stack"+str(midpoint).replace(".","p"));   list_sig.append(min_asig);   list_RVbounds1.append([min_RVbound1, max_RVbound1]);   list_RVbounds2.append([min_RVbound2, max_RVbound2])
                


    list_HJD_values, list_cut, list_norm, list_mod, list_res, list_refwl, list_input_files, list_sig, list_RVbounds1, list_RVbounds2 = np.asarray(list_HJD_values), np.asarray(list_cut), np.asarray(list_norm), np.asarray(list_mod), np.asarray(list_res), np.asarray(list_refwl), np.asarray(list_input_files), np.asarray(list_sig), np.asarray(list_RVbounds1), np.asarray(list_RVbounds2)


    # get a list of indices for sorting, then sort indices first by refwl then by the bin number
    indices = list(range(len(list_refwl)))
    args = sorted(indices, key=lambda i: (-list_refwl[i], list_HJD_values[i]))

    list_norm_wl_grids = [list_norm_wl_grids[i] for i in args];    list_normalised_flux = [list_normalised_flux[i] for i in args];    list_normalised_err = [list_normalised_err[i] for i in args]

    list_HJD_values, list_cut, list_norm, list_mod, list_res, list_refwl, list_input_files, list_sig, list_RVbounds1, list_RVbounds2 = list_HJD_values[args], list_cut[args], list_norm[args], list_mod[args], list_res[args], list_refwl[args], list_input_files[args], list_sig[args], list_RVbounds1[args], list_RVbounds2[args]


    max_wl=npamax(list_refwl)
    done_midpoint, final_shareRV = [], []
    for aHjd, aRefwl in zip(list_HJD_values, list_refwl):
        if aRefwl==max_wl:    final_shareRV.append(-1);   done_midpoint.append(aHjd)
        else:                 final_shareRV.append(npargwhere(aHjd==done_midpoint)[0][0])

    final_shareRV=np.asarray(final_shareRV)


    #for c, d,e,f,g,h,i,j,k,l in zip(list_HJD_values, list_mod, list_norm, list_cut, list_res, list_refwl, list_input_files, list_sig, list_RVbounds, final_shareRV): # list_norm_wl_grids, list_normalised_flux, list_normalised_err
    #    print(c,d,e,f,g,h,i,j,k,l)


    modelHa,  normaliseHa_all,  cut_Ha_all = list_mod,  list_norm,  list_cut
    resolutions=list_res;    reference_wl=list_refwl;    share_rv=final_shareRV;    RV_boundaries1=list_RVbounds1;    RV_boundaries2=list_RVbounds2;    HJD_values=list_HJD_values;    input_files=list_input_files;     sigma_clip=list_sig


    for rr, cc, nn in zip(reference_wl, cut_Ha_all, normaliseHa_all):
        cut_limits_min, cut_limits_max = rr+cc[0],  rr+cc[1]
        cut_limits_min_all.append(cut_limits_min);    cut_limits_max_all.append(cut_limits_max)

        mask_out_Ha_min, mask_out_Ha_max = rr+nn[0],  rr+nn[1]
        mask_out_Ha_min_all.append(mask_out_Ha_min);    mask_out_Ha_max_all.append(mask_out_Ha_max)
    
    
    np.savetxt("out/stacked_details.txt", np.asarray(binned_details),fmt="%s")








# now initiate the MCMC parameters
if starType1.startswith("D") and starType2.startswith("D"):
    #Uniform priors over which the variables will be allowed to vary
    p0T1=np.asarray(config_info["p0teff"])[0].astype(float)  # max teff is 14000 for DA
    p0T2=np.asarray(config_info["p0teff"])[1].astype(float)  # max teff is 14000 for DA
    p0logg1=np.asarray(config_info["p0logg"])[0].astype(float)
    p0logg2=np.asarray(config_info["p0logg"])[1].astype(float)
    p0HoverHe1=np.asarray(config_info["p0HoverHe"])[0].astype(float)
    p0HoverHe2=np.asarray(config_info["p0HoverHe"])[1].astype(float)
    if forced_Scaling==False:  p0scaling=np.asarray(config_info["p0scaling"])[0].astype(float)
    if found_parallax:
        p0parallax=[plax-10*plax_unc, plax+10*plax_unc]  # only used for plotting. gaussian prior
        if p0parallax[0]<0: p0parallax[0]=0.0001

    ndim, num_DBA = 0, 0
    p0range = np.array([]).reshape((0,2));   p0labels = np.array([])
    pos_min = np.array([]);   pos_max = np.array([])
    used_RV_boundaries, star_DBA = [], []



    if sys.argv[1] == "ATM":

        for cn, (fteff, flogg, fHHe, starType) in enumerate(zip([forced_teff1, forced_teff2], [forced_logg1, forced_logg2], [forced_HoverHe1, forced_HoverHe2], [starType1, starType2])):
            if fteff == 0 and flogg == 0:
                ndim+=2 # teff + logg
                if cn==0:
                    p0range = np.concatenate((p0range,[p0T1, p0logg1]));          p0labels = np.append(p0labels, ["T1", "logg1"])
                    pos_min = np.append(pos_min, [p0T1[0], p0logg1[0]]);   pos_max = np.append(pos_max, [p0T1[1], p0logg1[1]])
                elif cn==1:
                    p0range = np.concatenate((p0range,[p0T2, p0logg2]));          p0labels = np.append(p0labels, ["T2", "logg2"])
                    pos_min = np.append(pos_min, [p0T2[0], p0logg2[0]]);   pos_max = np.append(pos_max, [p0T2[1], p0logg2[1]])
            elif fteff == 0:
                ndim+=1 # teff
                if cn==0:
                    p0range = np.concatenate((p0range, [p0T1]));         p0labels = np.append(p0labels, ["T1"])
                    pos_min = np.append(pos_min, [p0T1[0]]);      pos_max = np.append(pos_max, [p0T1[1]])
                elif cn==1:
                    p0range = np.concatenate((p0range, [p0T2]));         p0labels = np.append(p0labels, ["T2"])
                    pos_min = np.append(pos_min, [p0T2[0]]);      pos_max = np.append(pos_max,[p0T2[1]])
            elif flogg == 0:
                ndim+=1 # logg
                if cn==0:
                    p0range = np.concatenate((p0range, [p0logg1]));         p0labels = np.append(p0labels, ["logg"])
                    pos_min = np.append(pos_min, [p0logg1[0]]);      pos_max = np.append(pos_max, [p0logg1[1]])
                if cn==1:
                    p0range = np.concatenate((p0range, [p0logg2]));         p0labels = np.append(p0labels, ["logg2"])
                    pos_min = np.append(pos_min, [p0logg2[0]]);      pos_max = np.append(pos_max, [p0logg2[1]])
            
            if starType=="DBA" and fHHe==0:
                num_DBA+=1
                star_DBA.append(str(cn+1))

        for xx in range(num_DBA):
            ndim+=1 # teff + logg + H/He
            if star_DBA[xx]=="1":
                p0range = np.concatenate((p0range, np.array([p0HoverHe1])));         p0labels = np.concatenate((p0labels, np.array(["H/He1"])))
                pos_min = np.concatenate((pos_min, np.array([p0HoverHe1[0]])));   pos_max = np.concatenate((pos_max, np.array([p0HoverHe1[1]])))
            else:
                p0range = np.concatenate((p0range, np.array([p0HoverHe2])));         p0labels = np.concatenate((p0labels, np.array(["H/He2"])))
                pos_min = np.concatenate((pos_min, np.array([p0HoverHe2[0]])));   pos_max = np.concatenate((pos_max, np.array([p0HoverHe2[1]])))

        

    if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
        if forced_K1==False and forced_K2==False:
            for i, (fi, arv, arvBounds1, arvBounds2) in enumerate(zip(input_files, share_rv, RV_boundaries1, RV_boundaries2)):
                if arv==-1:
                    # this is to configure the unique RVs in the mcmc
                    pos_min = np.append(pos_min, arvBounds1[0]);   pos_max = np.append(pos_max, arvBounds1[1])
                    pos_min = np.append(pos_min, arvBounds2[0]);   pos_max = np.append(pos_max, arvBounds2[1])
                    ndim+=2 #Number of individual RVs to vary in the mcmc run
                    
                    # this is for corner plotting
                    #try:
                    p0range = np.concatenate((p0range, np.array([arvBounds1])))
                    p0range = np.concatenate((p0range, np.array([arvBounds2])))
                    #except:
                    #    p0range = np.array([arvBounds1])
                    #    p0range = np.concatenate((p0range, np.array([arvBounds2])))
                    p0labels = np.append(p0labels,"RV1_"+str(i))
                    p0labels = np.append(p0labels,"RV2_"+str(i))
                    
                    
                    used_RV_boundaries.append(arvBounds1)
                    used_RV_boundaries.append(arvBounds2)
                else:   
                    if rank==1:   print("ignoring entry of index (zero-indexed) ", i, " as a separate set of RVs, it is paired with the RV of index ", arv)

        elif forced_K1==False:
            for i, (fi, arv, arvBounds1) in enumerate(zip(input_files, share_rv, RV_boundaries1)):
                if arv==-1:
                    # this is to configure the unique RVs in the mcmc
                    pos_min = np.append(pos_min, arvBounds1[0]);   pos_max = np.append(pos_max, arvBounds1[1])
                    ndim+=1 #Number of individual RVs to vary in the mcmc run
                    
                    # this is for corner plotting
                    #try:      
                    p0range = np.concatenate((p0range, np.array([arvBounds1])))
                    #except:   p0range = np.array([arvBounds1])
                    p0labels = np.append(p0labels,"RV1_"+str(i))
                    
                    
                    used_RV_boundaries.append(arvBounds1)
                else:   
                    if rank==1:   print("ignoring entry of index (zero-indexed) ", i, " as a separate set of RVs, it is paired with the RV of index ", arv)

        elif forced_K2==False:
            for i, (fi, arv, arvBounds2) in enumerate(zip(input_files, share_rv, RV_boundaries2)):
                if arv==-1:
                    # this is to configure the unique RVs in the mcmc
                    pos_min = np.append(pos_min, arvBounds2[0]);   pos_max = np.append(pos_max, arvBounds2[1])
                    ndim+=1 #Number of individual RVs to vary in the mcmc run
                    
                    # this is for corner plotting
                    #try:     
                    p0range = np.concatenate((p0range, np.array([arvBounds2])))
                    #except:  p0range = np.array([arvBounds2])
                    p0labels = np.append(p0labels,"RV2_"+str(i))
                    
                    
                    used_RV_boundaries.append(arvBounds2)
                else:   
                    if rank==1:   print("ignoring entry of index (zero-indexed) ", i, " as a separate set of RVs, it is paired with the RV of index ", arv)
        
    elif isinstance(forced_P0, float) and isinstance(forced_T0, float):
        if forced_K1=="Fit":
            pos_min = np.append(pos_min, p0K1[0]);   pos_max = np.append(pos_max, p0K1[1])
            ndim+=1
            
            # this is for corner plotting
            #try:     
            p0range = np.concatenate((p0range, np.array([p0K1])))
            #except:  p0range = np.array([p0K1])
            p0labels = np.append(p0labels,"K1")
        
        if forced_Vgamma1=="Fit":
            pos_min = np.append(pos_min, p0Vgamma1[0]);   pos_max = np.append(pos_max, p0Vgamma1[1])
            ndim+=1
            
            # this is for corner plotting
            #try:     
            p0range = np.concatenate((p0range, np.array([p0Vgamma1])))
            #except:  p0range = np.array([p0Vgamma1])
            p0labels = np.append(p0labels,"Vg1")

        if forced_K2=="Fit":
            pos_min = np.append(pos_min, p0K2[0]);   pos_max = np.append(pos_max, p0K2[1])
            ndim+=1
            
            # this is for corner plotting
            #try:     
            p0range = np.concatenate((p0range, np.array([p0K2])))
            #except:  p0range = np.array([p0K2])
            p0labels = np.append(p0labels,"K2")
        
        if forced_Vgamma2=="Fit":
            pos_min = np.append(pos_min, p0Vgamma2[0]);   pos_max = np.append(pos_max, p0Vgamma2[1])
            ndim+=1
            
            # this is for corner plotting
            #try:     
            p0range = np.concatenate((p0range, np.array([p0Vgamma2])))
            #except:  p0range = np.array([p0Vgamma2])
            p0labels = np.append(p0labels,"Vg2")
    else: raise ValueError



    if sys.argv[1] == "ATM":
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

elif (starType1.startswith("LL") or starType1.startswith("GG") or starType1.startswith("GL")) and (starType2.startswith("LL") or starType2.startswith("GG") or starType2.startswith("GL")):
    #Uniform priors over which the variables will be allowed to vary

    ndim = 0
    p0range = np.array([]);   p0labels = np.array([])
    pos_min = np.array([]);   pos_max = np.array([])
    used_RV_boundaries = []


    if sys.argv[1] == "ATM" and not "run_double_one_by_oneMCMC" in sys.argv:
        ndim+=8
        p0range = np.array([A1_1_boundaries, sigma1_1_boundaries, A1_2_boundaries, sigma1_2_boundaries, A2_1_boundaries, sigma2_1_boundaries, A2_2_boundaries, sigma2_2_boundaries]);   p0labels = np.array(["A1_1", "sig1_1", "A1_2", "sig1_2", "A2_1", "sig2_1", "A2_2", "sig2_2"])
        pos_min = np.array([A1_1_boundaries[0], sigma1_1_boundaries[0], A1_2_boundaries[0], sigma1_2_boundaries[0], A2_1_boundaries[0], sigma2_1_boundaries[0], A2_2_boundaries[0], sigma2_2_boundaries[0]])
        pos_max = np.array([A1_1_boundaries[1], sigma1_1_boundaries[1], A1_2_boundaries[1], sigma1_2_boundaries[1], A2_1_boundaries[1], sigma2_1_boundaries[1], A2_2_boundaries[1], sigma2_2_boundaries[1]])
            

        

    if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
        if forced_K1==False and forced_K2==False:
            for i, (fi, arv, arvBounds1, arvBounds2) in enumerate(zip(input_files, share_rv, RV_boundaries1, RV_boundaries2)):
                if arv==-1:
                    # this is to configure the unique RVs in the mcmc
                    pos_min = np.append(pos_min, arvBounds1[0]);   pos_max = np.append(pos_max, arvBounds1[1])
                    pos_min = np.append(pos_min, arvBounds2[0]);   pos_max = np.append(pos_max, arvBounds2[1])
                    ndim+=2 #Number of individual RVs to vary in the mcmc run
                    
                    # this is for corner plotting
                    #try:
                    p0range = np.concatenate((p0range, np.array([arvBounds1])))
                    p0range = np.concatenate((p0range, np.array([arvBounds2])))
                    #except:
                    #    p0range = np.array([arvBounds1])
                    #    p0range = np.concatenate((p0range, np.array([arvBounds2])))
                    p0labels = np.append(p0labels,"RV1_"+str(i))
                    p0labels = np.append(p0labels,"RV2_"+str(i))
                    
                    
                    used_RV_boundaries.append(arvBounds1)
                    used_RV_boundaries.append(arvBounds2)
                else:   
                    if rank==1:   print("ignoring entry of index (zero-indexed) ", i, " as a separate set of RVs, it is paired with the RV of index ", arv)

        elif forced_K1==False:
            for i, (fi, arv, arvBounds1) in enumerate(zip(input_files, share_rv, RV_boundaries1)):
                if arv==-1:
                    # this is to configure the unique RVs in the mcmc
                    pos_min = np.append(pos_min, arvBounds1[0]);   pos_max = np.append(pos_max, arvBounds1[1])
                    ndim+=1 #Number of individual RVs to vary in the mcmc run
                    
                    # this is for corner plotting
                    #try:      
                    p0range = np.concatenate((p0range, np.array([arvBounds1])))
                    #except:   p0range = np.array([arvBounds1])
                    p0labels = np.append(p0labels,"RV1_"+str(i))
                    
                    
                    used_RV_boundaries.append(arvBounds1)
                else:   
                    if rank==1:   print("ignoring entry of index (zero-indexed) ", i, " as a separate set of RVs, it is paired with the RV of index ", arv)

        elif forced_K2==False:
            for i, (fi, arv, arvBounds2) in enumerate(zip(input_files, share_rv, RV_boundaries2)):
                if arv==-1:
                    # this is to configure the unique RVs in the mcmc
                    pos_min = np.append(pos_min, arvBounds2[0]);   pos_max = np.append(pos_max, arvBounds2[1])
                    ndim+=1 #Number of individual RVs to vary in the mcmc run
                    
                    # this is for corner plotting
                    #try:     
                    p0range = np.concatenate((p0range, np.array([arvBounds2])))
                    #except:  p0range = np.array([arvBounds2])
                    p0labels = np.append(p0labels,"RV2_"+str(i))
                    
                    
                    used_RV_boundaries.append(arvBounds2)
                else:   
                    if rank==1:   print("ignoring entry of index (zero-indexed) ", i, " as a separate set of RVs, it is paired with the RV of index ", arv)
        
    elif isinstance(forced_P0, float) and isinstance(forced_T0, float):
        if forced_K1=="Fit":
            pos_min = np.append(pos_min, p0K1[0]);   pos_max = np.append(pos_max, p0K1[1])
            ndim+=1
            
            # this is for corner plotting
            #try:     
            p0range = np.concatenate((p0range, np.array([p0K1])))
            #except:  p0range = np.array([p0K1])
            p0labels = np.append(p0labels,"K1")
        
        if forced_Vgamma1=="Fit":
            pos_min = np.append(pos_min, p0Vgamma1[0]);   pos_max = np.append(pos_max, p0Vgamma1[1])
            ndim+=1
            
            # this is for corner plotting
            #try:     
            p0range = np.concatenate((p0range, np.array([p0Vgamma1])))
            #except:  p0range = np.array([p0Vgamma1])
            p0labels = np.append(p0labels,"Vg1")

        if forced_K2=="Fit":
            pos_min = np.append(pos_min, p0K2[0]);   pos_max = np.append(pos_max, p0K2[1])
            ndim+=1
            
            # this is for corner plotting
            #try:     
            p0range = np.concatenate((p0range, np.array([p0K2])))
            #except:  p0range = np.array([p0K2])
            p0labels = np.append(p0labels,"K2")
        
        if forced_Vgamma2=="Fit":
            pos_min = np.append(pos_min, p0Vgamma2[0]);   pos_max = np.append(pos_max, p0Vgamma2[1])
            ndim+=1
            
            # this is for corner plotting
            #try:     
            p0range = np.concatenate((p0range, np.array([p0Vgamma2])))
            #except:  p0range = np.array([p0Vgamma2])
            p0labels = np.append(p0labels,"Vg2")
    else: raise ValueError


elif starType1=="quadLorentz" and starType2=="quadLorentz":
    ndim = 0
    p0range = np.array([]).reshape((0,2));   p0labels = np.array([])
    pos_min = np.array([]);   pos_max = np.array([])
    used_RV_boundaries = []


    if sys.argv[1] == "ATM"  and not "run_double_one_by_oneMCMC" in sys.argv:
        ndim+=4
        p0range = np.array([A1_1_boundaries, sigma1_1_boundaries, A2_1_boundaries, sigma2_1_boundaries]);   p0labels = np.array(["A1_1", "sig1_1", "A2_1", "sig2_1"])
        pos_min = np.array([A1_1_boundaries[0], sigma1_1_boundaries[0], A2_1_boundaries[0], sigma2_1_boundaries[0]])
        pos_max = np.array([A1_1_boundaries[1], sigma1_1_boundaries[1], A2_1_boundaries[1], sigma2_1_boundaries[1]])
            

        

    if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
        if forced_K1==False and forced_K2==False:
            for i, (fi, arv, arvBounds1, arvBounds2) in enumerate(zip(input_files, share_rv, RV_boundaries1, RV_boundaries2)):
                if arv==-1:
                    # this is to configure the unique RVs in the mcmc
                    pos_min = np.append(pos_min, arvBounds1[0]);   pos_max = np.append(pos_max, arvBounds1[1])
                    pos_min = np.append(pos_min, arvBounds2[0]);   pos_max = np.append(pos_max, arvBounds2[1])
                    ndim+=2 #Number of individual RVs to vary in the mcmc run
                    
                    # this is for corner plotting
                    #try:
                    p0range = np.concatenate((p0range, np.array([arvBounds1])))
                    p0range = np.concatenate((p0range, np.array([arvBounds2])))
                    #except:
                    #    p0range = np.array([arvBounds1])
                    #    p0range = np.concatenate((p0range, np.array([arvBounds2])))
                    p0labels = np.append(p0labels,"RV1_"+str(i))
                    p0labels = np.append(p0labels,"RV2_"+str(i))
                    
                    
                    used_RV_boundaries.append(arvBounds1)
                    used_RV_boundaries.append(arvBounds2)
                else:   
                    if rank==1:   print("ignoring entry of index (zero-indexed) ", i, " as a separate set of RVs, it is paired with the RV of index ", arv)

        elif forced_K1==False:
            for i, (fi, arv, arvBounds1) in enumerate(zip(input_files, share_rv, RV_boundaries1)):
                if arv==-1:
                    # this is to configure the unique RVs in the mcmc
                    pos_min = np.append(pos_min, arvBounds1[0]);   pos_max = np.append(pos_max, arvBounds1[1])
                    ndim+=1 #Number of individual RVs to vary in the mcmc run
                    
                    # this is for corner plotting
                    #try:      
                    p0range = np.concatenate((p0range, np.array([arvBounds1])))
                    #except:   p0range = np.array([arvBounds1])
                    p0labels = np.append(p0labels,"RV1_"+str(i))
                    
                    
                    used_RV_boundaries.append(arvBounds1)
                else:   
                    if rank==1:   print("ignoring entry of index (zero-indexed) ", i, " as a separate set of RVs, it is paired with the RV of index ", arv)

        elif forced_K2==False:
            for i, (fi, arv, arvBounds2) in enumerate(zip(input_files, share_rv, RV_boundaries2)):
                if arv==-1:
                    # this is to configure the unique RVs in the mcmc
                    pos_min = np.append(pos_min, arvBounds2[0]);   pos_max = np.append(pos_max, arvBounds2[1])
                    ndim+=1 #Number of individual RVs to vary in the mcmc run
                    
                    # this is for corner plotting
                    #try:     
                    p0range = np.concatenate((p0range, np.array([arvBounds2])))
                    #except:  p0range = np.array([arvBounds2])
                    p0labels = np.append(p0labels,"RV2_"+str(i))
                    
                    
                    used_RV_boundaries.append(arvBounds2)
                else:   
                    if rank==1:   print("ignoring entry of index (zero-indexed) ", i, " as a separate set of RVs, it is paired with the RV of index ", arv)
        
    elif isinstance(forced_P0, float) and isinstance(forced_T0, float):
        if forced_K1=="Fit":
            pos_min = np.append(pos_min, p0K1[0]);   pos_max = np.append(pos_max, p0K1[1])
            ndim+=1
            
            # this is for corner plotting
            #try:     
            p0range = np.concatenate((p0range, np.array([p0K1])))
            #except:  p0range = np.array([p0K1])
            p0labels = np.append(p0labels,"K1")
        
        if forced_Vgamma1=="Fit":
            pos_min = np.append(pos_min, p0Vgamma1[0]);   pos_max = np.append(pos_max, p0Vgamma1[1])
            ndim+=1
            
            # this is for corner plotting
            #try:     
            p0range = np.concatenate((p0range, np.array([p0Vgamma1])))
            #except:  p0range = np.array([p0Vgamma1])
            p0labels = np.append(p0labels,"Vg1")

        if forced_K2=="Fit":
            pos_min = np.append(pos_min, p0K2[0]);   pos_max = np.append(pos_max, p0K2[1])
            ndim+=1
            
            # this is for corner plotting
            #try:     
            p0range = np.concatenate((p0range, np.array([p0K2])))
            #except:  p0range = np.array([p0K2])
            p0labels = np.append(p0labels,"K2")
        
        if forced_Vgamma2=="Fit":
            pos_min = np.append(pos_min, p0Vgamma2[0]);   pos_max = np.append(pos_max, p0Vgamma2[1])
            ndim+=1
            
            # this is for corner plotting
            #try:     
            p0range = np.concatenate((p0range, np.array([p0Vgamma2])))
            #except:  p0range = np.array([p0Vgamma2])
            p0labels = np.append(p0labels,"Vg2")
    else: raise ValueError

else: raise ValueError("WD-BASS Error: Star types must be two degenerate stars (Dx) or both be of: GG, LL")
    

psize = pos_max - pos_min
p0 = [pos_min + psize*np.random.rand(ndim) for i in range(nwalkers)]




if starType1=="quadLorentz" and starType2=="quadLorentz":  ## all of this is to do the polynomial fitting to the wings of spectral lines for when I do Lorentzian components
    # Polynomial degree (reduce if needed)
    poly_degree = 6
    
    def fast_polyval(coeffs, x):
        result = coeffs[0] * np.ones_like(x)
        for c in coeffs[1:]:
            result = result * x + c
        return result
    
    def poly_fit(x, coeffs):
        return fast_polyval(coeffs[::-1], x)
    
    # Define constraints
    def constraint_function(coeffs, x_scaled):
        x_left, x_right = x_scaled[0], x_scaled[-1]  # Edges
        y_left =  fast_polyval(coeffs[::-1], x_left) - 1  # Ensure y=1 at left
        y_right = fast_polyval(coeffs[::-1], x_right) - 1  # Ensure y=1 at right
        
        # First derivative constraints (flat edges)
        deriv_coeffs = np.arange(1, poly_degree + 1) * coeffs[1:]  # Skip the 0th coefficient
        
        deriv_left = fast_polyval(deriv_coeffs[::-1], x_left)
        deriv_right = fast_polyval(deriv_coeffs[::-1], x_right)
        
        return [y_left, y_right, deriv_left, deriv_right]
    	
    # Define loss function (least squares)
    def loss(coeffs, x_scaled, yquad, yquaderr):
        return npsum( np_square((fast_polyval(coeffs[::-1], x_scaled) - yquad) / yquaderr))
    
    import scipy.optimize as opt





    
    
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


try:
    if fit_phot_SED==False:
        mask = (wl_all_logg > min_wl_all_files) & (wl_all_logg < max_wl_all_files)
        Teff_all_logg = Teff_all_logg[mask];   wl_all_logg = wl_all_logg[mask];     flux_all_logg = flux_all_logg[mask];      logg_all_logg = logg_all_logg[mask]
except: None





@njit
def return_DAgrids(temperature_star, logg_star):
    pier_or_antoine="mixed"#"pier3Dphot_antoine1Dspec"
    minval=10000;   maxval=-10000
    #if temperature_star>25000: raise ValueError
    
    if pier_or_antoine == "mixed":
        if logg_star>=6.5:
            temdiff = temperature_star - unique_Teffs_synth
            Teff_min=unique_Teffs_synth[npargwhere(temdiff==npamin(temdiff[temdiff>0]))[0][0]]
            Teff_max=unique_Teffs_synth[npargwhere(temdiff==npamax(temdiff[temdiff<0]))[0][0]]
            
        else: raise ValueError
        
            
            # then find the nearest 2 loggs
        list_search = [6.5,7,7.5,8,8.5,9]
        
        for logg_opts in list_search:
            if logg_opts-logg_star<np.abs(logg_opts-maxval):    maxval = logg_opts
            if logg_star-logg_opts>0:   minval = logg_opts
        
        
        mask_logg = (logg_all_synth<=maxval) & (logg_all_synth>=minval) & (Teff_all_synth<=Teff_max) & (Teff_all_synth>=Teff_min)
        
        
        Grav1_N = logg_all_synth[mask_logg];    wl_all1_N=wl_all_synth[mask_logg];    flux1_N=flux_all_synth[mask_logg];    Teff1_N=Teff_all_synth[mask_logg]
        
        
        if len(npunique(Teff1_N))==3:
            un_teffs = npunique(Teff1_N)
            if temperature_star<un_teffs[1]:   newmask = (Teff1_N!=un_teffs[2])
            else:   newmask = (Teff1_N!=un_teffs[0])
            
            Grav1_N, wl_all1_N, flux1_N, Teff1_N = Grav1_N[newmask], wl_all1_N[newmask], flux1_N[newmask], Teff1_N[newmask]
        
        return Grav1_N, wl_all1_N, flux1_N, Teff1_N


@njit
def return_DBAgrids(temperature_star, logg_star, HoverHestar):
    # find nearest Teff in models in the DA/DBA/DC grid, then take the two for interpolation limits
    minval=10000;   maxval=-10000
    if logg_star>=7.5 and logg_star<8: minval=7.5; maxval=8
    elif logg_star>=8 and logg_star<8.5: minval=8; maxval=8.5
    elif logg_star>=8.5 and logg_star<=9: minval=8.5; maxval=9
    elif logg_star<7.5: raise ValueError("Minimum logg for grid is 7.5. You attempted with " + str(logg_star))
    else: print(logg_star); print(logg_star); print(logg_star);  raise ValueError
    
    
    
    if HoverHestar>=2 and HoverHestar<=5: HoverHe_min=2; HoverHe_max=5
    elif HoverHestar>=5 and HoverHestar<=8: HoverHe_min=5; HoverHe_max=8
    elif HoverHestar>=8 and HoverHestar<=30: HoverHe_min=8; HoverHe_max=30
    else: print(HoverHestar);  print(HoverHestar);  raise ValueError
    

    mask_logg = (logg_all_logg<=maxval) & (logg_all_logg>=minval) & (HoverHe_all_logg_DBA>=HoverHe_min) & (HoverHe_all_logg_DBA<=HoverHe_max)
    Grav_N = logg_all_logg[mask_logg];    wl_all_N=wl_all_logg[mask_logg];    flux_N=flux_all_logg[mask_logg];    Teff_N=Teff_all_logg[mask_logg];     HoverHe_N = HoverHe_all_logg_DBA[mask_logg]
    
    
    unique_Teffs_DBA=npunique(Teff_N)
    
    temdiff = temperature_star - unique_Teffs_DBA
    Teff_min=unique_Teffs_DBA[npargwhere(temdiff==npamin(temdiff[temdiff>0]))[0][0]]
    Teff_max=unique_Teffs_DBA[npargwhere(temdiff==npamax(temdiff[temdiff<0]))[0][0]]
    
    
    mask_logg_then_T = (Teff_N<=Teff_max) & (Teff_N>=Teff_min)
    
    
    if len(Grav_N[mask_logg_then_T])==0:
        raise ValueError("ERROR - Grid of DBA/DB/DC")
    
    Grav1_N = Grav_N[mask_logg_then_T];    wl_all1_N=wl_all_N[mask_logg_then_T];    flux1_N=flux_N[mask_logg_then_T];    Teff1_N=Teff_N[mask_logg_then_T];   HoverHe1_N=HoverHe_N[mask_logg_then_T]
    return Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N




@njit
def return_DBgrids(temperature_star, logg_star):
    # find nearest Teff in models in the DA/DBA/DC grid, then take the two for interpolation limits
    minval=10000;   maxval=-10000
    if logg_star>=7.5 and logg_star<8: minval=7.5; maxval=8
    elif logg_star>=8 and logg_star<8.5: minval=8; maxval=8.5
    elif logg_star>=8.5 and logg_star<=9: minval=8.5; maxval=9
    elif logg_star<7.5: raise ValueError("Minimum logg for grid is 7.5. You attempted with " + str(logg_star))
    else: print(logg_star); print(logg_star); print(logg_star);  raise ValueError
    
    
    
    mask_logg = (logg_all_logg<=maxval) & (logg_all_logg>=minval)
    Grav_N = logg_all_logg[mask_logg];    wl_all_N=wl_all_logg[mask_logg];    flux_N=flux_all_logg[mask_logg];    Teff_N=Teff_all_logg[mask_logg]
    
    
    unique_Teffs_DC=npunique(Teff_N)
    
    temdiff = temperature_star - unique_Teffs_DC
    Teff_min=unique_Teffs_DC[npargwhere(temdiff==npamin(temdiff[temdiff>0]))[0][0]]
    Teff_max=unique_Teffs_DC[npargwhere(temdiff==npamax(temdiff[temdiff<0]))[0][0]]
    
    
    mask_logg_then_T = (Teff_N<=Teff_max) & (Teff_N>=Teff_min)
    
    
    if len(Grav_N[mask_logg_then_T])==0:
        raise ValueError("ERROR - Grid of DBA/DB/DC")
    
    Grav1_N = Grav_N[mask_logg_then_T];    wl_all1_N=wl_all_N[mask_logg_then_T];    flux1_N=flux_N[mask_logg_then_T];    Teff1_N=Teff_N[mask_logg_then_T]
    return Grav1_N, wl_all1_N, flux1_N, Teff1_N


@njit
def return_DCgrids(temperature_star, logg_star):
    # find nearest Teff in models in the DA/DBA/DC grid, then take the two for interpolation limits
    minval=10000;   maxval=-10000
    if logg_star>=7.5 and logg_star<8: minval=7.5; maxval=8
    elif logg_star>=8 and logg_star<8.5: minval=8; maxval=8.5
    elif logg_star>=8.5 and logg_star<=9: minval=8.5; maxval=9
    elif logg_star<7.5: raise ValueError("Minimum logg for grid is 7.5. You attempted with " + str(logg_star))
    else: print(logg_star); print(logg_star); print(logg_star);  raise ValueError
    
    
    
    mask_logg = (logg_all_logg<=maxval) & (logg_all_logg>=minval)
    Grav_N = logg_all_logg[mask_logg];    wl_all_N=wl_all_logg[mask_logg];    flux_N=flux_all_logg[mask_logg];    Teff_N=Teff_all_logg[mask_logg]
    
    
    unique_Teffs_DC=npunique(Teff_N)
    
    temdiff = temperature_star - unique_Teffs_DC
    Teff_min=unique_Teffs_DC[npargwhere(temdiff==npamin(temdiff[temdiff>0]))[0][0]]
    Teff_max=unique_Teffs_DC[npargwhere(temdiff==npamax(temdiff[temdiff<0]))[0][0]]
    
    
    mask_logg_then_T = (Teff_N<=Teff_max) & (Teff_N>=Teff_min)
    
    
    if len(Grav_N[mask_logg_then_T])==0:
        raise ValueError("ERROR - Grid of DBA/DB/DC")
    
    Grav1_N = Grav_N[mask_logg_then_T];    wl_all1_N=wl_all_N[mask_logg_then_T];    flux1_N=flux_N[mask_logg_then_T];    Teff1_N=Teff_N[mask_logg_then_T]
    return Grav1_N, wl_all1_N, flux1_N, Teff1_N




@njit
def return_model_spectrum_DA(wl_all1_N, ref_wl, cut_limits_min, cut_limits_max, Grav1_N, flux1_N, Teff1_N, temperature_star, logg_star):
    if ref_wl>6000:    excess_slack=10 # angstroms to allow variation with RV of star. Here, Halpha goes +-450kms-1. Increase this if RV diff larger
    else:  excess_slack=7.5 # angstroms to allow variation with RV of star. Here, Hbeta goes +-500kms-1. Increase this if RV diff larger
    
    mask_logg_wl = ((wl_all1_N > ref_wl+cut_limits_min-excess_slack) & (wl_all1_N < ref_wl+cut_limits_max+excess_slack))
    Grav_N_N = Grav1_N[mask_logg_wl];   wl_all_N_N=wl_all1_N[mask_logg_wl];    flux_N_N=flux1_N[mask_logg_wl];    Teff_N_N=Teff1_N[mask_logg_wl]
    
    # interpolate for a model at the reference wavelength with this mcmc interation    
    wl_grid, unique_Ts, unique_Gs = npunique(wl_all_N_N), npunique(Teff_N_N), npunique(Grav_N_N)

    ffffsss1 = flux_N_N[(Teff_N_N==unique_Ts[0]) & (Grav_N_N==unique_Gs[0])]
    ffffsss2 = flux_N_N[(Teff_N_N==unique_Ts[0]) & (Grav_N_N==unique_Gs[1])]
    ffffsss3 = flux_N_N[(Teff_N_N==unique_Ts[1]) & (Grav_N_N==unique_Gs[0])]
    ffffsss4 = flux_N_N[(Teff_N_N==unique_Ts[1]) & (Grav_N_N==unique_Gs[1])]
    
    model_spectrum = (ffffsss1 * (unique_Ts[1] - temperature_star) * (unique_Gs[1] - logg_star) +           ffffsss3 * (temperature_star - unique_Ts[0]) * (unique_Gs[1] - logg_star) +            ffffsss2 * (unique_Ts[1] - temperature_star) * (logg_star - unique_Gs[0]) +            ffffsss4 * (temperature_star - unique_Ts[0]) * (logg_star - unique_Gs[0])           ) / ((unique_Ts[1] - unique_Ts[0]) * (unique_Gs[1] - unique_Gs[0]))
    
    #plt.plot(wl_grid, model_spectrum,c='r', ls='--');  plt.show()
    
    
    #plt.plot(wl_grid, model_spectrum,c='r'); plt.plot(wl_grid, model_spectrum_starx,c='k'); plt.show()
    
    if ref_wl>6500:     fine_grid=nplinspace(npamin(wl_grid),npamax(wl_grid),int((npamax(wl_grid) - npamin(wl_grid))*50)) # 0.02AA spacing
    elif ref_wl>4500:   fine_grid=nplinspace(npamin(wl_grid),npamax(wl_grid),int((npamax(wl_grid) - npamin(wl_grid))*33)) # 0.03AA spacing
    elif ref_wl>4200:   fine_grid=nplinspace(npamin(wl_grid),npamax(wl_grid),int((npamax(wl_grid) - npamin(wl_grid))*25)) # 0.04AA spacing
    else:               fine_grid=nplinspace(npamin(wl_grid),npamax(wl_grid),int((npamax(wl_grid) - npamin(wl_grid))*25)) # 0.04AA spacing
    
    
    #fine_grid=nplinspace(npamin(wl_grid),npamax(wl_grid),int((npamax(wl_grid) - npamin(wl_grid))*400)) # 0.0025AA spacing
    model_spectrum=interp(fine_grid, wl_grid, model_spectrum)
    
    return fine_grid, model_spectrum
    


#@njit
# temperature grid for DBs is irregular. Don't try to speed it up and accept the griddata way
def return_model_spectrum_DBA(wl_all1_N, ref_wl, cut_limits_min, cut_limits_max, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, temperature_star, logg_star, HoverHe_star):
    if ref_wl>6000:    excess_slack=10 # angstroms to allow variation with RV of star. Here, Halpha goes +-450kms-1. Increase this if RV diff larger
    else:  excess_slack=7.5 # angstroms to allow variation with RV of star. Here, Hbeta goes +-500kms-1. Increase this if RV diff larger
    
    mask_logg_wl = ((wl_all1_N > ref_wl+cut_limits_min-excess_slack) & (wl_all1_N < ref_wl+cut_limits_max+excess_slack))
    Grav_N_N = Grav1_N[mask_logg_wl];   wl_all_N_N=wl_all1_N[mask_logg_wl];    flux_N_N=flux1_N[mask_logg_wl];    Teff_N_N=Teff1_N[mask_logg_wl];    HoverHe_N_N = HoverHe1_N[mask_logg_wl]
    
    # interpolate for a model at the reference wavelength with this mcmc interation    
    wl_grid = npunique(wl_all_N_N)
    model_spectrum=griddata(np.array([Teff_N_N, wl_all_N_N, Grav_N_N, HoverHe_N_N]).T, flux_N_N, np.array([np.full((len(wl_grid),),temperature_star), wl_grid, np.full((len(wl_grid),),logg_star), np.full((len(wl_grid),),HoverHe_star)]).T, method="linear")
        
    
    if ref_wl>6500:     fine_grid=nplinspace(npamin(wl_grid),npamax(wl_grid),int((npamax(wl_grid) - npamin(wl_grid))*10)) # 0.1AA spacing
    elif ref_wl>4500:   fine_grid=nplinspace(npamin(wl_grid),npamax(wl_grid),int((npamax(wl_grid) - npamin(wl_grid))*10)) # 0.1AA spacing
    elif ref_wl>4200:   fine_grid=nplinspace(npamin(wl_grid),npamax(wl_grid),int((npamax(wl_grid) - npamin(wl_grid))*10)) # 0.1AA spacing
    else:               fine_grid=nplinspace(npamin(wl_grid),npamax(wl_grid),int((npamax(wl_grid) - npamin(wl_grid))*10)) # 0.1AA spacing
    
    model_spectrum=interp(fine_grid, wl_grid, model_spectrum)
    
    return fine_grid, model_spectrum



if forced_Scaling=="WD":
    install_path = os.environ['WD_BASS_INSTALL_DIR']
    loaded_Althaus = np.load(install_path + "/saved_MTR/Althaus_2013_full_nomasses.npy")
    loaded_Istrate = np.load(install_path + "/saved_MTR/Istrate_Z0p02_diffusion_nomasses.npy")
    loaded_CO = np.load(install_path + "/saved_MTR/table_valuesCO.npy")



checkT1, checklogg1, checkT2, checklogg2, checkHoverHe1, checkHoverHe2, checkscaling = False, False, False, False, False, False, False
if forced_teff1==0: checkT1=True
if forced_logg1==0: checklogg1=True
if forced_teff2==0: checkT2=True
if forced_logg2==0: checklogg2=True
if forced_HoverHe1==0 and starType1=="DBA": checkHoverHe1=True
if forced_HoverHe2==0 and starType2=="DBA": checkHoverHe2=True
if forced_Scaling==False: checkscaling=True




#Test that next step falls within boundaries allowed by the priors
def lnprior(theta, arguments):
    input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, normalised_wavelength, normalised_flux, normalised_err,  inp_resolution, used_RV_boundaries, HJD_values, sigma_clip, the_unique_wavelengths, spec1wl, spec1flux, spec2wl, spec2flux, inputScaling = arguments
    if "T1" in p0labels:            args = npargwhere(p0labels=="T1")[0][0];          T1 = theta[args]
    else: T1=forced_teff1
    if "logg1" in p0labels:         args = npargwhere(p0labels=="logg1")[0][0];       logg1 = theta[args]
    else: logg1=forced_logg1
    if "T2" in p0labels:            args = npargwhere(p0labels=="T2")[0][0];          T2 = theta[args]
    else: T2=forced_teff2
    if "logg2" in p0labels:         args = npargwhere(p0labels=="logg2")[0][0];       logg2 = theta[args]
    else: logg2=forced_logg2
    if "H/He1" in p0labels:         args = npargwhere(p0labels=="H/He1")[0][0];       HoverHe1 = theta[args]
    else:
        try: HoverHe1 = forced_HoverHe1
        except: None
    if "H/He2" in p0labels:         args = npargwhere(p0labels=="H/He2")[0][0];       HoverHe2 = theta[args]
    else:
        try: HoverHe2 = forced_HoverHe2
        except: None
    if "K1" in p0labels:            args = npargwhere(p0labels=="K1")[0][0];          mcmc_K1 = theta[args]
    if "K2" in p0labels:            args = npargwhere(p0labels=="K2")[0][0];          mcmc_K2 = theta[args]
    if "Vg1" in p0labels:           args = npargwhere(p0labels=="Vg1")[0][0];         mcmc_Vgamma1 = theta[args]
    if "Vg2" in p0labels:           args = npargwhere(p0labels=="Vg2")[0][0];         mcmc_Vgamma2 = theta[args]
    if "Scaling" in p0labels:       args = npargwhere(p0labels=="Scaling")[0][0];     Scaling = theta[args]
    if "RV1_0" in p0labels:         num_start_RVs = npargwhere(p0labels=="RV1_0")[0][0]
    elif "RV1_1" in p0labels:       num_start_RVs = npargwhere(p0labels=="RV1_1")[0][0]
    elif "RV1_2" in p0labels:       num_start_RVs = npargwhere(p0labels=="RV1_2")[0][0]
    elif "RV1_3" in p0labels:       num_start_RVs = npargwhere(p0labels=="RV1_3")[0][0]
    elif "RV1_4" in p0labels:       num_start_RVs = npargwhere(p0labels=="RV1_4")[0][0]
    if "Parallax" in p0labels:      args = npargwhere(p0labels=="Parallax")[0][0];     mcmc_parallax = theta[args]
    
    
    if T1<T2:
        return -np_inf
    
    if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
        if forced_Scaling==False: RV=theta[num_start_RVs:-1];  Scaling=theta[-1]
        else:                     RV=theta[num_start_RVs:]
        
        
        
    passed = True
    

    
    if checkT1:
        if not p0T1[0]<T1<p0T1[1]: passed = False
    if checklogg1:  
        if not p0logg1[0]<logg1<p0logg1[1]: passed = False
    if checkT2:  
        if not p0T2[0]<T2<p0T2[1]: passed = False
    if checklogg2:  
        if not p0logg2[0]<logg2<p0logg2[1]: passed = False
    if checkHoverHe1:  
        if not p0HoverHe1[0]<HoverHe1<p0HoverHe1[1]: passed = False
    if checkHoverHe2:  
        if not p0HoverHe2[0]<HoverHe2<p0HoverHe2[1]: passed = False
    if checkscaling:
        if not p0scaling[0]<Scaling<p0scaling[1]: passed = False
    if passed==False:  return -np_inf
    
    
    if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
        for anRV, anRVbound in zip(RV,used_RV_boundaries):  # go through all individual RVs and make sure they are in boundaries
            if not anRVbound[0] < anRV < anRVbound[1]:        return -np_inf
    else:
        if forced_K1=="Fit" and forced_K2=="Fit":
            if forced_Vgamma1=="Fit" and forced_Vgamma2=="Fit":
                if not (p0K1[0] < mcmc_K1 < p0K1[1] and p0K2[0] < mcmc_K2 < p0K2[1] and p0Vgamma1[0] < mcmc_Vgamma1 < p0Vgamma1[1] and p0Vgamma2[0] < mcmc_Vgamma2 < p0Vgamma2[1]):
                    return -np_inf
            else:
                if not (p0K1[0] < mcmc_K1 < p0K1[1] and p0K2[0] < mcmc_K2 < p0K2[1]):              return -np_inf
        elif forced_K1=="Fit":
            if forced_Vgamma1=="Fit":
                if not (p0K1[0] < mcmc_K1 < p0K1[1] and p0Vgamma1[0] < mcmc_Vgamma1 < p0Vgamma1[1]):    return -np_inf
            else:
                if not p0K1[0] < mcmc_K1 < p0K1[1]:       return -np_inf
        elif forced_K2=="Fit":
            if forced_Vgamma2=="Fit":
                if not (p0K2[0] < mcmc_K2 < p0K2[1] and p0Vgamma2[0] < mcmc_Vgamma2 < p0Vgamma2[1]):    return -np_inf
            else:
                if not p0K2[0] < mcmc_K2 < p0K2[1]:       return -np_inf
        
        if forced_K1!="Fit" and forced_Vgamma1=="Fit":
            if not p0Vgamma1[0] < mcmc_Vgamma1 < p0Vgamma1[1]:    return -np_inf
        if forced_K1!="Fit" and forced_Vgamma1=="Fit" and forced_K2!="Fit" and forced_Vgamma2=="Fit":
            if not (p0Vgamma1[0] < mcmc_Vgamma1 < p0Vgamma1[1] and p0Vgamma2[0] < mcmc_Vgamma2 < p0Vgamma2[1]):    return -np_inf
        if forced_K1!="Fit" and forced_Vgamma1!="Fit" and forced_K2!="Fit" and forced_Vgamma2=="Fit":
            if not p0Vgamma2[0] < mcmc_Vgamma2 < p0Vgamma2[1]:    return -np_inf
    
    if "Parallax" in p0labels: 
        parallax_prior = np.log(1.0/(npsqrt(2*np_pi)*plax_unc))-0.5*(mcmc_parallax-plax)**2/plax_unc**2 + np.log(1.0/(npsqrt(2*np_pi)*plax_unc))-0.5*(mcmc_parallax-plax)**2/plax_unc**2
        return parallax_prior
    return 0.0
    
    


def lnlike(theta, arguments):
    input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, normalised_wavelength, normalised_flux, normalised_err,  inp_resolution, used_RV_boundaries, HJD_values, sigma_clip, the_unique_wavelengths, spec1wl, spec1flux, spec2wl, spec2flux, inputScaling = arguments
    if type(input_files)==np.str_:
        input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, normalised_wavelength, normalised_flux, normalised_err,  inp_resolution, used_RV_boundaries, HJD_values, sigma_clip = [input_files], [share_rv], [reference_wl], [cut_Ha_all], [normaliseHa_all], [normalised_wavelength], [normalised_flux], [normalised_err],  [inp_resolution], [used_RV_boundaries], [HJD_values], [sigma_clip]
    
    if "T1" in p0labels:            args = npargwhere(p0labels=="T1")[0][0];          T1 = theta[args]
    else:                           T1 = forced_teff1
    
    if "logg1" in p0labels:         args = npargwhere(p0labels=="logg1")[0][0];       logg1 = theta[args]
    else:                           logg1 = forced_logg1
    
    if "T2" in p0labels:            args = npargwhere(p0labels=="T2")[0][0];          T2 = theta[args]
    else:                           T2 = forced_teff2
    
    if "logg2" in p0labels:         args = npargwhere(p0labels=="logg2")[0][0];       logg2 = theta[args]
    else:                           logg2 = forced_logg2
    
    if "H/He1" in p0labels:         args = npargwhere(p0labels=="H/He1")[0][0];       HoverHe1 = theta[args]
    else:
        try:    HoverHe1 = forced_HoverHe1
        except: None
    
    if "H/He2" in p0labels:         args = npargwhere(p0labels=="H/He2")[0][0];       HoverHe2 = theta[args]
    else:
        try:    HoverHe2 = forced_HoverHe2
        except: None
    
    if "K1" in p0labels:            args = npargwhere(p0labels=="K1")[0][0];          mcmc_K1 = theta[args]
    else:                           mcmc_K1 = forced_K1
    
    if "K2" in p0labels:            args = npargwhere(p0labels=="K2")[0][0];          mcmc_K2 = theta[args]
    else:                           mcmc_K2 = forced_K2
    
    if "Vg1" in p0labels:           args = npargwhere(p0labels=="Vg1")[0][0];         mcmc_Vgamma1 = theta[args]
    else:                           mcmc_Vgamma1 = forced_Vgamma1
    
    if "Vg2" in p0labels:           args = npargwhere(p0labels=="Vg2")[0][0];         mcmc_Vgamma2 = theta[args]
    else:                           mcmc_Vgamma2 = forced_Vgamma2
    
    if "Parallax" in p0labels:      args = npargwhere(p0labels=="Parallax")[0][0];    mcmc_parallax = theta[args]
    
    if "Scaling" in p0labels:       args = npargwhere(p0labels=="Scaling")[0][0];     Scaling = theta[args]
    else:                           Scaling = forced_Scaling
    
    if "RV1_0" in p0labels:         num_start_RVs = npargwhere(p0labels=="RV1_0")[0][0]

    
    for starcnt, starType in enumerate([starType1, starType2]):
        if starType=="DA":
            if not pier_or_antoine=="pier3Dphot_antoine1Dspec":
                if starcnt==0:     Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DAgrids(T1, logg1)
                elif starcnt==1:   Grav2_N, wl_all2_N, flux2_N, Teff2_N = return_DAgrids(T2, logg2)
            else:
                if starcnt==0:
                    Grav1_N, wl_all1_N, flux1_N, Teff1_N, photGrav1_N, photwl_all1_N, photflux1_N, photTeff1_N = return_DAgrids(T1, logg1)
                    if len(wl_all1_N)==0:   raise ValueError(len(wl_all1_N), len(Grav1_N), len(flux1_N), len(Teff1_N), T1, logg1)
                elif starcnt==1:
                    Grav2_N, wl_all2_N, flux2_N, Teff2_N, photGrav2_N, photwl_all2_N, photflux2_N, photTeff2_N = return_DAgrids(T2, logg2)
                    if len(wl_all2_N)==0:   raise ValueError(len(wl_all2_N), len(Grav2_N), len(flux2_N), len(Teff2_N), T2, logg2)
            

        elif starType=="DBA":
            if starcnt==0:    Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N = return_DBAgrids(T1, logg1, HoverHe1)
            elif starcnt==1:  Grav2_N, wl_all2_N, flux2_N, Teff2_N, HoverHe2_N = return_DBAgrids(T2, logg2, HoverHe2)
        
        elif starType1=="DB":
            if starcnt==0:    Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DBgrids(T1, logg1)
            elif starcnt==1:  Grav2_N, wl_all2_N, flux2_N, Teff2_N = return_DBgrids(T2, logg2)
        
        elif starType=="DC":
            if starcnt==0:    Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DCgrids(T1, logg1)
            elif starcnt==1:  Grav2_N, wl_all2_N, flux2_N, Teff2_N = return_DCgrids(T2, logg2)

    
    if forced_Scaling=="WD" and inputScaling==None:
        if starType1=="DA":     R1 = get_MTR(T1, logg=logg1, return_R_from_T_logg=True, loaded_Istrate=loaded_Istrate, loaded_CO=loaded_CO, loaded_Althaus=loaded_Althaus)
        elif starType1=="DBA" or starType1=="DB" or starType1=="DC":
            R1 = get_MTR_DB(T1, logg=logg1)
            #### This only applies if Camisassa grid is being used
            #if logg1<7.7529337 or logg1>8.6700679: raise ValueError("logg1=", logg2, "Min/max Bounds of logg for DB scaling is 7.7796246 and 8.6700679. If your logg is close to this, consider restricting p0logg further or turning off WD scaling")
        
        if starType2=="DA":     R2 = get_MTR(T2, logg=logg2, return_R_from_T_logg=True, loaded_Istrate=loaded_Istrate, loaded_CO=loaded_CO, loaded_Althaus=loaded_Althaus)
        elif starType2=="DBA"  or starType2=="DB" or starType2=="DC":
            R2 = get_MTR_DB(T2, logg=logg2)
            #### This only applies if Camisassa grid is being used
            #if logg2<7.7529337 or logg2>8.6700679: raise ValueError("logg2=", logg2, "Min/max Bounds of logg for DB scaling is 7.7796246 and 8.6700679. If your logg is close to this, consider restricting p0logg further or turning off WD scaling")
        Scaling= (R1**2/R2**2)
        
        if np.isnan(R1): raise ValueError(T1, logg1)
        if np.isnan(R2): raise ValueError(T2, logg2)
    elif isinstance(forced_Scaling, float): Scaling=forced_Scaling
    elif inputScaling!=None:
        Scaling=inputScaling
    
    
    if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
        if forced_Scaling==False: RV=theta[num_start_RVs:-1];  Scaling=theta[-1]
        else:                     RV=theta[num_start_RVs:]
    else:
        if forced_Scaling==False: Scaling=theta[-1]
    
    
    
    if fit_phot_SED:
        if starType1=="DA" and starType2=="DA":
            if pier_or_antoine=="pier3Dphot_antoine1Dspec":
                smeared_wl, smeared_flux, star1flux, star2flux = Fit_phot.fit_phot_SED_double(photGrav1_N, photwl_all1_N, photflux1_N, photTeff1_N, None, photGrav2_N, photwl_all2_N, photflux2_N, photTeff2_N, None, T1, logg1, None, T2, logg2, None, min_wl=theminww, max_wl=themaxww, starType1=starType1, starType2=starType2, R1=R1, R2=R2, parallax=mcmc_parallax, red=reddening_Ebv, return_indiv_stars="forSpectrum")
            else:
                try:
                    smeared_wl, smeared_flux, star1flux, star2flux = Fit_phot.fit_phot_SED_double(Grav1_N, wl_all1_N, flux1_N, Teff1_N, None, Grav2_N, wl_all2_N, flux2_N, Teff2_N, None, T1, logg1, None, T2, logg2, None, min_wl=theminww, max_wl=themaxww, starType1=starType1, starType2=starType2, R1=R1, R2=R2, parallax=mcmc_parallax, red=reddening_Ebv, return_indiv_stars="forSpectrum")
                except: raise ValueError(T1, logg1, T2, logg2, len(Teff1_N), len(Grav1_N), len(Teff2_N), len(Grav2_N))
        elif starType1=="DA" and starType2=="DBA":
            smeared_wl, smeared_flux, star1flux, star2flux = Fit_phot.fit_phot_SED_double(Grav1_N, wl_all1_N, flux1_N, Teff1_N, None, Grav2_N, wl_all2_N, flux2_N, Teff2_N, HoverHe2_N, T1, logg1, None, T2, logg2, HoverHe2, min_wl=theminww, max_wl=themaxww,starType1=starType1, starType2=starType2, R1=R1, R2=R2, parallax=mcmc_parallax, red=reddening_Ebv, return_indiv_stars="forSpectrum")
        elif starType1=="DBA" and starType2=="DA":
            smeared_wl, smeared_flux, star1flux, star2flux = Fit_phot.fit_phot_SED_double(Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N, Grav2_N, wl_all2_N, flux2_N, Teff2_N, None, T1, logg1, HoverHe1, T2, logg2, None, min_wl=theminww, max_wl=themaxww,starType1=starType1, starType2=starType2, R1=R1, R2=R2, parallax=mcmc_parallax, red=reddening_Ebv, return_indiv_stars="forSpectrum")
        elif starType1=="DBA" and starType2=="DBA":
            smeared_wl, smeared_flux, star1flux, star2flux = Fit_phot.fit_phot_SED_double(Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N, Grav2_N, wl_all2_N, flux2_N, Teff2_N, HoverHe2_N, T1, logg1, HoverHe1, T2, logg2, HoverHe2, min_wl=theminww, max_wl=themaxww,starType1=starType1, starType2=starType2, R1=R1, R2=R2, parallax=mcmc_parallax, red=reddening_Ebv, return_indiv_stars="forSpectrum")
        elif (starType1=="DA" and starType2=="DC") or (starType1=="DC" and starType2=="DA") or (starType1=="DA" and starType2=="DB") or (starType1=="DB" and starType2=="DA"):
            smeared_wl, smeared_flux, star1flux, star2flux = Fit_phot.fit_phot_SED_double(Grav1_N, wl_all1_N, flux1_N, Teff1_N, None, Grav2_N, wl_all2_N, flux2_N, Teff2_N, None, T1, logg1, None, T2, logg2, None, min_wl=theminww, max_wl=themaxww, starType1=starType1, starType2=starType2, R1=R1, R2=R2, parallax=mcmc_parallax, red=reddening_Ebv, return_indiv_stars="forSpectrum")
            
            
        
        rchisq_phot, chisq_phot = Fit_phot.process_photometry_in_each_pb(smeared_wl, smeared_flux, sedfilter, sed_wl, sedflux, sedfluxe, filter_dict=filter_dict, theminww_plot=theminww+50, themaxww_plot=themaxww+50, return_points_for_phot_model=False)
        
        #chisq_phot*=1000
        
        if np.round(chisq_phot,11)==0:
            plt.plot(smeared_wl, smeared_flux);   plt.title(str(T1) + "  "+ str(logg1) + "  "+ str(T2) + "  "+ str(logg2) + "   " + "chisq=" +str(chisq_phot));   plt.show()
            raise ValueError
    
    
    
    
    chisq_spec=0
    
    index = np.arange(0,len(input_files),1)
    
    if fit_phot_SED:
        no_need_to_recompute = npamin(smeared_wl) < npamin(the_unique_wavelengths)+npamin(cut_Ha_all.T) - 10 and npamax(smeared_wl) > npamax(the_unique_wavelengths)+npamax(cut_Ha_all.T) + 10
    else:   no_need_to_recompute = False
    
    for ref_wl in the_unique_wavelengths:
        mask_ref_wl = reference_wl==ref_wl
        ref_wl_cut_lims_min, ref_wl_cut_lims_max = cut_Ha_all[mask_ref_wl].T
        cut_wl_min, cut_wl_max = npamin(ref_wl_cut_lims_min), npamax(ref_wl_cut_lims_max)
        
        
        
        if no_need_to_recompute:
            if fit_phot_SED:
                mask =  (smeared_wl>ref_wl+cut_wl_min-10) & (smeared_wl<ref_wl+cut_wl_max+10)
                model_wl1, model_wl2, model_spectrum_star1, model_spectrum_star2 = smeared_wl[mask], smeared_wl[mask], star1flux[mask], star2flux[mask]
            elif spec1flux!=None and spec2flux!=None:
                mask1 =  (spec1wl>ref_wl+cut_wl_min-10) & (spec1wl<ref_wl+cut_wl_max+10)
                mask2 =  (spec2wl>ref_wl+cut_wl_min-10) & (spec2wl<ref_wl+cut_wl_max+10)
                model_wl1, model_wl2, model_spectrum_star1, model_spectrum_star2 = spec1wl[mask1], spec2wl[mask2], spec1flux[mask1], star2flux[mask2]

        
        else:
            if starType1=="DA" or starType1=="DC" or starType1=="DB":  model_wl1, model_spectrum_star1 = return_model_spectrum_DA(wl_all1_N, ref_wl, cut_wl_min, cut_wl_max, Grav1_N, flux1_N, Teff1_N, T1, logg1)
            elif starType1=="DBA": model_wl1, model_spectrum_star1 = return_model_spectrum_DBA(wl_all1_N, ref_wl, cut_wl_min, cut_wl_max, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, T1, logg1, HoverHe1)
            
            if starType2=="DA" or starType2=="DC" or starType2=="DB":  model_wl2, model_spectrum_star2 = return_model_spectrum_DA(wl_all2_N, ref_wl, cut_wl_min, cut_wl_max, Grav2_N, flux2_N, Teff2_N, T2, logg2)
            elif starType2=="DBA": model_wl2, model_spectrum_star2 = return_model_spectrum_DBA(wl_all2_N, ref_wl, cut_wl_min, cut_wl_max, Grav2_N, flux2_N, Teff2_N, HoverHe2_N, T2, logg2, HoverHe2)
        
        
        
        
        
        
        
        
        
        # go through all files one by one and do the same as above, but recalculating the interpolated spectrum, which makes things slower
        for (ii, in_fi, sh_rv, modha, cut_lim, norm_lim, aHJD, sigclipspec) in zip(index[mask_ref_wl], input_files[mask_ref_wl], share_rv[mask_ref_wl], modelHa[mask_ref_wl], cut_Ha_all[mask_ref_wl], normaliseHa_all[mask_ref_wl], HJD_values[mask_ref_wl], sigma_clip[mask_ref_wl]):  # the mcmc RVs are gone through 1 by one, uses the ii counter.
            
            
            if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
                if sh_rv == -1:  mcmc_rv1, mcmc_rv2 = RV[2*ii], RV[2*ii+1]
                else:            mcmc_rv1, mcmc_rv2 = RV[2*sh_rv], RV[2*sh_rv+1]
            else:
                phi = ((aHJD - forced_T0)%forced_P0)/forced_P0
                if forced_K1=="Fit":     
                    if forced_Vgamma1=="Fit":    mcmc_rv1 = mcmc_Vgamma1 + mcmc_K1*np_sin(2*np_pi*phi)
                    else:    mcmc_rv1 = forced_Vgamma1 + mcmc_K1*np_sin(2*np_pi*phi)
                else:                   
                    if forced_Vgamma1=="Fit":    mcmc_rv1 = mcmc_Vgamma1 + forced_K1*np_sin(2*np_pi*phi)
                    else: mcmc_rv1 = forced_Vgamma1 + forced_K1*np_sin(2*np_pi*phi)
                if forced_K2=="Fit":   
                    if forced_Vgamma2=="Fit":    mcmc_rv2 = mcmc_Vgamma2 + mcmc_K2*np_sin(2*np_pi*phi)
                    else:    mcmc_rv2 = forced_Vgamma2 + mcmc_K2*np_sin(2*np_pi*phi)
                else:   
                    if forced_Vgamma1=="Fit":    mcmc_rv2 = mcmc_Vgamma2 + forced_K2*np_sin(2*np_pi*phi)
                    else:    mcmc_rv2 = forced_Vgamma2 + forced_K2*np_sin(2*np_pi*phi)
                if forced_K1!="Fit" and forced_Vgamma1=="Fit":     mcmc_rv1 = mcmc_Vgamma1 + forced_K1*np_sin(2*np_pi*phi)
                if forced_K2!="Fit" and forced_Vgamma2=="Fit":     mcmc_rv2 = mcmc_Vgamma2 + forced_K2*np_sin(2*np_pi*phi)
                    
            
            
            
            normalised_wavelength = list_norm_wl_grids[ii];    normalised_flux = list_normalised_flux[ii];    normalised_err = list_normalised_err[ii];    inp_resolution = resolutions[ii]
            
            
            dlam1 = ref_wl*mcmc_rv1/speed_of_light
            dlam2 = ref_wl*mcmc_rv2/speed_of_light

            if high_RV_amp: modelHa_min+=dlam1; modelHa_max+=dlam1; cut_limits_min+=dlam1; cut_limits_max+=dlam1; norm_limits_min+=dlam1; norm_limits_max+=dlam1
            
            # interpolate the model onto the observation's wavelength grid
            interparr = Scaling*interp(model_wl2, model_wl1+dlam1, model_spectrum_star1)    + interp(model_wl2, model_wl2+dlam2, model_spectrum_star2)
            
            
            # Smear to the resolution desired
            resstd=0.5*(ref_wl/inp_resolution)/(model_wl2[10]-model_wl2[9])
            interparr = convolve(interparr, Gaussian1DKernel(stddev=resstd), boundary = 'extend')
            
            
            interparr = interp(normalised_wavelength, model_wl2, interparr)
                
            
            
            
            mask = ((normalised_wavelength>ref_wl+cut_lim[0]) & (normalised_wavelength<ref_wl+norm_lim[0])) | ((normalised_wavelength>ref_wl+norm_lim[1])   & (normalised_wavelength<ref_wl+cut_lim[1]))
            
            
            try:
                m,c =  polyfit(normalised_wavelength[mask], interparr[mask], deg=1)
                interparr /= (m*normalised_wavelength + c)
            except: return -np_inf
            
            
            if False:
                mask_rv=True
                lam_rv_remove = ref_wl*50/speed_of_light
                if mask_rv:   aaamask = (normalised_wavelength > ref_wl + dlam2 + lam_rv_remove) | (normalised_wavelength < ref_wl + dlam2 - lam_rv_remove)
                else:         aaamask = normalised_wavelength==normalised_wavelength
            else:   aaamask = normalised_wavelength==normalised_wavelength
            
            
            
            desired_range = (normalised_wavelength>ref_wl+modha[0])  &  (normalised_wavelength<ref_wl+modha[1]) & aaamask
            if sigclipspec!=-1:
                resid = normalised_flux-(interparr)
                resid_in_sig = resid/np.std(np.abs(resid[desired_range]))
                
                clip_mask = (np.abs(resid_in_sig[desired_range])>sigclipspec) & ((normalised_wavelength[desired_range] > ref_wl+5) | (normalised_wavelength[desired_range] < ref_wl-5))
            else: clip_mask=normalised_flux[desired_range]!=normalised_flux[desired_range]
            
            
            
            off = polyfit(normalised_wavelength[desired_range][~clip_mask], normalised_flux[desired_range][~clip_mask] - interparr[desired_range][~clip_mask], w=1/normalised_err[desired_range][~clip_mask], deg=0)[0]
            

            ## Now on to calculating chisq
            chisq_indiv = -0.5*npsum((np_square(normalised_flux[desired_range & aaamask][~clip_mask]-(off+interparr[desired_range & aaamask][~clip_mask])))/np_square(normalised_err[desired_range & aaamask][~clip_mask]))
            chisq_spec += chisq_indiv
            
            
    
    if fit_phot_SED:
        #if np.isnan(chisq_spec): plt.plot(normalised_wavelength, interparr1); plt.plot(normalised_wavelength, interparr2); plt.title(str(T1) + "  "+ str(logg1)+"  "+ str(T2)+"  "+ str(logg2));  plt.show()
        if np.isnan(chisq_spec) == True or np.isnan(chisq_phot) == True:    return -np_inf
        else:  
            if rank==1: print(chisq_spec,chisq_phot)
            return chisq_spec + chisq_phot
    else:
        if np.isnan(chisq_spec) == True :    return -np_inf
        else:   return chisq_spec
        
        


def lnprob(theta, arguments):
    lp = lnprior(theta, arguments)
    if not np.isfinite(lp):
        return -np_inf
    return lp + lnlike(theta, arguments)









#Test that next step falls within boundaries allowed by the priors
def lnprior_gauss_lorentz(theta, arguments):
    input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, normalised_wavelength, normalised_flux, normalised_err,  inp_resolution, used_RV_boundaries, HJD_values, sigma_clip, the_unique_wavelengths, spec1wl, spec1flux, spec2wl, spec2flux, inputScaling = arguments
    if "K1" in p0labels:            args = npargwhere(p0labels=="K1")[0][0];          mcmc_K1 = theta[args]
    if "K2" in p0labels:            args = npargwhere(p0labels=="K2")[0][0];          mcmc_K2 = theta[args]
    if "Vg1" in p0labels:           args = npargwhere(p0labels=="Vg1")[0][0];         mcmc_Vgamma1 = theta[args]
    if "Vg2" in p0labels:           args = npargwhere(p0labels=="Vg2")[0][0];         mcmc_Vgamma2 = theta[args]
    if "RV1_0" in p0labels:         num_start_RVs = npargwhere(p0labels=="RV1_0")[0][0]
    elif "RV1_1" in p0labels:       num_start_RVs = npargwhere(p0labels=="RV1_1")[0][0]
    elif "RV1_2" in p0labels:       num_start_RVs = npargwhere(p0labels=="RV1_2")[0][0]
    elif "RV1_3" in p0labels:       num_start_RVs = npargwhere(p0labels=="RV1_3")[0][0]
    elif "RV1_4" in p0labels:       num_start_RVs = npargwhere(p0labels=="RV1_4")[0][0]
    
    args = npargwhere(p0labels=="A1_1")[0][0];   gaussLorentz1_A1 = theta[args]
    args = npargwhere(p0labels=="A1_2")[0][0];   gaussLorentz1_A2 = theta[args]
    args = npargwhere(p0labels=="A2_1")[0][0];   gaussLorentz2_A1 = theta[args]
    args = npargwhere(p0labels=="A2_2")[0][0];   gaussLorentz2_A2 = theta[args]
    args = npargwhere(p0labels=="sig1_1")[0][0];   gaussLorentz1_sigma1 = theta[args]
    args = npargwhere(p0labels=="sig1_2")[0][0];   gaussLorentz1_sigma2 = theta[args]
    args = npargwhere(p0labels=="sig2_1")[0][0];   gaussLorentz2_sigma1 = theta[args]
    args = npargwhere(p0labels=="sig2_2")[0][0];   gaussLorentz2_sigma2 = theta[args]
    
    
    if not (A1_1_boundaries[0] < gaussLorentz1_A1 < A1_1_boundaries[1]):
        return -np_inf
    if not (A1_2_boundaries[0] < gaussLorentz1_A2 < A1_2_boundaries[1]):
        return -np_inf
    if not (A2_1_boundaries[0] < gaussLorentz2_A1 < A2_1_boundaries[1]):
        return -np_inf
    if not (A2_2_boundaries[0] < gaussLorentz2_A2 < A2_2_boundaries[1]):
        return -np_inf
    if not (sigma1_1_boundaries[0] < gaussLorentz1_sigma1 < sigma1_1_boundaries[1]):
        return -np_inf
    if not (sigma1_2_boundaries[0] < gaussLorentz1_sigma2 < sigma1_2_boundaries[1]):
        return -np_inf
    if not (sigma2_1_boundaries[0] < gaussLorentz2_sigma1 < sigma2_1_boundaries[1]):
        return -np_inf
    if not (sigma2_2_boundaries[0] < gaussLorentz2_sigma2 < sigma2_2_boundaries[1]):
        return -np_inf
    
    
    
    if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
        RV=theta[num_start_RVs:]
        
        
    
    if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
        for anRV, anRVbound in zip(RV,used_RV_boundaries):  # go through all individual RVs and make sure they are in boundaries
            if not anRVbound[0] < anRV < anRVbound[1]:        return -np_inf
    else:
        if forced_K1=="Fit" and forced_K2=="Fit":
            if forced_Vgamma1=="Fit" and forced_Vgamma2=="Fit":
                if not (p0K1[0] < mcmc_K1 < p0K1[1] and p0K2[0] < mcmc_K2 < p0K2[1] and p0Vgamma1[0] < mcmc_Vgamma1 < p0Vgamma1[1] and p0Vgamma2[0] < mcmc_Vgamma2 < p0Vgamma2[1]):
                    return -np_inf
            else:
                if not (p0K1[0] < mcmc_K1 < p0K1[1] and p0K2[0] < mcmc_K2 < p0K2[1]):              return -np_inf
        elif forced_K1=="Fit":
            if forced_Vgamma1=="Fit":
                if not (p0K1[0] < mcmc_K1 < p0K1[1] and p0Vgamma1[0] < mcmc_Vgamma1 < p0Vgamma1[1]):    return -np_inf
            else:
                if not p0K1[0] < mcmc_K1 < p0K1[1]:       return -np_inf
        elif forced_K2=="Fit":
            if forced_Vgamma2=="Fit":
                if not (p0K2[0] < mcmc_K2 < p0K2[1] and p0Vgamma2[0] < mcmc_Vgamma2 < p0Vgamma2[1]):    return -np_inf
            else:
                if not p0K2[0] < mcmc_K2 < p0K2[1]:       return -np_inf
        
        if forced_K1!="Fit" and forced_Vgamma1=="Fit":
            if not p0Vgamma1[0] < mcmc_Vgamma1 < p0Vgamma1[1]:    return -np_inf
        if forced_K1!="Fit" and forced_Vgamma1=="Fit" and forced_K2!="Fit" and forced_Vgamma2=="Fit":
            if not (p0Vgamma1[0] < mcmc_Vgamma1 < p0Vgamma1[1] and p0Vgamma2[0] < mcmc_Vgamma2 < p0Vgamma2[1]):    return -np_inf
        if forced_K1!="Fit" and forced_Vgamma1!="Fit" and forced_K2!="Fit" and forced_Vgamma2=="Fit":
            if not p0Vgamma2[0] < mcmc_Vgamma2 < p0Vgamma2[1]:    return -np_inf
    
    return 0.0
    
    

@njit
def gauss_2x_mcmc(x, a1, sigma1, a2, sigma2):#  no offset 'b' because I fit offset the observations to best match the profile anyway. Ignoring b1 b2 also means that b1=b2=0.  No need to include a mean here
    return a1*np.exp(-0.5 * np_square(x/sigma1))   +   a2*np.exp(-0.5 * np_square(x/sigma2))

@njit
def lorentz_2x_mcmc(x, a1, sigma1, a2, sigma2):#  no offset 'b' because I fit offset the observations to best match the profile anyway. Ignoring b1 b2 also means that b1=b2=0.  No need to include a mean here
    return (a1/np_pi)  *  (sigma1/(np_square(x) + np_square(sigma1)))  +  (a2/np_pi)  *  (sigma2/(np_square(x) + np_square(sigma2)))

@njit
def gauss_lorentz_mcmc(x, a1, sigma1, a2, sigma2):#  no offset 'b' because I fit offset the observations to best match the profile anyway. Ignoring b1 b2 also means that b1=b2=0.  No need to include a mean here
    return a1*np.exp(-0.5 * np_square(x/sigma1))  +  (a2/np_pi)  *  (sigma2/(np_square(x) + np_square(sigma2)))
    




def lnlike_gauss_lorentz(theta, arguments):
    input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, normalised_wavelength, normalised_flux, normalised_err,  inp_resolution, used_RV_boundaries, HJD_values, sigma_clip, the_unique_wavelengths, spec1wl, spec1flux, spec2wl, spec2flux, inputScaling = arguments
    if type(input_files)==np.str_:
        input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, normalised_wavelength, normalised_flux, normalised_err,  inp_resolution, used_RV_boundaries, HJD_values, sigma_clip = [input_files], [share_rv], [reference_wl], [cut_Ha_all], [normaliseHa_all], [normalised_wavelength], [normalised_flux], [normalised_err],  [inp_resolution], [used_RV_boundaries], [HJD_values], [sigma_clip]
    
    args = npargwhere(p0labels=="A1_1")[0][0];   gaussLorentz1_A1 = theta[args]
    args = npargwhere(p0labels=="A1_2")[0][0];   gaussLorentz1_A2 = theta[args]
    args = npargwhere(p0labels=="A2_1")[0][0];   gaussLorentz2_A1 = theta[args]
    args = npargwhere(p0labels=="A2_2")[0][0];   gaussLorentz2_A2 = theta[args]
    args = npargwhere(p0labels=="sig1_1")[0][0];   gaussLorentz1_sigma1 = theta[args]
    args = npargwhere(p0labels=="sig1_2")[0][0];   gaussLorentz1_sigma2 = theta[args]
    args = npargwhere(p0labels=="sig2_1")[0][0];   gaussLorentz2_sigma1 = theta[args]
    args = npargwhere(p0labels=="sig2_2")[0][0];   gaussLorentz2_sigma2 = theta[args]

    
    if "K1" in p0labels:            args = npargwhere(p0labels=="K1")[0][0];          mcmc_K1 = theta[args]
    else:                           mcmc_K1 = forced_K1
    
    if "K2" in p0labels:            args = npargwhere(p0labels=="K2")[0][0];          mcmc_K2 = theta[args]
    else:                           mcmc_K2 = forced_K2
    
    if "Vg1" in p0labels:           args = npargwhere(p0labels=="Vg1")[0][0];         mcmc_Vgamma1 = theta[args]
    else:                           mcmc_Vgamma1 = forced_Vgamma1
    
    if "Vg2" in p0labels:           args = npargwhere(p0labels=="Vg2")[0][0];         mcmc_Vgamma2 = theta[args]
    else:                           mcmc_Vgamma2 = forced_Vgamma2
    
    if "RV1_0" in p0labels:         num_start_RVs = npargwhere(p0labels=="RV1_0")[0][0]
    

    
    if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
        RV=theta[num_start_RVs:]
    
    
    
    chisq_spec=0
    
    index = np.arange(0,len(input_files),1)
    
    the_unique_wavelengths=[npamax(the_unique_wavelengths)]
    for ref_wl in the_unique_wavelengths:
        mask_ref_wl = reference_wl==ref_wl
        ref_wl_cut_lims_min, ref_wl_cut_lims_max = cut_Ha_all[mask_ref_wl].T
        cut_wl_min, cut_wl_max = npamin(ref_wl_cut_lims_min), npamax(ref_wl_cut_lims_max)
        model_wl1 = nplinspace(cut_wl_min, cut_wl_max, int(cut_wl_max-cut_wl_min)*10)
        model_wl2 = nplinspace(cut_wl_min, cut_wl_max, int(cut_wl_max-cut_wl_min)*10)
        
        
        model_spectrum_star1 = lorentz_2x_mcmc(model_wl1, gaussLorentz1_A1, gaussLorentz1_sigma1, 0, 1)
        model_spectrum_star2 = lorentz_2x_mcmc(model_wl2, gaussLorentz2_A1, gaussLorentz2_sigma1, 0, 1)
        
        
        if starType1=="GG":     model_spectrum_star1 = gauss_2x_mcmc(model_wl1, gaussLorentz1_A1, gaussLorentz1_sigma1, gaussLorentz1_A2, gaussLorentz1_sigma2)
        elif starType1=="LL":   model_spectrum_star1 = lorentz_2x_mcmc(model_wl1, gaussLorentz1_A1, gaussLorentz1_sigma1, gaussLorentz1_A2, gaussLorentz1_sigma2)
        elif starType1=="GL":   model_spectrum_star1 = gauss_lorentz_mcmc(model_wl1, gaussLorentz1_A1, gaussLorentz1_sigma1, gaussLorentz1_A2, gaussLorentz1_sigma2)
        if starType2=="GG":     model_spectrum_star2 = gauss_2x_mcmc(model_wl2, gaussLorentz2_A1, gaussLorentz2_sigma1, gaussLorentz2_A2, gaussLorentz2_sigma2)
        elif starType2=="LL":   model_spectrum_star2 = lorentz_2x_mcmc(model_wl2, gaussLorentz2_A1, gaussLorentz2_sigma1, gaussLorentz2_A2, gaussLorentz2_sigma2)
        elif starType2=="GL":   model_spectrum_star2 = gauss_lorentz_mcmc(model_wl2, gaussLorentz2_A1, gaussLorentz2_sigma1, gaussLorentz2_A2, gaussLorentz2_sigma2)
        
        
        model_wl1, model_wl2 = model_wl1 + ref_wl, model_wl2 + ref_wl
        
        
        
        
        # go through all files one by one and do the same as above, but recalculating the interpolated spectrum, which makes things slower
        for (ii, in_fi, sh_rv, modha, cut_lim, norm_lim, aHJD, sigclipspec) in zip(index[mask_ref_wl], input_files[mask_ref_wl], share_rv[mask_ref_wl], modelHa[mask_ref_wl], cut_Ha_all[mask_ref_wl], normaliseHa_all[mask_ref_wl], HJD_values[mask_ref_wl], sigma_clip[mask_ref_wl]):  # the mcmc RVs are gone through 1 by one, uses the ii counter.
            modelHa_min, modelHa_max = modha
            cut_limits_min, cut_limits_max = cut_lim
            norm_limits_min, norm_limits_max = norm_lim
            
            
            
            if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
                if sh_rv == -1:  mcmc_rv1, mcmc_rv2 = RV[2*ii], RV[2*ii+1]
                else:            mcmc_rv1, mcmc_rv2 = RV[2*sh_rv], RV[2*sh_rv+1]
            else:
                phi = ((aHJD - forced_T0)%forced_P0)/forced_P0
                if forced_K1=="Fit":     
                    if forced_Vgamma1=="Fit":    mcmc_rv1 = mcmc_Vgamma1 + mcmc_K1*np_sin(2*np_pi*phi)
                    else:    mcmc_rv1 = forced_Vgamma1 + mcmc_K1*np_sin(2*np_pi*phi)
                else:                   
                    if forced_Vgamma1=="Fit":    mcmc_rv1 = mcmc_Vgamma1 + forced_K1*np_sin(2*np_pi*phi)
                    else: mcmc_rv1 = forced_Vgamma1 + forced_K1*np_sin(2*np_pi*phi)
                if forced_K2=="Fit":   
                    if forced_Vgamma2=="Fit":    mcmc_rv2 = mcmc_Vgamma2 + mcmc_K2*np_sin(2*np_pi*phi)
                    else:    mcmc_rv2 = forced_Vgamma2 + mcmc_K2*np_sin(2*np_pi*phi)
                else:   
                    if forced_Vgamma1=="Fit":    mcmc_rv2 = mcmc_Vgamma2 + forced_K2*np_sin(2*np_pi*phi)
                    else:    mcmc_rv2 = forced_Vgamma2 + forced_K2*np_sin(2*np_pi*phi)
                if forced_K1!="Fit" and forced_Vgamma1=="Fit":     mcmc_rv1 = mcmc_Vgamma1 + forced_K1*np_sin(2*np_pi*phi)
                if forced_K2!="Fit" and forced_Vgamma2=="Fit":     mcmc_rv2 = mcmc_Vgamma2 + forced_K2*np_sin(2*np_pi*phi)
                    
            
            
            
            normalised_wavelength = list_norm_wl_grids[ii];    normalised_flux = list_normalised_flux[ii];    normalised_err = list_normalised_err[ii];    inp_resolution = resolutions[ii]
            
            
            dlam1 = ref_wl*mcmc_rv1/speed_of_light
            dlam2 = ref_wl*mcmc_rv2/speed_of_light
            
            if high_RV_amp: modelHa_min+=dlam1; modelHa_max+=dlam1; cut_limits_min+=dlam1; cut_limits_max+=dlam1; norm_limits_min+=dlam1; norm_limits_max+=dlam1

            
            
            # interpolate the model onto the observation's wavelength grid
            interparr = interp(model_wl2, model_wl1+dlam1, model_spectrum_star1)    +   interp(model_wl2, model_wl2+dlam2, model_spectrum_star2)   +   1
            
            
            # Smear to the resolution desired
            resstd=0.5*(ref_wl/inp_resolution)/(model_wl2[10]-model_wl2[9])
            interparr = convolve(interparr, Gaussian1DKernel(stddev=resstd), boundary = 'extend')
            interparr = interp(normalised_wavelength, model_wl2, interparr)
                
            
            
            mask = ((normalised_wavelength>ref_wl+cut_limits_min) & (normalised_wavelength<ref_wl+norm_limits_min)) | ((normalised_wavelength>ref_wl+norm_limits_max)   & (normalised_wavelength<ref_wl+cut_limits_max))
            
            
            try:
                m,c =  polyfit(normalised_wavelength[mask], interparr[mask], deg=1)
                interparr /= (m*normalised_wavelength + c)
            except: return -np_inf
            
            
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
            chisq_indiv = -0.5*npsum((np_square(normalised_flux[desired_range & aaamask][~clip_mask]-(off+interparr[desired_range & aaamask][~clip_mask])))/np_square(normalised_err[desired_range & aaamask][~clip_mask]))
            chisq_spec += chisq_indiv
            
            
    
    
    if np.isnan(chisq_spec) == True :    return -np_inf
    else:   return chisq_spec
        
        


def lnprob_gauss_lorentz(theta, arguments):
    lp = lnprior_gauss_lorentz(theta, arguments)
    if not np.isfinite(lp):
        return -np_inf
    return lp + lnlike_gauss_lorentz(theta, arguments)












#Test that next step falls within boundaries allowed by the priors
def lnprior_quad_lorentz(theta, arguments):
    input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, normalised_wavelength, normalised_flux, normalised_err,  inp_resolution, used_RV_boundaries, HJD_values, sigma_clip, the_unique_wavelengths, spec1wl, spec1flux, spec2wl, spec2flux, inputScaling = arguments
    if "K1" in p0labels:            args = npargwhere(p0labels=="K1")[0][0];          mcmc_K1 = theta[args]
    if "K2" in p0labels:            args = npargwhere(p0labels=="K2")[0][0];          mcmc_K2 = theta[args]
    if "Vg1" in p0labels:           args = npargwhere(p0labels=="Vg1")[0][0];         mcmc_Vgamma1 = theta[args]
    if "Vg2" in p0labels:           args = npargwhere(p0labels=="Vg2")[0][0];         mcmc_Vgamma2 = theta[args]
    if "RV1_0" in p0labels:         num_start_RVs = npargwhere(p0labels=="RV1_0")[0][0]
    elif "RV1_1" in p0labels:       num_start_RVs = npargwhere(p0labels=="RV1_1")[0][0]
    elif "RV1_2" in p0labels:       num_start_RVs = npargwhere(p0labels=="RV1_2")[0][0]
    elif "RV1_3" in p0labels:       num_start_RVs = npargwhere(p0labels=="RV1_3")[0][0]
    elif "RV1_4" in p0labels:       num_start_RVs = npargwhere(p0labels=="RV1_4")[0][0]
    
    
    if not "run_double_one_by_oneMCMC" in sys.argv:
        args = npargwhere(p0labels=="A1_1")[0][0];   gaussLorentz1_A1 = theta[args]
        args = npargwhere(p0labels=="A2_1")[0][0];   gaussLorentz2_A1 = theta[args]
        args = npargwhere(p0labels=="sig1_1")[0][0];   gaussLorentz1_sigma1 = theta[args]
        args = npargwhere(p0labels=="sig2_1")[0][0];   gaussLorentz2_sigma1 = theta[args]
        if not (A1_1_boundaries[0] < gaussLorentz1_A1 < A1_1_boundaries[1]):
            return -np_inf
        if not (A2_1_boundaries[0] < gaussLorentz2_A1 < A2_1_boundaries[1]):
            return -np_inf
        if not (sigma1_1_boundaries[0] < gaussLorentz1_sigma1 < sigma1_1_boundaries[1]):
            return -np_inf
        if not (sigma2_1_boundaries[0] < gaussLorentz2_sigma1 < sigma2_1_boundaries[1]):
            return -np_inf

    else:
        with open("out/result.out") as resultfile:
            result_lines = resultfile.readlines()
            for iii, lll in enumerate(result_lines):
                if "A1_1:" in lll:          gaussLorentz1_A1=float(result_lines[iii+1].split("\t")[0])
                elif "sig1_1:" in lll:      gaussLorentz1_sigma1=float(result_lines[iii+1].split("\t")[0])
                elif "A1_2:" in lll:        A1_2_med=float(result_lines[iii+1].split("\t")[0])
                elif "sig1_2:" in lll:      sig1_2_med=float(result_lines[iii+1].split("\t")[0])
                elif "A2_1:" in lll:        gaussLorentz2_A1=float(result_lines[iii+1].split("\t")[0])
                elif "sig2_1:" in lll:
                    gaussLorentz2_sigma1=float(result_lines[iii+1].split("\t")[0])
                    break
    
    
    
    
    if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
        RV=theta[num_start_RVs:]
        
        
    
    if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
        for anRV, anRVbound in zip(RV,used_RV_boundaries):  # go through all individual RVs and make sure they are in boundaries
            if not anRVbound[0] < anRV < anRVbound[1]:        return -np_inf
    else:
        if forced_K1=="Fit" and forced_K2=="Fit":
            if forced_Vgamma1=="Fit" and forced_Vgamma2=="Fit":
                if not (p0K1[0] < mcmc_K1 < p0K1[1] and p0K2[0] < mcmc_K2 < p0K2[1] and p0Vgamma1[0] < mcmc_Vgamma1 < p0Vgamma1[1] and p0Vgamma2[0] < mcmc_Vgamma2 < p0Vgamma2[1]):
                    return -np_inf
            else:
                if not (p0K1[0] < mcmc_K1 < p0K1[1] and p0K2[0] < mcmc_K2 < p0K2[1]):              return -np_inf
        elif forced_K1=="Fit":
            if forced_Vgamma1=="Fit":
                if not (p0K1[0] < mcmc_K1 < p0K1[1] and p0Vgamma1[0] < mcmc_Vgamma1 < p0Vgamma1[1]):    return -np_inf
            else:
                if not p0K1[0] < mcmc_K1 < p0K1[1]:       return -np_inf
        elif forced_K2=="Fit":
            if forced_Vgamma2=="Fit":
                if not (p0K2[0] < mcmc_K2 < p0K2[1] and p0Vgamma2[0] < mcmc_Vgamma2 < p0Vgamma2[1]):    return -np_inf
            else:
                if not p0K2[0] < mcmc_K2 < p0K2[1]:       return -np_inf
        
        if forced_K1!="Fit" and forced_Vgamma1=="Fit":
            if not p0Vgamma1[0] < mcmc_Vgamma1 < p0Vgamma1[1]:    return -np_inf
        if forced_K1!="Fit" and forced_Vgamma1=="Fit" and forced_K2!="Fit" and forced_Vgamma2=="Fit":
            if not (p0Vgamma1[0] < mcmc_Vgamma1 < p0Vgamma1[1] and p0Vgamma2[0] < mcmc_Vgamma2 < p0Vgamma2[1]):    return -np_inf
        if forced_K1!="Fit" and forced_Vgamma1!="Fit" and forced_K2!="Fit" and forced_Vgamma2=="Fit":
            if not p0Vgamma2[0] < mcmc_Vgamma2 < p0Vgamma2[1]:    return -np_inf
    
    return 0.0


def lnlike_quad_lorentz(theta, arguments):
    input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, normalised_wavelength, normalised_flux, normalised_err,  inp_resolution, used_RV_boundaries, HJD_values, sigma_clip, the_unique_wavelengths, spec1wl, spec1flux, spec2wl, spec2flux, inputScaling = arguments
    if type(input_files)==np.str_:
        input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, normalised_wavelength, normalised_flux, normalised_err,  inp_resolution, used_RV_boundaries, HJD_values, sigma_clip = [input_files], [share_rv], [reference_wl], [cut_Ha_all], [normaliseHa_all], [normalised_wavelength], [normalised_flux], [normalised_err],  [inp_resolution], [used_RV_boundaries], [HJD_values], [sigma_clip]
    
    
    if not "run_double_one_by_oneMCMC" in sys.argv:
        args = npargwhere(p0labels=="A1_1")[0][0];   gaussLorentz1_A1 = theta[args]
        args = npargwhere(p0labels=="A2_1")[0][0];   gaussLorentz2_A1 = theta[args]
        args = npargwhere(p0labels=="sig1_1")[0][0];   gaussLorentz1_sigma1 = theta[args]
        args = npargwhere(p0labels=="sig2_1")[0][0];   gaussLorentz2_sigma1 = theta[args]
    
    else:
        with open("out/result.out") as resultfile:
            result_lines = resultfile.readlines()
            for iii, lll in enumerate(result_lines):
                if "A1_1:" in lll:          gaussLorentz1_A1=float(result_lines[iii+1].split("\t")[0])
                elif "sig1_1:" in lll:      gaussLorentz1_sigma1=float(result_lines[iii+1].split("\t")[0])
                elif "A1_2:" in lll:        A1_2_med=float(result_lines[iii+1].split("\t")[0])
                elif "sig1_2:" in lll:      sig1_2_med=float(result_lines[iii+1].split("\t")[0])
                elif "A2_1:" in lll:        gaussLorentz2_A1=float(result_lines[iii+1].split("\t")[0])
                elif "sig2_1:" in lll:
                    gaussLorentz2_sigma1=float(result_lines[iii+1].split("\t")[0])
                    break
    
    

    
    if "K1" in p0labels:            args = npargwhere(p0labels=="K1")[0][0];          mcmc_K1 = theta[args]
    else:                           mcmc_K1 = forced_K1
    
    if "K2" in p0labels:            args = npargwhere(p0labels=="K2")[0][0];          mcmc_K2 = theta[args]
    else:                           mcmc_K2 = forced_K2
    
    if "Vg1" in p0labels:           args = npargwhere(p0labels=="Vg1")[0][0];         mcmc_Vgamma1 = theta[args]
    else:                           mcmc_Vgamma1 = forced_Vgamma1
    
    if "Vg2" in p0labels:           args = npargwhere(p0labels=="Vg2")[0][0];         mcmc_Vgamma2 = theta[args]
    else:                           mcmc_Vgamma2 = forced_Vgamma2
    
    if "RV1_0" in p0labels:         num_start_RVs = npargwhere(p0labels=="RV1_0")[0][0]
    

    
    if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
        RV=theta[num_start_RVs:]
    
    
    
    chisq_spec=0
    
    index = np.arange(0,len(input_files),1)
    
    the_unique_wavelengths=[npamax(the_unique_wavelengths)]
    for ref_wl in the_unique_wavelengths:
        mask_ref_wl = reference_wl==ref_wl
        ref_wl_cut_lims_min, ref_wl_cut_lims_max = cut_Ha_all[mask_ref_wl].T
        cut_wl_min, cut_wl_max = npamin(ref_wl_cut_lims_min), npamax(ref_wl_cut_lims_max)
        model_wl1 = nplinspace(cut_wl_min, cut_wl_max, int(cut_wl_max-cut_wl_min)*10)
        model_wl2 = nplinspace(cut_wl_min, cut_wl_max, int(cut_wl_max-cut_wl_min)*10)
        
        
        model_spectrum_star1 = lorentz_2x_mcmc(model_wl1, gaussLorentz1_A1, gaussLorentz1_sigma1, 0, 1)
        model_spectrum_star2 = lorentz_2x_mcmc(model_wl2, gaussLorentz2_A1, gaussLorentz2_sigma1, 0, 1)
        
        
        model_wl1, model_wl2 = model_wl1 + ref_wl, model_wl2 + ref_wl
        
        
        
        
        # go through all files one by one and do the same as above, but recalculating the interpolated spectrum, which makes things slower
        for (ii, in_fi, sh_rv, modha, cut_lim, norm_lim, aHJD, sigclipspec) in zip(index[mask_ref_wl], input_files[mask_ref_wl], share_rv[mask_ref_wl], modelHa[mask_ref_wl], cut_Ha_all[mask_ref_wl], normaliseHa_all[mask_ref_wl], HJD_values[mask_ref_wl], sigma_clip[mask_ref_wl]):  # the mcmc RVs are gone through 1 by one, uses the ii counter.
            modelHa_min, modelHa_max = modha
            cut_limits_min, cut_limits_max = cut_lim
            norm_limits_min, norm_limits_max = norm_lim
            
            
            
            if not (isinstance(forced_P0, float) and isinstance(forced_T0, float)):
                if sh_rv == -1:  mcmc_rv1, mcmc_rv2 = RV[2*ii], RV[2*ii+1]
                else:            mcmc_rv1, mcmc_rv2 = RV[2*sh_rv], RV[2*sh_rv+1]
            else:
                phi = ((aHJD - forced_T0)%forced_P0)/forced_P0
                if forced_K1=="Fit":     
                    if forced_Vgamma1=="Fit":    mcmc_rv1 = mcmc_Vgamma1 + mcmc_K1*np_sin(2*np_pi*phi)
                    else:    mcmc_rv1 = forced_Vgamma1 + mcmc_K1*np_sin(2*np_pi*phi)
                else:                   
                    if forced_Vgamma1=="Fit":    mcmc_rv1 = mcmc_Vgamma1 + forced_K1*np_sin(2*np_pi*phi)
                    else: mcmc_rv1 = forced_Vgamma1 + forced_K1*np_sin(2*np_pi*phi)
                if forced_K2=="Fit":   
                    if forced_Vgamma2=="Fit":    mcmc_rv2 = mcmc_Vgamma2 + mcmc_K2*np_sin(2*np_pi*phi)
                    else:    mcmc_rv2 = forced_Vgamma2 + mcmc_K2*np_sin(2*np_pi*phi)
                else:   
                    if forced_Vgamma1=="Fit":    mcmc_rv2 = mcmc_Vgamma2 + forced_K2*np_sin(2*np_pi*phi)
                    else:    mcmc_rv2 = forced_Vgamma2 + forced_K2*np_sin(2*np_pi*phi)
                if forced_K1!="Fit" and forced_Vgamma1=="Fit":     mcmc_rv1 = mcmc_Vgamma1 + forced_K1*np_sin(2*np_pi*phi)
                if forced_K2!="Fit" and forced_Vgamma2=="Fit":     mcmc_rv2 = mcmc_Vgamma2 + forced_K2*np_sin(2*np_pi*phi)
                    
            
            
            
            normalised_wavelength = list_norm_wl_grids[ii];    normalised_flux = list_normalised_flux[ii];    normalised_err = list_normalised_err[ii];    inp_resolution = resolutions[ii]
            
            
            dlam1 = ref_wl*mcmc_rv1/speed_of_light
            dlam2 = ref_wl*mcmc_rv2/speed_of_light
            
            if high_RV_amp: modelHa_min+=dlam1; modelHa_max+=dlam1; cut_limits_min+=dlam1; cut_limits_max+=dlam1; norm_limits_min+=dlam1; norm_limits_max+=dlam1

            
            
            # interpolate the model onto the observation's wavelength grid
            interparr = interp(model_wl2, model_wl1+dlam1, model_spectrum_star1)    +   interp(model_wl2, model_wl2+dlam2, model_spectrum_star2)
            
            
            # Smear to the resolution desired
            resstd=0.5*(ref_wl/inp_resolution)/(model_wl2[10]-model_wl2[9])
            interparr = convolve(interparr, Gaussian1DKernel(stddev=resstd), boundary = 'extend')
            interparr = interp(normalised_wavelength, model_wl2, interparr)
                
            
            mask=~((normalised_wavelength>ref_wl-4) & (normalised_wavelength<ref_wl+4))
            
            xquad = normalised_wavelength[mask];   yquad = (normalised_flux - interparr)[mask];   yquaderr = normalised_err[mask]
            
            # Normalize X values (this prevents numerical issues)
            x_mean = npmean(xquad)
            max_minus_min = npamax(xquad) - npamin(xquad)
            x_scaled = (xquad - x_mean) / max_minus_min  # Scale to [-1, 1]
            
            initial_guess = polyfit(x_scaled, yquad, deg=poly_degree)[::-1]  # Reverse order
            
            # Optimize with constraints
            constraints = [{'type': 'eq', 'fun': partial(constraint_function, x_scaled=x_scaled)}]
            result = opt.minimize(partial(loss, x_scaled=x_scaled, yquad=yquad, yquaderr=yquaderr), initial_guess, constraints=constraints, method='SLSQP', options={'ftol': 1e-4, 'gtol': 1e-4, 'maxiter': 200})
            
            # Get final polynomial fit
            quad = fast_polyval(result.x[::-1], (normalised_wavelength - x_mean) / max_minus_min)
            
            interparr += quad
            
            
            
            
            
            
            
            
            #mask = ((normalised_wavelength>ref_wl+cut_limits_min) & (normalised_wavelength<ref_wl+norm_limits_min)) | ((normalised_wavelength>ref_wl+norm_limits_max)   & (normalised_wavelength<ref_wl+cut_limits_max))
            
            
            #try:
            #    m,c =  polyfit(normalised_wavelength[mask], interparr[mask], deg=1)
            #    interparr /= (m*normalised_wavelength + c)
            #except: return -np_inf
            
            
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
            chisq_indiv = -0.5*npsum((np_square(normalised_flux[desired_range & aaamask][~clip_mask]-(off+interparr[desired_range & aaamask][~clip_mask])))/np_square(normalised_err[desired_range & aaamask][~clip_mask]))
            chisq_spec += chisq_indiv
            
            
    
    
    if np.isnan(chisq_spec) == True :    return -np_inf
    else:   return chisq_spec
        
        


def lnprob_quad_lorentz(theta, arguments):
    lp = lnprior_quad_lorentz(theta, arguments)
    if not np.isfinite(lp):
        return -np_inf
    return lp + lnlike_quad_lorentz(theta, arguments)








if sys.argv[1]=="RV" or sys.argv[1]=="RV_gauss" or sys.argv[1]=="RV_double_commonRV12":
    desired_wl=6562.81
    # load in atmospheric result
    found_out_scaling=False
    with open("out/result.out") as resultfile:
        result_lines = resultfile.readlines()
        for iii, lll in enumerate(result_lines):
            if starType1.startswith("D") or starType2.startswith("D"):
                if "Teff1:" in lll:
                    T1_med=float(result_lines[iii+1].split("\t")[0])
                elif "Teff2:" in lll:
                    T2_med=float(result_lines[iii+1].split("\t")[0])
                elif "Logg1:" in lll:
                    logg1_med=float(result_lines[iii+1].split("\t")[0])
                elif "Logg2:" in lll:
                    logg2_med=float(result_lines[iii+1].split("\t")[0])
                elif "H/He1:" in lll:
                    HoverHe1_med=float(result_lines[iii+1].split("\t")[0])
                elif "H/He2:" in lll:
                    HoverHe2_med=float(result_lines[iii+1].split("\t")[0])
                elif "Scaling:" in lll:
                    Scaling_med=float(result_lines[iii+1].split("\t")[0])
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
                elif "A2_1:" in lll:
                    A2_1_med=float(result_lines[iii+1].split("\t")[0])
                elif "sig2_1:" in lll:
                    sig2_1_med=float(result_lines[iii+1].split("\t")[0])
                elif "A2_2:" in lll:
                    A2_2_med=float(result_lines[iii+1].split("\t")[0])
                elif "sig2_2:" in lll:
                    sig2_2_med=float(result_lines[iii+1].split("\t")[0])
                
        resultfile.close()
        
    if starType1.startswith("D") or starType2.startswith("D"):
        try: T1_med
        except: T1_med = forced_teff1
        try: logg1_med
        except: logg1_med = forced_logg1
        try: T2_med
        except: T2_med = forced_teff2
        try: logg2_med
        except: logg2_med = forced_logg2
        if found_out_scaling==False:
            if forced_Scaling==False:    Scaling_med = Scaling_med
            elif forced_Scaling=="WD":
                if starType1=="DA":     R1 = get_MTR(T1_med, logg=logg1_med, return_R_from_T_logg=True, loaded_Istrate=loaded_Istrate, loaded_CO=loaded_CO, loaded_Althaus=loaded_Althaus)
                elif starType1=="DBA" or starType1=="DC" or starType1=="DB":  R1 = get_MTR_DB(T1_med, logg=logg1_med)
                
                if starType2=="DA":     R2 = get_MTR(T2_med, logg=logg2_med, return_R_from_T_logg=True, loaded_Istrate=loaded_Istrate, loaded_CO=loaded_CO, loaded_Althaus=loaded_Althaus)
                elif starType2=="DBA" or starType2=="DC" or starType2=="DB":  R2 = get_MTR_DB(T2_med, logg=logg2_med)
                Scaling_med=(R1**2/R2**2)
            else: Scaling_med=forced_Scaling
                
    
    
        pier_or_antoine="mixed"#"pier3Dphot_antoine1Dspec"
        for tempcnt, (temperature_star, logg_star, starType) in enumerate(zip([T1_med,T2_med],[logg1_med,logg2_med], [starType1,starType2])):
            minval=10000;   maxval=-10000
            if starType=="DA":
                if not pier_or_antoine=="pier3Dphot_antoine1Dspec":
                    if tempcnt==0:     Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DAgrids(temperature_star, logg_star)
                    elif tempcnt==1:   Grav2_N, wl_all2_N, flux2_N, Teff2_N = return_DAgrids(temperature_star, logg_star)
                else:
                    if tempcnt==0:     Grav1_N, wl_all1_N, flux1_N, Teff1_N, photGrav1_N, photwl_all1_N, photflux1_N, photTeff1_N = return_DAgrids(temperature_star, logg_star)
                    elif tempcnt==1:   Grav2_N, wl_all2_N, flux2_N, Teff2_N, photGrav2_N, photwl_all2_N, photflux2_N, photTeff2_N = return_DAgrids(temperature_star, logg_star)
                
            

            elif starType=="DBA":
                if tempcnt==0:    HoverHestar=HoverHe1_med;    Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N = return_DBAgrids(temperature_star, logg_star, HoverHestar)
                elif tempcnt==1:  HoverHestar=HoverHe2_med;    Grav2_N, wl_all2_N, flux2_N, Teff2_N, HoverHe2_N = return_DBAgrids(temperature_star, logg_star, HoverHestar)
            
            
            elif starType=="DB":
                if tempcnt==0:     Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DBgrids(temperature_star, logg_star)
                if tempcnt==1:     Grav2_N, wl_all2_N, flux2_N, Teff2_N = return_DBgrids(temperature_star, logg_star)
                
                
            elif starType=="DC":
                if tempcnt==0:    Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DCgrids(temperature_star, logg_star)
                elif tempcnt==1:  Grav2_N, wl_all2_N, flux2_N, Teff2_N = return_DCgrids(temperature_star, logg_star)
    
    else:   Scaling_med=1
    
    
    cut_limits_min, cut_limits_max = npamin(cut_Ha_all.T), npamax(cut_Ha_all.T)
    
    if starType1.startswith("D") or starType2.startswith("D"):
        if starType1=="DA": model_wl1, model_spectrum_star1 = return_model_spectrum_DA(wl_all1_N, desired_wl, cut_limits_min-20, cut_limits_max+20, Grav1_N, flux1_N, Teff1_N, T1_med, logg1_med)
        elif starType1=="DBA": model_wl1, model_spectrum_star1 = return_model_spectrum_DBA(wl_all1_N, desired_wl, cut_limits_min-20, cut_limits_max+20, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, T1_med, logg1_med, HoverHe1_med)
        elif starType1=="DC" or starType1=="DB": model_wl1, model_spectrum_star1 = return_model_spectrum_DA(wl_all1_N, desired_wl, cut_limits_min-20, cut_limits_max+20, Grav1_N, flux1_N, Teff1_N, T1_med, logg1_med)
        if starType2=="DA": model_wl2, model_spectrum_star2 = return_model_spectrum_DA(wl_all2_N, desired_wl, cut_limits_min-20, cut_limits_max+20, Grav2_N, flux2_N, Teff2_N, T2_med, logg2_med)
        elif starType2=="DBA": model_wl2, model_spectrum_star2 = return_model_spectrum_DBA(wl_all2_N, desired_wl, cut_limits_min-20, cut_limits_max+20, Grav2_N, flux2_N, Teff2_N, HoverHe2_N, T2_med, logg2_med, HoverHe2_med)
        elif starType2=="DC" or starType2=="DB": model_wl2, model_spectrum_star2 = return_model_spectrum_DA(wl_all2_N, desired_wl, cut_limits_min-20, cut_limits_max+20, Grav2_N, flux2_N, Teff2_N, T2_med, logg2_med)
        
        model_spectrum_star1*=Scaling_med
    
    
    else:
        model_wl1 = model_wl2 = nplinspace(cut_limits_min, cut_limits_max, int(cut_limits_max-cut_limits_min)*1000)
        if starType1=="GG":     model_spectrum_star1 = gauss_2x_mcmc(model_wl1, A1_1_med, sig1_1_med, A1_2_med, sig1_2_med) + 1
        elif starType1=="LL":   model_spectrum_star1 = lorentz_2x_mcmc(model_wl1, A1_1_med, sig1_1_med, A1_2_med, sig1_2_med) + 1
        elif starType1=="GL":   model_spectrum_star1 = gauss_lorentz_mcmc(model_wl1, A1_1_med, sig1_1_med, A1_2_med, sig1_2_med) + 1
        if starType2=="GG":     model_spectrum_star2 = gauss_2x_mcmc(model_wl2, A2_1_med, sig2_1_med, A2_2_med, sig2_2_med)
        elif starType2=="LL":   model_spectrum_star2 = lorentz_2x_mcmc(model_wl2, A2_1_med, sig2_1_med, A2_2_med, sig2_2_med)
        elif starType2=="GL":   model_spectrum_star2 = gauss_lorentz_mcmc(model_wl2, A2_1_med, sig2_1_med, A2_2_med, sig2_2_med)
        if starType1=="quadLorentz" and starType2=="quadLorentz":
            model_spectrum_star1 = lorentz_2x_mcmc(model_wl1, A1_1_med, sig1_1_med, 0, 1)
            model_spectrum_star2 = lorentz_2x_mcmc(model_wl2, A2_1_med, sig2_1_med, 0, 1)
        model_wl1, model_wl2 = model_wl1 + npamax(reference_wl), model_wl2 + npamax(reference_wl)
    
    
    from scipy.optimize import least_squares
    
    
    if sys.argv[1]=="RV_gauss":
        desired_refwl = npamax(reference_wl)
        mask_ref_wl = reference_wl==desired_refwl
        
        
        ###################################
        ##### ADD EXTRA GAUSSIAN IF DESIRED
        ###################################
        if "run_double_one_by_oneMCMC" in sys.argv:
            file_ignore=np.array([])
        else:
            try: file_ignore = np.asarray([config_info["file_ignore"]])[0]
            except: file_ignore = []
        
        
        
        
        ndim=4
        p0A1, p0A2 = [-0.05, 0], [-0.05, 0]
        p0std_dev1,  p0std_dev2 = [0, 1.5], [0, 1.5]
        
        pos_min = np.array([p0A1[0], p0std_dev1[0], p0A2[0], p0std_dev2[0]])
        pos_max = np.array([p0A1[1], p0std_dev1[1], p0A2[1], p0std_dev2[1]])
        p0range = np.array([p0A1, p0A2, p0std_dev1,  p0std_dev2])
        p0labels = np.array(["A1", "std1", "A2", "std2"])
        for ii, in_fi in enumerate(input_files[mask_ref_wl]):
            if not in_fi in file_ignore:
                rvbound1, rvbound2 = RV_boundaries1[ii], RV_boundaries2[ii]
                
                ndim+=1
                pos_min = np.append(pos_min, rvbound1[0]);   pos_max = np.append(pos_max, rvbound1[1])
                p0range = np.concatenate((p0range, np.array([rvbound1])))
                p0labels = np.append(p0labels,"RV1_"+str(ii))
                
                ndim+=1
                pos_min = np.append(pos_min, rvbound2[0]);   pos_max = np.append(pos_max, rvbound2[1])
                p0range = np.concatenate((p0range, np.array([rvbound2])))
                p0labels = np.append(p0labels,"RV2_"+str(ii))
        
        nwalkers=len(p0labels)+(len(p0labels)-4)*2 + 2
        burnin=100
        nsteps=50
        
        def lnprior_gauss(theta, arguments):
            A1, std_dev1, A2, std_dev2 = theta[0], theta[1], theta[2], theta[3]
            RVs = theta[4:]
            RVs_star1 = RVs[::2]
            RVs_star2 = RVs[1::2]
            
            
            if p0A1[0]<A1<p0A1[1] and p0std_dev1[0]<std_dev1<p0std_dev1[1] and p0A2[0]<A2<p0A2[1] and p0std_dev2[0]<std_dev2<p0std_dev2[1]:
                return 0.0
            else:   return np_inf
        
        def agauss(x, A, cen, std_dev):
            return A*np.exp(-np_square(x-cen)/(2*std_dev**2))

        def lnlike_gauss(theta, arguments):
            A1, std_dev1, A2, std_dev2 = theta[0], theta[1], theta[2], theta[3]
            RVs = theta[4:]
            RVs_star1 = RVs[::2]
            RVs_star2 = RVs[1::2]
            
            input_files, reference_wl, cut_Ha_all, normaliseHa_all, list_norm_wl_grids, list_normalised_flux, list_normalised_err,  resolutions, used_RV_boundaries, wl1, spec1, wl2, spec2, desired_refwl = arguments
            
            
            gauss1 = agauss(wl1-desired_refwl, A1, 0, std_dev1)
            gauss2 = agauss(wl1-desired_refwl, A2, 0, std_dev2)
            
            chisq=0
            for ii, (in_fi, modha, cut_lim, norm_lim, RV1, RV2) in enumerate(zip(input_files, modelHa, cut_Ha_all, normaliseHa_all, RVs_star1, RVs_star2)):
                wl_RV1 = RV1/speed_of_light * desired_refwl;    wl_RV2 = RV2/speed_of_light * desired_refwl
                normalised_wavelength, normalised_flux, normalised_err, inp_resolution = list_norm_wl_grids[ii], list_normalised_flux[ii], list_normalised_err[ii], resolutions[ii]
                rvbound1, rvbound2 = RV_boundaries1[ii], RV_boundaries2[ii]
                modelHa_min, modelHa_max = modha
                cut_limits_min, cut_limits_max = cut_lim
                norm_limits_min, norm_limits_max = norm_lim
                
                
                resAA = desired_refwl/inp_resolution
                
                
                smear_model_spectrum_star1 = convolve(spec1, Gaussian1DKernel(stddev=0.5*resAA/(wl1[10]-wl1[9])), boundary = 'extend')
                smear_model_spectrum_star2 = convolve(spec2, Gaussian1DKernel(stddev=0.5*resAA/(wl2[10]-wl2[9])), boundary = 'extend')
                
                smear_model_spectrum_star1 = interp(wl1, wl1+wl_RV1, smear_model_spectrum_star1)
                smear_model_spectrum_star2 = interp(wl1, wl2+wl_RV2, smear_model_spectrum_star2)
                
                
                interparr = smear_model_spectrum_star1 + smear_model_spectrum_star2
                
                
                mask = ((wl1>desired_refwl+cut_limits_min) & (wl1<desired_refwl+norm_limits_min)) | ((wl1>desired_refwl+norm_limits_max)   & (wl1<desired_refwl+cut_limits_max))
                
                m,c =  polyfit(wl1[mask], interparr[mask], deg=1)
                interparr /= (m*wl1 + c)
                interparr = interp(normalised_wavelength, wl1, interparr)
                
                
                smear_gauss1 = convolve(gauss1, Gaussian1DKernel(stddev=0.5*resAA/(wl1[10]-wl1[9])), boundary = 'extend')
                smear_gauss2 = convolve(gauss2, Gaussian1DKernel(stddev=0.5*resAA/(wl2[10]-wl2[9])), boundary = 'extend')
                
                smear_gauss1 = interp(normalised_wavelength, wl1+wl_RV1, smear_gauss1)
                smear_gauss2 = interp(normalised_wavelength, wl1+wl_RV2, smear_gauss2)
                
                
                #plt.plot(normalised_wavelength, interparr+smear_gauss1+smear_gauss2,c='g');  plt.plot(normalised_wavelength, interparr,c='r');  plt.plot(normalised_wavelength, normalised_flux, c='k')
                ##plt.plot(normalised_wavelength, smear_gauss1+1, c='orange');   plt.plot(normalised_wavelength, smear_gauss2+1, c='b')
                ##plt.plot(normalised_wavelength, smear_gauss1+smear_gauss2+1, c='purple')
                #plt.title(str(A1) + "   " + str(A2) + "   " + str(std_dev1) + "   " + str(std_dev2));   plt.show()
                
                
                desired_range = (normalised_wavelength>desired_wl-25)  &  (normalised_wavelength<desired_wl+25)
                chisq += -0.5*npsum((np_square(normalised_flux[desired_range]-interparr[desired_range]))/np_square(normalised_err[desired_range]))
                    
            return chisq
            
        
        
        def lnprob_gauss(theta, arguments):
            lp = lnprior_gauss(theta, arguments)
            if not np.isfinite(lp):
                return -np_inf
            return lp + lnlike_gauss(theta, arguments)
        
                    
        
                
                    
                    
        pool = MPIPool()
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
        print("Starting MCMC")
        


        psize = pos_max - pos_min
        p0 = [pos_min + psize*np.random.rand(ndim) for i in range(nwalkers)]
        
        
        
        
        if "run_double_one_by_oneMCMC" in sys.argv:   file_ignore=np.array([])
        else:   file_ignore=np.asarray([config_info["file_ignore"]])[0]
        
        process_gauss_inpf, process_gauss_rwl, process_gauss_cHa, process_gauss_nHa, process_gauss_normwl, process_gauss_normfl, process_gauss_normerr, process_gauss_res, process_gauss_rvbound = [], [], [], [], [], [], [], [], []
        
        if len(file_ignore)>0:
            for inpf, rwl, cHa, nHa, normwl, normfl, normerr, res, rvbound in zip(input_files, reference_wl, cut_Ha_all, normaliseHa_all, list_norm_wl_grids, list_normalised_flux, list_normalised_err,  resolutions, used_RV_boundaries):
                if not inpf in file_ignore and rwl==desired_refwl:
                    process_gauss_inpf.append(inpf);          process_gauss_rwl.append(rwl)
                    process_gauss_cHa.append(cHa);            process_gauss_nHa.append(nHa)
                    process_gauss_normwl.append(normwl);      process_gauss_normfl.append(normfl)
                    process_gauss_normerr.append(normerr);    process_gauss_res.append(res)
                    process_gauss_rvbound.append(rvbound)

        
        
        process_gauss_inpf = np.asarray(process_gauss_inpf);       process_gauss_rwl = np.asarray(process_gauss_rwl);         process_gauss_cHa = np.asarray(process_gauss_cHa)
        process_gauss_nHa = np.asarray(process_gauss_nHa);         process_gauss_normwl = np.asarray(process_gauss_normwl);   process_gauss_normfl = np.asarray(process_gauss_normfl)
        process_gauss_normerr = np.asarray(process_gauss_normerr); process_gauss_res = np.asarray(process_gauss_res);         process_gauss_rvbound = np.asarray(process_gauss_rvbound)
        
        
        #sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_gauss, pool=pool, args=[[input_files, reference_wl, cut_Ha_all, normaliseHa_all, list_norm_wl_grids, list_normalised_flux, list_normalised_err,  resolutions, used_RV_boundaries, model_wl1, model_spectrum_star1, model_wl2, model_spectrum_star2, desired_refwl]])
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_gauss, pool=pool, args=[[process_gauss_inpf, process_gauss_rwl, process_gauss_cHa, process_gauss_nHa, process_gauss_normwl, process_gauss_normfl, process_gauss_normerr,  process_gauss_res, process_gauss_rvbound, model_wl1, model_spectrum_star1, model_wl2, model_spectrum_star2, desired_refwl]])
        pos, prob, state = sampler.run_mcmc(p0,burnin,progress=True)
        samples = sampler.flatchain
        np.savetxt("RVfits/Gauss_MCMC_samples_burnin.dat",samples)
        np.savetxt("RVfits/Gauss_MCMC_burninpos.dat",pos)
        print("Finished burn-in")



        #Make a corner plot which shows values
        samples=nploadtxt("RVfits/Gauss_MCMC_samples_burnin.dat")

        sampler.reset()
        pos=nploadtxt("RVfits/Gauss_MCMC_burninpos.dat")
        

        pos, prob, state = sampler.run_mcmc(pos,nsteps,progress=True)
        print("MCMC completed")
        pool.close()

        samples = sampler.flatchain
        np.savetxt("RVfits/Gauss_MCMC_samples.dat",samples)

        #Make a corner plot which shows values
        samples=nploadtxt("RVfits/Gauss_MCMC_samples.dat")
        print(samples)
        
        A1_result = samples[::,0];    std_dev1_result = samples[::,1];    A2_result = samples[::,2];   std_dev2_result = samples[::,3]
        A1_med = np.percentile(A1_result, 50, axis=0);  A1_min = np.percentile(A1_result, 16, axis=0);   A1_max = np.percentile(A1_result, 84, axis=0)
        std_dev1_med = np.percentile(std_dev1_result, 50, axis=0);  std_dev1_min = np.percentile(std_dev1_result, 16, axis=0);   std_dev1_max = np.percentile(std_dev1_result, 84, axis=0)
        A2_med = np.percentile(A2_result, 50, axis=0);  A2_min = np.percentile(A2_result, 16, axis=0);   A2_max = np.percentile(A2_result, 84, axis=0)
        std_dev2_med = np.percentile(std_dev2_result, 50, axis=0);  std_dev2_min = np.percentile(std_dev2_result, 16, axis=0);   std_dev2_max = np.percentile(std_dev2_result, 84, axis=0)
        
        
        #if plot_corner==True:
        #    fig2 = corner.corner(samples, labels=p0labels, range=p0range, bins=30,smooth=True, quantiles=[0.16, 0.5, 0.84], show_titles=True, labels_args={"fontsize": 40}, title_fmt = '.3f')
        #    fig2.savefig("RVfits/Gauss_corner.pdf")
        #    del fig2
        #plt.close();  plt.clf()
        
        
        
        
        
        
    ############################
    ##### RUN RVS FOR A SOLUTION
    ############################

    def WDmodel(x, RV1, RV2, wl1, spec1, wl2, spec2,gauss1=np.array([1,1]),gauss2=np.array([1,1])):
        wl_RV1 = RV1/speed_of_light * desired_refwl
        wl_RV2 = RV2/speed_of_light * desired_refwl
        
        # interpolate the model onto the observation's wavelength grid
        interparr = interp(x, wl1+wl_RV1, spec1) + interp(x, wl2+wl_RV2, spec2)
        
        mask = ((x>desired_refwl+cut_limits_min) & (x<desired_refwl+norm_limits_min)) | ((x>desired_refwl+norm_limits_max)   & (x<desired_refwl+cut_limits_max))
        
        if not (starType1=="quadLorentz" and starType2=="quadLorentz"):
            m,c =  polyfit(x[mask], interparr[mask], deg=1)	
            interparr /= (m*x + c)

        
        if gauss1[0]!=1 and gauss2[0]!=1:
            interparr += interp(x,wl1+wl_RV1,gauss1) + interp(x,wl2+wl_RV2,gauss2)
        
        return interparr

    

    def residual(params, xxx, yyy, yyyerr,wl1,spec1,wl2,spec2,gauss1,gauss2,initial_guess=0, xdetails=[]):
        '''Calculates residuals for a given model'''
        
        if starType1=="quadLorentz" and starType2=="quadLorentz":
            RV1, RV2 = params
            x_mean, max_minus_min, x_scaled, constraints = xdetails
            model = WDmodel(xxx, RV1, RV2, wl1, spec1, wl2, spec2,gauss1,gauss2)# / (xxx*extra_norm + coff)
            
            args=np.argsort(xxx)
            yquad = (yyy - model)[args];   yquaderr = yyyerr[args]   # xquad = xxx[args]
            
            
            initial_guess = polyfit(x_scaled, yquad, deg=poly_degree)[::-1]  # Reverse order
            
            result = opt.minimize(partial(loss, x_scaled=x_scaled, yquad=yquad, yquaderr=yquaderr), initial_guess, constraints=constraints, method='SLSQP', options={'maxiter': 200})
            
            model+=fast_polyval(result.x[::-1], (xxx - x_mean) / max_minus_min)
            
        else:
            RV1, RV2, extra_norm, coff = params
            model = WDmodel(xxx, RV1, RV2, wl1, spec1, wl2, spec2,gauss1,gauss2)
            model/=(xxx*extra_norm + coff)
        
        return (yyy-model) / yyyerr
        
        
    
    def residual_commonRV(params, xxx, yyy, yyyerr,wl1,spec1,wl2,spec2,gauss1,gauss2):
        '''Calculates residuals for a given model'''
        RV1, extra_norm, coff = params
        model = (WDmodel(xxx, RV1, RV1, wl1, spec1, wl2, spec2,gauss1,gauss2))/(xxx*extra_norm + coff)
        return (yyy-model) / yyyerr


    def fit_bootstrap(params, xxx, yyy, yyyerr, bounds, num_its, wl1, spec1, wl2, spec2, gauss1=np.array([1,1]), gauss2=np.array([1,1]), Common_RV_both_stars=False):
        '''Samples the data a thousand times, fits each sample using least-squares,
        then reports median and standard deviation of the result'''
        if not (starType1=="quadLorentz" and starType2=="quadLorentz"):    RV1_guess, RV2_guess, extra_norm_guess, coff_guess = params
        else:
            RV1_guess, RV2_guess = params
            model = WDmodel(xxx, RV1_guess, RV2_guess, wl1, spec1, wl2, spec2,gauss1,gauss2)# / (xxx*extra_norm + coff)
            
            args=np.argsort(xxx)
            xquad = xxx[args];   yquad = (yyy - model)[args];   yquaderr = yyyerr[args]
            
            x_mean = npmean(xquad)
            max_minus_min = npamax(xquad) - npamin(xquad)
            x_scaled = (xquad - x_mean) / max_minus_min  # Scale to [-1, 1]
            
            initial_guess = polyfit(x_scaled, yquad, deg=poly_degree)[::-1]  # Reverse order
        
        
        RV1, RV2, extra_norm, coff = [], [], [], []
        data = np.array((xxx, yyy, yyyerr)).T
        if Common_RV_both_stars==False:
            for i in range(num_its):
                # Select 100% of our sample
                sample = data[np.random.choice(len(data),size=int(1.0*len(data)), replace=True)]
                # get the vectors back
                x = sample.T[0]
                y = sample.T[1]
                err = sample.T[2]
                
                
                if starType1=="quadLorentz" and starType2=="quadLorentz":
                    lsm_fit = least_squares(residual, [RV1_guess+np.random.normal(0, 1), RV2_guess+np.random.normal(0, 1)], method='trf',
                                args=(x,y,err,wl1,spec1,wl2,spec2,gauss1,gauss2,initial_guess,[x_mean, max_minus_min, x_scaled, [{'type': 'eq', 'fun': partial(constraint_function, x_scaled=x_scaled)}]]), ftol=1E-8, bounds=((bounds[0], bounds[2]), (bounds[1], bounds[3])))
                                
                    RV1.append(lsm_fit.x[0]);   RV2.append(lsm_fit.x[1])
                    
                else:
                    lsm_fit = least_squares(residual, [RV1_guess, RV2_guess, extra_norm_guess, coff_guess], method='trf',
                                args=(x,y,err,wl1,spec1,wl2,spec2,gauss1,gauss2), ftol=1E-8, bounds=((bounds[0], bounds[2], bounds[4], bounds[6]), (bounds[1], bounds[3], bounds[5], bounds[7])))
                                
                    RV1.append(lsm_fit.x[0]);   RV2.append(lsm_fit.x[1]);   extra_norm.append(lsm_fit.x[2]);   coff.append(lsm_fit.x[3])
                    
            if starType1=="quadLorentz" and starType2=="quadLorentz":
                return np.median(RV1), np.std(RV1), np.median(RV2), np.std(RV2)
            else:
                return np.median(RV1), np.std(RV1), np.median(RV2), np.std(RV2), np.median(extra_norm), np.std(extra_norm), np.median(coff), np.std(coff)
                
        else:
            for i in range(num_its):
                # Select 100% of our sample
                sample = data[np.random.choice(len(data),size=int(1.0*len(data)), replace=True)]
                # get the vectors back
                x = sample.T[0]
                y = sample.T[1]
                err = sample.T[2]
                lsm_fit = least_squares(residual_commonRV, [RV1_guess, extra_norm_guess, coff_guess], method='trf',
                                args=(x,y,err,wl1,spec1,wl2,spec2,gauss1,gauss2), ftol=1E-8, bounds=((bounds[0], bounds[4], bounds[6]), (bounds[1], bounds[5], bounds[7])))
                
                RV1.append(lsm_fit.x[0]);    extra_norm.append(lsm_fit.x[1]);   coff.append(lsm_fit.x[2])
            
            return np.median(RV1), np.std(RV1), np.median(extra_norm), np.std(extra_norm), np.median(coff), np.std(coff)
    
    

    def residual2(params, xxx, yyy, yyyerr,wl1,spec1,wl2,spec2,gauss1,gauss2,model_range):
        '''Calculates residuals for a given model'''

        RV1, RV2, offset = params
        model = WDmodel(xxx, RV1, RV2, wl1, spec1, wl2, spec2,gauss1,gauss2) + offset
        mask_narrow = (xxx > desired_refwl-model_range) & (xxx < desired_refwl+model_range)
        return (yyy[mask_narrow]-model[mask_narrow]) / yyyerr[mask_narrow]
    
    
    
    def residual2_quadlorentz(params, xxx, yyy, yyyerr,wl1,spec1,wl2,spec2,gauss1,gauss2,model_range):
        '''Calculates residuals for a given model'''
        RV1, RV2 = params
        lorentzComponents=interp(xxx, wl1+RV1/speed_of_light * desired_refwl, spec1)  + interp(xxx, wl2+RV2/speed_of_light * desired_refwl, spec2)
        
        
        xquad = xxx;   yquad = yyy - lorentzComponents;   yquaderr=yyyerr
        
        # Normalize X values (this prevents numerical issues)
        x_mean = npmean(xquad)
        max_minus_min = npamax(xquad) - npamin(xquad)
        x_scaled = (xquad - x_mean) / max_minus_min  # Scale to [-1, 1]
        
        initial_guess = polyfit(x_scaled, yquad, deg=poly_degree)[::-1]  # Reverse order
        
        # Optimize with constraints
        result = opt.minimize(partial(loss, x_scaled=x_scaled, yquad=yquad, yquaderr=yquaderr), initial_guess, constraints=constraints, method='SLSQP', options={'ftol': 1e-4, 'gtol': 1e-4, 'maxiter': 200})
        
        # Get final polynomial fit
        quad = fast_polyval(result.x[::-1], (normalised_wavelength - x_mean) / max_minus_min)
        
        
        model=quad+lorentzComponents
        mask_narrow = (xxx > desired_refwl-model_range) & (xxx < desired_refwl+model_range)
        
        
        plt.plot(xxx[mask_narrow], yyy[mask_narrow])
        plt.plot(xxx[mask_narrow], model[mask_narrow])
        plt.show()
        
        
        return (yyy[mask_narrow]-model[mask_narrow]) / yyyerr[mask_narrow]
    
    
    
    def residual2_commonRV(params, xxx, yyy, yyyerr,wl1,spec1,wl2,spec2,gauss1,gauss2,model_range):
        '''Calculates residuals for a given model'''
        RV1, offset = params
        model = WDmodel(xxx, RV1, RV1, wl1, spec1, wl2, spec2,gauss1,gauss2) + offset

        mask_narrow = (xxx > desired_refwl-model_range) & (xxx < desired_refwl+model_range)
        return (yyy[mask_narrow]-model[mask_narrow]) / yyyerr[mask_narrow]


    def fit_bootstrap2(params, xxx, yyy, yyyerr, bounds, num_its, wl1, spec1, wl2, spec2, gauss1=np.array([1,1]), gauss2=np.array([1,1]), Common_RV_both_stars=False, model_range=40, expected_RV1err=0,expected_RV2err=0):
        '''Samples the data a thousand times, fits each sample using least-squares,
        then reports median and standard deviation of the result'''
        RV1_guess, RV2_guess, off_guess = params
        RV1, RV2, offset = [], [], []
        data = np.array((xxx, yyy, yyyerr)).T
        if Common_RV_both_stars==False:
            for i in range(num_its):
                sample = data[np.random.choice(len(data),size=int(1.0*len(data)), replace=True)]
                x = sample.T[0];    y = sample.T[1];    err = sample.T[2]
                if starType1=="quadLorentz" and starType2=="quadLorentz":
                    noiseRV1=np.random.normal(0, expected_RV1err)
                    noiseRV2=np.random.normal(0, expected_RV2err)
                    if RV1_guess+noiseRV1<bounds[0]:    noiseRV1=bounds[0]-RV1_guess
                    elif RV1_guess+noiseRV1>bounds[1]:  noiseRV1=bounds[1]-RV1_guess
                    
                    if RV2_guess+noiseRV2<bounds[2]:    noiseRV2=bounds[2]-RV2_guess
                    elif RV2_guess+noiseRV2>bounds[3]:  noiseRV2=bounds[3]-RV2_guess
                    
                    
                    lsm_fit = least_squares(residual2_quadlorentz, [RV1_guess+noiseRV1, RV2_guess+noiseRV2], method='trf',
                                args=(x,y,err,wl1,spec1,wl2,spec2,gauss1,gauss2, model_range), ftol=1E-8, bounds=((bounds[0], bounds[2]), (bounds[1], bounds[3])))
                
                    RV1.append(lsm_fit.x[0]);    RV2.append(lsm_fit.x[1])
                else:
                    lsm_fit = least_squares(residual2, [RV1_guess, RV2_guess, off_guess], method='trf',
                                args=(x,y,err,wl1,spec1,wl2,spec2,gauss1,gauss2, model_range), ftol=1E-8, bounds=((bounds[0], bounds[2], bounds[4]), (bounds[1], bounds[3], bounds[5])))
                
                    RV1.append(lsm_fit.x[0]);    RV2.append(lsm_fit.x[1]);    offset.append(lsm_fit.x[2])
            
            if starType1=="quadLorentz" and starType2=="quadLorentz":
                return np.median(RV1), np.std(RV1), np.median(RV2), np.std(RV2), 0, 0
            else:
                return np.median(RV1), np.std(RV1), np.median(RV2), np.std(RV2), np.median(offset), np.std(offset)
        else:
            
            for i in range(num_its):
                sample = data[np.random.choice(len(data),size=int(1.0*len(data)), replace=True)]
                x = sample.T[0];    y = sample.T[1];    err = sample.T[2]
                lsm_fit = least_squares(residual2_commonRV, [RV1_guess, off_guess], method='trf',
                                args=(x,y,err,wl1,spec1,wl2,spec2,gauss1,gauss2), ftol=1E-8, bounds=((bounds[0], bounds[4]), (bounds[1], bounds[5])))
                
                RV1.append(lsm_fit.x[0]);    offset.append(lsm_fit.x[1])
            
            return np.median(RV1), np.std(RV1), np.median(offset), np.std(offset)
            
    




    try: os.mkdir("RVfits")
    except: None
    
    
    if sys.argv[1]=="RV_double_commonRV12":  commonRV_1_2=True
    else: commonRV_1_2=False
    
        
    rvfilename=[]
    desired_refwl = npamax(reference_wl)
    mask_ref_wl = reference_wl==desired_refwl
    ref_wl = desired_refwl
    
    
    
    for ii, (in_fi, modha, cut_lim, norm_lim, aHJD, sigclipspec) in enumerate(zip(input_files[mask_ref_wl], modelHa[mask_ref_wl], cut_Ha_all[mask_ref_wl], normaliseHa_all[mask_ref_wl], HJD_values[mask_ref_wl], sigma_clip[mask_ref_wl])):
        issue_here=False
        
        if npamin(modha)>-40 or npamax(modha)<40:   model_range = npamin(np.abs(modha))
        else:   model_range=npamin(np.abs(modha))
        
        
        try:
            
            normalised_wavelength, normalised_flux, normalised_err, inp_resolution = list_norm_wl_grids[ii], list_normalised_flux[ii], list_normalised_err[ii], resolutions[ii]
            rvbound1, rvbound2 = RV_boundaries1[ii], RV_boundaries2[ii]
            modelHa_min, modelHa_max = modha
            cut_limits_min, cut_limits_max = cut_lim
            norm_limits_min, norm_limits_max = norm_lim
            
            resAA = desired_refwl/inp_resolution
            
            
            smear_model_spectrum_star1 = convolve(model_spectrum_star1, Gaussian1DKernel(stddev=0.5*resAA/(model_wl1[10]-model_wl1[9])), boundary = 'extend')
            smear_model_spectrum_star2 = convolve(model_spectrum_star2, Gaussian1DKernel(stddev=0.5*resAA/(model_wl2[10]-model_wl2[9])), boundary = 'extend')
            smear_model_spectrum_star2 = interp(model_wl1, model_wl2, smear_model_spectrum_star2)
            
            
            
            if sys.argv[1]=="RV_gauss":
                gauss1 = agauss(model_wl1, A1_med, desired_wl, std_dev1_med)
                gauss2 = agauss(model_wl2, A2_med, desired_wl, std_dev2_med)
                
                smear_gauss1 = convolve(gauss1, Gaussian1DKernel(stddev=0.5*resAA/(model_wl1[10]-model_wl1[9])), boundary = 'extend')
                smear_gauss2 = convolve(gauss2, Gaussian1DKernel(stddev=0.5*resAA/(model_wl2[10]-model_wl2[9])), boundary = 'extend')
            
            
            
            
            Low_SNR=True  # note that this improvement to normalisation only makes a noticeable difference when the SNR is low (because the normalisation of the data is not good). For speed, it can absolutely be ignored, but if you're analysing a small sample or on a case by case basis I recommend to leave it on. I noticed an accuracy of a couple 100ms-1 improved in SNR=25 in the wings for Halpha fits and a big improvement in precision for SNR=10-15 data.
            
            #if npamin(modha)>-40 or npamax(modha)<40: Low_SNR=False
            RV1, RV2 = npmean(rvbound1), npmean(rvbound2)   # this is the initial guess for the bootstrap
            if Low_SNR:
                # First try to improve the accuracy of the normalisation using the model and the full wavelength range supplied
                if sys.argv[1]=="RV_gauss":
                    if not (starType1=="quadLorentz" and starType2=="quadLorentz"):
                        returnvars = fit_bootstrap([RV1, RV2, 1, 1],  normalised_wavelength,  normalised_flux,  normalised_err, bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1], -1,1, -20, 20], num_its=10, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, gauss1=smear_gauss1, gauss2=smear_gauss2, Common_RV_both_stars=commonRV_1_2)
                    else:
                        returnvars = fit_bootstrap([RV1, RV2],  normalised_wavelength,  normalised_flux,  normalised_err, bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1]], num_its=10, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, gauss1=smear_gauss1, gauss2=smear_gauss2, Common_RV_both_stars=commonRV_1_2)
                else:
                    if not (starType1=="quadLorentz" and starType2=="quadLorentz"):
                        returnvars = fit_bootstrap([RV1, RV2, 1, 1],  normalised_wavelength,  normalised_flux,  normalised_err, bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1], -1,1, -20, 20], num_its=10, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, Common_RV_both_stars=commonRV_1_2)
                    else:
                        returnvars = fit_bootstrap([RV1, RV2],  normalised_wavelength,  normalised_flux,  normalised_err, bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1]], num_its=10, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, Common_RV_both_stars=commonRV_1_2)
                try: 
                    if not (starType1=="quadLorentz" and starType2=="quadLorentz"):  RV1, RV1err, RV2, RV2err, extra_norm, extra_norm_err, coff, coff_err = returnvars
                    else:  RV1, RV1err, RV2, RV2err = returnvars
                except: RV1, RV1err, extra_norm, extra_norm_err, coff, coff_err = returnvars
                
                
                ######===============================
                # clip out the bad points
                if sigclipspec!=-1:
                    if not (starType1=="quadLorentz" and starType2=="quadLorentz"):
                        resid = normalised_flux - (WDmodel(normalised_wavelength, RV1, RV2, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2))/(normalised_wavelength*extra_norm + coff)
                    else:
                        resid = normalised_flux - WDmodel(normalised_wavelength, RV1, RV2, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2)
                    resid_in_sig = resid/np.std(np.abs(resid))
                    clip_mask = (np.abs(resid_in_sig)<sigclipspec) | ((normalised_wavelength < desired_refwl+5 + np.amax(np.array([desired_refwl*RV1/speed_of_light, desired_refwl*RV2/speed_of_light]))) & (normalised_wavelength > desired_refwl-5  + np.amin(np.array([desired_refwl*RV1/speed_of_light, desired_refwl*RV2/speed_of_light]))  ))
                else: clip_mask=normalised_flux==normalised_flux
                
                
                
                ######===============================
                # Try to improve the accuracy of the normalisation again with the full wavelength range, this time with the bad points sigma clipped
                if sys.argv[1]=="RV_gauss":
                    returnvars = fit_bootstrap([RV1, RV2, extra_norm, coff],  normalised_wavelength[clip_mask],  normalised_flux[clip_mask],  normalised_err[clip_mask], bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1], -1,1, -20, 20], num_its=100, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, gauss1=smear_gauss1, gauss2=smear_gauss2, Common_RV_both_stars=commonRV_1_2)
                else:
                    try:
                        if not (starType1=="quadLorentz" and starType2=="quadLorentz"):
                            returnvars = fit_bootstrap([RV1, RV2, extra_norm, coff],  normalised_wavelength[clip_mask],  normalised_flux[clip_mask],  normalised_err[clip_mask], bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1], -1,1, -20, 20], num_its=50, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, Common_RV_both_stars=commonRV_1_2)
                        else:
                            returnvars = fit_bootstrap([RV1, RV2],  normalised_wavelength[clip_mask],  normalised_flux[clip_mask],  normalised_err[clip_mask], bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1], -1,1, -20, 20], num_its=50, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, Common_RV_both_stars=commonRV_1_2)
                    except:
                        issue_here=True
                        if not (starType1=="quadLorentz" and starType2=="quadLorentz"):
                            returnvars = fit_bootstrap([RV1, RV2, extra_norm, coff],  normalised_wavelength,  normalised_flux,  normalised_err, bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1], -1,1, -20, 20], num_its=50, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, Common_RV_both_stars=commonRV_1_2)
                        else:
                            returnvars = fit_bootstrap([RV1, RV2],  normalised_wavelength,  normalised_flux,  normalised_err, bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1], -1,1, -20, 20], num_its=50, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, Common_RV_both_stars=commonRV_1_2)
                
                try:
                    if not (starType1=="quadLorentz" and starType2=="quadLorentz"):  RV1, RV1err, RV2, RV2err, extra_norm, extra_norm_err, coff, coff_err = returnvars
                    else:  RV1, RV1err, RV2, RV2err = returnvars
                except: RV1, RV1err, extra_norm, extra_norm_err, coff, coff_err = returnvars
                
                
                ######===============================
                # Apply the tweaked normalisation to the data
                #plt.plot(normalised_wavelength, normalised_flux, c='k')
                if not (starType1=="quadLorentz" and starType2=="quadLorentz"):
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
                resid = normalised_flux - (WDmodel(normalised_wavelength, RV1, RV2, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2))
                
                
                if (starType1=="quadLorentz" and starType2=="quadLorentz"):
                    
                    mask=~((normalised_wavelength>desired_refwl-4) & (normalised_wavelength<desired_refwl+4))
                    
                    xquad = normalised_wavelength[mask];   yquad = normalised_flux[mask];   yquaderr=normalised_err[mask]
                    
                    # Normalize X values (this prevents numerical issues)
                    x_mean = npmean(xquad)
                    max_minus_min = npamax(xquad) - npamin(xquad)
                    x_scaled = (xquad - x_mean) / max_minus_min  # Scale to [-1, 1]
                    
                    initial_guess = polyfit(x_scaled, yquad, deg=poly_degree)[::-1]  # Reverse order
                    
                    # Optimize with constraints
                    constraints = [{'type': 'eq', 'fun': partial(constraint_function, x_scaled=x_scaled)}]
                    result = opt.minimize(partial(loss, x_scaled=x_scaled, yquad=yquad, yquaderr=yquaderr), initial_guess, constraints=constraints, method='SLSQP', options={'ftol': 1e-4, 'gtol': 1e-4, 'maxiter': 200})
                    
                    # Get final polynomial fit
                    quad = fast_polyval(result.x[::-1], (normalised_wavelength - x_mean) / max_minus_min)
                    

                    resid = normalised_flux - (WDmodel(normalised_wavelength, RV1, RV2, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2)) - quad
                    
                    
                    resid_in_sig = resid/np.std(np.abs(resid))
                    clip_mask = (np.abs(resid_in_sig)<sigclipspec) | ((normalised_wavelength < desired_refwl+5 + np.amax(np.array([desired_refwl*RV1/speed_of_light, desired_refwl*RV2/speed_of_light]))   ) & (normalised_wavelength > desired_refwl-5  +  np.amin(np.array([desired_refwl*RV1/speed_of_light, desired_refwl*RV2/speed_of_light]))  ))
                    
                    
                else:
                    resid = normalised_flux - (WDmodel(normalised_wavelength, RV1, RV2, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2) )#+ off)
                    
                    resid_in_sig = resid/np.std(np.abs(resid))
                    clip_mask = (np.abs(resid_in_sig)<sigclipspec) | ((normalised_wavelength < desired_refwl+5 + np.amax(np.array([desired_refwl*RV1/speed_of_light, desired_refwl*RV2/speed_of_light]))   ) & (normalised_wavelength > desired_refwl-5  +  np.amin(np.array([desired_refwl*RV1/speed_of_light, desired_refwl*RV2/speed_of_light]))  ))
                    
                
                
            
                
            else: clip_mask=normalised_flux==normalised_flux
            
            
            
            
            if sigclipspec!=-1:
                # fit for RV1 RV2 with a narrowed range around the line cores
                if not (starType1=="quadLorentz" and starType2=="quadLorentz"):
                    if not issue_here:
                        returnvars = fit_bootstrap2([RV1, RV2, 0],  normalised_wavelength[clip_mask],  normalised_flux[clip_mask],  normalised_err[clip_mask], bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1], -0.2, 0.2], num_its=10, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, model_range=model_range, Common_RV_both_stars=commonRV_1_2, expected_RV1err=RV1err,expected_RV2err=RV2err)
                    else:
                        returnvars = fit_bootstrap2([RV1, RV2, 0],  normalised_wavelength,  normalised_flux,  normalised_err, bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1], -20.0, 20.0], num_its=10, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, model_range=model_range, Common_RV_both_stars=commonRV_1_2, expected_RV1err=RV1err,expected_RV2err=RV2err)
                
                    try: RV1, RV1err, RV2, RV2err, off, off_err = returnvars
                    except: RV1, RV1err, off, off_err = returnvars
                
                
                else:
                    try:
                        returnvars = fit_bootstrap([RV1, RV2, extra_norm, coff],  normalised_wavelength[clip_mask],  normalised_flux[clip_mask],  normalised_err[clip_mask], bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1], -1,1, -20, 20], num_its=10, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, Common_RV_both_stars=commonRV_1_2)
                    except:
                        issue_here=True
                        if not (starType1=="quadLorentz" and starType2=="quadLorentz"):
                            returnvars = fit_bootstrap([RV1, RV2, extra_norm, coff],  normalised_wavelength,  normalised_flux,  normalised_err, bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1], -1,1, -20, 20], num_its=10, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, Common_RV_both_stars=commonRV_1_2)
                        else:
                            returnvars = fit_bootstrap([RV1, RV2],  normalised_wavelength,  normalised_flux,  normalised_err, bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1]], num_its=10, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, Common_RV_both_stars=commonRV_1_2)
            
                    try: 
                        if not (starType1=="quadLorentz" and starType2=="quadLorentz"):  RV1, RV1err, RV2, RV2err, extra_norm, extra_norm_err, coff, coff_err = returnvars
                        else:  RV1, RV1err, RV2, RV2err = returnvars
                    except: RV1, RV1err, extra_norm, extra_norm_err, coff, coff_err = returnvars
                    
                    off=0
                
            
                
                # clip the narrowed area
                mask_narrow = ((normalised_wavelength > desired_refwl-40) & (normalised_wavelength < desired_refwl+40)) | ((normalised_wavelength<desired_refwl+norm_limits_min+5  + np.amax(np.array([desired_refwl*RV1/speed_of_light, desired_refwl*RV2/speed_of_light])) ) | (normalised_wavelength>desired_refwl+norm_limits_max-5   +  np.amin(np.array([desired_refwl*RV1/speed_of_light, desired_refwl*RV2/speed_of_light]))  ))
                
                
                
                if (starType1=="quadLorentz" and starType2=="quadLorentz"):
                    
                    
                    mask=~((normalised_wavelength>desired_refwl-4) & (normalised_wavelength<desired_refwl+4))
                    
                    xquad = normalised_wavelength[mask];   yquad = normalised_flux[mask];   yquaderr = normalised_err[mask]
                    
                    # Normalize X values (this prevents numerical issues)
                    x_mean = npmean(xquad)
                    max_minus_min = npamax(xquad) - npamin(xquad)
                    x_scaled = (xquad - x_mean) / max_minus_min  # Scale to [-1, 1]
                    
                    initial_guess = polyfit(x_scaled, yquad, deg=poly_degree)[::-1]  # Reverse order
                    
                    # Optimize with constraints
                    constraints = [{'type': 'eq', 'fun': partial(constraint_function, x_scaled=x_scaled)}]
                    result = opt.minimize(partial(loss, x_scaled=x_scaled, yquad=yquad, yquaderr=yquaderr), initial_guess, constraints=constraints, method='SLSQP', options={'ftol': 1e-4, 'gtol': 1e-4, 'maxiter': 200})
                    
                    # Get final polynomial fit
                    quad = fast_polyval(result.x[::-1], (normalised_wavelength - x_mean) / max_minus_min)
                    
                    
                    
                    resid = normalised_flux[mask_narrow] - (WDmodel(normalised_wavelength, RV1, RV2, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2) + off)[mask_narrow] - quad[mask_narrow]
                    
                else:
                    resid = normalised_flux[mask_narrow] - (WDmodel(normalised_wavelength, RV1, RV2, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2) + off)[mask_narrow]
                
                
                
                resid_in_sig = resid/np.std(np.abs(resid))
                clip_mask2 = ((np.abs(resid_in_sig)<sigclipspec) | ((normalised_wavelength[mask_narrow] < desired_refwl+5 + np.amax(np.array([desired_refwl*RV1/speed_of_light, desired_refwl*RV2/speed_of_light]))) & (normalised_wavelength[mask_narrow] > desired_refwl-5    +  np.amin(np.array([desired_refwl*RV1/speed_of_light, desired_refwl*RV2/speed_of_light]))  ))) | ((normalised_wavelength[mask_narrow]>desired_refwl+norm_limits_max-5 + np.amin(np.array([desired_refwl*RV1/speed_of_light, desired_refwl*RV2/speed_of_light])) ) | (normalised_wavelength[mask_narrow]<desired_refwl+norm_limits_min+5 + np.amax(np.array([desired_refwl*RV1/speed_of_light, desired_refwl*RV2/speed_of_light])) ))
                
            else: clip_mask2=normalised_flux==normalised_flux; off=0; mask_narrow=normalised_wavelength==normalised_wavelength
            
            
            if not (starType1=="quadLorentz" and starType2=="quadLorentz"):
            
                # fit for RV1 RV2 with a narrowed range around the line cores for a final result
                if not issue_here:
                    returnvars = fit_bootstrap2([RV1, RV2, off],  normalised_wavelength[mask_narrow][clip_mask2],  normalised_flux[mask_narrow][clip_mask2],  normalised_err[mask_narrow][clip_mask2], bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1], -0.2, 0.2], num_its=200, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, model_range=model_range, Common_RV_both_stars=commonRV_1_2, expected_RV1err=RV1err,expected_RV2err=RV2err)
                else:
                    returnvars = fit_bootstrap2([RV1, RV2, off],  normalised_wavelength[mask_narrow],  normalised_flux[mask_narrow],  normalised_err[mask_narrow], bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1], -20.0, 20.0], num_its=200, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, model_range=model_range, Common_RV_both_stars=commonRV_1_2, expected_RV1err=RV1err,expected_RV2err=RV2err)
            

                try: RV1, RV1err, RV2, RV2err, off, off_err = returnvars
                except: 
                    RV1, RV1err, off, off_err = returnvars
                    RV2=RV1
                    RV2err=RV1err
            
            
            else:
                try:
                    returnvars = fit_bootstrap([RV1, RV2, extra_norm, coff],  normalised_wavelength[mask_narrow][clip_mask2],  normalised_flux[mask_narrow][clip_mask2],  normalised_err[mask_narrow][clip_mask2], bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1], -1,1, -20, 20], num_its=100, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, Common_RV_both_stars=commonRV_1_2)
                    
                except:
                    issue_here=True
                    if not (starType1=="quadLorentz" and starType2=="quadLorentz"):
                        returnvars = fit_bootstrap([RV1, RV2, extra_norm, coff],  normalised_wavelength,  normalised_flux,  normalised_err, bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1], -1,1, -20, 20], num_its=100, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, Common_RV_both_stars=commonRV_1_2)
                    else:
                        returnvars = fit_bootstrap([RV1, RV2],  normalised_wavelength,  normalised_flux,  normalised_err, bounds=[rvbound1[0], rvbound1[1], rvbound2[0], rvbound2[1]], num_its=100, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, Common_RV_both_stars=commonRV_1_2)
            
                try: 
                    if not (starType1=="quadLorentz" and starType2=="quadLorentz"):  RV1, RV1err, RV2, RV2err, extra_norm, extra_norm_err, coff, coff_err = returnvars
                    else:  RV1, RV1err, RV2, RV2err = returnvars
                except: RV1, RV1err, extra_norm, extra_norm_err, coff, coff_err = returnvars
                
                off=0
            
            
            
            if starType1=="quadLorentz" and starType2=="quadLorentz":
                
                mask=~((normalised_wavelength>desired_refwl-4) & (normalised_wavelength<desired_refwl+4))
                    
                xquad = normalised_wavelength[mask];   yquad = normalised_flux[mask];   yquaderr = normalised_err[mask]
                    
                # Normalize X values (this prevents numerical issues)
                x_mean = npmean(xquad)
                max_minus_min = npamax(xquad) - npamin(xquad)
                x_scaled = (xquad - x_mean) / max_minus_min  # Scale to [-1, 1]
                    
                initial_guess = polyfit(x_scaled, yquad, deg=poly_degree)[::-1]  # Reverse order
                    
                # Optimize with constraints
                constraints = [{'type': 'eq', 'fun': partial(constraint_function, x_scaled=x_scaled)}]
                result = opt.minimize(partial(loss, x_scaled=x_scaled, yquad=yquad, yquaderr=yquaderr), initial_guess, constraints=constraints, method='SLSQP', options={'ftol': 1e-4, 'gtol': 1e-4, 'maxiter': 200})
                    
                # Get final polynomial fit
                quad = fast_polyval(result.x[::-1], (normalised_wavelength - x_mean) / max_minus_min)
                    
                
                smear_model_spectrum_star1 += 0.5*interp(model_wl1, normalised_wavelength, quad)
                smear_model_spectrum_star2 += 0.5*interp(model_wl2, normalised_wavelength, quad)
                
                
            
            if True:
                print(in_fi)
                fig, (ax, ax2) = plt.subplots(1,2,gridspec_kw={'width_ratios': [2, 2]})
                plt.subplots_adjust(wspace=.0)
                #fine_grid = nplinspace(npamin(normalised_wavelength[clip_mask]), npamax(normalised_wavelength[clip_mask]), 100000)
                fine_grid = normalised_wavelength
                if sys.argv[1]=="RV_gauss":
                    ax.plot(normalised_wavelength,  normalised_flux,c='k');    ax.plot(fine_grid, WDmodel(fine_grid, RV1, RV2, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, gauss1=smear_gauss1, gauss2=smear_gauss2) + off, c='orange'); ax.scatter(normalised_wavelength[~clip_mask],  normalised_flux[~clip_mask],c='pink'); ax.scatter(normalised_wavelength[mask_narrow][~clip_mask2],  normalised_flux[mask_narrow][~clip_mask2],c='r')
                else:
                    ax.plot(normalised_wavelength,  normalised_flux,c='k');    ax.plot(fine_grid, WDmodel(fine_grid, RV1, RV2, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2) + off, c='orange'); ax.scatter(normalised_wavelength[~clip_mask],  normalised_flux[~clip_mask],c='pink'); ax.scatter(normalised_wavelength[mask_narrow][~clip_mask2],  normalised_flux[mask_narrow][~clip_mask2],c='r')
                mask_small_region = (normalised_wavelength>desired_refwl-17) & (normalised_wavelength<desired_refwl+17)
                if sys.argv[1]=="RV_gauss":
                    ax2.plot(speed_of_light*(normalised_wavelength[mask_small_region]-desired_refwl)/desired_refwl,  normalised_flux[mask_small_region],c='k');    ax2.plot(speed_of_light*(fine_grid[mask_small_region]-desired_refwl)/desired_refwl, WDmodel(fine_grid, RV1, RV2, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2, gauss1=smear_gauss1, gauss2=smear_gauss2)[mask_small_region] + off, c='orange')
                else:
                    ax2.plot(speed_of_light*(normalised_wavelength[mask_small_region]-desired_refwl)/desired_refwl,  normalised_flux[mask_small_region],c='k');    ax2.plot(speed_of_light*(fine_grid[mask_small_region]-desired_refwl)/desired_refwl, WDmodel(fine_grid, RV1, RV2, wl1=model_wl1, spec1=smear_model_spectrum_star1, wl2=model_wl2, spec2=smear_model_spectrum_star2)[mask_small_region]  + off, c='orange')
                ax2.set_yticks([])
                ax2.axvline(-100,c='grey', ls='dotted')
                ax2.axvline(100,c='grey', ls='dotted')
                

                if sys.argv[1]=="RV_double_commonRV12":
                    rvfilename.append([in_fi, aHJD, RV1, RV1err])
                    ax.set_title("V1=" + str(np.round(RV1,1)) + " " + str(np.round(RV1err,1)));  ax2.set_title("R=" + str(inp_resolution) + "  " + str(rvbound1))
                else:
                    rvfilename.append([in_fi, aHJD, RV1, RV1err, RV2, RV2err])
                    ax.set_title("V1=" + str(np.round(RV1,1)) + " " + str(np.round(RV1err,1)) + "  V2=" + str(np.round(RV2,1)) + " " + str(np.round(RV2err,1)));  ax2.set_title("R=" + str(inp_resolution) + "  " + str(rvbound1) + "|" + str(rvbound2))
                
                
                #np.savetxt("RVfits/plot/"+str(in_fi)+"data.dat", np.array([normalised_wavelength,  normalised_flux]).T)
                #np.savetxt("RVfits/plot/"+str(in_fi)+"fit1.dat", np.array([model_wl1, smear_model_spectrum_star1]).T, header=str(off))
                #np.savetxt("RVfits/plot/"+str(in_fi)+"fit2.dat", np.array([model_wl2, smear_model_spectrum_star2]).T, header=str(off))
                
                plt.savefig("RVfits/"+str(in_fi)+".png", dpi=300) # plt.show()
                plt.clf(); plt.close()
                
        except Exception as e: 
            print_exception_details(e)
            rvfilename.append([in_fi, aHJD, "ERROR", "ERROR", "ERROR", "ERROR"])
    
    np.savetxt("RVfits/RVfit_results.dat", rvfilename, fmt="%s", delimiter="\t",header="Observed velocities. For relativistic correction to source velocities (~0.14kms-1 for 200kms-1),  vsource  =  vobs/np.sqrt((1+vobs/c) / (1-vobs/c) )")
    
    
    

elif sys.argv[1]=="ATM":
    pool = MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
    print("Starting MCMC")
    
    if starType1.startswith("D") or starType2.startswith("D"):
        if forced_teff1!=0 and forced_logg1!=0 and forced_teff2!=0 and forced_logg2!=0 and starType1=="DA" and starType2=="DA":
            Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DAgrids(forced_teff1, forced_logg1)
            Grav2_N, wl_all2_N, flux2_N, Teff2_N = return_DAgrids(forced_teff2, forced_logg2)
            model_wl1, model_spectrum_star1 = return_model_spectrum_DA(wl_all1_N, 5000, -5000, 5000, Grav1_N, flux1_N, Teff1_N, forced_teff1, forced_logg1)
            model_wl2, model_spectrum_star2 = return_model_spectrum_DA(wl_all2_N, 5000, -5000, 5000, Grav2_N, flux2_N, Teff2_N, forced_teff2, forced_logg2)
            
            R1 = get_MTR(forced_teff1, logg=forced_logg1, return_R_from_T_logg=True, loaded_Istrate=loaded_Istrate, loaded_CO=loaded_CO, loaded_Althaus=loaded_Althaus)
            R2 = get_MTR(forced_teff2, logg=forced_logg2, return_R_from_T_logg=True, loaded_Istrate=loaded_Istrate, loaded_CO=loaded_CO, loaded_Althaus=loaded_Althaus)
            inputScaling= (R1**2/R2**2)
            
            
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool, args=[[input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, list_norm_wl_grids, list_normalised_flux, list_normalised_err,  resolutions, used_RV_boundaries, HJD_values, sigma_clip, npunique(reference_wl), model_wl1, model_spectrum_star1, model_wl2, model_spectrum_star2, inputScaling]])
            
            
        else:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool, args=[[input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, list_norm_wl_grids, list_normalised_flux, list_normalised_err,  resolutions, used_RV_boundaries, HJD_values, sigma_clip, npunique(reference_wl), None, None, None, None, None]])
    
    elif starType1=="quadLorentz" and starType2=="quadLorentz":
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_quad_lorentz, pool=pool, args=[[input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, list_norm_wl_grids, list_normalised_flux, list_normalised_err,  resolutions, used_RV_boundaries, HJD_values, sigma_clip, npunique(reference_wl), None, None, None, None, None]])
    
    else:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_gauss_lorentz, pool=pool, args=[[input_files, share_rv, reference_wl, cut_Ha_all, normaliseHa_all, list_norm_wl_grids, list_normalised_flux, list_normalised_err,  resolutions, used_RV_boundaries, HJD_values, sigma_clip, npunique(reference_wl), None, None, None, None, None]])
        
        
    if nwalkers/ndim<2:   raise ValueError("Your number of walkers is less than twice the size of the number of dimensions. Num dimensions = ", ndim, " nwalkers = ", nwalkers)
    
    pos, prob, state = sampler.run_mcmc(p0,burnin,progress=True)
    samples = sampler.flatchain
    np.savetxt("out/MCMC_samples_burnin.dat",samples)
    np.savetxt("out/MCMC_burninpos.dat",pos)
    print("Finished burn-in")



    #Make a corner plot which shows values
    samples=nploadtxt("out/MCMC_samples_burnin.dat")
    
    if not "run_double_one_by_oneMCMC" in sys.argv:

        if plot_burnin:
            fig = corner.corner(samples, labels=p0labels, range=p0range, bins=30,smooth=True, quantiles=[0.16, 0.5, 0.84], show_titles=True, labels_args={"fontsize": 40}, title_fmt = '.3f')
            fig.savefig("out/corner_burnin.pdf")
            del fig
    plt.close();  plt.clf()

    sampler.reset()
    pos=nploadtxt("out/MCMC_burninpos.dat")
    

    pos, prob, state = sampler.run_mcmc(pos,nsteps,progress=True)
    print("MCMC completed")
    pool.close()

    if not "run_double_one_by_oneMCMC" in sys.argv:
        samples = sampler.flatchain
        np.savetxt("out/MCMC_samples.dat",samples)

        #Make a corner plot which shows values
        samples=nploadtxt("out/MCMC_samples.dat")
        if plot_corner:
            fig2 = corner.corner(samples, labels=p0labels, range=p0range, bins=30,smooth=True, quantiles=[0.16, 0.5, 0.84], show_titles=True, labels_args={"fontsize": 40}, title_fmt = '.3f')
            fig2.savefig("out/corner.pdf")
            del fig2
    plt.close();  plt.clf()
    
    
    
    #try: os.mkdir("out/corner")
    #except: None
    
    #try:
    #    for ll, (fi, arv, arvBounds1, arvBounds2) in enumerate(zip(input_files, share_rv, RV_boundaries1, RV_boundaries2)):
    #    
    #        fig3 = corner.corner(samples.T[ll*2:ll*2+2].T, labels=["RV1", "RV2"], range=[arvBounds1, arvBounds2], bins=30, smooth=True, quantiles=[0.16, 0.5, 0.84], show_titles=True, labels_args={"fontsize": 40}, title_fmt = '.3f')
    #        fig3.savefig("out/corner/corner_"+str(fi).split("/")[-1]+".png")
    #        del fig3
    #except: None



    lines_to_write = []
    print(p0labels)
    if "T1" in p0labels:
        args = npargwhere(p0labels=="T1")[0][0];          Teff1_result = samples[::,args]
        T1_med = np.percentile(Teff1_result, 50, axis=0);  T1_min = np.percentile(Teff1_result, 16, axis=0);   T1_max = np.percentile(Teff1_result, 84, axis=0)
        lines_to_write.append("Teff1:\n")
        lines_to_write.append(str(T1_med) + "\t" + str(np.percentile(Teff1_result, 16, axis=0) - T1_med) + "\t" + str(np.percentile(Teff1_result, 84, axis=0) - T1_med) + "\n")
    else:   T1_med = forced_teff1;   T1_min = 0;   T1_max=0
        
    if "T2" in p0labels:
        args = npargwhere(p0labels=="T2")[0][0];          Teff2_result = samples[::,args]
        T2_med = np.percentile(Teff2_result, 50, axis=0);  T2_min = np.percentile(Teff2_result, 16, axis=0);   T2_max = np.percentile(Teff2_result, 84, axis=0)
        lines_to_write.append("Teff2:\n")
        lines_to_write.append(str(T2_med) + "\t" + str(np.percentile(Teff2_result, 16, axis=0) - T2_med) + "\t" + str(np.percentile(Teff2_result, 84, axis=0) - T2_med) + "\n")
    else:   T2_med = forced_teff2;   T2_min = 0;   T2_max=0

    if "logg1" in p0labels:
        args = npargwhere(p0labels=="logg1")[0][0];       logg1_result = samples[::,args]
        logg1_med = np.percentile(logg1_result, 50, axis=0);   logg1_min = np.percentile(logg1_result, 16, axis=0);   logg1_max = np.percentile(logg1_result, 84, axis=0)
        lines_to_write.append("Logg1:\n")
        lines_to_write.append(str(logg1_med) + "\t" + str(np.percentile(logg1_result, 16, axis=0) - logg1_med) + "\t" + str(np.percentile(logg1_result, 84, axis=0) - logg1_med) + "\n")
    else:   logg1_med = forced_logg1;   logg1_min = 0;   logg1_max=0

    if "logg2" in p0labels:
        args = npargwhere(p0labels=="logg2")[0][0];       logg2_result = samples[::,args]
        logg2_med = np.percentile(logg2_result, 50, axis=0);   logg2_min = np.percentile(logg2_result, 16, axis=0);   logg2_max = np.percentile(logg2_result, 84, axis=0)
        lines_to_write.append("Logg2:\n")
        lines_to_write.append(str(logg2_med) + "\t" + str(np.percentile(logg2_result, 16, axis=0) - logg2_med) + "\t" + str(np.percentile(logg2_result, 84, axis=0) - logg2_med) + "\n")
    else:   logg2_med = forced_logg2;   logg2_min = 0;   logg2_max=0

    if "H/He1" in p0labels:
        args = npargwhere(p0labels=="H/He1")[0][0];       HoverHe1_result = samples[::,args]
        HoverHe1_med = np.percentile(HoverHe1_result, 50, axis=0);   HoverHe1_min = np.percentile(HoverHe1_result, 16, axis=0);   HoverHe1_max = np.percentile(HoverHe1_result, 84, axis=0)
        lines_to_write.append("H/He1:\n")
        lines_to_write.append(str(HoverHe1_med) + "\t" + str(np.percentile(HoverHe1_result, 16, axis=0) - HoverHe1_med) + "\t" + str(np.percentile(HoverHe1_result, 84, axis=0) - HoverHe1_med) + "\n")
    else:   HoverHe1_med = forced_HoverHe1;   HoverHe1_min = 0;   HoverHe1_max=0

    if "H/He2" in p0labels:
        args = npargwhere(p0labels=="H/He2")[0][0];       HoverHe2_result = samples[::,args]
        HoverHe2_med = np.percentile(HoverHe2_result, 50, axis=0);   HoverHe2_min = np.percentile(HoverHe2_result, 16, axis=0);   HoverHe2_max = np.percentile(HoverHe2_result, 84, axis=0)
        lines_to_write.append("H/He2:\n")
        lines_to_write.append(str(HoverHe2_med) + "\t" + str(np.percentile(HoverHe2_result, 16, axis=0) - HoverHe2_med) + "\t" + str(np.percentile(HoverHe2_result, 84, axis=0) - HoverHe2_med) + "\n")
    else:   HoverHe2_med = forced_HoverHe2;   HoverHe2_min = 0;   HoverHe2_max=0


    if "K1" in p0labels:
        args = npargwhere(p0labels=="K1")[0][0];          K1_result = samples[::,args]
        K1_med = np.percentile(K1_result, 50, axis=0);     K1_min = np.percentile(K1_result, 16, axis=0);   K1_max = np.percentile(K1_result, 84, axis=0)
        lines_to_write.append("K1:\n")
        lines_to_write.append(str(K1_med) + "\t" + str(np.percentile(K1_result, 16, axis=0) - K1_med) + "\t" + str(np.percentile(K1_result, 84, axis=0) - K1_med) + "\n")
    else:   K1_med = forced_K1;   K1_min = 0;   K1_max=0

    if "K2" in p0labels:
        args = npargwhere(p0labels=="K2")[0][0];          K2_result = samples[::,args]
        K2_med = np.percentile(K2_result, 50, axis=0);     K2_min = np.percentile(K2_result, 16, axis=0);   K2_max = np.percentile(K2_result, 84, axis=0)
        lines_to_write.append("K2:\n")
        lines_to_write.append(str(K2_med) + "\t" + str(np.percentile(K2_result, 16, axis=0) - K2_med) + "\t" + str(np.percentile(K2_result, 84, axis=0) - K2_med) + "\n")
    else:   K2_med = forced_K2;   K2_min = 0;   K2_max=0

    if "Vg1" in p0labels:
        args = npargwhere(p0labels=="Vg1")[0][0];         Vgamma1_result = samples[::,args]
        Vgamma1_med = np.percentile(Vgamma1_result, 50, axis=0);   Vgamma1_min = np.percentile(Vgamma1_result, 16, axis=0);   Vgamma1_max = np.percentile(Vgamma1_result, 84, axis=0)
        lines_to_write.append("Vg1:\n")
        lines_to_write.append(str(Vgamma1_med) + "\t" + str(np.percentile(Vgamma1_result, 16, axis=0) - Vgamma1_med) + "\t" + str(np.percentile(Vgamma1_result, 84, axis=0) - Vgamma1_med) + "\n")
    else:   Vgamma1_med = forced_Vgamma1;   Vgamma1_min = 0;   Vgamma1_max=0

    if "Vg2" in p0labels:
        args = npargwhere(p0labels=="Vg2")[0][0];         Vgamma2_result = samples[::,args]
        Vgamma2_med = np.percentile(Vgamma2_result, 50, axis=0);   Vgamma2_min = np.percentile(Vgamma2_result, 16, axis=0);   Vgamma2_max = np.percentile(Vgamma2_result, 84, axis=0)
        lines_to_write.append("Vg2:\n")
        lines_to_write.append(str(Vgamma2_med) + "\t" + str(np.percentile(Vgamma2_result, 16, axis=0) - Vgamma2_med) + "\t" + str(np.percentile(Vgamma2_result, 84, axis=0) - Vgamma2_med) + "\n")
    else:   Vgamma2_med = forced_Vgamma2;   Vgamma2_min = 0;   Vgamma2_max=0

    if "Scaling" in p0labels:
        args = npargwhere(p0labels=="Scaling")[0][0];     Scaling_result = samples[::,args]
        Scaling_med = np.percentile(Scaling_result, 50, axis=0);   Scaling_min = np.percentile(Scaling_result, 16, axis=0);   Scaling_max = np.percentile(Scaling_result, 84, axis=0)
        lines_to_write.append("Scaling:\n")
        lines_to_write.append(str(Scaling_med) + "\t" + str(np.percentile(Scaling_result, 16, axis=0) - Scaling_med) + "\t" + str(np.percentile(Scaling_result, 84, axis=0) - Scaling_med) + "\n")
    else:   Scaling_med = forced_Scaling;   Scaling_min = 0;   Scaling_max=0
    
    if "Parallax" in p0labels:
        args = npargwhere(p0labels=="Parallax")[0][0];     parallax_result = samples[::,args]
        parallax_med = np.percentile(parallax_result, 50, axis=0);   parallax_min = np.percentile(parallax_result, 16, axis=0);   parallax_max = np.percentile(parallax_result, 84, axis=0)
        lines_to_write.append("Parallax:\n")
        lines_to_write.append(str(parallax_med) + "\t" + str(np.percentile(parallax_result, 16, axis=0) - parallax_med) + "\t" + str(np.percentile(parallax_result, 84, axis=0) - parallax_med) + "\n")
    
    
    
    if "A1_1" in p0labels:
        args = npargwhere(p0labels=="A1_1")[0][0];         A1_1_result = samples[::,args]
        A1_1_med = np.percentile(A1_1_result, 50, axis=0);   A1_1_min = np.percentile(A1_1_result, 16, axis=0);   A1_1_max = np.percentile(A1_1_result, 84, axis=0)
        lines_to_write.append("A1_1:\n")
        lines_to_write.append(str(A1_1_med) + "\t" + str(np.percentile(A1_1_result, 16, axis=0) - A1_1_med) + "\t" + str(np.percentile(A1_1_result, 84, axis=0) - A1_1_med) + "\n")
    
    if "sig1_1" in p0labels:
        args = npargwhere(p0labels=="sig1_1")[0][0];         sig1_1_result = samples[::,args]
        sig1_1_med = np.percentile(sig1_1_result, 50, axis=0);   sig1_1_min = np.percentile(sig1_1_result, 16, axis=0);   sig1_1_max = np.percentile(sig1_1_result, 84, axis=0)
        lines_to_write.append("sig1_1:\n")
        lines_to_write.append(str(sig1_1_med) + "\t" + str(np.percentile(sig1_1_result, 16, axis=0) - sig1_1_med) + "\t" + str(np.percentile(sig1_1_result, 84, axis=0) - sig1_1_med) + "\n")
    
    if "A1_2" in p0labels:
        args = npargwhere(p0labels=="A1_2")[0][0];         A1_2_result = samples[::,args]
        A1_2_med = np.percentile(A1_2_result, 50, axis=0);   A1_2_min = np.percentile(A1_2_result, 16, axis=0);   A1_2_max = np.percentile(A1_2_result, 84, axis=0)
        lines_to_write.append("A1_2:\n")
        lines_to_write.append(str(A1_2_med) + "\t" + str(np.percentile(A1_2_result, 16, axis=0) - A1_2_med) + "\t" + str(np.percentile(A1_2_result, 84, axis=0) - A1_2_med) + "\n")
    
    if "sig1_2" in p0labels:
        args = npargwhere(p0labels=="sig1_2")[0][0];         sig1_2_result = samples[::,args]
        sig1_2_med = np.percentile(sig1_2_result, 50, axis=0);   sig1_2_min = np.percentile(sig1_2_result, 16, axis=0);   sig1_2_max = np.percentile(sig1_2_result, 84, axis=0)
        lines_to_write.append("sig1_2:\n")
        lines_to_write.append(str(sig1_2_med) + "\t" + str(np.percentile(sig1_2_result, 16, axis=0) - sig1_2_med) + "\t" + str(np.percentile(sig1_2_result, 84, axis=0) - sig1_2_med) + "\n")
    
    if "A2_1" in p0labels:
        args = npargwhere(p0labels=="A2_1")[0][0];         A2_1_result = samples[::,args]
        A2_1_med = np.percentile(A2_1_result, 50, axis=0);   A2_1_min = np.percentile(A2_1_result, 16, axis=0);   A2_1_max = np.percentile(A2_1_result, 84, axis=0)
        lines_to_write.append("A2_1:\n")
        lines_to_write.append(str(A2_1_med) + "\t" + str(np.percentile(A2_1_result, 16, axis=0) - A2_1_med) + "\t" + str(np.percentile(A2_1_result, 84, axis=0) - A2_1_med) + "\n")
    
    if "sig2_1" in p0labels:
        args = npargwhere(p0labels=="sig2_1")[0][0];         sig2_1_result = samples[::,args]
        sig2_1_med = np.percentile(sig2_1_result, 50, axis=0);   sig2_1_min = np.percentile(sig2_1_result, 16, axis=0);   sig2_1_max = np.percentile(sig2_1_result, 84, axis=0)
        lines_to_write.append("sig2_1:\n")
        lines_to_write.append(str(sig2_1_med) + "\t" + str(np.percentile(sig2_1_result, 16, axis=0) - sig2_1_med) + "\t" + str(np.percentile(sig2_1_result, 84, axis=0) - sig2_1_med) + "\n")
    
    if "A2_2" in p0labels:
        args = npargwhere(p0labels=="A2_2")[0][0];         A2_2_result = samples[::,args]
        A2_2_med = np.percentile(A2_2_result, 50, axis=0);   A2_2_min = np.percentile(A2_2_result, 16, axis=0);   A2_2_max = np.percentile(A2_2_result, 84, axis=0)
        lines_to_write.append("A2_2:\n")
        lines_to_write.append(str(A2_2_med) + "\t" + str(np.percentile(A2_2_result, 16, axis=0) - A2_2_med) + "\t" + str(np.percentile(A2_2_result, 84, axis=0) - A2_2_med) + "\n")
    
    if "sig2_2" in p0labels:
        args = npargwhere(p0labels=="sig2_2")[0][0];         sig2_2_result = samples[::,args]
        sig2_2_med = np.percentile(sig2_2_result, 50, axis=0);   sig2_2_min = np.percentile(sig2_2_result, 16, axis=0);   sig2_2_max = np.percentile(sig2_2_result, 84, axis=0)
        lines_to_write.append("sig2_2:\n")
        lines_to_write.append(str(sig2_2_med) + "\t" + str(np.percentile(sig2_2_result, 16, axis=0) - sig2_2_med) + "\t" + str(np.percentile(sig2_2_result, 84, axis=0) - sig2_2_med) + "\n")





    # save just the median for individual epoch pngs with fit and observations
    if forced_Scaling==False:    Scaling_med = Scaling_med
    elif forced_Scaling=="WD": 
        if starType1=="DA":      R1 = get_MTR(T1_med, logg=logg1_med, return_R_from_T_logg=True, loaded_Istrate=loaded_Istrate, loaded_CO=loaded_CO, loaded_Althaus=loaded_Althaus)
        elif starType1=="DBA" or starType1=="DC" or starType1=="DB":   R1 = get_MTR_DB(T1_med, logg=logg1_med)
        if starType2=="DA":      R2 = get_MTR(T2_med, logg=logg2_med, return_R_from_T_logg=True, loaded_Istrate=loaded_Istrate, loaded_CO=loaded_CO, loaded_Althaus=loaded_Althaus)
        elif starType2=="DBA" or starType2=="DC" or starType2=="DB":   R2 = get_MTR_DB(T2_med, logg=logg2_med)
        Scaling_med=(R1**2/R2**2)
    else: Scaling_med=forced_Scaling




    allRV1s, allRV2s = [], []
    allRV1s_minerr, allRV1s_maxerr, allRV2s_minerr, allRV2s_maxerr = [], [], [], []
    if "RV1_0" in p0labels:
        num_start_RVs = npargwhere(p0labels=="RV1_0")[0][0]
        if "Scaling" in p0labels or "Parallax" in p0labels:   RVvals_result = samples[::,num_start_RVs:-1]
        else: RVvals_result = samples[::,num_start_RVs:]
        for cnt, item in enumerate(RVvals_result.T):
            if cnt%2==0:   allRV1s.append(np.percentile(item, 50, axis=0));   allRV1s_minerr.append(np.percentile(item, 16, axis=0));   allRV1s_maxerr.append(np.percentile(item, 84, axis=0))
            else:          allRV2s.append(np.percentile(item, 50, axis=0));   allRV2s_minerr.append(np.percentile(item, 16, axis=0));   allRV2s_maxerr.append(np.percentile(item, 84, axis=0))
    else:
        for aHJD in HJD_values:
            phi = ((aHJD - forced_T0)%forced_P0)/forced_P0
            allRV1s.append(Vgamma1_med + K1_med*np_sin(2*np_pi*phi));   allRV1s_minerr.append(0);   allRV1s_maxerr.append(0)
            allRV2s.append(Vgamma2_med + K2_med*np_sin(2*np_pi*phi));   allRV2s_minerr.append(0);   allRV2s_maxerr.append(0)
            


    lines_to_write.append("RV1:\n")
    for arvresult1, arvresult1_min, arvresult1_max, ahjd in zip(allRV1s, allRV1s_minerr, allRV1s_maxerr, HJD_values):
        if forced_K1==False:
            lines_to_write.append(str(ahjd) + "\t" + str(arvresult1) + "\t" + str(arvresult1_min) + "\t" + str(arvresult1_max) + "\n")
        else:
            lines_to_write.append(str(ahjd) + "\t" + str(arvresult1) + "\n")

    lines_to_write.append("RV2:\n")
    for arvresult2, arvresult2_min, arvresult2_max, ahjd in zip(allRV2s, allRV2s_minerr, allRV2s_maxerr, HJD_values):
        if forced_K2==False:     lines_to_write.append(str(ahjd) + "\t" + str(arvresult2) + "\t" + str(arvresult2_min) + "\t" + str(arvresult2_max) + "\n")
        else:                    lines_to_write.append(str(ahjd) + "\t" + str(arvresult2) + "\n")
    
    
    try: os.mkdir("out")
    except: None
    
    
    if "run_double_one_by_oneMCMC" in sys.argv:  
        with open("out/result"+str(int(sys.argv[-1])+1)+".out", 'w') as output_file:
            for aline in lines_to_write:
                output_file.write(aline)
    else:
        with open("out/result.out", 'w') as output_file:
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



        if starType1.startswith("D") or starType2.startswith("D"):
            if "Teff1:" in lll:
                T1_med=float(result[iii+1].split("\t")[0])
                T1_minus=float(result[iii+1].split("\t")[1])
                T1_plus=float(result[iii+1].split("\t")[2])
            elif "Teff2:" in lll:
                T2_med=float(result[iii+1].split("\t")[0])
                T2_minus=float(result[iii+1].split("\t")[1])
                T2_plus=float(result[iii+1].split("\t")[2])
            elif "Logg1:" in lll:
                logg1_med=float(result[iii+1].split("\t")[0])
                logg1_minus=float(result[iii+1].split("\t")[1])
                logg1_plus=float(result[iii+1].split("\t")[2])
            elif "Logg2:" in lll:
                logg2_med=float(result[iii+1].split("\t")[0])
                logg2_minus=float(result[iii+1].split("\t")[1])
                logg2_plus=float(result[iii+1].split("\t")[2])
            elif "H/He1:" in lll:
                HoverHe1_med=float(result[iii+1].split("\t")[0])
                HoverHe1_minus=float(result[iii+1].split("\t")[1])
                HoverHe1_plus=float(result[iii+1].split("\t")[2])
            elif "H/He2:" in lll:
                HoverHe2_med=float(result[iii+1].split("\t")[0])
                HoverHe2_minus=float(result[iii+1].split("\t")[1])
                HoverHe2_plus=float(result[iii+1].split("\t")[2])
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
            elif "A2_1:" in lll:
                A2_1_med=float(result[iii+1].split("\t")[0])
                A2_1_minus=float(result[iii+1].split("\t")[1])
                A2_1_plus=float(result[iii+1].split("\t")[2])
            elif "sig2_1:" in lll:
                sig2_1_med=float(result[iii+1].split("\t")[0])
                sig2_1_minus=float(result[iii+1].split("\t")[1])
                sig2_1_plus=float(result[iii+1].split("\t")[2])
            elif "A2_2:" in lll:
                A2_2_med=float(result[iii+1].split("\t")[0])
                A2_2_minus=float(result[iii+1].split("\t")[1])
                A2_2_plus=float(result[iii+1].split("\t")[2])
            elif "sig2_2:" in lll:
                sig2_2_med=float(result[iii+1].split("\t")[0])
                sig2_2_minus=float(result[iii+1].split("\t")[1])
                sig2_2_plus=float(result[iii+1].split("\t")[2])
            
    
    if forced_teff1!=0: T1_med=forced_teff1
    if forced_teff2!=0: T2_med=forced_teff2
    if forced_logg1!=0: logg1_med=forced_logg1
    if forced_logg2!=0: logg2_med=forced_logg2
    if forced_HoverHe1!=0: HoverHe1_med=forced_HoverHe1
    if forced_HoverHe2!=0: HoverHe2_med=forced_HoverHe2
    
    if "run_double_one_by_oneMCMC" in sys.argv:
        filesout = os.listdir(os.getcwd()+"/out")
        files_one_by_one = []
        for ff in filesout:
            if "result" in ff and ".out" in ff:
                files_one_by_one.append(ff)
        
        try: files_one_by_one.remove("result_one_by_one.out")
        except: None
        
        try: files_one_by_one.remove("result.out")
        except: None
        
        biggest_num=-1
        for ff in files_one_by_one:
            num = int(ff.split("result")[-1].split(".out")[0])
            if num>biggest_num:
                biggest_num=num
        
        res = "result"+str(biggest_num)+".out"
        
        result = open(os.getcwd()+"/out/"+res).readlines()
        
        allRV1s, allRV2s = [], []
        allRV1sminus, allRV2sminus =[], []
        allRV1smax, allRV2smax =[], []
        
        for iii, lll in enumerate(result):
            if "RV1:" in lll:
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
        
        
        
    
    
    
    if starType1.startswith("D"):
        if forced_Scaling==False:    Scaling_med = Scaling_med
        elif forced_Scaling=="WD": 
            if starType1=="DA":      R1 = get_MTR(T1_med, logg=logg1_med, return_R_from_T_logg=True, loaded_Istrate=loaded_Istrate, loaded_CO=loaded_CO, loaded_Althaus=loaded_Althaus)
            elif starType1=="DBA" or starType1=="DC"  or starType1=="DB":   R1 = get_MTR_DB(T1_med, logg=logg1_med)
            if starType2=="DA":      R2 = get_MTR(T2_med, logg=logg2_med, return_R_from_T_logg=True, loaded_Istrate=loaded_Istrate, loaded_CO=loaded_CO, loaded_Althaus=loaded_Althaus)
            elif starType2=="DBA" or starType2=="DC"  or starType2=="DB":   R2 = get_MTR_DB(T2_med, logg=logg2_med)
            Scaling_med=(R1**2/R2**2)  #(T1_med/T2_med)**4 * 
        else: Scaling_med=forced_Scaling
    
    else:
        Scaling_med = 1
    
    
    
    
    
    if starType1.startswith("D") or starType2.startswith("D"):
        if T2_med > T1_med:  raise ValueError("ERROR: T2 > T1. The fit is likely useless. Inspect the config file and the corner")

            
        for tempcnt, (temperature_star, logg_star, starType, HoverHestar) in enumerate(zip([T1_med,T2_med],[logg1_med,logg2_med], [starType1,starType2], [HoverHe1_med, HoverHe2_med])):
            minval=10000;   maxval=-10000
            if starType=="DA":
                if not pier_or_antoine=="pier3Dphot_antoine1Dspec":
                    if tempcnt==0:     Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DAgrids(temperature_star, logg_star)
                    elif tempcnt==1:   Grav2_N, wl_all2_N, flux2_N, Teff2_N = return_DAgrids(temperature_star, logg_star)
                else:
                    if tempcnt==0:     Grav1_N, wl_all1_N, flux1_N, Teff1_N, photGrav1_N, photwl_all1_N, photflux1_N, photTeff1_N = return_DAgrids(temperature_star, logg_star)
                    elif tempcnt==1:   Grav2_N, wl_all2_N, flux2_N, Teff2_N, photGrav2_N, photwl_all2_N, photflux2_N, photTeff2_N = return_DAgrids(temperature_star, logg_star)
                    
            

            elif starType=="DBA":
                if tempcnt==0:    HoverHestar=HoverHe1_med;    Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N = return_DBAgrids(temperature_star, logg_star, HoverHestar)
                elif tempcnt==1:  HoverHestar=HoverHe2_med;    Grav2_N, wl_all2_N, flux2_N, Teff2_N, HoverHe2_N = return_DBAgrids(temperature_star, logg_star, HoverHestar)
            
            elif starType=="DB":
                if tempcnt==0:     Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DBgrids(temperature_star, logg_star)
                elif tempcnt==1:   Grav2_N, wl_all2_N, flux2_N, Teff2_N = return_DBgrids(temperature_star, logg_star)
            
            elif starType=="DC":
                if tempcnt==0:     Grav1_N, wl_all1_N, flux1_N, Teff1_N = return_DCgrids(temperature_star, logg_star)
                elif tempcnt==1:   Grav2_N, wl_all2_N, flux2_N, Teff2_N = return_DCgrids(temperature_star, logg_star)
            
                




        # if relevant, plot the photometric solution
        if fit_phot_SED:
            if starType1=="DA" and starType2=="DA":
                if pier_or_antoine=="pier3Dphot_antoine1Dspec":
                    smeared_wl, smeared_flux, smeared_flux1, smeared_flux2 = Fit_phot.fit_phot_SED_double(photGrav1_N, photwl_all1_N, photflux1_N, photTeff1_N, None, photGrav2_N, photwl_all2_N, photflux2_N, photTeff2_N, None, T1_med, logg1_med, None, T2_med, logg2_med, None, min_wl=theminww, max_wl=themaxww, starType1=starType1, starType2=starType2, R1=R1, R2=R2, parallax=parallax_med, red=reddening_Ebv, return_indiv_stars=True)
                else:
                    smeared_wl, smeared_flux, smeared_flux1, smeared_flux2 = Fit_phot.fit_phot_SED_double(Grav1_N, wl_all1_N, flux1_N, Teff1_N, None, Grav2_N, wl_all2_N, flux2_N, Teff2_N, None, T1_med, logg1_med, None, T2_med, logg2_med, None, min_wl=theminww, max_wl=themaxww, starType1=starType1, starType2="DA", R1=R1, R2=R2, parallax=parallax_med, red=reddening_Ebv, return_indiv_stars=True)
            
            elif starType1=="DA" and starType2=="DBA":
                smeared_wl, smeared_flux, smeared_flux1, smeared_flux2 = Fit_phot.fit_phot_SED_double(Grav1_N, wl_all1_N, flux1_N, Teff1_N, None, Grav2_N, wl_all2_N, flux2_N, Teff2_N, HoverHe2_N, T1_med, logg1_med, None, T2_med, logg2_med, HoverHe2_med, min_wl=theminww, max_wl=themaxww,starType1=starType1, starType2="DA", R1=R1, R2=R2, parallax=parallax_med, red=reddening_Ebv, return_indiv_stars=True)
            elif starType1=="DBA" and starType2=="DA":
                smeared_wl, smeared_flux, smeared_flux1, smeared_flux2 = Fit_phot.fit_phot_SED_double(Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N, Grav2_N, wl_all2_N, flux2_N, Teff2_N, None, T1_med, logg1_med, HoverHe1_med, T2_med, logg2_med, None, min_wl=theminww, max_wl=themaxww,starType1=starType1, starType2=starType2, scaling=Scaling_med, return_indiv_stars=True)
            elif starType1=="DBA" and starType2=="DBA":
                smeared_wl, smeared_flux, smeared_flux1, smeared_flux2 = Fit_phot.fit_phot_SED_double(Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N, Grav2_N, wl_all2_N, flux2_N, Teff2_N, HoverHe2_N, T1_med, logg1_med, HoverHe1_med, T2_med, logg2_med, HoverHe2_med, min_wl=theminww, max_wl=themaxww,starType1=starType1, starType2=starType2, scaling=Scaling_med, return_indiv_stars=True)
                
            elif (starType1=="DA" and starType2=="DC") or (starType1=="DC" and starType2=="DA") or (starType1=="DA" and starType2=="DB") or (starType1=="DB" and starType2=="DA"):
                smeared_wl, smeared_flux, smeared_flux1, smeared_flux2 = Fit_phot.fit_phot_SED_double(Grav1_N, wl_all1_N, flux1_N, Teff1_N, None, Grav2_N, wl_all2_N, flux2_N, Teff2_N, None, T1_med, logg1_med, None, T2_med, logg2_med, None, min_wl=theminww, max_wl=themaxww, starType1=starType1, starType2=starType2, R1=R1, R2=R2, parallax=parallax_med, red=reddening_Ebv, return_indiv_stars=True)
            
            
            rchisq_phot, chisq_phot = Fit_phot.process_photometry_in_each_pb(smeared_wl, smeared_flux, sedfilter, sed_wl, sedflux, sedfluxe, filter_dict=filter_dict, plot_solution=True, theminww_plot=theminww+50, themaxww_plot=themaxww+50, specStar1=smeared_flux1, specStar2=smeared_flux2, return_points_for_phot_model=False)









    # go and plot the spectra with the fits, don't plot anything with a shared wavelength, unless it's the gaussian/lorentzian fits
    for ii, (in_fi, sh_rv, ref_wl, modha, cut_lim, norm_lim, mask_out_Ha_min, mask_out_Ha_max, aHJD, inp_resolution, sigclipspec) in enumerate(zip(input_files, share_rv, reference_wl, modelHa, cut_Ha_all, normaliseHa_all, mask_out_Ha_min_all, mask_out_Ha_max_all, HJD_values, resolutions, sigma_clip)):
        if ((starType1.startswith("GG") or starType1.startswith("LL") or starType1.startswith("GL")) and sh_rv==-1) or (sh_rv==-1 and not ii in share_rv):
            norm_limits_min, norm_limits_max = norm_lim
            cut_limits_min, cut_limits_max = cut_lim
            modelHa_min, modelHa_max = modha
            
            if high_RV_amp: 
                if sh_wl == -1:  rv_med1 = allRV1s[ii]
                else:            rv_med1 = allRV1s[sh_wl]
                dlam1 = ref_wl*rv_med1/speed_of_light
                modelHa_min+=dlam1; modelHa_max+=dlam1; cut_limits_min+=dlam1; cut_limits_max+=dlam1; norm_limits_min+=dlam1; norm_limits_max+=dlam1
                
                mask_input_spectrum = ((normalised_wavelength>ref_wl+cut_limits_min)   & (normalised_wavelength<ref_wl+cut_limits_max))
                normalised_wavelength, normalised_flux, normalised_err = normalised_wavelength[mask_input_spectrum], normalised_flux[mask_input_spectrum], normalised_err[mask_input_spectrum]
            
            fig, (ax, ax2) = plt.subplots(2)
            
            if starType1.startswith("D") and starType2.startswith("D"):
                if starType1=="DA" or starType1=="DC" or starType1=="DB": model_wl1, model_spectrum_star1 = return_model_spectrum_DA(wl_all1_N, ref_wl, cut_limits_min, cut_limits_max, Grav1_N, flux1_N, Teff1_N, T1_med, logg1_med)
                elif starType1=="DBA": model_wl1, model_spectrum_star1 = return_model_spectrum_DBA(wl_all1_N, ref_wl, cut_limits_min, cut_limits_max, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, T1_med, logg1_med, HoverHe1_med)
                
                if starType2=="DA" or starType2=="DC" or starType2=="DB": model_wl2, model_spectrum_star2 = return_model_spectrum_DA(wl_all2_N, ref_wl, cut_limits_min, cut_limits_max, Grav2_N, flux2_N, Teff2_N, T2_med, logg2_med)
                elif starType2=="DBA": model_wl2, model_spectrum_star2 = return_model_spectrum_DBA(wl_all2_N, ref_wl, cut_limits_min, cut_limits_max, Grav2_N, flux2_N, Teff2_N, HoverHe2_N, T2_med, logg2_med, HoverHe2_med)
            else:
                model_wl1 = nplinspace(cut_limits_min, cut_limits_max, int(cut_limits_max-cut_limits_min)*10) #+ ref_wl
                model_wl2 = nplinspace(cut_limits_min, cut_limits_max, int(cut_limits_max-cut_limits_min)*10) #+ ref_wl
                if starType1=="GG":     model_spectrum_star1 = gauss_2x_mcmc(model_wl1, A1_1_med, sig1_1_med, A1_2_med, sig1_2_med)
                elif starType1=="LL":   model_spectrum_star1 = lorenz_2x_mcmc(model_wl1, A1_1_med, sig1_1_med, A1_2_med, sig1_2_med)
                elif starType1=="GL":   model_spectrum_star1 = gauss_lorentz_mcmc(model_wl1, A1_1_med, sig1_1_med, A1_2_med, sig1_2_med)
                if starType2=="GG":     model_spectrum_star2 = gauss_2x_mcmc(model_wl2, A2_1_med, sig2_1_med, A2_2_med, sig2_2_med)
                elif starType2=="LL":   model_spectrum_star2 = lorentz_2x_mcmc(model_wl2, A2_1_med, sig2_1_med, A2_2_med, sig2_2_med)
                elif starType2=="GL":   model_spectrum_star2 = gauss_lorentz_mcmc(model_wl2, A2_1_med, sig2_1_med, A2_2_med, sig2_2_med)
                if starType1=="quadLorentz" and starType2=="quadLorentz":
                    model_spectrum_star1 = lorentz_2x_mcmc(model_wl1, A1_1_med, sig1_1_med, 0, 1)
                    model_spectrum_star2 = lorentz_2x_mcmc(model_wl2, A2_1_med, sig2_1_med, 0, 1)
                model_wl1, model_wl2 = model_wl1 + ref_wl, model_wl2 + ref_wl
            
            
            if sh_rv == -1:  rv_med1 = allRV1s[ii];  rv_med2 = allRV2s[ii]
            else:            rv_med1 = allRV1s[sh_rv];  rv_med2 =allRV2s[sh_rv]
            
            
            filename = input_files[ii]
            if "/" in filename:    filename=filename.split("/")[-1]
            normalised_wavelength = list_norm_wl_grids[ii];    normalised_flux = list_normalised_flux[ii];    normalised_err = list_normalised_err[ii]

            
            
            
            dlam1 = ref_wl*rv_med1/speed_of_light
            dlam2 = ref_wl*rv_med2/speed_of_light
            desired_range = (normalised_wavelength>ref_wl+modelHa_min)  &  (normalised_wavelength<ref_wl+modelHa_max)
            
            
            
            
            # interpolate the model onto the observation's wavelength grid
            f1 = interp(model_wl2, model_wl1+dlam1, model_spectrum_star1)
            f2 = interp(model_wl2, model_wl2+dlam2, model_spectrum_star2)
            
            interparr = Scaling_med*f1 + f2
            
            if starType1=="GG" or starType1=="LL" or starType1=="GL" or starType2=="GG" or starType2=="LL" or starType2=="GL":   interparr += 1
            
            # Smear to the resolution desired
            interparr = convolve(interparr, Gaussian1DKernel(stddev=0.5*(ref_wl/inp_resolution)/(model_wl2[10]-model_wl2[9])), boundary = 'extend')
            
            interparr = interp(normalised_wavelength, model_wl2, interparr)
        
            if starType1=="quadLorentz" and starType2=="quadLorentz":
                
                mask=~((normalised_wavelength>ref_wl-4) & (normalised_wavelength<ref_wl+4))
                    
                xquad = normalised_wavelength[mask];   yquad = normalised_flux[mask];   yquaderr = normalised_err[mask]
                    
                # Normalize X values (this prevents numerical issues)
                x_mean = npmean(xquad)
                max_minus_min = npamax(xquad) - npamin(xquad)
                x_scaled = (xquad - x_mean) / max_minus_min  # Scale to [-1, 1]
                    
                initial_guess = polyfit(x_scaled, yquad, deg=poly_degree)[::-1]  # Reverse order
                    
                # Optimize with constraints
                constraints = [{'type': 'eq', 'fun': partial(constraint_function, x_scaled=x_scaled)}]
                result = opt.minimize(partial(loss, x_scaled=x_scaled, yquad=yquad, yquaderr=yquaderr), initial_guess, constraints=constraints, method='SLSQP', options={'ftol': 1e-4, 'gtol': 1e-4, 'maxiter': 200})
                    
                # Get final polynomial fit
                quad = fast_polyval(result.x[::-1], (normalised_wavelength - x_mean) / max_minus_min)
                
                interparr += quad
        
            
            mask = ((normalised_wavelength>ref_wl+cut_limits_min) & (normalised_wavelength<ref_wl+norm_limits_min)) | ((normalised_wavelength>ref_wl+norm_limits_max)   & (normalised_wavelength<ref_wl+cut_limits_max))
                
            try:   m,c =  polyfit(normalised_wavelength[mask], interparr[mask], deg=1)
            except Exception as eee:
                print_exception_details(eee)
                plt.clf()
                plt.plot(normalised_wavelength, normalised_flux);    plt.plot(normalised_wavelength, interparr)
                plt.axvline(ref_wl+cut_limits_min, c='k');     plt.axvline(ref_wl+cut_limits_max, c='k')
                plt.axvline(ref_wl+norm_limits_min, c='r');    plt.axvline(ref_wl+norm_limits_max, c='r')
                plt.title(str() + in_fi + " " + str(ref_wl) + "    "+str(cut_limits_min) + "  " + str(norm_limits_min) + "  " + str(norm_limits_max) + "  " + str(cut_limits_max))
                plt.show();    plt.clf()
                m,c =  polyfit(normalised_wavelength[mask], interparr[mask], deg=1)
            
            
            #plt.plot(normalised_wavelength, interparr);   plt.plot(nplinspace(6000, 7000, 1000), m*nplinspace(6000, 7000, 1000) + c);  plt.show()
            interparr /= (m*normalised_wavelength + c)
            
            
            
            
            if sigclipspec!=-1:
                resid = normalised_flux-interparr
                resid_in_sig = resid/np.std(np.abs(resid[desired_range]))
                clip_mask = (np.abs(resid_in_sig)>sigclipspec) & desired_range & ((normalised_wavelength > ref_wl+5) | (normalised_wavelength < ref_wl-5))
            else: clip_mask=normalised_wavelength!=normalised_wavelength
            
            off =  polyfit(normalised_wavelength[desired_range & ~clip_mask], normalised_flux[desired_range & ~clip_mask] - interparr[desired_range & ~clip_mask], w=1/normalised_err[desired_range & ~clip_mask], deg=0)[0]
            interparr+=off
            
            
            
            
            
            mask1 = (normalised_wavelength>ref_wl+cut_limits_min) & (normalised_wavelength<ref_wl+norm_limits_min)
            mask2 = (normalised_wavelength>ref_wl+norm_limits_max) & (normalised_wavelength<ref_wl+cut_limits_max)
            
            if False:
                plt.clf();   plt.plot(normalised_wavelength[mask1], interparr[mask1], c='g');   plt.plot(normalised_wavelength[mask2], interparr[mask2], c='r');   plt.show()
            
            
            
            normwl_to_rv = speed_of_light * (normalised_wavelength-ref_wl)/ref_wl 

            mask_norm = ((normalised_wavelength>ref_wl+cut_limits_min) & (normalised_wavelength<ref_wl+norm_limits_min) ) | (normalised_wavelength<ref_wl+cut_limits_max) & ((normalised_wavelength>ref_wl+norm_limits_max))
            
            try:  os.mkdir("out/forPlotting")
            except: None
            np.savetxt("out/forPlotting/fit_"+in_fi+"_"+str(ref_wl)+".dat", np.array([normalised_wavelength[desired_range] - ref_wl, interparr[desired_range]]).T)
            np.savetxt("out/forPlotting/data_for_fit_"+in_fi+"_"+str(ref_wl)+".dat", np.array([normalised_wavelength - ref_wl, normalised_flux]).T)
                
            
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

            if starType1=="DBA":    ax.scatter(0,1,alpha=0, label="H/He1 = " + str(HoverHe1_med)); ax2.scatter(0,1,alpha=0, label="H/He1 = " + str(HoverHe1_med))
            if starType2=="DBA":    ax.scatter(0,1,alpha=0, label="H/He2 = " + str(HoverHe2_med)); ax2.scatter(0,1,alpha=0, label="H/He2 = " + str(HoverHe2_med))
            
            
            ax.scatter(0,1,alpha=0, label="V1 = " + str(np.round(rv_med1,3)));        ax2.scatter(0,1,alpha=0, label="V1 = " + str(np.round(rv_med1,3)))
            ax.scatter(0,1,alpha=0, label="V2 = " + str(np.round(rv_med2,3)));        ax2.scatter(0,1,alpha=0, label="V2 = " + str(np.round(rv_med2,3)))
            ax.scatter(0,1,alpha=0, label="wl = " + str(ref_wl));                 ax2.scatter(0,1,alpha=0, label="wl = " + str(ref_wl))
            if starType1.startswith("D"):
                ax.scatter(0,1,alpha=0, label="T1 = " + str(np.round(T1_med,0)));          ax2.scatter(0,1,alpha=0, label="T1 = " + str(np.round(T1_med,0)))
                ax.scatter(0,1,alpha=0, label="Logg1 = " + str(np.round(logg1_med,3)));    ax2.scatter(0,1,alpha=0, label="Logg1 = " + str(np.round(logg1_med,3)))
            if starType2.startswith("D"):
                ax.scatter(0,1,alpha=0, label="T2 = " + str(np.round(T2_med,0)));          ax2.scatter(0,1,alpha=0, label="T2 = " + str(np.round(T2_med,0)))
                ax.scatter(0,1,alpha=0, label="Logg2 = " + str(np.round(logg2_med,3)));    ax2.scatter(0,1,alpha=0, label="Logg2 = " + str(np.round(logg2_med,3)))
            ax2.set_xlabel("RV [kms-1]");
            ax.legend(loc="lower right", handlelength=0)
            #ax2.legend(loc="lower right", handlelength=0)
            trim_plot2_xvals_min, trim_plot2_xvals_max = -2500, 2500
            ax2.set_xlim(trim_plot2_xvals_min, trim_plot2_xvals_max)
            ax.set_ylim(0.3,1.2);    ax2.set_ylim(0.35,1.1)
            plt.savefig("out/"+filename.split(".dat")[0]+"_"+str(ref_wl)+".png", dpi=300)
            if ii in plot_fit and not starType1=="quadLorentz":   plt.show()
            plt.clf();    plt.close()





    plot_fancy=False
    if plot_fancy==True:
        minstored=np_inf
        maxstored=-np_inf
        import matplotlib
        import matplotlib.pyplot as plt
        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['mathtext.fontset'] = 'stix'
        matplotlib.rcParams['font.family'] = 'STIXGeneral'
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42
        plt.rcParams["xtick.direction"] = "in"
        plt.rcParams["ytick.direction"] = "in"
        plt.rcParams['xtick.top'] = True
        plt.rcParams['ytick.right'] = True
        matplotlib.rcParams["savefig.dpi"] = 100
        matplotlib.rcParams["font.size"] = 14
        plt.rcParams['xtick.minor.visible'] = True
        plt.rcParams['ytick.minor.visible'] = True
    
    #### Now plot just the shared rv combinations
    where_share_rv = (share_rv != -1)
    if not starType1.startswith("D"):
        where_share_rv[where_share_rv==True] = False
    plt.clf()
    for unique_shared_rv in npunique(share_rv[where_share_rv]):
        norm_wl_combinations, normalised_flux_combinations, normalised_err_combinations = [], [], []
        if plot_fancy==False:
            fig = plt.figure()
            gs = fig.add_gridspec(2, 2)
            ax = fig.add_subplot(gs[0, :])
            ax2 = fig.add_subplot(gs[1, 0])
            ax3 = fig.add_subplot(gs[1, 1])
        else:
            fig, ax = plt.subplots()
        
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

        
        toappendtext=[]
        for cn_shared, (in_fi, ref_wl, modha, cut_lim, norm_lim, normalised_wavelength, normalised_flux, normalised_err, inp_resolution, mask_out_Ha_min, mask_out_Ha_max, sigclipspec) in enumerate(zip(input_files_combinations, reference_wl_combinations, modelHa_combinations, cut_Ha_combinations, normaliseHa_combinations, norm_wl_combinations, normalised_flux_combinations, normalised_err_combinations, resolutions_combinations, mask_out_Ha_min_all_combinations, mask_out_Ha_max_all_combinations, sigma_clip_combinations)):
            cut_limits_min, cut_limits_max = cut_lim
            norm_limits_min, norm_limits_max = norm_lim
            modelHa_min, modelHa_max = modha
            if high_RV_amp: 
                rv_med1 = allRV1s[unique_shared_rv]
                dlam1 = ref_wl*rv_med1/speed_of_light
                modelHa_min+=dlam1; modelHa_max+=dlam1; cut_limits_min+=dlam1; cut_limits_max+=dlam1; norm_limits_min+=dlam1; norm_limits_max+=dlam1
                
                mask_input_spectrum = ((normalised_wavelength>ref_wl+cut_limits_min)   & (normalised_wavelength<ref_wl+cut_limits_max))
                normalised_wavelength, normalised_flux, normalised_err = normalised_wavelength[mask_input_spectrum], normalised_flux[mask_input_spectrum], normalised_err[mask_input_spectrum]
            
            
            if starType1=="DA" or starType1=="DC" or starType1=="DB": model_wl1, model_spectrum_star1 = return_model_spectrum_DA(wl_all1_N, ref_wl, cut_limits_min, cut_limits_max, Grav1_N, flux1_N, Teff1_N, T1_med, logg1_med)
            elif starType1=="DBA": model_wl1, model_spectrum_star1 = return_model_spectrum_DBA(wl_all1_N, ref_wl, cut_limits_min, cut_limits_max, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, T1_med, logg1_med, HoverHe1_med)
            
            if starType2=="DA" or starType2=="DC" or starType2=="DB": model_wl2, model_spectrum_star2 = return_model_spectrum_DA(wl_all2_N, ref_wl, cut_limits_min, cut_limits_max, Grav2_N, flux2_N, Teff2_N, T2_med, logg2_med)
            elif starType2=="DBA": model_wl2, model_spectrum_star2 = return_model_spectrum_DBA(wl_all2_N, ref_wl, cut_limits_min, cut_limits_max, Grav2_N, flux2_N, Teff2_N, HoverHe2_N, T2_med, logg2_med, HoverHe2_med)
                
            
            
            
            rv_med1 = allRV1s[unique_shared_rv];  rv_med2 =allRV2s[unique_shared_rv]
            rv_med1_min = allRV1sminus[unique_shared_rv];  rv_med2_min =allRV2sminus[unique_shared_rv]
            rv_med1_max = allRV1smax[unique_shared_rv];  rv_med2_max =allRV2smax[unique_shared_rv]
            
            
            
            filename = in_fi
            if "/" in filename:    filename=filename.split("/")[-1]




            dlam1 = ref_wl*rv_med1/speed_of_light
            dlam2 = ref_wl*rv_med2/speed_of_light
            desired_range = (normalised_wavelength>ref_wl+modelHa_min)  &  (normalised_wavelength<ref_wl+modelHa_max)
            
            
            if False:


                # interpolate the model onto the observation's wavelength grid
                interparr1 = interp(normalised_wavelength, model_wl1+dlam1, model_spectrum_star1)
                interparr2 = interp(normalised_wavelength, model_wl2+dlam2, model_spectrum_star2)
                interparr = Scaling_med*interparr1 + interparr2
                
                interparr = convolve(interparr, Gaussian1DKernel(stddev=ref_wl/inp_resolution), boundary = 'extend')
            
            
            else:
                # interpolate the model onto the observation's wavelength grid
                f1 = interp(model_wl2, model_wl1+dlam1, model_spectrum_star1)
                f2 = interp(model_wl2, model_wl2+dlam2, model_spectrum_star2)
                interparr = Scaling_med*f1 + f2
                
                
                # Smear to the resolution desired
                resstd=0.5*(ref_wl/inp_resolution)/(model_wl2[10]-model_wl2[9])
                interparr = convolve(interparr, Gaussian1DKernel(stddev=resstd), boundary = 'extend')
                    
                interparr = interp(normalised_wavelength, model_wl2, interparr)
            
            
            
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
                print_exception_details(eee)
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
            
            
            if sigclipspec!=-1:
                resid = normalised_flux-interparr
                resid_in_sig = resid/np.std(np.abs(resid[desired_range]))
                clip_mask = (np.abs(resid_in_sig)>sigclipspec) & desired_range & ((normalised_wavelength > ref_wl+5) | (normalised_wavelength < ref_wl-5))
            else: clip_mask=normalised_wavelength!=normalised_wavelength
            
            off =  polyfit(normalised_wavelength[desired_range & ~clip_mask], normalised_flux[desired_range & ~clip_mask] - interparr[desired_range & ~clip_mask], w=1/normalised_err[desired_range & ~clip_mask], deg=0)[0]
            interparr+=off
            

            mask_norm = ((normalised_wavelength>ref_wl+cut_limits_min) & (normalised_wavelength<ref_wl+norm_limits_min) ) | (normalised_wavelength<ref_wl+cut_limits_max) & ((normalised_wavelength>ref_wl+norm_limits_max))
            
            
            
            if plot_fancy==False:
                offset=cn_shared*-0.3 * 1.05**cn_shared
            else:
                offset=cn_shared*0.27 * 0.91**cn_shared
                if ref_wl<4000:      offset-=(0.01+0.01*cn_shared)
                #if cn_shared==6:    offset-=0.09
                
                mask=input_files_combinations==in_fi
                mins, maxes = modelHa_combinations[mask].T
                textvalmin, textvalmax = npamin(mins), npamax(maxes)
                difference = textvalmax - textvalmin
                
                
                if np.round(ref_wl,0)==6563:   toappendtext.append([npamin(mins)+difference/25, npamax(interparr+offset)-0.15, r"H$\alpha$", npamin(mins), npamax(maxes)])
                elif np.round(ref_wl,0)==4861: toappendtext.append([npamin(mins)+difference/25,npamax(interparr+offset)-0.1,r"H$\beta$", npamin(mins), npamax(maxes)])
                elif np.round(ref_wl,0)==4340: toappendtext.append([npamin(mins)+difference/25,npamax(interparr+offset)-0.075,r"H$\gamma$", npamin(mins), npamax(maxes)])
                elif np.round(ref_wl,0)==4102: toappendtext.append([npamin(mins)+difference/25,npamax(interparr+offset)-0.05-cn_shared*0.01,r"H$\delta$", npamin(mins), npamax(maxes)])
                elif np.round(ref_wl,0)==3970: toappendtext.append([npamin(mins)+difference/25,npamax(interparr+offset)-0.045-cn_shared*0.01,r"H$\epsilon$", npamin(mins), npamax(maxes)])
                elif np.round(ref_wl,0)==3889: toappendtext.append([npamin(mins)+difference/25,npamax(interparr+offset)-0.025-cn_shared*0.01,r"H$\zeta$", npamin(mins), npamax(maxes)])
                elif np.round(ref_wl,0)==3835: toappendtext.append([npamin(mins)+difference/25,npamax(interparr+offset)+0.02-cn_shared*0.01,r"H$\eta$", npamin(mins), npamax(mins)])
                else: name="UNDEF"
            
            
            if True:#ref_wl>6000:
                try:  os.mkdir("out/forPlotting")
                except: None
                np.savetxt("out/forPlotting/fit_"+in_fi+"_"+str(ref_wl)+".dat", np.array([normalised_wavelength[desired_range] - ref_wl, interparr[desired_range]]).T)
                np.savetxt("out/forPlotting/data_for_fit_"+in_fi+"_"+str(ref_wl)+".dat", np.array([normalised_wavelength - ref_wl, normalised_flux]).T)
            
                
            if plot_fancy==False:
                ax.plot(normalised_wavelength - ref_wl,normalised_flux + offset, c='k')
                ax2.plot(speed_of_light * (normalised_wavelength - ref_wl)/ref_wl,normalised_flux + offset, c='k')
                if ref_wl==npamax(reference_wl_combinations):
                    an_x = speed_of_light * (normalised_wavelength - ref_wl)/ref_wl
                    amask = (an_x>-1750) & (an_x<1750)
                    ax3.plot(an_x[amask], normalised_flux[amask] + offset, c='k')
            else:
                mask_model = (normalised_wavelength>ref_wl+norm_limits_min) & (normalised_wavelength<ref_wl+norm_limits_max)
                ax.plot(normalised_wavelength[mask_model] - ref_wl,normalised_flux[mask_model] + offset, c='k')
            
            x, y = normalised_wavelength[mask_norm], normalised_flux[mask_norm]
            x,y =normalised_wavelength[~desired_range],  interparr[~desired_range]
            if plot_fancy==False:
                ax.plot(x[x<ref_wl] - ref_wl, y[x<ref_wl] + offset, c='b', label="_blue = Continued model fit")
                ax2.plot(speed_of_light * (x[x<ref_wl] - ref_wl)/ref_wl, y[x<ref_wl] + offset, c='b', label="_blue = Continued model fit")
                if ref_wl==npamax(reference_wl_combinations):
                    an_x = speed_of_light * (x[x<ref_wl] - ref_wl)/ref_wl
                    a_y = y[x<ref_wl] + offset
                    amask = (an_x>-1750) & (an_x<1750)
                    ax3.plot(an_x[amask], a_y[amask], c='b', label="_blue = Continued model fit")
                ax.plot(x[x>ref_wl] - ref_wl, y[x>ref_wl] + offset, c='b')
                ax2.plot(speed_of_light * (x[x>ref_wl] - ref_wl)/ref_wl, y[x>ref_wl] + offset, c='b')
                if ref_wl==npamax(reference_wl_combinations):
                    an_x = speed_of_light * (x[x>ref_wl] - ref_wl)/ref_wl
                    amask = (an_x>-1750) & (an_x<1750)
                    a_y = y[x>ref_wl] + offset
                    ax3.plot(an_x[amask], a_y[amask], c='b')
            

                ax.plot(x[x<ref_wl] - ref_wl, y[x<ref_wl] + offset, c='r', label="_red = For normalising")
                ax2.plot(speed_of_light * (x[x<ref_wl] - ref_wl)/ref_wl, y[x<ref_wl] + offset, c='r', label="_red = For normalising")
                if ref_wl==npamax(reference_wl_combinations):
                    an_x = speed_of_light * (x[x<ref_wl] - ref_wl)/ref_wl
                    amask = (an_x>-1750) & (an_x<1750)
                    a_y = y[x<ref_wl] + offset
                    ax2.plot(an_x[amask], a_y[amask], c='r', label="_red = For normalising")
                ax.plot(x[x>ref_wl] - ref_wl, y[x>ref_wl] + offset, c='r')
                ax2.plot(speed_of_light * (x[x>ref_wl] - ref_wl)/ref_wl, y[x>ref_wl] + offset, c='r')
                if ref_wl==npamax(reference_wl_combinations):
                    an_x = speed_of_light * (x[x>ref_wl] - ref_wl)/ref_wl
                    amask = (an_x>-1750) & (an_x<1750)
                    a_y = y[x>ref_wl] + offset
                    ax3.plot(an_x[amask], a_y[amask], c='r')
                ax.plot(normalised_wavelength[desired_range] - ref_wl, interparr[desired_range] + offset, c='orange')
                
                
                
                ax2.plot(speed_of_light * (normalised_wavelength[desired_range] - ref_wl)/ref_wl, interparr[desired_range] + offset, c='orange')
                if ref_wl==npamax(reference_wl_combinations):
                    an_x = speed_of_light * (normalised_wavelength[desired_range] - ref_wl)/ref_wl
                    amask = (an_x>-1750) & (an_x<1750)
                    a_y = interparr[desired_range] + offset
                    ax3.plot(an_x[amask], a_y[amask], c='orange');  ax3.text(-400,1, str(ref_wl))
                
                
            else:
                ax.plot(normalised_wavelength[desired_range] - ref_wl, interparr[desired_range] + offset, c='r',linewidth=1.25)
                
                # record minimum of entire plot
                if npamin(interparr[desired_range] + offset) < minstored:
                    minstored=npamin(interparr[desired_range] + offset)
                if npamin(normalised_flux[mask_model] + offset) < minstored:
                    minstored=npamin(normalised_flux[mask_model] + offset)
                if npamax(interparr[desired_range] + offset) > maxstored:
                    maxstored=npamax(interparr[desired_range] + offset)
                if npamax(normalised_flux[mask_model] + offset) > maxstored:
                    maxstored=npamax(normalised_flux[mask_model] + offset)
                
                
            
            
            
            if sigclipspec!=-1:
                resid = normalised_flux-interparr
                resid_in_sig = resid/np.std(np.abs(resid[desired_range]))
                clip_mask = (np.abs(resid_in_sig)>sigclipspec) & desired_range & ((normalised_wavelength > ref_wl+5) | (normalised_wavelength < ref_wl-5))
                if plot_fancy==False:
                    ax.scatter(normalised_wavelength[clip_mask] - ref_wl, normalised_flux[clip_mask] + offset, c='r', s=4)
                    ax2.scatter(speed_of_light * (normalised_wavelength[clip_mask] - ref_wl)/ref_wl, normalised_flux[clip_mask] + offset, c='r', s=4)
                    if ref_wl==npamax(reference_wl_combinations):
                        xx=speed_of_light * (normalised_wavelength[clip_mask] - ref_wl)/ref_wl
                        amask = (xx>-1750) & (xx<1750)
                        yy=normalised_flux[clip_mask] + offset
                        ax3.scatter(xx[amask], yy[amask], c='r', s=8)
            
            
        
        if starType1=="DBA":    plt.scatter(0,1,alpha=0, label="H/He1 = " + str(HoverHe1_med))
        if starType2=="DBA":    plt.scatter(0,1,alpha=0, label="H/He2 = " + str(HoverHe2_med))
        
        if plot_fancy==False:
            ax.scatter(0,1,alpha=0, label="V1 = " + str(np.round(rv_med1,3)))
            ax.scatter(0,1,alpha=0, label="V2 = " + str(np.round(rv_med2,3)))
            ax.scatter(0,1,alpha=0, label="wl = " + str(ref_wl))
            ax.scatter(0,1,alpha=0, label="T1 = " + str(np.round(T1_med,0)))
            ax.scatter(0,1,alpha=0, label="Logg1 = " + str(np.round(logg1_med,3)))
            ax.scatter(0,1,alpha=0, label="T2 = " + str(np.round(T2_med,0)))
            ax.scatter(0,1,alpha=0, label="Logg2 = " + str(np.round(logg2_med,3)))
            ax.set_title(filename)
            
            ax2.scatter(0,1,alpha=0, label="V1 = " + str(np.round(rv_med1,3)))
            ax2.scatter(0,1,alpha=0, label="V2 = " + str(np.round(rv_med2,3)))
            ax2.scatter(0,1,alpha=0, label="wl = " + str(ref_wl))
            ax2.scatter(0,1,alpha=0, label="T1 = " + str(np.round(T1_med,0)))
            ax2.scatter(0,1,alpha=0, label="Logg1 = " + str(np.round(logg1_med,3)))
            ax2.scatter(0,1,alpha=0, label="T2 = " + str(np.round(T2_med,0)))
            ax2.scatter(0,1,alpha=0, label="Logg2 = " + str(np.round(logg2_med,3)))
        
            ax3.set_xlim(-1750,1750)
            ax2.set_xlabel("RV [kms-1]");   ax3.set_xlabel("RV [kms-1]")
            
            ax.legend(loc='lower right', handlelength=0, fontsize=10, frameon=False)
            plt.savefig("out/"+filename.split(".dat")[0]+"_all_shared_rvs"+".png", dpi=300)
            
            
        else:
            #ax.scatter(0,1,alpha=0, label="wl = " + str(ref_wl))
            #ax.scatter(0,1,alpha=0, label="T$_1$ = " + str(int(np.round(T1_med,-1))) +"\,$\pm$\,"+ str(int(np.round(max(T1_plus,-1*T1_minus),-1))) + "\,K")
            #ax.scatter(0,1,alpha=0, label="$\log$(g$_1$) = " + str(np.round(logg1_med,2)) +"\,$\pm$\,"+ str(np.round(max(logg1_plus,-1*logg1_minus),2)))
            #ax.scatter(0,1,alpha=0, label="V$_1$ = " + str(np.round(rv_med1,1)) +"\,$\pm$\,"+ str(np.round(max(rv_med1_max,-1*rv_med1_min),1)) + "\,km\,s$^{-1}$")
            #ax.scatter(0,1,alpha=0, label="T$_2$ = " + str(int(np.round(T2_med,-1))) +"\,$\pm$\,"+ str(int(np.round(max(T2_plus,-1*T2_minus),-1))) + "\,K")
            #ax.scatter(0,1,alpha=0, label="$\log$(g$_2$) = " + str(np.round(logg2_med,2)) +"\,$\pm$\,"+ str(np.round(max(logg2_plus,-1*logg2_minus),2)))
            #ax.scatter(0,1,alpha=0, label="V$_2$ = " + str(np.round(rv_med2,1)) +"\,$\pm$\,"+ str(np.round(max(rv_med2_max,-1*rv_med2_min),1)) + "\,km\,s$^{-1}$")
            
            
            ax.scatter(0,1,alpha=0, label="T$_1$ = " + str(int(np.round(T1_med,-1))) + "\,K")
            ax.scatter(0,1,alpha=0, label="$\log$(g$_1$) = " + str(np.round(logg1_med,2)))
            if rv_med1<0:
                ax.scatter(0,1,alpha=0, label="V$_1$ = $-$" + str(np.abs(np.round(rv_med1,1)))+"\,km\,s$^{-1}$")
            else:
                ax.scatter(0,1,alpha=0, label="V$_1$ = +" + str(np.round(rv_med1,1))+"\,km\,s$^{-1}$")
            ax.scatter(0,1,alpha=0, label="T$_2$ = " + str(int(np.round(T2_med,-1))) + "\,K")
            ax.scatter(0,1,alpha=0, label="$\log$(g$_2$) = " + str(np.round(logg2_med,2)))
            if rv_med2<0:
                ax.scatter(0,1,alpha=0, label="V$_2$ = $-$" + str(np.abs(np.round(rv_med2,1))) + "\,km\,s$^{-1}$")
            else:
                ax.scatter(0,1,alpha=0, label="V$_2$ = +" + str(np.round(rv_med2,1)) + "\,km\,s$^{-1}$")
            
        
            ax.set_xlabel(r"$\lambda$ [\AA]");   ax.set_ylabel("Relative flux")
            
            
            x,y,lab,xmin,xmax=np.asarray(toappendtext).T
            x,xmin,xmax=x.astype(float), xmin.astype(float), xmax.astype(float)
            ax.set_xlim(npamin(xmin), npamax(xmax))
            
            
            ax.set_ylim(minstored-0.25,maxstored+0.05)
            h, l = ax.get_legend_handles_labels()
            print(l)
            
            
            if float(l[0].split("=")[-1].split("\\")[0])<10000 and float(l[3].split("=")[-1].split("\\")[0])>10000:
                l[0] = l[0].split("=")[0] + "=\,\,\,  " + l[0].split("=")[1]
            elif float(l[3].split("=")[-1].split("\\")[0])<10000 and float(l[0].split("=")[-1].split("\\")[0])>10000:
                l[3] = l[3].split("=")[0] + "=\,\,\,  " + l[3].split("=")[1]
            
            if l[1][-2]==".":
                l[1]=l[1]+"0"
            if l[4][-2]==".":
                l[4]=l[4]+"0"


            textvalmin, textvalmax = npamin(xmin), npamax(xmax)
            start_for_plot = textvalmin + 0.15*(textvalmax - textvalmin)
            

            ax.text(start_for_plot, minstored-0.085, l[0] + "      " + l[1] + " dex   " + l[2], fontsize=11)
            ax.text(start_for_plot, minstored-0.185, l[3] + "      " + l[4] + " dex   " + l[5], fontsize=11)
            
            
            
            minBalmerText = npamin(x)
            for x,y,lab,xmin,xmax in toappendtext:
                ax.text(minBalmerText,y,lab, fontsize=11)
            
            plt.savefig("out/"+filename.split(".dat")[0]+"_all_shared_rvs_fancy"+".pdf", dpi=300)
        if plot_fit[0]!=False and unique_shared_rv in plot_fit and not starType1=="quadLorentz":   raise ValueError(unique_shared_rv,plot_fit); plt.show()
        plt.clf();   plt.close()





