import numpy as np
import os
from load_All_Models import load_models
from scipy.interpolate import griddata
import matplotlib
import matplotlib.pyplot as plt
import warnings
from importable_MTR_function import get_MTR
warnings.filterwarnings('ignore')
from scipy import integrate
from astropy.constants import R_sun,pc
from numba import njit, jit
from dust_extinction.parameter_averages import G23
import astropy.units as u

class Fit_phot(object):
	def air_wl_to_vacuum_wl(lambda_air):
		""" convert air to vacuum wavelengths"""
		
		s = 10**4 / lambda_air
		n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - np.square(s)) + 0.0001599740894897 / (38.92568793293 - np.square(s))
		lambda_vac = lambda_air *  n
		return lambda_vac


	

	@njit
	def return_model_spectrum_DA(wl_all1_N, ref_wl, cut_limits_min, cut_limits_max, Grav1_N, flux1_N, Teff1_N, temperature_star, logg_star):
		""" return model spectrum to fit photometry with a DA """
	
		Grav_N_N = Grav1_N;   wl_all_N_N=wl_all1_N;    flux_N_N=flux1_N;    Teff_N_N=Teff1_N
		
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

		
		model_spectrum = (ffffsss1 * (unique_Ts[1] - temperature_star) * (unique_Gs[1] - logg_star) +           ffffsss3 * (temperature_star - unique_Ts[0]) * (unique_Gs[1] - logg_star) +            ffffsss2 * (unique_Ts[1] - temperature_star) * (logg_star - unique_Gs[0]) +            ffffsss4 * (temperature_star - unique_Ts[0]) * (logg_star - unique_Gs[0])           ) / ((unique_Ts[1] - unique_Ts[0]) * (unique_Gs[1] - unique_Gs[0]))
		
		#plt.plot(wl_grid, model_spectrum,c='r', ls='--');   plt.show()
			
		#plt.plot(wl_grid, model_spectrum,c='r'); plt.plot(wl_grid, model_spectrum_starx,c='k'); plt.show()
		
		
		fine_grid=np.linspace(np.amin(wl_grid),np.amax(wl_grid),int((np.amax(wl_grid) - np.amin(wl_grid))*10)) # 0.1AA spacing
		
		
		#fine_grid=np.linspace(np.amin(wl_grid),np.amax(wl_grid),int((np.amax(wl_grid) - np.amin(wl_grid))*400)) # 0.0025AA spacing
		model_spectrum=np.interp(fine_grid, wl_grid, model_spectrum)
		
		return fine_grid, model_spectrum
		


	#@njit
	# temperature grid for DBs is irregular. Don't try to speed it up and accept the griddata way
	def return_model_spectrum_DBA(wl_all1_N, ref_wl, cut_limits_min, cut_limits_max, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, temperature_star, logg_star, HoverHe_star):
		""" return model spectrum to fit photometry with a DB/DBA/DC with varied H/He """
		
		Grav_N_N = Grav1_N;   wl_all_N_N=wl_all1_N;    flux_N_N=flux1_N;    Teff_N_N=Teff1_N;    HoverHe_N_N = HoverHe1_N
		
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
	# temperature grid for DBs is irregular. Don't try to speed it up and accept the griddata way
	def return_model_spectrum_subdwarf(wl_all1_N, min_wl, max_wl, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, temperature_star, logg_star, HoverHe_star):
		mask_logg_wl = ((wl_all1_N >= min_wl) & (wl_all1_N <= max_wl))
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



	
		
	def fit_phot_SED_double(Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N, Grav2_N, wl_all2_N, flux2_N, Teff2_N, HoverHe2_N, T1, logg1, HoverHe1, T2, logg2, HoverHe2, min_wl=2000, max_wl=10000,starType1="DA", starType2="DA", R1=None, R2=None, parallax=None, red=None,return_indiv_stars=False):
		""" combine two spectra with scaling and reden the model to fit photometry """
		
		mask_logg_wl_1 = (wl_all1_N > min_wl) & (wl_all1_N < max_wl)
		if starType1=="DA" or starType1=="DB" or starType1=="DC":
			Grav1_N, wl_all1_N, flux1_N, Teff1_N = Grav1_N[mask_logg_wl_1], wl_all1_N[mask_logg_wl_1], flux1_N[mask_logg_wl_1], Teff1_N[mask_logg_wl_1]
			model_wl1, model_spectrum_star1 = Fit_phot.return_model_spectrum_DA(wl_all1_N, 0, 0, 0, Grav1_N, flux1_N, Teff1_N, T1, logg1)
			
		elif starType1=="DBA":
			Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N = Grav1_N[mask_logg_wl_1], wl_all1_N[mask_logg_wl_1], flux1_N[mask_logg_wl_1], Teff1_N[mask_logg_wl_1], HoverHe1_N[mask_logg_wl_1]
			model_wl1, model_spectrum_star1 = Fit_phot.return_model_spectrum_DBA(wl_all1_N, 0, 0, 0, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, T1, logg1, HoverHe_star=HoverHe1)
		
			
		
		mask_logg_wl_2 = (wl_all2_N > min_wl) & (wl_all2_N < max_wl)
		if starType2=="DA" or starType2=="DB" or starType2=="DC":
			Grav2_N, wl_all2_N, flux2_N, Teff2_N = Grav2_N[mask_logg_wl_2], wl_all2_N[mask_logg_wl_2], flux2_N[mask_logg_wl_2], Teff2_N[mask_logg_wl_2]
			model_wl2, model_spectrum_star2 = Fit_phot.return_model_spectrum_DA(wl_all2_N, 0, 0, 0, Grav2_N, flux2_N, Teff2_N, T2, logg2)
			
			
		elif starType2=="DBA":
			Grav2_N, wl_all2_N, flux2_N, Teff2_N, HoverHe2_N = Grav2_N[mask_logg_wl_2], wl_all2_N[mask_logg_wl_2], flux2_N[mask_logg_wl_2], Teff2_N[mask_logg_wl_2], HoverHe2_N[mask_logg_wl_2]
			model_wl2, model_spectrum_star2 = Fit_phot.return_model_spectrum_DBA(wl_all2_N, 0, 0, 0, Grav2_N, flux2_N, Teff2_N, HoverHe2_N, T2, logg2, HoverHe_star=HoverHe2)
		
		
		
		if not np.array_equal(model_wl1,model_wl2):
			# i'm allowing slack on the models so they extend the full wl range of the transmission curve, so here I choose the biggest min and smallest max to ensure that both stars have a valid interpolation range
			min_wl1 = np.amin(model_wl1);  min_wl2 = np.amin(model_wl2);   min_wl = np.amax(np.array([min_wl1, min_wl2]))
			max_wl1 = np.amax(model_wl1);  max_wl2 = np.amax(model_wl2);   max_wl = np.amin(np.array([max_wl1, max_wl2]))
		
			
			model_wl = np.linspace(min_wl, max_wl, int((max_wl - min_wl)*10))
			model_spectrum_star1=np.interp(model_wl, model_wl1, model_spectrum_star1)
			model_spectrum_star2=np.interp(model_wl, model_wl2, model_spectrum_star2)
		else: model_wl=model_wl1
		
		
		D=(1000/parallax) * pc.value
		
		
		factor1=4*(R1*R_sun.value)**2;  factor2=4*(R2*R_sun.value)**2
		
		
		if False:
			plt.plot(model_wl, factor1 * model_spectrum_star1,c='r');   plt.plot(model_wl, factor2 * model_spectrum_star2,c='g')
			plt.title(str(R1) + "   " + 	str(R2));   plt.show();   plt.close()
		
		spec = (factor1 * model_spectrum_star1  +  factor2 * model_spectrum_star2) *  1E23 * np.pi/D**2 # was in erg/cm^2/s/Hz, putting into Jy
		
		ext = G23(Rv=3.1)
		spec *= ext.extinguish(model_wl*u.AA, Ebv=red)
		
		
		
		#plt.plot(model_wl, spec,c='r');   plt.plot(model_wl, spectrum_ext,c='g');   plt.show()
		
		if return_indiv_stars==False:  return model_wl, spec
		elif return_indiv_stars=="forSpectrum":
			return model_wl, spec, model_spectrum_star1, model_spectrum_star2
		else:
			spec1 = factor1 * model_spectrum_star1 *  1E23 * np.pi/D**2 # was in erg/cm^2/s/Hz, putting into Jy
			spec2 = factor2 * model_spectrum_star2 *  1E23 * np.pi/D**2 # was in erg/cm^2/s/Hz, putting into Jy
			spec1 *= ext.extinguish(model_wl*u.AA, Ebv=red)
			spec2 *= ext.extinguish(model_wl*u.AA, Ebv=red)
			return model_wl, spec, spec1, spec2
	
	
	
	
	def fit_phot_SED_single(Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N, T1, logg1, HoverHe1, min_wl=2000, max_wl=10000, starType1="DA", R1=None, parallax=None, red=None):
		""" scale one spectrum and reden the model to fit photometry """
		if not isinstance(R1, list):
			mask_logg_wl_1 = (wl_all1_N > min_wl) & (wl_all1_N < max_wl)
			if starType1=="DA" or starType1=="DB" or starType1=="DC":
				Grav1_N, wl_all1_N, flux1_N, Teff1_N = Grav1_N[mask_logg_wl_1], wl_all1_N[mask_logg_wl_1], flux1_N[mask_logg_wl_1], Teff1_N[mask_logg_wl_1]
				model_wl1, model_spectrum_star1 = Fit_phot.return_model_spectrum_DA(wl_all1_N, 0, 0, 0, Grav1_N, flux1_N, Teff1_N, T1, logg1)

			elif starType1=="DBA":
				Grav1_N, wl_all1_N, flux1_N, Teff1_N, HoverHe1_N = Grav1_N[mask_logg_wl_1], wl_all1_N[mask_logg_wl_1], flux1_N[mask_logg_wl_1], Teff1_N[mask_logg_wl_1], HoverHe1_N[mask_logg_wl_1]
				model_wl1, model_spectrum_star1 = Fit_phot.return_model_spectrum_DBA(wl_all1_N, 0, 0, 0, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, T1, logg1, HoverHe_star=HoverHe1)
				
			elif starType1.startswith("sd"):
				model_wl1, model_spectrum_star1 = Fit_phot.return_model_spectrum_subdwarf(wl_all1_N, min_wl, max_wl, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, T1, logg1, HoverHe1)
				
			
			D=(1000/parallax) *pc.value
			
			
			factor1=4*(R1*R_sun.value)**2  # here R1 is in solR
			
			
			if False:
				plt.plot(model_wl1, factor1 * model_spectrum_star1,c='r')
				plt.title(str(R1));   plt.show();   plt.close()
			

			spec = factor1 * model_spectrum_star1 *  1E23 * np.pi/D**2 # was in erg/cm^2/s/Hz, putting into Jy
			ext = G23(Rv=3.1);   spec *= ext.extinguish(model_wl1*u.AA, Ebv=red)
			
			return model_wl1, spec
			
		
			
		else:  #  this is a place holder in case I want to do fitting R^2/D^2 one day with no parallax information at all
			model_wl1, model_spectrum_star1 = Fit_phot.return_model_spectrum_subdwarf(wl_all1_N, min_wl, max_wl, Grav1_N, flux1_N, Teff1_N, HoverHe1_N, T1, logg1, HoverHe1)
			
			atype, mcmc_Rsquared_over_Dsquared = R1
			
			D=(1000/parallax) *pc.value
			
			R1 = np.sqrt(mcmc_Rsquared_over_Dsquared * D**2)  # here R1 is in m
			
			factor1=4*R1**2
			
			
			if False:
				plt.plot(model_wl1, factor1 * model_spectrum_star1,c='r')
				plt.title(str(R1));   plt.show();   plt.close()
			

			spec = factor1 * model_spectrum_star1 *  1E23 * np.pi/D**2 # was in erg/cm^2/s/Hz, putting into Jy
			ext = G23(Rv=3.1);   spec *= ext.extinguish(model_wl1*u.AA, Ebv=red)
			
			return model_wl1, spec
	
	
	
	
	def get_url_CDS(RADec, rad):
		""" lookup photometric data from CDS """
		
		RADec_split = RADec.split(":")
		RADec_new_split=RADec_split[2].split(" ")
		url="http://vizier.u-strasbg.fr/viz-bin/sed?-c="
		url+=RADec_split[0]  +  "%20"  +  RADec_split[1]  +  "%20"  +  RADec_new_split[0]
		if "+" in RADec:
		    url+="%20%2B"
		    url+=RADec_new_split[1][0:]
		elif "-" in RADec:
		    url+="%20"
		    url+=RADec_new_split[1]
		url+="%20"  +  RADec_split[3]  +  "%20"  +  RADec_split[4]
		url+="&-c.rs="  +  str(rad)
		return url



		
	def get_data(url, searchOnline):
		""" get CDS url to lookup and download photometry """
		
		#raise ValueError(url)
		import urllib.request
		print(url);   print(url);   print(url);   print(url)
		if searchOnline:
			if not "photSED.vot" in os.listdir(os.getcwd()+"/out"):    urllib.request.urlretrieve(url, "out/photSED.vot")
		from astropy.io.votable import parse
		votab = parse("out/photSED.vot")
		
		
		for table in votab.iter_tables():
			data = table.array
			#print(data)
			catalogue=(data["_tabname"])
			ra=(data["_RAJ2000"])
			dec=(data["_DEJ2000"])
			sedfreq=data["sed_freq"]
			sedflux=data["sed_flux"]
			sedfluxe=data["sed_eflux"]
			sedfilter=data["sed_filter"]
			
		
		sed_wl = 299792458/sedfreq  # c/f to get wl
		
		#plt.scatter(sed_wl, sedflux, c='r');  plt.show();  plt.clf()
		
		return sedfilter, sed_wl, sedflux, sedfluxe, catalogue
		
	def process_photometry_in_each_pb(model_wl, model_flux, filters, sed_wl, sedflux, sedfluxe, filter_dict, plot_solution=False, theminww_plot=1000, themaxww_plot=10000, specStar1=None, specStar2=None, single_or_double="double", return_points_for_phot_model=False):
		""" integrate the model spectrum in each passband and compute chisq compared with the observed data """
		
		#### integrate the model over the transmission filter. If I find photometry from a space based satellite, convert air to vacuum wavelengths.		
		chisq_indiv=np.array([])
		list_wl_bpass, list_flux_bpass, list_flux_sed, list_fluxe = [], [], [], []
		for cnt, filt in enumerate(filters):
			try:
				filter_wl, filter_transmission = filter_dict[filt][0]
				
				
				maskfilt_gr_0=(filter_transmission>=0)
				filter_wl=filter_wl[maskfilt_gr_0];    filter_transmission=filter_transmission[maskfilt_gr_0]
				
				mask_filt = (model_wl>=np.amin(filter_wl))  &  (model_wl<=np.amax(filter_wl))
				input_wl = model_wl[mask_filt];    input_flux = model_flux[mask_filt]
				
				if "gaia" in filt.lower() or "galex" in filt.lower() or "wise" in filt.lower():
					input_wl = Fit_phot.air_wl_to_vacuum_wl(input_wl)

				filter_transmission_on_model_grid = np.interp(input_wl, filter_wl, filter_transmission)
				
				#plt.plot(input_wl, filter_transmission_on_model_grid);   plt.plot(input_wl, input_flux)
				#plt.show();   plt.clf()
				
				if filter_dict[filt][2] == "ENERGY":
					flux_in_the_bandpass = integrate.simpson(filter_transmission_on_model_grid * input_flux, x=input_wl, dx=0.1)   /   integrate.simpson(filter_transmission_on_model_grid, x=input_wl, dx=0.1)
				elif filter_dict[filt][2] == "PHOTON":
					flux_in_the_bandpass = integrate.simpson(filter_transmission_on_model_grid * input_flux * input_wl, x=input_wl, dx=0.1)   /   integrate.simpson(filter_transmission_on_model_grid * input_wl, x=input_wl, dx=0.1)
				else: raise ValueError
				
				
				list_wl_bpass.append(sed_wl[cnt]);   list_flux_bpass.append(flux_in_the_bandpass);    list_flux_sed.append(sedflux[cnt]);    list_fluxe.append(sedfluxe[cnt])
				
				
			except Exception as e: None

		
		list_wl_bpass, list_flux_bpass, list_flux_sed, list_fluxe = np.asarray(list_wl_bpass), np.asarray(list_flux_bpass), np.asarray(list_flux_sed), np.asarray(list_fluxe)

		if False:  #### if parallax not given. make this an option in the future simply to fit the shape of the SED. Could be good for sources with no/crappy parallax
			def func_norm(x,norm):
				return x/norm
			guessed_multiplier  =  np.amax(model_flux) / np.amax(sedflux)
			from scipy.optimize import curve_fit
			norm, pcov = curve_fit(func_norm, list_flux_bpass, list_flux_sed, sigma=list_fluxe, p0=[guessed_multiplier])
			
			
			if plot_solution==True:
				plt.clf()
				fig, (ax, ax2) = plt.subplots(2, gridspec_kw={'height_ratios': [3, 1]})
				ax.plot(model_wl, model_flux,c='g', label="model fit")
				if single_or_double=="double":
					ax.plot(model_wl, specStar1,c='grey', ls='--', label="star1 fit")
					ax.plot(model_wl, specStar2,c='k', ls='--', label="star2 fit")
				ax.scatter(list_wl_bpass, list_flux_bpass, c='r', label='integrated flux')
				ax.errorbar(list_wl_bpass, list_flux_sed*norm, yerr=list_fluxe*norm, fmt='.b', label='cds data')
				ax2.errorbar(list_wl_bpass, list_flux_bpass/norm - list_flux_sed, yerr=list_fluxe, fmt='.b')
				ax2.axhline(0,ls='--',c='grey')
				ax.set_ylabel("Flux")
				ax2.set_ylabel("Flux")
				ax.legend()#;  ax2.legend()
				ax.set_xlim(theminww_plot, themaxww_plot)
				ax2.set_xlim(theminww_plot, themaxww_plot)
				if themaxww_plot<20000: ax.set_xscale("linear"); ax2.set_xscale("linear"); ax2.set_xlabel("Wavelength [AA]")
				else: ax.set_xscale("log"); ax2.set_xscale("log"); ax2.set_xlabel("log10( Wavelength [AA] )")
				plt.savefig(os.getcwd()+"/out/photometric_fit.png")
				plt.clf();  plt.close()
			
			chisq=-0.5*np.sum(np.square((list_flux_bpass/norm - list_flux_sed)   /    list_fluxe))
		
		
		
		chisq=-0.5*np.sum(np.square((list_flux_bpass - list_flux_sed)   /    list_fluxe))
		
		if np.isnan(chisq):
			for z,zz,zzz in zip(filters,list_flux_sed, list_fluxe):
				print(z,zz, zzz)
			raise ValueError("Issue with chisq calculation for photometric SED fitting. Likely that one of the errors are nans/0. I have printed all filters with the flux and fluxerr above.")
		
		
		plot_fancy=False
		if plot_solution==True and plot_fancy==False:
			plt.clf()
			fig, (ax, ax2) = plt.subplots(2, gridspec_kw={'height_ratios': [3, 1]})
			ax.plot(model_wl, model_flux,c='k', label="model fit")
			if single_or_double=="double":
				ax.plot(model_wl, specStar1,c='g', ls='--', label="star1 fit");
				ax.plot(model_wl, specStar2,c='r', ls='--', label="star2 fit")
			ax.scatter(list_wl_bpass, list_flux_bpass, c='orange', label='integrated flux')
			ax.errorbar(list_wl_bpass, list_flux_sed, yerr=list_fluxe, fmt='.b', label='cds data')
			ax2.errorbar(list_wl_bpass, list_flux_bpass - list_flux_sed, yerr=list_fluxe, fmt='.b')
			ax2.axhline(0,ls='--',c='grey')
			ax2.set_xlabel("Wavelength [AA]")
			ax.set_ylabel("Flux");   ax2.set_ylabel("Flux")
			ax.legend(loc='lower right')#;  ax2.legend()
			ax.set_xlim(theminww_plot, themaxww_plot);   ax2.set_xlim(theminww_plot, themaxww_plot)
			ax.set_title(chisq)
			plt.savefig(os.getcwd()+"/out/photometric_fit.png")
			#plt.show()
			plt.clf();   plt.close()
		elif plot_solution==True and plot_fancy==True:
			
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
 

			
			plt.clf()
			fig, (ax, ax2) = plt.subplots(2, gridspec_kw={'height_ratios':[4.5,1]}, sharex=True)
			fig.subplots_adjust(hspace=0.08)
			#ax.plot(model_wl, model_flux/1E-15,c='k', label="model fit")
			ax.plot(model_wl, model_flux*1000,c='k', label="model fit")
			if single_or_double=="double":
				#ax.plot(model_wl, specStar1/1E-15,c='g', ls='--', label="star1 fit");   ax.plot(model_wl, specStar2/1E-15,c='r', ls='--', label="star2 fit")
				ax.plot(model_wl, specStar1*1000,c='g', ls='--', label="star1 fit");
				ax.plot(model_wl, specStar2*1000,c='r', ls='--', label="star2 fit")
			#ax.scatter(list_wl_bpass, list_flux_bpass/1E-15, marker="x", c='orange', label='integrated flux')
			#ax.errorbar(list_wl_bpass, list_flux_sed/1E-15, yerr=list_fluxe/1E-15, fmt='.k', label='cds data')
			ax.scatter(list_wl_bpass, list_flux_bpass*1000, marker="x", c='orange', label='integrated flux')
			ax.errorbar(list_wl_bpass, list_flux_sed*1000, yerr=list_fluxe*1000, fmt='.k', label='cds data')
			
			
			ax2.errorbar(list_wl_bpass, 100*(list_flux_bpass-list_flux_sed)/list_flux_sed, yerr=list_fluxe/list_flux_sed, fmt='.k', label='cds data')
			ax2.axhline(0,ls='--',c='grey')
			ax2.set_xlim(np.amin(model_wl), np.amax(model_wl))
			ax2.set_ylabel(r"\%\,F", labelpad=5)
			ax2val=np.amax(np.abs(100*(list_flux_bpass-list_flux_sed)/list_flux_sed))
			from matplotlib.ticker import MaxNLocator
			ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
			ax2.set_ylim(np.round(-1.2*ax2val,0), np.round(1.2*ax2val,))
			ax2.set_yticks([-2,0,2])
			
			
			
			#ax.set_ylabel("Flux [10$^{-15}$ erg\,cm$^{-2}$\,s${-1}$ $A^{-1}$]")
			ax.set_ylabel(r"F$_{\nu}$ [mJy]");   ax2.set_xlabel(r"Wavelength [\AA]")
			#ax.legend(loc='lower right')#;  ax2.legend()
			ax.set_xlim(np.amin(model_wl), np.amax(model_wl));   ax.set_ylim(0,1.05*np.amax(model_flux*1000))
			#ax.set_title(chisq)
			plt.savefig(os.getcwd()+"/out/photometric_fit_fancy.pdf");   #plt.show()
			plt.clf();   plt.close()
			
		
		
		if return_points_for_phot_model:  return chisq/(len(list_flux_bpass)-1), chisq, list_wl_bpass, list_flux_bpass, list_flux_sed, list_fluxe
		else:                             return chisq/(len(list_flux_bpass)-1), chisq  # reduced chisq
	
	
