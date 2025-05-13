import os, sys
install_path = os.environ['WD_BASS_INSTALL_DIR']
sys.path.append(install_path)
import numpy as np
from scipy import interpolate
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from numpy import load
from numpy import polyfit

def save_tables_output_DA():
	""" save cooling tracks into .dat files for DA He and DA CO core"""
	
	dir_str='/home/james/Music/tables_for_MTR/'
	directory = os.fsencode("/home/james/Music/tables_for_MTR")
	all_temps=[];  all_logg=[];  all_mass=[]

	for files in os.listdir(directory):
		filename = os.fsdecode(files)
		if filename.startswith("Table_Mass") and not ("0.2" in filename or "0.3" in filename):
			try: fi = np.loadtxt(dir_str+str(filename), unpack=True, dtype=np.str_, skiprows=2, max_rows=61)
			except: fi = np.loadtxt(dir_str+str(filename), unpack=True, dtype=np.str_, skiprows=2, max_rows=43)

			aTeff=fi[0].astype(float);         alogg=fi[1].astype(float)
			
			for i in range(len(aTeff)):
			    all_temps.append(aTeff[i]);    all_logg.append(alogg[i]);     all_mass.append(float(filename[-3:])) # mass here is found by the file name, could just as easily be a column in the table
			    
	
	#for tt, mm in zip(all_temps, all_mass):
	#	print(tt, mm)
	#raise ValueError
	all_tempsCO=np.asarray(all_temps);  all_loggCO=np.asarray(all_logg);  all_massCO=np.asarray(all_mass)

	all_gCO=10**all_loggCO
	all_radiusCO=np.sqrt(all_massCO / (all_gCO/(9.319541*2942.206218)))
	
	


	# ELM below
	all_temps=[]; all_logg=[]; all_mass=[]
	G=6.67430E-11;    one_solM = 1.98847E30;    one_solR = 6.957E8

	d = np.loadtxt('/home/james/Music/tables_for_MTR/elm.dat', unpack=True, dtype=np.str_, delimiter="|")
	TeffELM = 10**d[0].astype(float); loggELM = d[1].astype(float); MassELM = d[2].astype(float); MassELM_err = d[3].astype(float); AgeELM = d[4].astype(float); AgeELM_err = d[5].astype(float)

	for i in range(len(TeffELM)):
		if MassELM[i] <= 0.45 and MassELM[i] > 0  and MassELM_err[i] > 0 and AgeELM[i] >= 0 and AgeELM_err[i] > 0: # and TeffELM[i]<=30000
			th_g = 10**loggELM[i] * 0.01
			th_radius=np.sqrt(G*MassELM[i]*one_solM/th_g)/one_solR
			if th_radius < 0.1:
				all_temps.append(TeffELM[i]);     all_logg.append(loggELM[i]);     all_mass.append(MassELM[i])


	all_tempsHe=np.asarray(all_temps); all_loggHe=np.asarray(all_logg); all_massHe=np.asarray(all_mass)


	all_gHe=10**all_loggHe
	# [cm/s2]
	all_gHe*=0.01
	all_radiusHe=np.sqrt(G*all_massHe*one_solM/all_gHe)/one_solR
	
	np.savetxt("../saved_MTR/table_valuesCO.dat", np.array([all_tempsCO, all_loggCO, all_massCO, all_radiusCO]).T)
	np.save("../saved_MTR/table_valuesCO.npy", np.array([all_tempsCO, all_loggCO, all_massCO, all_radiusCO]))
	
	np.savetxt("../saved_MTR/table_valuesHe.dat", np.array([all_tempsHe, all_loggHe, all_massHe, all_radiusHe]).T)
	np.save("../saved_MTR/table_valuesHe.npy", np.array([all_tempsHe, all_loggHe, all_massHe, all_radiusHe]))

def save_Althaus_2013_full():  #  I depend on Istrate for lower masses, so only need the higher mass ones here
	all_T, all_logg, all_R, all_M = np.array([]), np.array([]), np.array([]), np.array([])
	os.chdir("/home/james/Downloads/ELMtracks")
	for f in os.listdir(os.getcwd()):
		if f.endswith(".trk"):# and ("03207" in f or "03630" in f or "04352" in f) : # I don't care too much about ELMs here. That was improved in https://ui.adsabs.harvard.edu/abs/2018A%26A...614A..49C/abstract  and could be added in the future to WD-BASS. Going to sample the masses at a slightly lower resolution for M<0.2 as this is the speed bottleneck otherwise
			print(f)
			LOG_L, LOG_TEFF, T_c, Ro_c, Hc, Hec, percent_Con_s, percent_Con_c, Log_age_Myr, Masa, M_dot, modelo, Log_Lpp, Log_Lcno, Log_LHe, Log_LCC, int_dS_dt, Log_Lnu,  Log_MHtot, Log_HeBuf, Masa_HFC, Masa_HeFC, Log_grav, R, L_H, Sep_Fase, periodo_orb, masa_secun = np.loadtxt(f, unpack=True, comments="#", dtype=str)
			
			
			
			LOG_L, LOG_TEFF, T_c, Ro_c, Hc, Hec, percent_Con_s, percent_Con_c, Log_age_Myr, Masa, M_dot, modelo, Log_Lpp, Log_Lcno, Log_LHe, Log_LCC, int_dS_dt, Log_Lnu,  Log_MHtot, Log_HeBuf, Masa_HFC, Masa_HeFC, Log_grav, R, L_H, Sep_Fase, periodo_orb, masa_secun   =    LOG_L.astype(float), LOG_TEFF.astype(float), T_c.astype(float), Ro_c.astype(float), Hc.astype(float), Hec.astype(float), percent_Con_s.astype(float), percent_Con_c.astype(float), Log_age_Myr.astype(float), Masa.astype(float), M_dot.astype(float), modelo.astype(float), Log_Lpp.astype(float), Log_Lcno.astype(float), Log_LHe.astype(float), Log_LCC.astype(float), int_dS_dt.astype(float), Log_Lnu.astype(float),  Log_MHtot.astype(float), Log_HeBuf.astype(float), Masa_HFC.astype(float), Masa_HeFC.astype(float), Log_grav.astype(float), R.astype(float), L_H.astype(float), Sep_Fase.astype(float), periodo_orb.astype(float), masa_secun.astype(float)
			
			
			Teff = 10**LOG_TEFF
			
			
			#ignore_starting_Ages
			
			dR, dT= np.array([]), np.array([])
			for cnt, (RR, TT) in enumerate(zip(R,Teff)):
				try:    dT = np.append(dT, Teff[cnt+1]-TT);    dR = np.append(dR, R[cnt+1]-RR)
				except: None
			
			
			mask = ((dR<0) & (dT<0) & (R[:-1]>0.05) & (Teff[:-1]==np.amax(Teff)))
			#print(f, mass)
			
			turn_around = np.argwhere(mask==True)[-1][0]  + 50
			
			
			mask  =  (Teff[turn_around:] > 4300)  &  ((Teff[turn_around:]>9000) | (np.in1d(Teff[turn_around:],  Teff[::2])))    &   ((Teff[turn_around:]>6000)  |  (np.in1d(Teff[turn_around:],  Teff[::5])))   &   ((Teff[turn_around:]<22500)  |  (np.in1d(Teff[turn_around:],  Teff[::2])))
			
			
			mass_star = float("0."+str(f[1:5]))
			
			if True:# Masa[0]>0.3:
				#plt.plot(10**LOG_TEFF, Log_grav)
				if len(R[turn_around:][mask]) < 1000:
					Teff = Teff[turn_around:][mask][::2]
					Log_grav=Log_grav[turn_around:][mask][::2]
					R=R[turn_around:][mask][::2]
					
					
					plt.scatter(Teff, R)
					plt.title(str(mass_star) + "   "  + str(len(Teff)) + "  " + str(np.amin(Teff)))
					
				else:
					Teff = Teff[turn_around:][mask]
					Log_grav=Log_grav[turn_around:][mask]
					R=R[turn_around:][mask]
					plt.scatter(Teff[::4], R[::4])
					plt.title(str(mass_star) + "   "  + str(len(Teff[::3])) + "  " + str(np.amin(Teff[::3])))
				
				#plt.show()
			
			
			mass_star = np.full((len(Teff),), mass_star)
			
			all_T, all_logg, all_R, all_M  =  np.append(all_T, Teff),  np.append(all_logg, Log_grav),  np.append(all_R, R),  np.append(all_M, mass_star)
			
	
	
	np.save(install_path+ "/saved_MTR/Althaus_2013_full.npy", np.array([all_T, all_logg, all_M, all_R]))
	np.save(install_path+ "/saved_MTR/Althaus_2013_full_nomasses.npy", np.array([all_T, all_logg, all_R]))
	
	
			
			

def save_tables_output_DB():
	""" save cooling tracks for DB/DBA/DC into numpy files for H/He=10^-5"""

	dir_str='/home/james/Music/tables_for_MTR/'
	directory = os.fsencode("/home/james/Music/tables_for_MTR")
	all_temps=[];  all_logg=[];  all_mass=[]

	for files in os.listdir(directory):
		filename = os.fsdecode(files)
		if filename.startswith("Table_Mass"):
			checkpassed=False
			testf=open(dir_str+str(filename)).readlines()
			for cn, l in enumerate(testf):
				if "Pure-helium grid" in l:    checkpassed=True
				if checkpassed==True:
					if not "Teff" in l and not "grid" in l:
						line=l.split(" ")
						line = list(filter(None, line))
						all_temps.append(float(line[0]))
						all_logg.append(float(line[1]))
						all_mass.append(float(filename[-3:]))
			
			
	all_temps = np.asarray(all_temps).astype(float)
	all_logg = np.asarray(all_logg).astype(float)
	all_massDB = np.asarray(all_mass).astype(float)


	all_tempsDB=np.asarray(all_temps);  all_loggDB=np.asarray(all_logg);  all_massDB=np.asarray(all_mass)

	all_gDB=10**all_loggDB
	all_radiusDB=np.sqrt(all_massDB / (all_gDB/(9.319541*2942.206218)))

	
	#np.savetxt("../saved_grids_npy/tableBedardDB.dat", np.array([all_tempsDB, all_loggDB, all_radiusDB]).T)  # K, logg, solR
	np.save(install_path + "/saved_grids_npy/tableBedardDB", np.array([all_tempsDB, all_loggDB, all_radiusDB]))  # K, logg, solR


def get_MTR(T, M=None, R=None, logg=None, compute_logg=False, return_R=False, return_M=False, return_R_from_T_logg=False, Althaus_or_Istrate="Istrate", loaded_Istrate=[], loaded_CO=[], loaded_Althaus=[]):
	
	#print(np.amin(all_radiusHe), np.amax(all_radiusHe), np.amin(all_tempsHe), np.amax(all_tempsHe))

	if return_R==True:
		#all_tempsCO, all_loggCO, all_massCO, all_radiusCO = np.loadtxt("/home/james/python_scripts_path/dwd_fit_package/saved_MTR/table_valuesCO.dat", unpack=True)
		all_tempsCO, all_loggCO, all_massCO, all_radiusCO = load(install_path + "/saved_MTR/table_valuesCO.npy")
		#all_tempsHe, all_loggHe, all_massHe, all_radiusHe = np.loadtxt(install_path + "/saved_MTR/table_valuesHe.dat", unpack=True)
		all_tempsHe, all_loggHe, all_massHe, all_radiusHe = load(install_path + "/saved_MTR/table_valuesHe.npy")
		
		
		rad_He = float(griddata(np.array([all_tempsHe,all_massHe]).T,all_radiusHe,np.array([T, M]).T, method='linear')[0])
		rad_CO = float(griddata(np.array([all_tempsCO,all_massCO]).T,all_radiusCO,np.array([T, M]).T, method='linear')[0])
		
		G=6.67430E-11;  one_solM = 1.989E30;   one_solR = 695700000
		if compute_logg:
			g_He=np.log10(   1000*G*M*one_solM*1000/np.square(rad_He*one_solR*100) )
			g_CO=np.log10(   1000*G*M*one_solM*1000/np.square(rad_CO*one_solR*100) )
			return rad_He, rad_CO, g_He, g_CO
		else:
			return rad_He, rad_CO
	elif return_R_from_T_logg==True:
		try:
			if logg>7.71: raise ValueError  #  Althaus has a maximum logg of 7.7, Istrate 7.62
			
			if Althaus_or_Istrate=="Althaus":# or corresponding_mass_CO>0.392:
				#all_tempsHe, all_loggHe, all_massHe, all_radiusHe = np.loadtxt("/home/james/python_scripts_path/dwd_fit_package/saved_MTR/table_valuesHe.dat", unpack=True)
				if False:   all_tempsHe, all_loggHe, all_massHe, all_radiusHe = load(install_path + "/saved_MTR/table_valuesHe.npy")
				else: 
					if len(loaded_Althaus)==0:   all_tempsHe, all_loggHe, all_radiusHe  =  loaded_Althaus
					else:  all_tempsHe, all_loggHe, all_radiusHe = load(install_path + "/saved_MTR/Althaus_2013_full_nomasses.npy")
				
				
				radius_He = float(griddata(np.array([all_tempsHe,all_loggHe]).T,all_radiusHe,np.array([T, logg]).T, method='linear')[0])
			else:
				
				try:
					if logg>7.625: raise ValueError  # istrate grid max logg is 7.625  (I checked)
					if len(loaded_Istrate)==0:  all_tempsHe_Ist, all_radiusHe_Ist, all_loggHe_Ist=load(install_path + "/saved_MTR/Istrate_Z0p02_diffusion_nomasses.npy")
					else:  all_tempsHe_Ist, all_radiusHe_Ist, all_loggHe_Ist  =  loaded_Istrate
					
					if logg>=7:
						maskIst = (all_loggHe_Ist>logg-0.25) & (all_loggHe_Ist<logg+0.25)
						if T>20000: maskIst=maskIst & (all_tempsHe_Ist>T-5000)
						elif T>15000: maskIst=maskIst & (all_tempsHe_Ist>T-3500)   &  (all_tempsHe_Ist<22500)
						elif T>10000: maskIst=maskIst & (all_tempsHe_Ist>T-2000)   &  (all_tempsHe_Ist<17500)
						elif T>8500: maskIst=maskIst & (all_tempsHe_Ist>T-1500)   &  (all_tempsHe_Ist<12500)
						elif T>7000: maskIst=maskIst & (all_tempsHe_Ist>T-1000)   &  (all_tempsHe_Ist<10000)
						elif T>6000: maskIst=maskIst & (all_tempsHe_Ist>T-800)   &  (all_tempsHe_Ist<8500)
						else: maskIst=maskIst & (all_tempsHe_Ist>T-800)   &  (all_tempsHe_Ist<8000)
					else:
						if T>15000: maskIst= all_tempsHe_Ist>11000
						elif T>10000: maskIst= all_tempsHe_Ist>7000
						elif T>8000: maskIst= all_tempsHe_Ist>6000
					
					
					#plt.scatter(all_tempsHe, all_loggHe);   plt.show()
					
					#all_massHe, all_tempsHe, all_radiusHe, all_loggHe = all_massHe[mask], all_tempsHe[mask], all_radiusHe[mask], all_loggHe[mask]
					radius_He = float(griddata(np.array([all_tempsHe_Ist[maskIst],all_loggHe_Ist[maskIst]]).T,all_radiusHe_Ist[maskIst],np.array([T, logg]).T, method='linear')[0])
					if np.isnan(radius_He): raise ValueError
				except:
					#all_tempsHe, all_loggHe, all_massHe, all_radiusHe = np.loadtxt(install_path + "/saved_MTR/table_valuesHe.dat", unpack=True)
					
					if False:   all_tempsHe, all_loggHe, all_massHe, all_radiusHe = load(install_path + "/saved_MTR/table_valuesHe.npy")
					else:       
						if len(loaded_Althaus)==0:   all_tempsHe, all_loggHe, all_radiusHe = load(install_path + "/saved_MTR/Althaus_2013_full_nomasses.npy")
						else:  all_tempsHe, all_loggHe, all_radiusHe  =  loaded_Althaus
					
					mask = (all_loggHe>logg-0.35)    &  (all_loggHe<logg+0.35)
					if T>20000: mask=mask & (all_tempsHe>T-5000)
					elif T>15000: mask=mask & (all_tempsHe>T-3500)  &  (all_tempsHe<22500)
					elif T>10000: mask=mask & (all_tempsHe>T-2000)  &  (all_tempsHe<17500)
					elif T>8500: mask=mask & (all_tempsHe>T-1500)   &  (all_tempsHe<12500)
					elif T>7000: mask=mask & (all_tempsHe>T-1000)  &  (all_tempsHe<10000)
					elif T>6000: mask=mask & (all_tempsHe>T-800)  &  (all_tempsHe<8500)
					else: mask=mask & (all_tempsHe>T-800)  &  (all_tempsHe<8000)
	
	
					radius_He = float(griddata(np.array([all_tempsHe[mask],all_loggHe[mask]]).T,all_radiusHe[mask],np.array([T, logg]).T, method='linear')[0])
					
					
					
					if np.isnan(radius_He): 
						if logg<7.6:
							radius_He = float(griddata(np.array([all_tempsHe_Ist[maskIst],all_loggHe_Ist[maskIst]]).T,all_radiusHe_Ist[maskIst],np.array([T, logg]).T, method='nearest')[0])
							
						else:  raise ValueError
			return radius_He	
		except:
			#all_tempsCO, all_loggCO, all_massCO, all_radiusCO = np.loadtxt("/home/james/python_scripts_path/dwd_fit_package/saved_MTR/table_valuesCO.dat", unpack=True)
			if len(loaded_CO)==0: all_tempsCO, all_loggCO, all_massCO, all_radiusCO = load(install_path + "/saved_MTR/table_valuesCO.npy")
			else:  all_tempsCO, all_loggCO, all_massCO, all_radiusCO  =  loaded_CO
			#mask = all_tempsCO==all_tempsCO
			mask=(all_tempsCO>=4500)  &  (all_tempsCO>=T-5000)  &  (all_tempsCO<=T+5000)  &  (all_loggCO<logg+0.5)  &  (all_loggCO>logg-0.5)
			radius_CO = float(griddata(np.array([all_tempsCO[mask],all_loggCO[mask]]).T,all_radiusCO[mask],np.array([T, logg]).T, method='linear')[0])
			#print("CO")
			return radius_CO
			
		if False:
			#all_tempsCO, all_loggCO, all_massCO, all_radiusCO = np.loadtxt("/home/james/python_scripts_path/dwd_fit_package/saved_MTR/table_valuesCO.dat", unpack=True)
			all_tempsCO, all_loggCO, all_massCO, all_radiusCO = load(install_path + "/saved_MTR/table_valuesCO.npy")
			import matplotlib.pyplot as plt
			plt.clf()
			for i in np.linspace(7,9,1000):
				plt.scatter(i, float(griddata(np.array([all_tempsHe,all_loggHe]).T,all_radiusHe,np.array([T, i]).T, method='linear')[0]), c='k')
				plt.scatter(i, float(griddata(np.array([all_tempsCO,all_loggCO]).T,all_radiusCO,np.array([T, i]).T, method='linear')[0]), c='r')
			plt.show()
			plt.clf()
			for i in np.linspace(7,9,1000):
				plt.scatter(i, float(griddata(np.array([all_tempsHe,all_loggHe]).T,all_massHe,np.array([T, i]).T, method='linear')[0]), c='k')
				plt.scatter(i, float(griddata(np.array([all_tempsCO,all_loggCO]).T,all_massCO,np.array([T, i]).T, method='linear')[0]), c='r')
			plt.show()
			# Conclude from this that if logg<=7.6 (M=0.448), interpolate on ELM grid. Otherwise, CO grid
		
	elif return_M==True:
	
		if Althaus_or_Istrate=="Althaus":# or corresponding_mass_CO>0.392:
			#all_tempsHe, all_loggHe, all_massHe, all_radiusHe = np.loadtxt("/home/james/python_scripts_path/dwd_fit_package/saved_MTR/table_valuesHe.dat", unpack=True)
			
			if False:   all_tempsHe, all_loggHe, all_massHe, all_radiusHe = load(install_path + "/saved_MTR/table_valuesHe.npy")
			else:
				all_tempsHe, all_loggHe, all_massHe, all_radiusHe = load(install_path + "/saved_MTR/Althaus_2013_full.npy")
			
			radius_He = float(griddata(np.array([all_tempsHe,all_loggHe]).T,all_radiusHe,np.array([T, logg]).T, method='linear')[0])
		else:
			#if logg>7.55 and logg<7.7  and T>7100: logg=7.545
			try:
				if logg>7.625: raise ValueError  # istrate grid max logg is 7.5493763049546985  (I checked, exact number)
				
				if len(loaded_Istrate)==0:  all_massHe, all_tempsHe, all_radiusHe, all_loggHe=load(install_path + "/saved_MTR/Istrate_Z0p02_diffusion_nomasses.npy")
				else:  all_tempsHe, all_radiusHe, all_loggHe  =  loaded_Istrate
				
				if logg>=7:
					maskIst = (all_loggHe>logg-0.25) & (all_loggHe<logg+0.25)
					if T>20000: maskIst=maskIst & (all_tempsHe>T-5000)
					elif T>15000: maskIst=maskIst & (all_tempsHe>T-3500)   &  (all_tempsHe<22500)
					elif T>10000: maskIst=maskIst & (all_tempsHe>T-2000)   &  (all_tempsHe<17500)
					elif T>8500: maskIst=maskIst & (all_tempsHe>T-1500)   &  (all_tempsHe<12500)
					elif T>7000: maskIst=maskIst & (all_tempsHe>T-1000)   &  (all_tempsHe<10000)
					elif T>6000: maskIst=maskIst & (all_tempsHe>T-800)   &  (all_tempsHe<8500)
					else: maskIst=maskIst & (all_tempsHe>T-800)   &  (all_tempsHe<8000)
				else:
					if T>15000: maskIst= all_tempsHe>11000
					elif T>10000: maskIst= all_tempsHe>7000
					elif T>8000: maskIst= all_tempsHe>6000

				
				radius_He = float(griddata(np.array([all_tempsHe[maskIst],all_loggHe[maskIst]]).T,all_radiusHe[maskIst],np.array([T, logg]).T, method='linear')[0])
			except:
				if False:   all_tempsHe, all_loggHe, all_massHe, all_radiusHe = load(install_path + "/saved_MTR/table_valuesHe.npy")
				else: 
					all_tempsHe, all_loggHe, all_massHe, all_radiusHe = load(install_path + "/saved_MTR/Althaus_2013_full.npy")
				radius_He = float(griddata(np.array([all_tempsHe,all_loggHe]).T,all_radiusHe,np.array([T, logg]).T, method='linear')[0])
				


		all_tempsCO, all_loggCO, all_massCO, all_radiusCO = load(install_path + "/saved_MTR/table_valuesCO.npy")
		mask=all_tempsCO>=4500
		radius_CO = float(griddata(np.array([all_tempsCO[mask],all_loggCO[mask]]).T,all_radiusCO[mask],np.array([T, logg]).T, method='linear')[0])
	
	
	
	
		G=6.67430E-11;    one_solM = 1.98847E30;    one_solR = 6.957E8

		if logg!=None:
			mass_He = float(griddata(np.array([all_tempsHe,all_loggHe]).T,all_massHe,np.array([T, logg]).T, method='linear')[0])
			mass_CO = float(griddata(np.array([all_tempsCO,all_loggCO]).T,all_massCO,np.array([T, logg]).T, method='linear')[0])
		if compute_logg:
			g_He=np.log10(   1000*G*mass_He*one_solM*1000/np.square(radius_He*one_solR*100) )
			g_CO=np.log10(   1000*G*mass_CO*one_solM*1000/np.square(radius_CO*one_solR*100) )
			return mass_He, mass_CO, g_He, g_CO
		else:
			return mass_He, mass_CO



def get_MR(T, M=None, R=None, logg=None, compute_logg=False, return_R=False, return_M=False, return_R_from_T_logg=False, Althaus_or_Istrate="Istrate", loaded_Istrate=[], loaded_CO=[], loaded_Althaus=[]):
	#  the code is not good here for He cores. Haven't needed to do it, so I'm leaving it half finished. CO only
	try:
		raise ValueError
		if logg>7.71: raise ValueError  #  Althaus has a maximum logg of 7.7, Istrate 7.62
		
		if Althaus_or_Istrate=="Althaus":# or corresponding_mass_CO>0.392:
			#all_tempsHe, all_loggHe, all_massHe, all_radiusHe = np.loadtxt("/home/james/python_scripts_path/dwd_fit_package/saved_MTR/table_valuesHe.dat", unpack=True)
			if False:   all_tempsHe, all_loggHe, all_massHe, all_radiusHe = load(install_path + "/saved_MTR/table_valuesHe.npy")
			else: 
				if len(loaded_Althaus)==0:   all_tempsHe, all_loggHe, all_radiusHe  =  loaded_Althaus
				else:  all_tempsHe, all_loggHe, all_radiusHe = load(install_path + "/saved_MTR/Althaus_2013_full_nomasses.npy")
			
			
			masked=(all_tempsHe>20000)
			all_tempsHe, all_loggHe, all_massHe, all_radiusHe = all_tempsHe[masked], all_loggHe[masked], all_massHe[masked], all_radiusHe[masked]
			
			sorted_loggs = np.argsort(all_loggHe)
			
			m4, m3, m2, m, c = polyfit(all_loggHe[sorted_loggs], all_radiusHe[sorted_loggs], deg=4)
			
			radius_He = logg**4 * m4  +  logg**3 * m3  +  logg**2 * m2  +  logg * m  +  c
			
			#plt.plot(all_loggHe[sorted_loggs], all_loggHe[sorted_loggs]**4 * m4  +  all_loggHe[sorted_loggs]**3 * m3  +  all_loggHe[sorted_loggs]**2 * m2  +  all_loggHe[sorted_loggs] * m  +  c)
			#plt.scatter(all_loggHe, all_radiusHe)
			#plt.show()
			
			
			
		else:
			
			try:
				if logg>7.625: raise ValueError  # istrate grid max logg is 7.625  (I checked)
				if len(loaded_Istrate)==0:  all_tempsHe_Ist, all_radiusHe_Ist, all_loggHe_Ist=load(install_path + "/saved_MTR/Istrate_Z0p02_diffusion_nomasses.npy")
				else:  all_tempsHe_Ist, all_radiusHe_Ist, all_loggHe_Ist  =  loaded_Istrate
				
				
				masked=(all_tempsHe_Ist>10000) & (all_loggHe_Ist>7.05)
				all_tempsHe_Ist, all_loggHe_Ist, all_radiusHe_Ist = all_tempsHe_Ist[masked], all_loggHe_Ist[masked], all_radiusHe_Ist[masked]
				
				sorted_loggs = np.argsort(all_loggHe_Ist)
				
				m4, m3, m2, m, c = polyfit(all_loggHe_Ist[sorted_loggs], all_radiusHe_Ist[sorted_loggs], deg=4)
				
				radius_He = logg**4 * m4  +  logg**3 * m3  +  logg**2 * m2  +  logg * m  +  c
				
				plt.scatter(all_loggHe_Ist[sorted_loggs], all_radiusHe_Ist[sorted_loggs])
				plt.plot(all_loggHe_Ist[sorted_loggs], all_loggHe_Ist[sorted_loggs]**4 * m4  +  all_loggHe_Ist[sorted_loggs]**3 * m3  +  all_loggHe_Ist[sorted_loggs]**2 * m2  +  all_loggHe_Ist[sorted_loggs] * m  +  c)
				plt.show()
				
				
				if np.isnan(radius_He): raise ValueError
			except:
				#all_tempsHe, all_loggHe, all_massHe, all_radiusHe = np.loadtxt(install_path + "/saved_MTR/table_valuesHe.dat", unpack=True)
				
				if False:   all_tempsHe, all_loggHe, all_massHe, all_radiusHe = load(install_path + "/saved_MTR/table_valuesHe.npy")
				else:       
					if len(loaded_Althaus)==0:   all_tempsHe, all_loggHe, all_radiusHe = load(install_path + "/saved_MTR/Althaus_2013_full_nomasses.npy")
					else:  all_tempsHe, all_loggHe, all_radiusHe  =  loaded_Althaus
				

				masked=(all_tempsHe>25000)
				all_tempsHe, all_loggHe, all_massHe, all_radiusHe = all_tempsHe[masked], all_loggHe[masked], all_massHe[masked], all_radiusHe[masked]
				
				sorted_loggs = np.argsort(all_loggHe)
				
				m4, m3, m2, m, c = polyfit(all_loggHe[sorted_loggs], all_radiusHe[sorted_loggs], deg=4)
				
				radius_He = logg**4 * m4  +  logg**3 * m3  +  logg**2 * m2  +  logg * m  +  c
				
				
		return radius_He	
	except:
		if logg<7.6: raise ValueError("I have not coded this properly. Do not use logg < 7.6 as prob He core. MR relationship is only correct for CO core")
		if len(loaded_CO)==0: all_tempsCO, all_loggCO, all_massCO, all_radiusCO = load(install_path + "/saved_MTR/table_valuesCO.npy")
		else:  all_tempsCO, all_loggCO, all_massCO, all_radiusCO  =  loaded_CO
		
		
		masked=(all_loggCO>7.55) & (all_tempsCO>26000)
		all_tempsCO, all_loggCO, all_massCO, all_radiusCO = all_tempsCO[masked], all_loggCO[masked], all_massCO[masked], all_radiusCO[masked]
		
		sorted_loggs = np.argsort(all_loggCO)
		radius_CO = np.interp(logg, all_loggCO[sorted_loggs], all_radiusCO[sorted_loggs])
		
		m4, m3, m2, m, c = polyfit(all_loggCO[sorted_loggs], all_radiusCO[sorted_loggs], deg=4)
		
		radius_CO = logg**4 * m4  +  logg**3 * m3  +  logg**2 * m2  +  logg * m  +  c
		
		#plt.plot(all_loggCO[sorted_loggs], all_loggCO[sorted_loggs]**4 * m4  +  all_loggCO[sorted_loggs]**3 * m3  +  all_loggCO[sorted_loggs]**2 * m2  +  all_loggCO[sorted_loggs] * m  +  c)
		#plt.scatter(all_loggCO, all_radiusCO)
		#plt.show()
		
		return radius_CO
			
	
		
def get_MTR_DB(T, M=None, R=None, logg=None, compute_logg=True, return_R=False, return_M=False, return_R_from_T_logg=False, track="Bedard"):
	if track=="Bedard":
		try: all_tempDB, all_loggDB, all_radiusDB = load(install_path + "/saved_grids_npy/tableBedardDB.npy")
		except: 
			save_tables_output_DB()
			all_tempDB, all_loggDB, all_radiusDB = load(install_path + "/saved_grids_npy/tableBedardDB.npy")
		
		mask=(all_tempDB<T+4000) & (all_tempDB>T-4000)
		
		rad = float(griddata(np.array([all_tempDB[mask],all_loggDB[mask]]).T,all_radiusDB[mask],np.array([T, logg]).T, method='linear')[0])
		return rad
		
		
	elif track=="Camisassa":
		try: all_tempDB, all_loggDB, all_radiusDB = load(install_path + "/saved_grids_npy/DB_Camisassa.npy")
		except:
			all_massesDB, all_tempDB, all_loggDB, all_radiusDB = [], [], [], []
			dir_str = install_path + "/DBmassradius/DB/Z002"
			for fi in os.listdir(dir_str):
				if ".trk" in fi:
					logL, logTeff, logTc, logRho_c, logAbundanceHCentre, logAbundanceHeCentre, OuterMassFraction, InnerMassConv, logAge, Mass, LogMdot, modelnum, LogLpp, LogLcno, LogLHe, LogLC, LogL_g, LogLnu, LogMHtot, LogMHeBuffer, MassHfreeCore, MassHefreeCore, Logg, RoverRsun, _, _ = np.loadtxt(dir_str+"/"+fi, comments="#", unpack=True)
					
					mask=(10**logTeff<35000) & (10**logTeff>4000)
					
					#plt.scatter(10**logTeff[mask], RoverRsun[mask])
					#plt.title(str(fi))
					#plt.show()
					
					for mm, tt, lglg, rr in zip(Mass[mask], 10**logTeff[mask], Logg[mask], RoverRsun[mask]):
						all_tempDB.append(tt);   all_loggDB.append(lglg);   all_radiusDB.append(rr)
						#all_massesDB.append(mm)

			all_tempDB, all_loggDB, all_radiusDB = np.asarray(all_tempDB), np.asarray(all_loggDB), np.asarray(all_radiusDB)
		
			np.save(install_path + "/saved_grids_npy/DB_Camisassa.npy", np.array([all_tempDB, all_loggDB, all_radiusDB]))
		
			#all_massesDB=np.asarray(all_massesDB)
			for gg in np.unique(np.round(all_loggDB,1)):
				mask=np.round(all_loggDB,1)==gg
				print(gg, np.unique(all_tempDB[mask]))
		
		mask=(all_tempDB<T+4000) & (all_tempDB>T-4000)
		
		rad = float(griddata(np.array([all_tempDB[mask],all_loggDB[mask]]).T,all_radiusDB[mask],np.array([T, logg]).T, method='linear')[0])
		return rad


def get_age(T, M):
	""" get total age of CO WD """
	
	print("Calculating age... Output will be age for He core models followed by age for CO core models")
	
	dir_str='/home/james/Music/tables_for_MTR/'
	directory = os.fsencode("/home/james/Music/tables_for_MTR")
	all_temps=[];  all_logg=[];  all_mass=[];  all_ages=[]

	for files in os.listdir(directory):
		filename = os.fsdecode(files)
		if filename.startswith("Table_Mass"):
			fi = np.loadtxt(dir_str+str(filename), unpack=True, dtype=np.str_, skiprows=2, max_rows=43) 

			aTeff=fi[0].astype(float);         alogg=fi[1].astype(float);    a_age=fi[-1].astype(float)
			
			for i in range(len(aTeff)):
			    all_temps.append(aTeff[i]);    all_logg.append(alogg[i]);     all_mass.append(float(filename[-3:]));  all_ages.append(float(a_age[i]))  # mass here is found by the file name, could just as easily be a column in the table


	all_tempsCO=np.asarray(all_temps);  all_loggCO=np.asarray(all_logg);  all_massCO=np.asarray(all_mass);    all_ages=np.asarray(all_ages)

	all_gCO=10**all_loggCO
	all_radiusCO=np.sqrt(all_massCO / (all_gCO/(9.319541*2942.206218)))
	
	age_CO = float(griddata(np.array([all_tempsCO,all_massCO]).T,all_ages,np.array([T, M]).T, method='linear')[0])
	
	
	
	
	
	if False:
		
		### get ages for low mass WD
		
		all_T, all_logg, all_R, all_M = np.array([]), np.array([]), np.array([]), np.array([])
		all_age_Myr = np.array([])
		os.chdir("/home/james/Downloads/ELMtracks")
		for f in os.listdir(os.getcwd()):
			if f.endswith(".trk"):# and ("03207" in f or "03630" in f or "04352" in f) : # I don't care too much about ELMs here. That was improved in https://ui.adsabs.harvard.edu/abs/2018A%26A...614A..49C/abstract  and could be added in the future to WD-BASS. Going to sample the masses at a slightly lower resolution for M<0.2 as this is the speed bottleneck otherwise
				LOG_L, LOG_TEFF, T_c, Ro_c, Hc, Hec, percent_Con_s, percent_Con_c, Log_age_Myr, Masa, M_dot, modelo, Log_Lpp, Log_Lcno, Log_LHe, Log_LCC, int_dS_dt, Log_Lnu,  Log_MHtot, Log_HeBuf, Masa_HFC, Masa_HeFC, Log_grav, R, L_H, Sep_Fase, periodo_orb, masa_secun = np.loadtxt(f, unpack=True, comments="#", dtype=str)
				
				
				
				LOG_L, LOG_TEFF, T_c, Ro_c, Hc, Hec, percent_Con_s, percent_Con_c, Log_age_Myr, Masa, M_dot, modelo, Log_Lpp, Log_Lcno, Log_LHe, Log_LCC, int_dS_dt, Log_Lnu,  Log_MHtot, Log_HeBuf, Masa_HFC, Masa_HeFC, Log_grav, R, L_H, Sep_Fase, periodo_orb, masa_secun   =    LOG_L.astype(float), LOG_TEFF.astype(float), T_c.astype(float), Ro_c.astype(float), Hc.astype(float), Hec.astype(float), percent_Con_s.astype(float), percent_Con_c.astype(float), Log_age_Myr.astype(float), Masa.astype(float), M_dot.astype(float), modelo.astype(float), Log_Lpp.astype(float), Log_Lcno.astype(float), Log_LHe.astype(float), Log_LCC.astype(float), int_dS_dt.astype(float), Log_Lnu.astype(float),  Log_MHtot.astype(float), Log_HeBuf.astype(float), Masa_HFC.astype(float), Masa_HeFC.astype(float), Log_grav.astype(float), R.astype(float), L_H.astype(float), Sep_Fase.astype(float), periodo_orb.astype(float), masa_secun.astype(float)
				
				
				Teff = 10**LOG_TEFF
				
				
				#ignore_starting_Ages
				
				dR, dT= np.array([]), np.array([])
				for cnt, (RR, TT) in enumerate(zip(R,Teff)):
					try:    dT = np.append(dT, Teff[cnt+1]-TT);    dR = np.append(dR, R[cnt+1]-RR)
					except: None
				
				
				mask = ((dR<0) & (dT<0) & (R[:-1]>0.05) & (Teff[:-1]==np.amax(Teff)))
				#print(f, mass)
				
				turn_around = np.argwhere(mask==True)[-1][0]  + 50
				
				
				mask  =  (Teff[turn_around:] > 4400)  &  ((Teff[turn_around:]>9000) | (np.in1d(Teff[turn_around:],  Teff[::2])))    &   ((Teff[turn_around:]>6000)  |  (np.in1d(Teff[turn_around:],  Teff[::5])))   &   ((Teff[turn_around:]<22500)  |  (np.in1d(Teff[turn_around:],  Teff[::2])))
				
				
				mass_star = float("0."+str(f[1:5]))
				
				
				if True:# Masa[0]>0.3:
					#plt.plot(10**LOG_TEFF, Log_grav)
					if len(R[turn_around:][mask]) < 1000:
						Teff = Teff[turn_around:][mask][::2]
						Log_grav=Log_grav[turn_around:][mask][::2]
						R=R[turn_around:][mask][::2]
						Log_age_Myr=Log_age_Myr[turn_around:][mask][::2]
						
						
						plt.scatter(Teff, R)
						plt.title(str(mass_star) + "   "  + str(len(Teff)) + "  " + str(np.amin(Teff)))
						
					else:
						Teff = Teff[turn_around:][mask]
						Log_grav=Log_grav[turn_around:][mask]
						R=R[turn_around:][mask]
						Log_age_Myr=Log_age_Myr[turn_around:][mask]
						plt.scatter(Teff[::4], R[::4])
						plt.title(str(mass_star) + "   "  + str(len(Teff[::3])) + "  " + str(np.amin(Teff[::3])))
					
					#plt.show()
				
				
				mass_star = np.full((len(Teff),), mass_star)
				
				all_T, all_logg, all_R, all_M, all_age_Myr  =  np.append(all_T, Teff),  np.append(all_logg, Log_grav),  np.append(all_R, R),  np.append(all_M, mass_star),  np.append(all_age_Myr, Log_age_Myr)
				
				
				
				
		age_He = float(griddata(np.array([all_T,all_M]).T,all_age_Myr,np.array([T, M]).T, method='linear')[0])  *  1E9







	os.chdir("/home/james/Music/Istrate_cooling/diff")
	all_z=[]
	all_mass=[]
	fnames=[]
	for filename in os.listdir(os.getcwd()):
		if filename.startswith("z_0.02"):
		    fnames.append(filename)
		    spl=filename.split("_")
		    all_z.append(spl[1])
		    all_mass.append(spl[3])

	all_z=np.asarray(all_z).astype(float)
	all_mass=np.asarray(all_mass).astype(float)


	ages_star, masses_star, temps_star, Rs_star, loggs_star=[], [], [], [], []
	for ff, mass in zip(fnames, all_mass):
		if True:  # trimmed some points that barely will change the output because another mass entry is very similar and doing this speeds up the runtime
			a=np.loadtxt(ff,unpack=True)
			AgeSt=a[1]
			MassSt=a[2] # Msun
			logTeff=a[9] # K
			Teff=10**logTeff
			logR=a[11] # Rsun
			R=10**logR
			logg=a[12] # cm/s2
			
			
			
			
			#ignore_starting_Ages
			
			dR, dT= np.array([]), np.array([])
			for cnt, (RR, TT) in enumerate(zip(R,Teff)):
				try:    dT = np.append(dT, Teff[cnt+1]-TT);    dR = np.append(dR, R[cnt+1]-RR)
				except: None
			
			
			mask = ((dR<0) & (dT<0) & (R[:-1]>0.05) & (Teff[:-1]==np.amax(Teff)))
			
			turn_around = np.argwhere(mask==True)[-1][0]
			
			
			
			
			amass=np.full((len(logg[turn_around:]),),mass)
			ateff=Teff[turn_around:]
			aR=R[turn_around:]
			alogg=logg[turn_around:]
			
			aages=AgeSt[turn_around:]
			
			
			
			masses_star = np.append(masses_star, amass)
			temps_star = np.append(temps_star, ateff)
			Rs_star = np.append(Rs_star, aR)
			loggs_star = np.append(loggs_star, alogg)
			ages_star = np.append(ages_star, aages)
	
	
	age_He = float(griddata(np.array([temps_star,masses_star]).T,ages_star,np.array([T, M]).T, method='linear')[0])

	
	return age_He/1E9, age_CO/1E9


#print(get_age(16800, 0.82)/1E9, "Gyr")
	
#save_tables_output_DB()
#save_tables_output_DA()


#save_Althaus_2013_full()


#save_tables_output_DA()

if False:
	teff1=21300;  teff1err=300
	teff2=11200;  teff2err=500
	logg1=7.86;   logg1err=0.04
	logg2=8.19;   logg2err=0.05

	medmass1=np.asarray(get_MTR(teff1, logg=logg1, return_M=True, Althaus_or_Istrate="Istrate", loaded_Istrate=[], loaded_CO=[], loaded_Althaus=[]))
	mask=np.isnan(medmass1)
	medmass1=medmass1[~mask][0]
	medmass2=np.asarray(get_MTR(teff2, logg=logg2, return_M=True, Althaus_or_Istrate="Istrate", loaded_Istrate=[], loaded_CO=[], loaded_Althaus=[]))
	mask=np.isnan(medmass2)
	medmass2=medmass2[~mask][0]


	mass1err1=np.asarray(get_MTR(teff1+teff1err, logg=logg1+logg1err, return_M=True, Althaus_or_Istrate="Istrate", loaded_Istrate=[], loaded_CO=[], loaded_Althaus=[]))
	mask=np.isnan(mass1err1)
	mass1err1=mass1err1[~mask][0]
	mass1err2=np.asarray(get_MTR(teff1+teff1err, logg=logg1-logg1err, return_M=True, Althaus_or_Istrate="Istrate", loaded_Istrate=[], loaded_CO=[], loaded_Althaus=[]))
	mask=np.isnan(mass1err2)
	mass1err2=mass1err2[~mask][0]
	mass1err3=np.asarray(get_MTR(teff1-teff1err, logg=logg1+logg1err, return_M=True, Althaus_or_Istrate="Istrate", loaded_Istrate=[], loaded_CO=[], loaded_Althaus=[]))
	mask=np.isnan(mass1err3)
	mass1err3=mass1err3[~mask][0]
	mass1err4=np.asarray(get_MTR(teff1-teff1err, logg=logg1-logg1err, return_M=True, Althaus_or_Istrate="Istrate", loaded_Istrate=[], loaded_CO=[], loaded_Althaus=[]))
	mask=np.isnan(mass1err4)
	mass1err4=mass1err4[~mask][0]


	mass1arr=np.array([mass1err1-medmass1, mass1err2-medmass1, mass1err3-medmass1, mass1err4-medmass1])
	mass1max, mass1min=np.amax(mass1arr), np.amin(mass1arr)

	mass2err1=np.asarray(get_MTR(teff2+teff2err, logg=logg2+logg2err, return_M=True, Althaus_or_Istrate="Istrate", loaded_Istrate=[], loaded_CO=[], loaded_Althaus=[]))
	mask=np.isnan(mass2err1)
	mass2err1=mass2err1[~mask][0]
	mass2err2=np.asarray(get_MTR(teff2+teff2err, logg=logg2-logg2err, return_M=True, Althaus_or_Istrate="Istrate", loaded_Istrate=[], loaded_CO=[], loaded_Althaus=[]))
	mask=np.isnan(mass2err2)
	mass2err2=mass2err2[~mask][0]
	mass2err3=np.asarray(get_MTR(teff2-teff2err, logg=logg2+logg2err, return_M=True, Althaus_or_Istrate="Istrate", loaded_Istrate=[], loaded_CO=[], loaded_Althaus=[]))
	mask=np.isnan(mass2err3)
	mass2err3=mass2err3[~mask][0]
	mass2err4=np.asarray(get_MTR(teff2-teff2err, logg=logg2-logg2err, return_M=True, Althaus_or_Istrate="Istrate", loaded_Istrate=[], loaded_CO=[], loaded_Althaus=[]))
	mask=np.isnan(mass2err4)
	mass2err4=mass2err4[~mask][0]


	mass2arr=np.array([mass2err1-medmass2, mass2err2-medmass2, mass2err3-medmass2, mass2err4-medmass2])
	mass2max, mass2min=np.amax(mass2arr), np.amin(mass2arr)


	print("M1:", medmass1, "+", np.round(mass1max,3), np.round(mass1min,3))
	print("M2:", medmass2, "+", np.round(mass2max,3), np.round(mass2min,3))




