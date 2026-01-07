import numpy as np
import os

class load_models(object):
	# old load DA Pier 3D
	#def load_models_DA(minwl=3000,maxwl=8000):
	#	cwd="/home/james/python_scripts_path/dwd_fit_package/models"

	#	# get wavelengths
	#	big_flux_list=[]
	#	flux, Teff, Grav=[], [], []

	#	wl_all, wl_3747, wl_2691= [], [], []

	#	case_2691, case_3747 = True, True
	#	found_3747, found_2691 = False, False
		
	#	for filename in os.listdir(cwd):
	#		Effective_counter=0
	#		if "3D" in filename:# and "800_3D" in filename:
	#			trigger_wl=True
	#			wl=[]
	#			f = open(cwd+"/"+filename).readlines()
	#			for count, line in enumerate(f):
	#				if not count == 0:
	#					if "Effective" in f[count]:
	#						Effective_counter+=1
	#						trigger_wl=False
	#						teff=line[24:35]
	#						grav=line[47:58]
	#					    
	#					    
	#					else:
	#						if trigger_wl==True:
	#							try:
	#								if float(f[0])==3747 and case_3747==True:
	#									found_3747=True
	#									wl_3747.append(float(line[1:10]))
	#									wl_3747.append(float(line[11:20]))
	#									wl_3747.append(float(line[21:30]))
	#									wl_3747.append(float(line[32:40]))
	#									wl_3747.append(float(line[42:50]))
	#									wl_3747.append(float(line[52:60]))
	#									wl_3747.append(float(line[62:70]))
	#									wl_3747.append(float(line[72:80]))
	#									wl_3747.append(float(line[82:90]))
	#									wl_3747.append(float(line[92:100]))
	#								elif float(f[0])==2691 and case_2691==True:
	#									found_2691=True
	#									wl_2691.append(float(line[:8]))
	#									wl_2691.append(float(line[8:16]))
	#									wl_2691.append(float(line[16:24]))
	#									wl_2691.append(float(line[24:32]))
	#									wl_2691.append(float(line[32:40]))
	#									wl_2691.append(float(line[40:48]))
	#									wl_2691.append(float(line[48:56]))
	#									wl_2691.append(float(line[56:64]))
	#									wl_2691.append(float(line[64:72]))
	#									wl_2691.append(float(line[72:80]))
							    
	#							except:
	#								None
					
	#			if found_3747==True:
	#				case_3747=False
	#			if found_2691==True:
	#				case_2691=False




	#	#print(len(wl_2691))
	#	#print(len(wl_3747))






	#	for filename in os.listdir(cwd):
	#		Effective_counter=0
	#		Teff_just_this_file=[]
	#		if "3D" in filename:# and "800_3D" in filename:
	#			trigger_wl=True
	#			wl=[]
	#			#print(filename)
	#			f = open(cwd+"/"+filename).readlines()
	#			for count, line in enumerate(f):
	#				if not count == 0:
	#					if "Effective" in f[count]:
	#						Effective_counter+=1
	#						trigger_wl=False
	#						teff=line[24:35]
	#						grav=line[47:58]
	#			    
	#						Teff_just_this_file.append(float(teff))
	#		
	#					else:
	#						if trigger_wl==True:
	#							try:
	#								if float(f[0])==3747:
	#									wl.append(float(line[1:10]))
	#									wl.append(float(line[11:20]))
	#									wl.append(float(line[21:30]))
	#									wl.append(float(line[32:40]))
	#									wl.append(float(line[42:50]))
	#									wl.append(float(line[52:60]))
	#									wl.append(float(line[62:70]))
	#									wl.append(float(line[72:80]))
	#									wl.append(float(line[82:90]))
	#									wl.append(float(line[92:100]))
	#								elif float(f[0])==2691:
	#									wl.append(float(line[:8]))
	#									wl.append(float(line[8:16]))
	#									wl.append(float(line[16:24]))
	#									wl.append(float(line[24:32]))
	#									wl.append(float(line[32:40]))
	#									wl.append(float(line[40:48]))
	#									wl.append(float(line[48:56]))
	#									wl.append(float(line[56:64]))
	#									wl.append(float(line[64:72]))
	#									wl.append(float(line[72:80]))
	#			            
	#							except Exception as e:
	#								None
	#						else:
	#							try:
	#								#flux.append(float(line[:13]))
	#								#print(line[61:72])
	#			            
	#								flux.append(float(line[:13]))
	#								Teff.append(float(teff))
	#								Grav.append(float(grav))
	#								flux.append(float(line[13:25]))
	#								Teff.append(float(teff))
	#								Grav.append(float(grav))
	#								flux.append(float(line[25:37]))
	#								Teff.append(float(teff))
	#								Grav.append(float(grav))
	#								flux.append(float(line[37:49]))
	#								Teff.append(float(teff))
	#								Grav.append(float(grav))
	#								flux.append(float(line[49:61]))
	#								Teff.append(float(teff))
	#								Grav.append(float(grav))
	#								flux.append(float(line[61:72]))
	#								Teff.append(float(teff))
	#								Grav.append(float(grav))
	#			            
	#			            
	#			            
	#			            
	#			            
	#							except Exception as e:
	#								None
	#			repeated_wls = np.tile(wl, Effective_counter)
	#			for l in range(len(repeated_wls)):
	#				wl_all.append(repeated_wls[l])
	#		
	#		
	#		
	#		
	#		
	#			
	#	wl_all=np.asarray(wl_all).astype("float")
	#	flux=np.asarray(flux).astype("float")
	#	Teff=np.asarray(Teff).astype("float")
	#	Grav=np.asarray(Grav).astype("float")
	#	mask_wls = ((wl_all<=maxwl) & (wl_all>=minwl))
	#	
	#	
	#	return wl_all[mask_wls], flux[mask_wls], Teff[mask_wls], Grav[mask_wls]


	def load_models_DA(minwl=3000,maxwl=8000, toConvolve=True):
		install_path = os.environ['WD_BASS_INSTALL_DIR']
		cwd=install_path + "/models"
		

		# get wavelengths
		big_flux_list, wl_all=[], []
		flux, Teff, Grav=[], [], []


		for filename in os.listdir(cwd):
			Effective_counter=0
			Teff_just_this_file=[]
			if "3D" in filename and not " " in filename and int(filename[0])<7:# and "800_3D" in filename:
				trigger_wl=True
				wl=[]
				#print(filename)
				f = open(cwd+"/"+filename).readlines()
				for count, line in enumerate(f):
					if not count == 0:
						if "Effective" in f[count]:
							Effective_counter+=1
							trigger_wl=False
							teff=line[24:35]
							grav=line[47:58]
				    
							Teff_just_this_file.append(float(teff))
			
						else:
							if trigger_wl==True:
								try:
									if float(f[0])==3747:
										wl.append(float(line[:10]))
										wl.append(float(line[11:20]))
										wl.append(float(line[20:30]))
										wl.append(float(line[30:40]))
										wl.append(float(line[40:50]))
										wl.append(float(line[50:60]))
										wl.append(float(line[60:70]))
										wl.append(float(line[70:80]))
										wl.append(float(line[80:90]))
										wl.append(float(line[90:100]))
									elif float(f[0])==2691:
										wl.append(float(line[:10]))
										wl.append(float(line[11:20]))
										wl.append(float(line[20:30]))
										wl.append(float(line[30:40]))
										wl.append(float(line[40:50]))
										wl.append(float(line[50:60]))
										wl.append(float(line[60:70]))
										wl.append(float(line[70:80]))
										wl.append(float(line[80:90]))
										wl.append(float(line[90:100]))
				            
								except Exception as e:
									None
							else:
								try:
									#flux.append(float(line[:13]))
									#print(line[61:72])
				            
									flux.append(float(line[:13]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									flux.append(float(line[13:25]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									flux.append(float(line[25:37]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									flux.append(float(line[37:49]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									flux.append(float(line[49:61]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									flux.append(float(line[61:72]))
									Teff.append(float(teff))
									Grav.append(float(grav))
				            
				            
				            
				            
				            
								except Exception as e:
									None
				repeated_wls = np.tile(wl, Effective_counter)
				for l in range(len(repeated_wls)):
					wl_all.append(repeated_wls[l])
			
			
			
			
				
		wl_all=np.asarray(wl_all).astype("float")
		flux=np.asarray(flux).astype("float")
		Teff=np.asarray(Teff).astype("float")
		Grav=np.asarray(Grav).astype("float")
		
		
		if True:
			new_wl, new_flux, new_teff, new_grav=[], [], [], []
			from scipy.interpolate import interp1d
			import matplotlib.pyplot as plt
			for ggg in np.unique(Grav):
				for TTT in np.unique(Teff):
					mask = (Teff==TTT) & (Grav==ggg)
					awl, aflux, ateff, agrav  =  wl_all[mask], flux[mask], Teff[mask], Grav[mask]
					if len(awl)>1:# and np.round(np.log10(ggg),0) >=8 and TTT>8000:
						minleft, maxleft, minright, maxright = 6562.8-150, 6561, 6564.7, 6562.8+150
						maskhalpha=(awl<maxright) & (awl>maxleft)
						maskhalphawings=((awl>minright) & (awl<maxright)) | ((awl<maxleft) & (awl>minleft))
						
						
						maskhalphawings1=((awl<maxleft) & (awl>minleft))
						maskhalphawings2=((awl>minright) & (awl<maxright))
						
						
						
						#plt.plot(awl[maskhalpha], aflux[maskhalpha])
						#plt.plot(awl[maskhalphawings], aflux[maskhalphawings])
						#plt.scatter(awl[maskhalpha], aflux[maskhalpha])
						#plt.scatter(awl[maskhalphawings], aflux[maskhalphawings])
						
						
						
						leftwing = interp1d(awl[maskhalphawings1], aflux[maskhalphawings1],kind="cubic")
						rightwing = interp1d(awl[maskhalphawings2], aflux[maskhalphawings2],kind="cubic")
						
						arrleft=np.linspace(np.amin(awl[maskhalphawings1]), np.amax(awl[maskhalphawings1]), 100)
						arrright=np.linspace(np.amin(awl[maskhalphawings2]), np.amax(awl[maskhalphawings2]), 100)
						#plt.plot(arrleft, leftwing(arrleft), c='r')
						#plt.plot(arrright, rightwing(arrright), c='r')
						
						
						
						
						new_wl_ax = np.concatenate((awl[awl<minleft], arrleft, awl[(awl<minright) & (awl>maxleft)], arrright, awl[awl>maxright]))
						new_flux_ax = np.concatenate((aflux[awl<minleft], leftwing(arrleft), aflux[(awl<minright) & (awl>maxleft)], rightwing(arrright), aflux[awl>maxright]))
						newtemp=np.full((len(new_wl_ax),),TTT)
						newG=np.full((len(new_wl_ax),),ggg)
						
						if False:
							#plt.title(str(TTT) + "  " + str(ggg))
							#plt.show()
							plt.clf()
							plt.plot(new_wl_ax, new_flux_ax)
							plt.scatter(new_wl_ax, new_flux_ax)
							plt.show()
							print(TTT,ggg)
						
						new_wl = np.append(new_wl, new_wl_ax)
						new_flux = np.append(new_flux, new_flux_ax)
						new_teff = np.append(new_teff, newtemp)
						new_grav = np.append(new_grav, newG)
		
		
		
		
			wl_all, flux, Teff, Grav = new_wl, new_flux, new_teff, new_grav
		
		
		
		
		
		
		
		
		if toConvolve:
			from scipy.interpolate import interp1d
		
			# https://arxiv.org/pdf/1302.2013.pdf  approximate convection broadening applied to models based on Fig. 11, squaring the y axis.
			xs = np.array([1000, 5900, 5997, 7012, 8032, 9035, 9520, 10018, 10530, 11004, 11531, 12022, 12505, 12999, 14000, 14100, 50000]);   ys = np.square(np.array([0, 0, 0.3,0.6,1.35,1.5,2,2,2.5,3.25,3.8,4.6,4.8,5.2,5.2,0,0])) / 3E5  #  approximate the kms-1 smear with a wavelength of 6000AA. underpredicted at Ha, overpredicted at the others. 5000, 5900, 14000 and 50000 are put in to make sure it never applies to non-3D atmospheres. Here, it's v^2/c^2 not times the wavelength since I change insert this down below
			convection_broadening_interpolator = interp1d(xs, ys, fill_value="extrapolate")
			
			#import matplotlib.pyplot as plt
			#plt.plot(xs, ys);  plt.xlim(5000,15000);  plt.show()
		
		
		
		
			#import matplotlib.pyplot as plt
			from astropy.convolution import Gaussian1DKernel,convolve
			
			
			
			
			
			
			
			unique_Grav=np.unique(Grav)
			unique_Grav=unique_Grav[(np.round(np.log10(unique_Grav),0)>=7)]
			for zteff in np.unique(Teff):
				if zteff<=14000 and zteff>=5000:
					for zgrav in unique_Grav:
						try:
							grab_model_mask = ((Teff==zteff) & (Grav==zgrav))
							
							if len(wl_all[grab_model_mask])>1:
								
								for refwl in [6562.79, 4861.35, 4340.472, 4101.734, 3970.075, 3889.064, 9546, 10050, 10940, 12820, 18750]:
									grab_mask2 = (wl_all>refwl-25) & (wl_all<refwl+25)
									fine_wl = np.linspace(np.amin(wl_all[grab_model_mask & grab_mask2]), np.amax(wl_all[grab_model_mask & grab_mask2]), 12000)
									fine_flux = np.interp(fine_wl, wl_all[grab_model_mask & grab_mask2], flux[grab_model_mask & grab_mask2])
								
									smeared_flux1 = convolve(fine_flux, Gaussian1DKernel(stddev=refwl*convection_broadening_interpolator(zteff)/(fine_wl[10]-fine_wl[9])), boundary = 'extend')
									
									mask3=(wl_all>refwl-20) & (wl_all<refwl+20)

									
									smeared_flux = np.interp(wl_all[grab_model_mask & grab_mask2 & mask3], fine_wl, smeared_flux1)
									
									
									
									flux[grab_model_mask & grab_mask2 & mask3] = smeared_flux
								
						except Exception as e: print(e); raise ValueError
		
		mask_wls = ((wl_all<=maxwl) & (wl_all>=minwl))

		print(Teff[mask_wls], Grav[mask_wls])
		
		wl_all, flux, Teff, Grav = wl_all[mask_wls], flux[mask_wls], Teff[mask_wls], Grav[mask_wls]
		
		
		
		
		
		
		return wl_all, flux, Teff, Grav

	def zzload_models_DBA(minwl=3000,maxwl=8000):
		cwd="/home/james/Desktop/FitDA_properly/3d_dba/3d_dba/3D_DBA"

		# get wavelengths
		big_flux_list=[]
		wl_all, flux=[], []
		Teff, Grav, H_over_He=[], [], []


		list_found=0
		for filename in os.listdir(cwd):
			Effective_counter=0
			Teff_just_this_file=[]
			if "DBA" in filename and not ".py" in filename:# and "800_3D" in filename:
				trigger_wl=True
				wl=[]
				#print(filename)
				f = open(cwd+"/"+filename).readlines()
				for count, line in enumerate(f):
					if not count == 0:
						if "Effective" in f[count]:
							#print(list_found/6,count)
							Effective_counter+=1
							trigger_wl=False
							teff=line[28:35]
							grav=line[48:57]
							h_over_he = line[64:-1]
							#print(teff, "|", grav, "|",  h_over_he)
				    
							Teff_just_this_file.append(float(teff))
							list_found=0
			
						else:
							if trigger_wl==True:
								try:
									wl.append(float(line[1:10]))
									wl.append(float(line[11:20]))
									wl.append(float(line[21:30]))
									wl.append(float(line[32:40]))
									wl.append(float(line[42:50]))
									wl.append(float(line[52:60]))
									wl.append(float(line[62:70]))
									wl.append(float(line[72:80]))
									wl.append(float(line[82:90]))
									wl.append(float(line[92:100]))
									#print(line[1:10], line[11:20], line[21:30], line[32:40], line[42:50], line[52:60], line[62:70], line[72:80], line[82:90], line[92:100])
				            
								except Exception as e:
									None
							else:
								try:
									#flux.append(float(line[:13]))
									#print(line[61:72])
									try:    flux.append(float(line[1:12]))
									except: flux.append(float(line[1:12].split("-")[0]+"E-"+line[1:12].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									H_over_He.append(float(h_over_he))
									try:    flux.append(float(line[13:24]))
									except: flux.append(float(line[13:24].split("-")[0]+"E-"+line[13:24].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									H_over_He.append(float(h_over_he))
									2#print("01,12|", line[1:12], line[13:24])
									
									try:    flux.append(float(line[25:36]))
									except: flux.append(float(line[25:36].split("-")[0]+"E-"+line[25:36].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									H_over_He.append(float(h_over_he))
									try:    flux.append(float(line[37:48]))
									except: flux.append(float(line[37:48].split("-")[0]+"E-"+line[37:48].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									H_over_He.append(float(h_over_he))
									#print("23,26|", line[25:36], line[37:48])
									
									try:     flux.append(float(line[49:60]))
									except:  flux.append(float(line[49:60].split("-")[0]+"E-"+line[49:60].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									H_over_He.append(float(h_over_he))
									try:     flux.append(float(line[61:75]))
									except:  flux.append(float(line[61:75].split("-")[0]+"E-"+line[61:71].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									H_over_He.append(float(h_over_he))
									#print("48,60|", line[49:60], line[61:72])
									list_found+=6
				            
								except Exception as e:
									print(e)
				
				repeated_wls = np.tile(wl, Effective_counter)
				for l in range(len(repeated_wls)):
					wl_all.append(repeated_wls[l])
		
		
		
		
		wl_all=np.asarray(wl_all).astype("float")
		flux=np.asarray(flux).astype("float")
		Teff=np.asarray(Teff).astype("float")
		Grav=np.asarray(Grav).astype("float")
		H_over_He=np.asarray(H_over_He).astype("float")
		
		
		mask_wl_all = ((wl_all > minwl) & (wl_all <= maxwl))
		
		
		return wl_all[mask_wl_all], flux[mask_wl_all], Teff[mask_wl_all], Grav[mask_wl_all], H_over_He[mask_wl_all]


	def load_model_1D_DA(minwl=3000,maxwl=8000): # use this when outside the bounds of 3Ds. Even grid
		# get wavelengths
		trigger_wl=True
		flux, Teff, Grav, wl_all=[], [], [], []
		workingdir="/home/james/Desktop/FitDA_properly/1d_DA_NLTE/"
		for filename in os.listdir(workingdir):
			Effective_counter=0
			if not filename.endswith(".py") and not filename.endswith(".tar"):
				trigger_wl=True
				wl=[]
				f = open(workingdir+filename).readlines()
				for count, line in enumerate(f):
					if not count == 0:
						a = (f[count].split(" "))
						if "Effective" in f[count]:
							Effective_counter+=1
							trigger_wl=False
							thing=True
							while thing==True:
								thing=False
								for j in range(len(a)):
									if a[j] == "":
										a.pop(j)
										thing=True
										break
							for k in range(len(a)):
								if a[k] == "temperature":
									teff=float(a[k+2])
								if a[k] == "gravity":
									grav = float(a[k+2])
				    
				    
				    
						else:
							for i in range(len(a)):
								if not a[i] == "":
									if "\n" in a[i] :
										a[i]=a[i][:-1]
										#print(a[i], a[0])
									if trigger_wl==True:
										wl.append(float(a[i]))
									else:
										if "-1" in a[i] and not "E-1" in a[i]:
											a[i] = (a[i][:-4]+"E"+a[i][-4:])
										flux.append(float(a[i]))
										Teff.append(teff)
										Grav.append(grav)
				                
				repeated_wls = np.tile(wl, Effective_counter)
				for l in range(len(repeated_wls)):
					wl_all.append(repeated_wls[l])
				
			
			
		# grav is always the same, teff changes
		mask_wl = (np.asarray(wl_all) > minwl) & (np.asarray(wl_all) <= maxwl)
			
		ret_Teff=np.asarray(Teff)[mask_wl]
		ret_grav=np.asarray(Grav)[mask_wl]
		ret_wl_all=np.asarray(wl_all)[mask_wl]
		ret_flux=np.asarray(flux)[mask_wl]

		return ret_wl_all, ret_flux, ret_Teff, np.round(np.log10(ret_grav),1)
	
	
	
	def load_models_DA_AntoineNLTE(minwl=3000,maxwl=8000):
		install_path = os.environ['WD_BASS_INSTALL_DIR']
		cwd=install_path+"/Antoine_NLTE"

		# get wavelengths
		big_flux_list=[]
		wl_all, flux, Teff, Grav = [], [], [], []


		for filename in os.listdir(cwd):
			Effective_counter=0
			if "NLTE" in filename:# and "800_3D" in filename:
				trigger_wl=True
				wl=[]
				f = open(cwd+"/"+filename).readlines()
				for count, line in enumerate(f):
					if not count == 0:
						if "Effective" in f[count]:
							Effective_counter+=1
							trigger_wl=False
							teff=line[28:35]
							grav=line[47:57].split(",")[0]
							print(float(teff), float(grav))
							
						    
						    
						else:
							if trigger_wl==True:
								try:
								
									wl.append(float(line[3:11]))
									wl.append(float(line[11:21]))
									wl.append(float(line[23:32]))
									wl.append(float(line[32:41]))
									wl.append(float(line[43:53]))
									wl.append(float(line[53:63]))
									wl.append(float(line[63:73]))
									wl.append(float(line[73:83]))
									wl.append(float(line[83:92]))
									wl.append(float(line[93:102]))
								except:  None
				break
					




		for filename in os.listdir(cwd):
			Effective_counter=0
			Teff_just_this_file=[]
			if "NLTE" in filename:# and "800_3D" in filename:
				trigger_wl=True
				#print(filename)
				f = open(cwd+"/"+filename).readlines()
				for count, line in enumerate(f):
					if not count == 0:
						if "Effective" in f[count]:
							Effective_counter+=1
							trigger_wl=False
							teff=line[28:35]
							grav=line[47:57].split(",")[0]
				    
							Teff_just_this_file.append(float(teff))
			
						else:
							if trigger_wl==False:
								#try:
								if "-" in line[:13] and not "E" in line[:13]: split=line[:13].split("-");  flux.append(float(split[0]+"E-"+split[1]))
								else: flux.append(float(line[:13]))
								
								Teff.append(float(teff));  Grav.append(float(grav))
								
								if "-" in line[13:25] and not "E" in line[13:25]: split=line[13:25].split("-");  flux.append(float(split[0]+"E-"+split[1]))
								else: flux.append(float(line[13:25]))
								
								Teff.append(float(teff));  Grav.append(float(grav))
								
								if "-" in line[25:37] and not "E" in line[25:37]: split=line[25:37].split("-");  flux.append(float(split[0]+"E-"+split[1]))
								else: flux.append(float(line[25:37]))
								
								Teff.append(float(teff));   Grav.append(float(grav))
								
								if "-" in line[37:49] and not "E" in line[37:49]: split=line[37:49].split("-");  flux.append(float(split[0]+"E-"+split[1]))
								else: flux.append(float(line[37:49]))
								
								Teff.append(float(teff));   Grav.append(float(grav))
								
								if "-" in line[49:61] and not "E" in line[49:61]: split=line[49:61].split("-");  flux.append(float(split[0]+"E-"+split[1]))
								else: flux.append(float(line[49:61]))
								
								Teff.append(float(teff));   Grav.append(float(grav))
								
								if "-" in line[61:72] and not "E" in line[61:72]: split=line[61:72].split("-");  flux.append(float(split[0]+"E-"+split[1]))
								else: flux.append(float(line[61:72]))
								
								Teff.append(float(teff));   Grav.append(float(grav))
					    
					    
								#except Exception as e:
								#	None
				repeated_wls = np.tile(wl, Effective_counter)
				for l in range(len(repeated_wls)):
					wl_all.append(repeated_wls[l])
			
				
				print(len(Teff), len(Grav), len(wl_all))
			
			
			
				
		wl_all=np.asarray(wl_all).astype("float")
		flux=np.asarray(flux).astype("float")
		Teff=np.asarray(Teff).astype("float")
		Grav=np.asarray(Grav).astype("float")
		
		mask_wls = ((wl_all>=minwl) & (wl_all<=maxwl))
		
		return wl_all[mask_wls], flux[mask_wls], Teff[mask_wls], Grav[mask_wls]
	
	
	
	def load_models_DA_AntoineLTE(minwl=3000,maxwl=8000):
		install_path = os.environ['WD_BASS_INSTALL_DIR']
		cwd=install_path + "/Antoine_NLTE"

		# get wavelengths
		big_flux_list=[]
		wl_all, flux, Teff, Grav=[], [], [], []


		for filename in os.listdir(cwd):
			Effective_counter=0
			if "LTE" in filename and not "NLTE" in filename:# and "800_3D" in filename:
				trigger_wl=True
				wl=[]
				f = open(cwd+"/"+filename).readlines()
				for count, line in enumerate(f):
					if not count == 0:
						if "Effective" in f[count]:
							Effective_counter+=1
							trigger_wl=False
							teff=line[28:35]
							grav=line[47:57].split(",")[0]
							print(float(teff), float(grav))
							
						    
						    
						else:
							if trigger_wl==True:
								try:
								
									wl.append(float(line[3:11]))
									wl.append(float(line[11:21]))
									wl.append(float(line[23:32]))
									wl.append(float(line[32:41]))
									wl.append(float(line[43:53]))
									wl.append(float(line[53:63]))
									wl.append(float(line[63:73]))
									wl.append(float(line[73:83]))
									wl.append(float(line[83:92]))
									wl.append(float(line[93:102]))
								except:  None
				break
					




		for filename in os.listdir(cwd):
			Effective_counter=0
			Teff_just_this_file=[]
			if "LTE" in filename and not "NLTE" in filename:# and "800_3D" in filename:
				trigger_wl=True
				#print(filename)
				f = open(cwd+"/"+filename).readlines()
				for count, line in enumerate(f):
					if not count == 0:
						if "Effective" in f[count]:
							Effective_counter+=1
							trigger_wl=False
							teff=line[28:35]
							grav=line[47:57].split(",")[0]
				    
							Teff_just_this_file.append(float(teff))
			
						else:
							if trigger_wl==False:
								#try:
								if "-" in line[:13] and not "E" in line[:13]: split=line[:13].split("-");  flux.append(float(split[0]+"E-"+split[1]))
								else: flux.append(float(line[:13]))
								
								Teff.append(float(teff));  Grav.append(float(grav))
								
								if "-" in line[13:25] and not "E" in line[13:25]: split=line[13:25].split("-");  flux.append(float(split[0]+"E-"+split[1]))
								else: flux.append(float(line[13:25]))
								
								Teff.append(float(teff));  Grav.append(float(grav))
								
								if "-" in line[25:37] and not "E" in line[25:37]: split=line[25:37].split("-");  flux.append(float(split[0]+"E-"+split[1]))
								else: flux.append(float(line[25:37]))
								
								Teff.append(float(teff));   Grav.append(float(grav))
								
								if "-" in line[37:49] and not "E" in line[37:49]: split=line[37:49].split("-");  flux.append(float(split[0]+"E-"+split[1]))
								else: flux.append(float(line[37:49]))
								
								Teff.append(float(teff));   Grav.append(float(grav))
								
								if "-" in line[49:61] and not "E" in line[49:61]: split=line[49:61].split("-");  flux.append(float(split[0]+"E-"+split[1]))
								else: flux.append(float(line[49:61]))
								
								Teff.append(float(teff));   Grav.append(float(grav))
								
								if "-" in line[61:72] and not "E" in line[61:72]: split=line[61:72].split("-");  flux.append(float(split[0]+"E-"+split[1]))
								else: flux.append(float(line[61:72]))
								
								Teff.append(float(teff));   Grav.append(float(grav))
					    
					    
								#except Exception as e:
								#	None
				repeated_wls = np.tile(wl, Effective_counter)
				for l in range(len(repeated_wls)):
					wl_all.append(repeated_wls[l])
			
				
				print(len(Teff), len(Grav), len(wl_all))
			
			
			
				
		wl_all=np.asarray(wl_all).astype("float")
		flux=np.asarray(flux).astype("float")
		Teff=np.asarray(Teff).astype("float")
		Grav=np.asarray(Grav).astype("float")
		
		mask_wls = ((wl_all>=minwl) & (wl_all<=maxwl))
		
		return wl_all[mask_wls], flux[mask_wls], Teff[mask_wls], Grav[mask_wls]
	
	
	def load_models_DA_3D_NLTE(minwl=3000,maxwl=8000, toConvolve=True):
		install_path = os.environ['WD_BASS_INSTALL_DIR']
		try: 
			wl_all, flux, Teff, Grav = np.load(install_path + "/saved_grids_npy/DA_3D_NLTE.npy")
			maskwl = (wl_all<maxwl) & (wl_all>minwl)
			return wl_all[maskwl], flux[maskwl], Teff[maskwl], Grav[maskwl]
		except:
			saveit=True
			minwl, maxwl = 1, 9E99
			
			
			wl_all_3Dp, flux_3Dp, Teff_3Dp, Grav_3Dp = load_models.load_models_DA_pier_3D_new(minwl=minwl,maxwl=maxwl,toConvolve=False)
			
			
			amask = (Teff_3Dp>=5000) & (Teff_3Dp%500 == 0)
			wl_all_3Dp, flux_3Dp, Teff_3Dp, Grav_3Dp = wl_all_3Dp[amask], flux_3Dp[amask], Teff_3Dp[amask], Grav_3Dp[amask]
			
			
			
			wl_all_NLTE, flux_NLTE, Teff_NLTE, Grav_NLTE = load_models.load_models_DA_AntoineNLTE(minwl=minwl,maxwl=maxwl)
			wl_all_LTE, flux_LTE, Teff_LTE, Grav_LTE = load_models.load_models_DA_AntoineLTE(minwl=minwl,maxwl=maxwl)


			Grav_3Dp = np.round(np.log10(Grav_3Dp),1)
			Grav_NLTE = np.round(np.log10(Grav_NLTE),2)
			Grav_LTE = np.round(np.log10(Grav_LTE),2)
			
			
			
			
			wl_all_1Dp, flux_1Dp, Teff_1Dp, Grav_1Dp = load_models.load_model_1D_DA(minwl=minwl,maxwl=maxwl)
			
			
			amask=  (Grav_1Dp<7) & ((Teff_1Dp==14500) | (Teff_1Dp==15000) | (Teff_1Dp==15500) | (Teff_1Dp==16000) | (Teff_1Dp==16500) | (Teff_1Dp==17000) | (Teff_1Dp==20000) | (Teff_1Dp==25000) | (Teff_1Dp==30000) | (Teff_1Dp==35000) | (Teff_1Dp==40000) )
       
			wl_all_1Dp, flux_1Dp, Teff_1Dp, Grav_1Dp = wl_all_1Dp[amask], flux_1Dp[amask], Teff_1Dp[amask], Grav_1Dp[amask]
			
			
			
			
			newwl=np.unique(wl_all_3Dp)
			for ateff in np.unique(Teff_1Dp): # unique 3Ds and 1Ds are the same
				mask=(Teff_1Dp==ateff) & (Grav_1Dp==6.5)
				
				newflux=np.interp(newwl, wl_all_1Dp[mask], flux_1Dp[mask])
				
				
				wl_all_3Dp = np.concatenate((wl_all_3Dp, newwl))
				flux_3Dp = np.concatenate((flux_3Dp, newflux))
				
				Teff_3Dp, Grav_3Dp = np.concatenate((Teff_3Dp, np.full((len(newwl),), ateff))), np.concatenate((Grav_3Dp, np.full((len(newwl),), 6.5)))
				
					
		
			
			for ateff in np.unique(Teff_NLTE): # unique 3Ds and 1Ds are the same
				for alogg in [6.5, 7, 7.5, 8, 8.5, 9]:
					
					mask_3D = (Teff_3Dp==ateff) & (Grav_3Dp==alogg)
					
					
					mask_NLTE = (Teff_NLTE==ateff) & (Grav_NLTE==alogg)
					mask_LTE = (Teff_LTE==ateff) & (Grav_LTE==alogg)
					
					NLTE_over_LTE = flux_NLTE[mask_NLTE]/flux_LTE[mask_LTE]
					
					
					wl_new, flux_new = wl_all_3Dp[mask_3D], flux_3Dp[mask_3D]
					
					
					new_NLTE_over_LTE = np.interp(wl_new, wl_all_NLTE[mask_NLTE], NLTE_over_LTE)
					
					flux_3Dp[mask_3D] *= new_NLTE_over_LTE
					
			
			
			

			if toConvolve:
				from scipy.interpolate import interp1d
			
				# https://arxiv.org/pdf/1302.2013.pdf  approximate convection broadening applied to models based on Fig. 11, squaring the y axis.
				xs = np.array([1000, 5900, 5997, 7012, 8032, 9035, 9520, 10018, 10530, 11004, 11531, 12022, 12505, 12999, 14000, 14100, 50000]);   ys = np.square(np.array([0, 0, 0.3,0.6,1.35,1.5,2,2,2.5,3.25,3.8,4.6,4.8,5.2,5.2,0,0])) / 3E5  #  approximate the kms-1 smear with a wavelength of 6000AA. underpredicted at Ha, overpredicted at the others. 5000, 5900, 14000 and 50000 are put in to make sure it never applies to non-3D atmospheres. Here, it's v^2/c^2 not times the wavelength since I change insert this down below
				convection_broadening_interpolator = interp1d(xs, ys, fill_value="extrapolate")
				
				#import matplotlib.pyplot as plt
				#plt.plot(xs, ys);  plt.xlim(5000,15000);  plt.show()
			
			
			
			
				#import matplotlib.pyplot as plt
				from astropy.convolution import Gaussian1DKernel,convolve
				
				
				
				
				
				unique_Grav=np.unique(Grav_3Dp)
				unique_Grav=unique_Grav[unique_Grav>=7]
				for zteff in np.unique(Teff_3Dp):
					if zteff<=14000 and zteff>=5000:
						for zgrav in unique_Grav:
							try:
								grab_model_mask = ((Teff_3Dp==zteff) & (Grav_3Dp==zgrav))
								
								if len(wl_all_3Dp[grab_model_mask])>1:
									#if zteff>10500:
									#	plt.plot(wl_all_3Dp[grab_model_mask], flux_3Dp[grab_model_mask],c='r')
									
									
									
									for refwl in [6562.79, 4861.35, 4340.472, 4101.734, 3970.075, 3889.064, 9546, 10050, 10940, 12820, 18750]:
										if refwl<11000:
											grab_mask2 = (wl_all_3Dp>refwl-25) & (wl_all_3Dp<refwl+25)
										else:
											grab_mask2 = (wl_all_3Dp>refwl-50) & (wl_all_3Dp<refwl+50)
										fine_wl = np.linspace(np.amin(wl_all_3Dp[grab_model_mask & grab_mask2]), np.amax(wl_all_3Dp[grab_model_mask & grab_mask2]), 12000)
										fine_flux = np.interp(fine_wl, wl_all_3Dp[grab_model_mask & grab_mask2], flux_3Dp[grab_model_mask & grab_mask2])
										
										
										smeared_flux1 = convolve(fine_flux, Gaussian1DKernel(stddev=refwl*convection_broadening_interpolator(zteff)/(fine_wl[10]-fine_wl[9])), boundary = 'extend')
										
										
										mask3=(wl_all_3Dp>refwl-20) & (wl_all_3Dp<refwl+20)

										smeared_flux = np.interp(wl_all_3Dp[grab_model_mask & grab_mask2 & mask3], fine_wl, smeared_flux1)
										
										
										if not np.isnan(smeared_flux[0]):
											flux_3Dp[grab_model_mask & grab_mask2 & mask3] = smeared_flux
										
									#if zteff>10500:
									#	plt.plot(wl_all_3Dp[grab_model_mask], flux_3Dp[grab_model_mask],c='g',ls='--')
									#	plt.title(str(zteff) + "   " + str(np.log10(zgrav)));    plt.show(); plt.clf()
							except Exception as e: print(e); raise ValueError(zteff, np.log10(zgrav), refwl, np.amin(flux_3Dp[grab_model_mask]), np.amax(flux_3Dp[grab_model_mask]))
			
			
			
			
			
			wl_all_1Dp, flux_1Dp, Teff_1Dp, Grav_1Dp = load_models.load_model_1D_DA(minwl=minwl,maxwl=maxwl)
			mask_under_5000K = (Teff_1Dp<5000) & (Teff_1Dp>=4000)
			
			new1Dwl, new1Dflux, new1Dteff, new1Dlogg = np.array([]), np.array([]), np.array([]), np.array([])
			wl_grid_3D = np.unique(wl_all_3Dp)
			for uniqueT_1D in np.unique(Teff_1Dp[mask_under_5000K]):
				for uniqueG_1D in np.unique(Grav_1Dp[mask_under_5000K]):
					mask = (Teff_1Dp==uniqueT_1D) & (Grav_1Dp==uniqueG_1D)
					getwls, getfluxes = wl_all_1Dp[mask], flux_1Dp[mask]
					new1Dwl = np.concatenate((new1Dwl, wl_grid_3D))
					new1Dflux = np.concatenate((new1Dflux, np.interp(wl_grid_3D, getwls, getfluxes)))
					new1Dteff = np.concatenate((new1Dteff, np.full((len(wl_grid_3D),), uniqueT_1D)))
					new1Dlogg = np.concatenate((new1Dlogg, np.full((len(wl_grid_3D),), uniqueG_1D)))
			
			
			
			
			wl_all_3Dp, flux_3Dp, Teff_3Dp, Grav_3Dp = np.concatenate((new1Dwl, wl_all_3Dp)), np.concatenate((new1Dflux, flux_3Dp)), np.concatenate((new1Dteff, Teff_3Dp)), np.concatenate((new1Dlogg, Grav_3Dp))
			

			if saveit==True:
				np.save(install_path + "/saved_grids_npy/DA_3D_NLTE.npy", np.array([wl_all_3Dp, flux_3Dp, Teff_3Dp, Grav_3Dp]))
			
			
			#raise ValueError(wl_all_3Dp.shape, flux_3Dp.shape, Teff_3Dp.shape, Grav_3Dp.shape)
			
			
			return wl_all_3Dp, flux_3Dp, Teff_3Dp, Grav_3Dp
		
	def load_models_DA_pier_3D_new(minwl=3000,maxwl=8000, toConvolve=True):
		install_path = os.environ['WD_BASS_INSTALL_DIR']
		cwd=install_path + "/Pier_new_3d"
		
		### THIS IS IN VACUUM WAVELENGTHS

		# get wavelengths
		big_flux_list, wl_all=[], []
		flux, Teff, Grav=[], [], []
		try: 
			wl_all, flux, Teff, Grav = np.load(install_path + "/saved_grids_npy/load_models_DA_pier_3D_new.npy")
			maskwl = (wl_all<maxwl) & (wl_all>minwl)
			return wl_all[maskwl], flux[maskwl], Teff[maskwl], Grav[maskwl]
		except:
			saveit=True
			minwl, maxwl = 1, 9E99
			
			if toConvolve==False: saveit=False
			
			for filename in os.listdir(cwd):
				Effective_counter=0
				Teff_just_this_file=[]
				if "3D" in filename and not ".tar" in filename and not " " in filename:# and "800_3D" in filename:
					trigger_wl=True
					wl=[]
					#print(filename)
					f = open(cwd+"/"+filename).readlines()
					for count, line in enumerate(f):
						if not count == 0:
							if "Effective" in f[count]:
								Effective_counter+=1
								trigger_wl=False
								teff=line[24:35]
								grav=line[47:58]
					    
								Teff_just_this_file.append(float(teff))
				
							else:
								if trigger_wl==True:
									try:
										if float(f[0])==6787:
											wl.append(float(line[:10]))
											wl.append(float(line[11:20]))
											wl.append(float(line[20:30]))
											wl.append(float(line[30:40]))
											wl.append(float(line[40:50]))
											wl.append(float(line[50:60]))
											wl.append(float(line[60:70]))
											wl.append(float(line[70:80]))
											wl.append(float(line[80:90]))
											wl.append(float(line[90:100]))
											
						    
									except Exception as e:
										None
								
								else:
									try:
										#flux.append(float(line[:13]))
										#print(line[61:72])
										
										if "-" in line[:13] and not "E" in line[:13]:
											thing=line[:13].split("-")
											flux.append(float(thing[0]+"E-"+thing[-1]))
										else:
											flux.append(float(line[:13]))
										Teff.append(float(teff))
										Grav.append(float(grav))
										if "-" in line[13:25] and not "E" in line[13:25]:
											thing=line[13:25].split("-")
											flux.append(float(thing[0]+"E-"+thing[-1]))
										else:
											flux.append(float(line[13:25]))
										Teff.append(float(teff))
										Grav.append(float(grav))
										if "-" in line[25:37] and not "E" in line[25:37]:
											thing=line[25:37].split("-")
											flux.append(float(thing[0]+"E-"+thing[-1]))
										else:
											flux.append(float(line[25:37]))
										Teff.append(float(teff))
										Grav.append(float(grav))
										if "-" in line[37:49] and not "E" in line[37:49]:
											thing=line[37:49].split("-")
											flux.append(float(thing[0]+"E-"+thing[-1]))
										else:
											flux.append(float(line[37:49]))
										Teff.append(float(teff))
										Grav.append(float(grav))
										if "-" in line[49:61] and not "E" in line[49:61]:
											thing=line[49:61].split("-")
											flux.append(float(thing[0]+"E-"+thing[-1]))
										else:
											flux.append(float(line[49:61]))
										Teff.append(float(teff))
										Grav.append(float(grav))
										if "-" in line[61:72] and not "E" in line[61:72]:
											thing=line[61:72].split("-")
											flux.append(float(thing[0]+"E-"+thing[-1]))
										else:
											flux.append(float(line[61:72]))
										Teff.append(float(teff))
										Grav.append(float(grav))
						    
									except Exception as e:
										None#raise ValueError(e)
					#print(filename, Effective_counter, len(flux), len(wl)*Effective_counter, len(wl))
					#print(np.amin(wl), np.amax(wl))
					repeated_wls = np.tile(wl, Effective_counter)
					for l in range(len(repeated_wls)):
						wl_all.append(repeated_wls[l])
				
				
				
				
			
					
			wl_all=np.asarray(wl_all).astype("float")
			flux=np.asarray(flux).astype("float")
			Teff=np.asarray(Teff).astype("float")
			Grav=np.asarray(Grav).astype("float")
			
			
			
			##https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
			s = 10**4 / wl_all
			n = 1 + 0.0000834254 + 0.02406147 / (130 - s**2) + 0.00015998 / (38.9 - s**2)
			wl_all = wl_all / n
			
			#raise ValueError(wl_all.shape, flux.shape, Teff.shape, Grav.shape)
			
			if False:
				new_wl, new_flux, new_teff, new_grav=[], [], [], []
				from scipy.interpolate import interp1d
				import matplotlib.pyplot as plt
				for ggg in np.unique(Grav):
					for TTT in np.unique(Teff):
						mask = (Teff==TTT) & (Grav==ggg)
						awl, aflux, ateff, agrav  =  wl_all[mask], flux[mask], Teff[mask], Grav[mask]
						if len(awl)>1:# and np.round(np.log10(ggg),0) >=8 and TTT>8000:
							minleft, maxleft, minright, maxright = 6562.8-150, 6561, 6564.7, 6562.8+150
							maskhalpha=(awl<maxright) & (awl>maxleft)
							maskhalphawings=((awl>minright) & (awl<maxright)) | ((awl<maxleft) & (awl>minleft))
							
							
							maskhalphawings1=((awl<maxleft) & (awl>minleft))
							maskhalphawings2=((awl>minright) & (awl<maxright))
							
							
							
							#plt.plot(awl[maskhalpha], aflux[maskhalpha])
							#plt.plot(awl[maskhalphawings], aflux[maskhalphawings])
							#plt.scatter(awl[maskhalpha], aflux[maskhalpha])
							#plt.scatter(awl[maskhalphawings], aflux[maskhalphawings])
							
							
							
							leftwing = interp1d(awl[maskhalphawings1], aflux[maskhalphawings1],kind="cubic")
							rightwing = interp1d(awl[maskhalphawings2], aflux[maskhalphawings2],kind="cubic")
							
							arrleft=np.linspace(np.amin(awl[maskhalphawings1]), np.amax(awl[maskhalphawings1]), 100)
							arrright=np.linspace(np.amin(awl[maskhalphawings2]), np.amax(awl[maskhalphawings2]), 100)
							#plt.plot(arrleft, leftwing(arrleft), c='r')
							#plt.plot(arrright, rightwing(arrright), c='r')
							
							
							
							
							new_wl_ax = np.concatenate((awl[awl<minleft], arrleft, awl[(awl<minright) & (awl>maxleft)], arrright, awl[awl>maxright]))
							new_flux_ax = np.concatenate((aflux[awl<minleft], leftwing(arrleft), aflux[(awl<minright) & (awl>maxleft)], rightwing(arrright), aflux[awl>maxright]))
							newtemp=np.full((len(new_wl_ax),),TTT)
							newG=np.full((len(new_wl_ax),),ggg)
							
							if False:
								#plt.title(str(TTT) + "  " + str(ggg))
								#plt.show()
								plt.clf()
								plt.plot(new_wl_ax, new_flux_ax)
								plt.scatter(new_wl_ax, new_flux_ax)
								plt.show()
								print(TTT,ggg)
							
							new_wl = np.append(new_wl, new_wl_ax)
							new_flux = np.append(new_flux, new_flux_ax)
							new_teff = np.append(new_teff, newtemp)
							new_grav = np.append(new_grav, newG)
			
			
				wl_all, flux, Teff, Grav = new_wl, new_flux, new_teff, new_grav
			
			
			
			
			
			
			
			if toConvolve:
				from scipy.interpolate import interp1d
			
				# https://arxiv.org/pdf/1302.2013.pdf  approximate convection broadening applied to models based on Fig. 11, squaring the y axis.
				xs = np.array([1000, 5900, 5997, 7012, 8032, 9035, 9520, 10018, 10530, 11004, 11531, 12022, 12505, 12999, 14000, 14100, 50000]);   ys = np.square(np.array([0, 0, 0.3,0.6,1.35,1.5,2,2,2.5,3.25,3.8,4.6,4.8,5.2,5.2,0,0])) / 3E5  #  approximate the kms-1 smear with a wavelength of 6000AA. underpredicted at Ha, overpredicted at the others. 5000, 5900, 14000 and 50000 are put in to make sure it never applies to non-3D atmospheres. Here, it's v^2/c^2 not times the wavelength since I change insert this down below
				convection_broadening_interpolator = interp1d(xs, ys, fill_value="extrapolate")
				
				#import matplotlib.pyplot as plt
				#plt.plot(xs, ys);  plt.xlim(5000,15000);  plt.show()
			
			
			
			
				#import matplotlib.pyplot as plt
				from astropy.convolution import Gaussian1DKernel,convolve
				
				
				
				
				
				
				
				unique_Grav=np.unique(Grav)
				unique_Grav=unique_Grav[(np.round(np.log10(unique_Grav),0)>=7)]
				for zteff in np.unique(Teff):
					if zteff<=14000 and zteff>=5000:
						for zgrav in unique_Grav:
							try:
								grab_model_mask = ((Teff==zteff) & (Grav==zgrav))
								
								if len(wl_all[grab_model_mask])>1:
									
									for refwl in [6562.79, 4861.35, 4340.472, 4101.734, 3970.075, 3889.064, 9546, 10050, 10940, 12820, 18750]:
										if refwl<11000:
											grab_mask2 = (wl_all>refwl-25) & (wl_all<refwl+25)
										else:
											grab_mask2 = (wl_all>refwl-50) & (wl_all<refwl+50)
										fine_wl = np.linspace(np.amin(wl_all[grab_model_mask & grab_mask2]), np.amax(wl_all[grab_model_mask & grab_mask2]), 12000)
										fine_flux = np.interp(fine_wl, wl_all[grab_model_mask & grab_mask2], flux[grab_model_mask & grab_mask2])
										
										
										smeared_flux1 = convolve(fine_flux, Gaussian1DKernel(stddev=refwl*convection_broadening_interpolator(zteff)/(fine_wl[10]-fine_wl[9])), boundary = 'extend')
										
										mask3=(wl_all>refwl-20) & (wl_all<refwl+20)
										
										smeared_flux = np.interp(wl_all[grab_model_mask & grab_mask2 & mask3], fine_wl, smeared_flux1)
										
										flux[grab_model_mask & grab_mask2 & mask3] = smeared_flux
							except Exception as e: print(e); raise ValueError(zteff, np.log10(zgrav), refwl, np.amin(flux[grab_model_mask]), np.amax(flux[grab_model_mask]))
			
			
			
			
			
			includeELM=True
			if includeELM:
				
				wl_all_logg6p5, flux_logg6p5, Teff_logg6p5, Grav_logg6p5 = load_models.load_models_DA(minwl=minwl, maxwl=maxwl)
				
				mmmm = np.round(np.log10(Grav_logg6p5),1)==6.5
				wl_all_logg6p5, flux_logg6p5, Teff_logg6p5, Grav_logg6p5 = wl_all_logg6p5[mmmm], flux_logg6p5[mmmm], Teff_logg6p5[mmmm], Grav_logg6p5[mmmm]
				
				
				mask5000 = (Teff>=5000) & (Teff!=5250)
				wl_all, flux, Teff, Grav = wl_all[mask5000], flux[mask5000], Teff[mask5000], Grav[mask5000]
				
				
				for gg in np.unique(Grav_logg6p5):
					for tt in np.unique(Teff_logg6p5):
						masktt = (Teff_logg6p5==tt) & (Grav_logg6p5==gg)
						wl_temp, flux_temp, teff_temp, grav_temp = wl_all_logg6p5[masktt], flux_logg6p5[masktt], Teff_logg6p5[masktt], Grav_logg6p5[masktt]
						
						flux_newgrid = np.interp(np.unique(wl_all), wl_temp, flux_temp)
						
						
						wl_all=np.append(wl_all, np.unique(wl_all))
						flux=np.append(flux,flux_newgrid)
						Teff=np.append(Teff, np.full((len(np.unique(wl_all)),),tt))
						Grav=np.append(Grav,np.full((len(np.unique(wl_all)),),gg))
				
			
			
			mask_wls = ((wl_all<=maxwl) & (wl_all>=minwl))

			
			wl_all, flux, Teff, Grav = wl_all[mask_wls], flux[mask_wls], Teff[mask_wls], Grav[mask_wls]
			
		
		#mask=(np.log10(Grav)==6) & (Teff==6000)
		#plt.plot(wl_all[mask], flux[mask])
		#plt.title(str(len(wl_all[mask])))
		#plt.show()
		
		
		if saveit==True:
			np.save(install_path + "/saved_grids_npy/load_models_DA_pier_3D_new.npy", np.array([wl_all, flux, Teff, Grav]))
		return wl_all, flux, Teff, Grav
				

	def load_models_DBA_1eMinus5(minwl=3000,maxwl=8000):
		# get wavelengths
		big_flux_list=[]
		wl_all, flux = [], []
		Teff, Grav = [], []
		
		install_path = os.environ['WD_BASS_INSTALL_DIR']
		cwd=install_path+"/saved_grids_npy/he-grid_v3"


		list_found=0
		for filename in os.listdir(cwd):
			Effective_counter=0
			Teff_just_this_file=[]
			if not ".py" in filename and "050" in filename:
				trigger_wl=True
				wl=[]
				f = open(cwd+"/"+filename).readlines()
				for count, line in enumerate(f):
					if not count == 0:
						if "Effective" in f[count]:
							#print(list_found/6,count)
							Effective_counter+=1
							trigger_wl=False
							teff=line[28:35]
							grav=line[48:57]
				    
							Teff_just_this_file.append(float(teff))
							list_found=0
			
						else:
							if trigger_wl==True:
								try:
									#print(line[1:10], line[11:20], line[21:30], line[32:40], line[42:50], line[52:60], line[62:70], line[72:80], line[82:90], line[92:100])
									wl.append(float(line[1:10]))
									wl.append(float(line[11:20]))
									wl.append(float(line[21:30]))
									wl.append(float(line[31:40]))
									wl.append(float(line[41:50]))
									wl.append(float(line[51:60]))
									wl.append(float(line[61:70]))
									wl.append(float(line[71:80]))
									wl.append(float(line[81:90]))
									wl.append(float(line[91:100]))
									#print(line[1:10], line[11:20], line[21:30], line[32:40], line[42:50], line[52:60], line[62:70], line[72:80], line[82:90], line[92:100])
				            
								except Exception as e:
									None
							else:
								try:
									#flux.append(float(line[:13]))
									#print(line[61:72])
									try:    flux.append(float(line[1:12]))
									except: flux.append(float(line[1:12].split("-")[0]+"E-"+line[1:12].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									try:    flux.append(float(line[13:24]))
									except: flux.append(float(line[13:24].split("-")[0]+"E-"+line[13:24].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									2#print("01,12|", line[1:12], line[13:24])
									
									try:    flux.append(float(line[25:36]))
									except: flux.append(float(line[25:36].split("-")[0]+"E-"+line[25:36].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									try:    flux.append(float(line[37:48]))
									except: flux.append(float(line[37:48].split("-")[0]+"E-"+line[37:48].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									#print("23,26|", line[25:36], line[37:48])
									
									try:     flux.append(float(line[49:60]))
									except:  flux.append(float(line[49:60].split("-")[0]+"E-"+line[49:60].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									try:     flux.append(float(line[61:75]))
									except:  flux.append(float(line[61:75].split("-")[0]+"E-"+line[61:71].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									#print("48,60|", line[49:60], line[61:72])
									list_found+=6
				            
								except Exception as e: None
									#print(e)
				
				repeated_wls = np.tile(wl, Effective_counter)
				for l in range(len(repeated_wls)):
					wl_all.append(repeated_wls[l])
		
		
		
		
		
		
		wl_all=np.asarray(wl_all).astype("float")
		flux=np.asarray(flux).astype("float")
		Teff=np.asarray(Teff).astype("float")
		Grav=np.asarray(Grav).astype("float")
		
		
		
		maskTeffs=Teff<11000
		
		wl_all,  flux,  Teff,  Grav  =  wl_all[maskTeffs],  flux[maskTeffs],  Teff[maskTeffs],  Grav[maskTeffs]
		
		wl_all,  flux,  Teff,  Grav  =  wl_all.tolist(),  flux.tolist(),  Teff.tolist(),  Grav.tolist()
		
		
		cwd=install_path+"/saved_grids_npy/3D_DBA"
		
		list_found=0
		for filename in os.listdir(cwd):
			Effective_counter=0
			Teff_just_this_file=[]
			if "DBA" in filename and not ".py" in filename and "HtoHe-5" in filename:
				trigger_wl=True
				wl=[]
				#print(filename)
				f = open(cwd+"/"+filename).readlines()
				for count, line in enumerate(f):
					if not count == 0:
						if "Effective" in f[count]:
							#print(list_found/6,count)
							Effective_counter+=1
							trigger_wl=False
							teff=line[28:35]
							grav=line[48:57]
				    
							Teff_just_this_file.append(float(teff))
							list_found=0
			
						else:
							if trigger_wl==True:
								try:
									wl.append(float(line[1:10]))
									wl.append(float(line[11:20]))
									wl.append(float(line[21:30]))
									wl.append(float(line[32:40]))
									wl.append(float(line[42:50]))
									wl.append(float(line[52:60]))
									wl.append(float(line[62:70]))
									wl.append(float(line[72:80]))
									wl.append(float(line[82:90]))
									wl.append(float(line[92:100]))
									#print(line[1:10], line[11:20], line[21:30], line[32:40], line[42:50], line[52:60], line[62:70], line[72:80], line[82:90], line[92:100])
				            
								except Exception as e:
									None
							else:
								try:
									#flux.append(float(line[:13]))
									#print(line[61:72])
									try:    flux.append(float(line[1:12]))
									except: flux.append(float(line[1:12].split("-")[0]+"E-"+line[1:12].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									try:    flux.append(float(line[13:24]))
									except: flux.append(float(line[13:24].split("-")[0]+"E-"+line[13:24].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									2#print("01,12|", line[1:12], line[13:24])
									
									try:    flux.append(float(line[25:36]))
									except: flux.append(float(line[25:36].split("-")[0]+"E-"+line[25:36].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									try:    flux.append(float(line[37:48]))
									except: flux.append(float(line[37:48].split("-")[0]+"E-"+line[37:48].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									#print("23,26|", line[25:36], line[37:48])
									
									try:     flux.append(float(line[49:60]))
									except:  flux.append(float(line[49:60].split("-")[0]+"E-"+line[49:60].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									try:     flux.append(float(line[61:75]))
									except:  flux.append(float(line[61:75].split("-")[0]+"E-"+line[61:71].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									#print("48,60|", line[49:60], line[61:72])
									list_found+=6
				            
								except Exception as e:
									print(e)
				
				repeated_wls = np.tile(wl, Effective_counter)
				for l in range(len(repeated_wls)):
					wl_all.append(repeated_wls[l])
		
		
		
		
		
		
		
		
		wl_all=np.asarray(wl_all).astype("float")
		flux=np.asarray(flux).astype("float")
		Teff=np.asarray(Teff).astype("float")
		Grav=np.asarray(Grav).astype("float")
		
		maskMod500 = Teff%250!=0
		Teff[maskMod500]=np.round(Teff[maskMod500],-3)
		
		
		Grav[Grav>100]=np.log10(Grav[Grav>100])
		Grav=np.round(Grav,2)
		
		
		mask_wl_all = ((wl_all > minwl) & (wl_all <= maxwl))
		
		
		return wl_all[mask_wl_all], flux[mask_wl_all], Teff[mask_wl_all], Grav[mask_wl_all]
	
	
	def load_models_Pier_DA_ELM(minwl=3000,maxwl=8000):
		# get wavelengths
		trigger_wl=True
		flux, Teff, Grav, wl_all=[], [], [], []
		install_path = os.environ['WD_BASS_INSTALL_DIR']
		workingdir=install_path+"/saved_grids_npy/grid_elm_new/"
		for filename in os.listdir(workingdir):
			Effective_counter=0
			if not filename.endswith(".py") and not filename.endswith(".tar"):
				trigger_wl=True
				wl=[]
				f = open(workingdir+filename).readlines()
				for count, line in enumerate(f):
					if not count == 0:
						a = (f[count].split(" "))
						if "Effective" in f[count]:
							Effective_counter+=1
							trigger_wl=False
							thing=True
							while thing==True:
								thing=False
								for j in range(len(a)):
									if a[j] == "":
										a.pop(j)
										thing=True
										break
							for k in range(len(a)):
								if a[k] == "temperature":
									teff=float(a[k+2])
								if a[k] == "gravity":
									grav = float(a[k+2])
				    
				    
				    
						else:
							for i in range(len(a)):
								if not a[i] == "":
									if "\n" in a[i] :
										a[i]=a[i][:-1]
										#print(a[i], a[0])
									if trigger_wl==True:
										wl.append(float(a[i]))
									else:
										if "-1" in a[i] and not "E-1" in a[i]:
											a[i] = (a[i][:-4]+"E"+a[i][-4:])
										flux.append(float(a[i]))
										Teff.append(teff)
										Grav.append(grav)
				                
				repeated_wls = np.tile(wl, Effective_counter)
				for l in range(len(repeated_wls)):
					wl_all.append(repeated_wls[l])
				
			
			
		# grav is always the same, teff changes
		mask_wl = (np.asarray(wl_all) > minwl) & (np.asarray(wl_all) <= maxwl)
			
		ret_Teff=np.asarray(Teff)[mask_wl]
		ret_grav=np.asarray(Grav)[mask_wl]
		ret_wl_all=np.asarray(wl_all)[mask_wl]
		ret_flux=np.asarray(flux)[mask_wl]

		return ret_wl_all, ret_flux, ret_Teff, np.round(np.log10(ret_grav),3)
	
	
	
	def load_models_DBA(minwl=3000,maxwl=8000):
		# get wavelengths
		big_flux_list=[]
		wl_all, flux = [], []
		Teff, Grav = [], []
		
		install_path = os.environ['WD_BASS_INSTALL_DIR']
		cwd=install_path+"/saved_grids_npy/he-grid_v3"

		# get wavelengths
		big_flux_list=[]
		wl_all, flux = [], []
		Teff, Grav, HoverHe = [], [], []
		

		list_found=0
		for filename in os.listdir(cwd):
			Effective_counter=0
			Teff_just_this_file=[]
			if not ".py" in filename:
				HHe = filename[4:7]
				trigger_wl=True
				wl=[]
				f = open(cwd+"/"+filename).readlines()
				for count, line in enumerate(f):
					if not count == 0:
						if "Effective" in f[count]:
							#print(list_found/6,count)
							Effective_counter+=1
							trigger_wl=False
							teff=line[28:35]
							grav=line[48:57]
				    
							Teff_just_this_file.append(float(teff))
							list_found=0
			
						else:
							if trigger_wl==True:
								try:
									#print(line[1:10], line[11:20], line[21:30], line[32:40], line[42:50], line[52:60], line[62:70], line[72:80], line[82:90], line[92:100])
									wl.append(float(line[1:10]))
									wl.append(float(line[11:20]))
									wl.append(float(line[21:30]))
									wl.append(float(line[31:40]))
									wl.append(float(line[41:50]))
									wl.append(float(line[51:60]))
									wl.append(float(line[61:70]))
									wl.append(float(line[71:80]))
									wl.append(float(line[81:90]))
									wl.append(float(line[91:100]))
									#print(line[1:10], line[11:20], line[21:30], line[32:40], line[42:50], line[52:60], line[62:70], line[72:80], line[82:90], line[92:100])
				            
								except Exception as e:
									None
							else:
								try:
									#flux.append(float(line[:13]))
									#print(line[61:72])
									try:    flux.append(float(line[1:12]))
									except: flux.append(float(line[1:12].split("-")[0]+"E-"+line[1:12].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									HoverHe.append(HHe)
									try:    flux.append(float(line[13:24]))
									except: flux.append(float(line[13:24].split("-")[0]+"E-"+line[13:24].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									HoverHe.append(HHe)
									2#print("01,12|", line[1:12], line[13:24])
									
									try:    flux.append(float(line[25:36]))
									except: flux.append(float(line[25:36].split("-")[0]+"E-"+line[25:36].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									HoverHe.append(HHe)
									try:    flux.append(float(line[37:48]))
									except: flux.append(float(line[37:48].split("-")[0]+"E-"+line[37:48].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									HoverHe.append(HHe)
									#print("23,26|", line[25:36], line[37:48])
									
									try:     flux.append(float(line[49:60]))
									except:  flux.append(float(line[49:60].split("-")[0]+"E-"+line[49:60].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									HoverHe.append(HHe)
									try:     flux.append(float(line[61:75]))
									except:  flux.append(float(line[61:75].split("-")[0]+"E-"+line[61:71].split("-")[1]))
									Teff.append(float(teff))
									Grav.append(float(grav))
									HoverHe.append(HHe)
									#print("48,60|", line[49:60], line[61:72])
									list_found+=6
				            
								except Exception as e: None
									#print(e)
				
				repeated_wls = np.tile(wl, Effective_counter)
				for l in range(len(repeated_wls)):
					wl_all.append(repeated_wls[l])
		
		
		
		
		
		
		wl_all=np.asarray(wl_all).astype("float")
		flux=np.asarray(flux).astype("float")
		Teff=np.asarray(Teff).astype("float")
		Grav=np.asarray(Grav).astype("float")
		HoverHe=np.asarray(HoverHe)
		
		wavelength_grid = np.unique(wl_all)
		
		
		if False:  maskTeffs=(Teff<11000) & (HoverHe=="020" | HoverHe=="050" | HoverHe=="300")
		else:      maskTeffs=Teff==Teff
		
		wl_all,  flux,  Teff,  Grav,  HoverHe  =  wl_all[maskTeffs],  flux[maskTeffs],  Teff[maskTeffs],  Grav[maskTeffs],  HoverHe[maskTeffs]
		
		wl_all,  flux,  Teff,  Grav,  HoverHe  =  wl_all.tolist(),  flux.tolist(),  Teff.tolist(),  Grav.tolist(),  HoverHe.tolist()
		
		
		
		
		
		wl_all_3D,  flux_3D,  Teff_3D,  Grav_3D,  HoverHe_3D = [], [], [], [], []
		
		cwd=install_path+"/saved_grids_npy/3D_DBA"
		
		list_found=0
		for filename in os.listdir(cwd):
			Effective_counter=0
			Teff_just_this_file=[]
			if "DBA" in filename and not ".py" in filename and ("HtoHe-30" in filename or "HtoHe-5" in filename or "HtoHe-2" in filename):
				if "e-30" in filename:  HHe="300"
				elif "e-5" in filename:  HHe="050"
				elif "e-2" in filename:  HHe="020"
				
				trigger_wl=True
				wl=[]
				#print(filename)
				f = open(cwd+"/"+filename).readlines()
				for count, line in enumerate(f):
					if not count == 0:
						if "Effective" in f[count]:
							#print(list_found/6,count)
							Effective_counter+=1
							trigger_wl=False
							teff=line[28:35]
							grav=line[48:57]
				    
							Teff_just_this_file.append(float(teff))
							list_found=0
			
						else:
							if trigger_wl==True:
								try:
									wl.append(float(line[1:10]))
									wl.append(float(line[11:20]))
									wl.append(float(line[21:30]))
									wl.append(float(line[32:40]))
									wl.append(float(line[42:50]))
									wl.append(float(line[52:60]))
									wl.append(float(line[62:70]))
									wl.append(float(line[72:80]))
									wl.append(float(line[82:90]))
									wl.append(float(line[92:100]))
									#print(line[1:10], line[11:20], line[21:30], line[32:40], line[42:50], line[52:60], line[62:70], line[72:80], line[82:90], line[92:100])
				            
								except Exception as e:
									None
							else:
								try:
									#flux.append(float(line[:13]))
									#print(line[61:72])
									try:    flux_3D.append(float(line[1:12]))
									except: flux_3D.append(float(line[1:12].split("-")[0]+"E-"+line[1:12].split("-")[1]))
									Teff_3D.append(float(teff))
									Grav_3D.append(float(grav))
									HoverHe_3D.append(HHe)
									try:    flux_3D.append(float(line[13:24]))
									except: flux_3D.append(float(line[13:24].split("-")[0]+"E-"+line[13:24].split("-")[1]))
									Teff_3D.append(float(teff))
									Grav_3D.append(float(grav))
									HoverHe_3D.append(HHe)
									2#print("01,12|", line[1:12], line[13:24])
									
									try:    flux_3D.append(float(line[25:36]))
									except: flux_3D.append(float(line[25:36].split("-")[0]+"E-"+line[25:36].split("-")[1]))
									Teff_3D.append(float(teff))
									Grav_3D.append(float(grav))
									HoverHe_3D.append(HHe)
									try:    flux_3D.append(float(line[37:48]))
									except: flux_3D.append(float(line[37:48].split("-")[0]+"E-"+line[37:48].split("-")[1]))
									Teff_3D.append(float(teff))
									Grav_3D.append(float(grav))
									HoverHe_3D.append(HHe)
									#print("23,26|", line[25:36], line[37:48])
									
									try:     flux_3D.append(float(line[49:60]))
									except:  flux_3D.append(float(line[49:60].split("-")[0]+"E-"+line[49:60].split("-")[1]))
									Teff_3D.append(float(teff))
									Grav_3D.append(float(grav))
									HoverHe_3D.append(HHe)
									try:     flux_3D.append(float(line[61:75]))
									except:  flux_3D.append(float(line[61:75].split("-")[0]+"E-"+line[61:71].split("-")[1]))
									Teff_3D.append(float(teff))
									Grav_3D.append(float(grav))
									HoverHe_3D.append(HHe)
									#print("48,60|", line[49:60], line[61:72])
									list_found+=6
				            
								except Exception as e:
									print(e)
				
				repeated_wls = np.tile(wl, Effective_counter)
				for l in range(len(repeated_wls)):
					wl_all_3D.append(repeated_wls[l])
		
		
		
		if False:
			Teff_3D=np.round(np.asarray(Teff_3D),-3)
			Grav_3D=np.asarray(Grav_3D)
			HoverHe_3D=np.asarray(HoverHe_3D)
			wl_all_3D=np.asarray(wl_all_3D)
			flux_3D=np.asarray(flux_3D)
			
			
			wl_all,  flux,  Teff,  Grav,  HoverHe  =  np.asarray(wl_all),  np.asarray(flux),  np.asarray(Teff),  np.asarray(Grav),  np.asarray(HoverHe)
			
			raise ValueError(np.unique(Teff_3D))
			
			for Ts in np.unique(Teff_3D):
				for Gs in np.unique(Grav_3D):
					for HHes in np.unique(HoverHe_3D):
						mask=(Teff_3D==Ts) & (Grav_3D==Gs) & (HoverHe_3D==HHes)
						wltemp=wl_all_3D[mask]
						fluxtemp=flux_3D[mask]
						
						fluxinterp_newgrid=np.interp(wavelength_grid, wltemp, fluxtemp)
						
						wl_all = np.append(wl_all, np.full((len(wavelength_grid),),wavelength_grid))
						flux = np.append(flux, np.full((len(wavelength_grid),),fluxinterp_newgrid))
						Teff = np.append(Teff, np.full((len(wavelength_grid),),Ts))
						Grav = np.append(Grav, np.full((len(wavelength_grid),),Gs))
						HoverHe = np.append(HoverHe, np.full((len(wavelength_grid),),HHes))
					
					
		
		
		
		HoverHe=np.asarray(HoverHe).astype(float)/10
		wl_all=np.asarray(wl_all).astype("float")
		flux=np.asarray(flux).astype("float")
		Teff=np.asarray(Teff).astype("float")
		Grav=np.asarray(Grav).astype("float")
		
		if False:
			maskMod500 = Teff%250!=0
			Teff[maskMod500]=np.round(Teff[maskMod500],-3)
		
		
		Grav[Grav>100]=np.log10(Grav[Grav>100])
		Grav=np.round(Grav,2)
		
		
		mask_wl_all = ((wl_all > minwl) & (wl_all <= maxwl))  &  (Teff>=4000)
		
		
		return wl_all[mask_wl_all], flux[mask_wl_all], Teff[mask_wl_all], Grav[mask_wl_all], HoverHe[mask_wl_all]



	def load_models_Subdwarfs(minwl=3000,maxwl=8000, Hemin=-10, Hemax=10, loggmin=2, loggmax=10, teffmin=1, teffmax=100000, half_Teff_res=True):
		if True: #try: 
			#raise ValueError
			root="/home/james/HarryTest/"
			wl_all, flux, Teff, Grav, He_all = [], [], [], [], []
			for fil in ["Subdwarf_Grid_Hemin-1.0_Hemax0.0_Z0.00.npy", "Subdwarf_Grid_Hemin-2.0_Hemax-1.0_Z0.00.npy", "Subdwarf_Grid_Hemin-3.0_Hemax-2.0_Z0.00.npy", "Subdwarf_Grid_Hemin-4.0_Hemax-3.0_Z0.00.npy", "Subdwarf_Grid_Hemin-5.0_Hemax-4.0_Z0.00.npy", "Subdwarf_Grid_Hemin-6.0_Hemax-5.0_Z0.00.npy"]:
				Hemin_indiv = float(fil.split("Hemin")[-1].split("_")[0])
				Hemax_indiv = float(fil.split("Hemax")[-1].split("_")[0])
				if Hemin_indiv>=Hemin and Hemax_indiv<=Hemax:
					w, f, T, G, He = np.load(root+fil)
					maskwl = (w<maxwl) & (w>minwl) & (T<=teffmax) & (T>=teffmin) & (G<=loggmax) & (G>=loggmin)
					if half_Teff_res:   maskwl = maskwl & (T%2000==0)
					
					wl_all.extend(w[maskwl]); flux.extend(f[maskwl]); Teff.extend(T[maskwl]); Grav.extend(G[maskwl]);  He_all.extend(He[maskwl])
			
			wl_all, flux, Teff, Grav, He_all = np.asarray(wl_all), np.asarray(flux), np.asarray(Teff), np.asarray(Grav), np.asarray(He_all)
			return wl_all, flux, Teff, Grav, He_all
		
		#except Exception as qq:
		#	raise ValueError("Not liked input: subdwarfs", qq)
	
	
	
	


## https://www.aanda.org/articles/aa/olm/2013/11/aa22318-13/aa22318-13.html
#def ML18_to_3D_dTeff(Teff,logg):
#	# JM: NOT TESTED
#	A1=1.0947335E-03
#	A2=-1.8716231E-01
#	A3=1.9350009E-02
#	A4=6.4821613E-01
#	A5=-2.2863187E-01
#	A6=5.8699232E-01
#	A7=-1.0729871E-01
#	A8=1.1009070E-01
#	
#	Teff0=(Teff-10000.0)/1000.00
#	logg0=(logg-8.00000)
#	
#	
#	Shift=A1+(A2+A7*Teff0+A8*logg0)*np.exp(-(A3+A5*Teff0+A6*logg0)**2*((Teff0-A4)**2))
#	
#	return Shift*1000.00
#
## https://www.aanda.org/articles/aa/olm/2013/11/aa22318-13/aa22318-13.html
#def ML18_to_3D_dlogg(Teff,logg):
#	# JM: NOT TESTED
#	A1=7.5209868E-04
#	A2=-9.2086619E-01
#	A3=3.1253746E-01
#	A4=-1.0348176E+01
#	A5=6.5854716E-01
#	A6=4.2849862E-01
#	A7=-8.8982873E-02
#	A8=1.0199718E+01
#	A9=4.9277883E-02
#	A10=-8.6543477E-01
#	A11=3.6232756E-03
#	A12=-5.8729354E-02
#	
#	Teff0=(Teff-10000.0)/1000.00
#	logg0=(logg-8.00000)
#	ML18_to_3D_dlogg=(A1+A5*np.exp(-A6*((Teff0-A7)**2)))+A2*np.exp(-A3*((Teff0-(A4+A8*np.exp(-(A9+A11*Teff0+A12*logg0)**2*((Teff0-A10)**2))))**2))
#	
#	return ML18_to_3D_dlogg







# check convolution works 
#=============#=============#=============#=============#=============#=============#=============#=============#=============#=============#=============
#import matplotlib.pyplot as plt
#wl_all_3Dp1, flux_3Dp1, Teff_3Dp1, Grav_3Dp1 = load_models.load_models_DA_3D_NLTE(minwl=1, maxwl=9E9,toConvolve=False)
#wl_all_3Dp2, flux_3Dp2, Teff_3Dp2, Grav_3Dp2 = load_models.load_models_DA_3D_NLTE(minwl=1, maxwl=9E9,toConvolve=True)

#maskG8 = Grav_3Dp1==8

#wl_all_3Dp1, flux_3Dp1, Teff_3Dp1, Grav_3Dp1 = wl_all_3Dp1[maskG8], flux_3Dp1[maskG8], Teff_3Dp1[maskG8], Grav_3Dp1[maskG8]
#wl_all_3Dp2, flux_3Dp2, Teff_3Dp2, Grav_3Dp2 = wl_all_3Dp2[maskG8], flux_3Dp2[maskG8], Teff_3Dp2[maskG8], Grav_3Dp2[maskG8]

#for TT in np.unique(Teff_3Dp1):
#	mask=(Teff_3Dp1==TT) & ((wl_all_3Dp1>3500) & (wl_all_3Dp1<7000))
#	plt.plot(wl_all_3Dp2[mask], flux_3Dp2[mask],c='r');  plt.plot(wl_all_3Dp1[mask], flux_3Dp1[mask], c='g')
#	plt.title(str(TT));  plt.show()
#=============#=============#=============#=============#=============#=============#=============#=============#=============#=============#=============

#load_models.load_models_DA_pier_3D_new(minwl=1,maxwl=9E9)
#wl_all, flux, Teff, Grav  =  load_models.load_models_DBA_1eMinus5(minwl=1,maxwl=10000)





#wl_all, flux, Teff, Grav, HoverHe  =  load_models.load_models_DBA(minwl=3000,maxwl=8000)
#print(np.unique(Grav))



if False:
	wl, fl, teff, logg = load_models.load_models_Pier_DA_ELM(minwl=3500,maxwl=7000)

	import matplotlib.pyplot as plt
	for i in np.unique(teff):
		for j in np.unique(logg):
			plt.figure(figsize=(16,10))
			mask = (teff==i) & (logg==j)
			plt.plot(wl[mask], fl[mask],c='k')
			plt.xlabel("wl")
			plt.ylabel("fl")
			plt.title(str(i) + "   " + str(j))
			plt.show()

	print(np.unique(teff), np.unique(logg))

if False:
	import matplotlib.pyplot as plt
	import sys
	sys.path.append("/home/james/PostDocPierMachineLearning")
	from normalise_Spectra import normalise_spectra


	wl_all, flux, Teff, Grav, HoverHe = load_models.load_models_DBA(minwl=3700,maxwl=4300)
	for i in np.unique(Teff):
		if i<10000 or i%1000!=0: continue
		for j in np.unique(Grav):
			for k in np.unique(HoverHe):
				if float(k) == 2.0: continue
				plt.figure(figsize=(16,10))
				mask = (Teff==i) & (Grav==j) & (HoverHe==k)
				
				#flux_norm = normalise_spectra.normalise_spectra(np.array([wl_all[mask]]), np.array([0]), np.array([flux[mask]]), np.array([np.full((len(flux[mask]),),1)]), np.array(["DB"]))
				
				plt.plot(wl_all[mask], flux[mask],c='k')
				plt.title("Teff="+str(i) + "   logg=" + str(j) + "   HoverHe" + str(k))
				plt.show()
		
		
		
		
		
		
		
		
		

