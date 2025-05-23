import sys, os
import numpy as np
from subprocess import call
import subprocess, yaml
install_path = os.environ['WD_BASS_INSTALL_DIR']
from natsort import natsorted


filesout = os.listdir(os.getcwd()+"/out")
files_one_by_one = []
for ff in filesout:
	if "result" in ff and ".out" in ff:
		files_one_by_one.append(ff)

try: files_one_by_one.remove("result_one_by_one.out")
except: None

try: files_one_by_one.remove("result.out")
except: None

for ff in files_one_by_one:
	os.remove(os.getcwd()+"/out/"+ff)



single_double=sys.argv[1]

if single_double=="single":
	with open('example_Config_sgl.yaml') as file:
		config_info = yaml.load(file, Loader=yaml.FullLoader)
		input_files = np.asarray(config_info["spectraHa"])
		reference_wl = np.asarray(config_info["refwave"])
		modelHa = np.asarray(config_info["modelHa"])
		normaliseHa_all = np.asarray(config_info["normaliseHa"])
		cut_Ha_all = np.asarray(config_info["cut_Ha_max"])
		resolutions = np.asarray(config_info["resolutions"])
		share_rv = np.asarray(config_info["share_rv"])
		starType1 = np.asarray(config_info["starType"])
		RV_boundaries1 = np.asarray(config_info["RV_boundaries1"])
		HJD_values = np.asarray(config_info["HJD_Values"])

		file_ignore=np.asarray([config_info["file_ignore"]])[0]
	    
elif single_double=="double":
	with open('example_Config_dbl.yaml') as file:
		config_info = yaml.load(file, Loader=yaml.FullLoader)
		input_files = np.asarray(config_info["spectraHa"])
		reference_wl = np.asarray(config_info["refwave"])
		modelHa = np.asarray(config_info["modelHa"])
		normaliseHa_all = np.asarray(config_info["normaliseHa"])
		cut_Ha_all = np.asarray(config_info["cut_Ha_max"])
		resolutions = np.asarray(config_info["resolutions"])
		share_rv = np.asarray(config_info["share_rv"])
		starType1, starType2 = np.asarray(config_info["starType"])
		RV_boundaries1 = np.asarray(config_info["RV_boundaries1"])
		RV_boundaries2 = np.asarray(config_info["RV_boundaries2"])
		HJD_values = np.asarray(config_info["HJD_Values"])

		file_ignore=np.asarray([config_info["file_ignore"]])[0]
else:  raise ValueError("Unclear if single or double for file Initiate_spectral_fitting_one_by_one.py")





if len(file_ignore) > 0:  dowhile=True
else: dowhile=False

while dowhile:
	stopwhile=True
	original_input_files = input_files
	original_share_rv = share_rv
	if False:
		for aaaa in file_ignore:
			try:
				original_file_index = np.argwhere(original_input_files==aaaa)[0][0]
			except: raise ValueError(aaaa)
			for zcnt, aaa_sharerv in enumerate(original_share_rv):
				if aaa_sharerv>=original_file_index:
					share_rv[zcnt]-=1
			
			#print(len(normaliseHa_all),len(cut_Ha_all),len(resolutions),len(reference_wl),len(RV_boundaries1),len(RV_boundaries2), len(HJD_values), len(input_files))
			mask_ignore_files = input_files!=aaaa
			modelHa = modelHa[mask_ignore_files]
			normaliseHa_all=normaliseHa_all[mask_ignore_files]
			cut_Ha_all=cut_Ha_all[mask_ignore_files]
			resolutions=resolutions[mask_ignore_files]
			reference_wl=reference_wl[mask_ignore_files]
			share_rv=share_rv[mask_ignore_files]
			RV_boundaries1=RV_boundaries1[mask_ignore_files]
			try: RV_boundaries2=RV_boundaries2[mask_ignore_files]
			except: None
			HJD_values=HJD_values[mask_ignore_files]
			input_files=input_files[mask_ignore_files]
			
			file_ignore_mask=file_ignore!=aaaa
			file_ignore=file_ignore[file_ignore_mask]
			
			
			stopwhile=False
			break
	if stopwhile: dowhile=False


done_wanted_index=[]
finished_index=[]
for i in np.unique(share_rv[share_rv!=-1]):
	if True: #i==3:
		wanted_index = (share_rv==i) | (share_rv==share_rv[i])
		
		if len(share_rv)==0:
			continue

		print(i, len(share_rv))


		print("Handling share_rv=", i)
		fp = open("output.txt", "w")
		if single_double=="single":
			p = subprocess.Popen(["mpiexec", "-np", "5", "python3", install_path+"single_FitSpectrum_MCMC_multfiles_faster2.py", "ATM", "run_single_one_by_oneMCMC", str(i)], stdout=fp)
		elif single_double=="double":
			p = subprocess.Popen(["mpiexec", "-np", "5", "python3", install_path+"double_FitSpectrum_MCMC_multfiles_faster.py", "ATM", "run_double_one_by_oneMCMC", str(i)], stdout=fp)
		p.wait()
		fp.close()
		
		done_wanted_index.append(wanted_index)
		
		finished_index.append(i)


for cn, shrv in enumerate(share_rv):
	if not cn in finished_index and shrv==-1: #i==3:
		wanted_index = (share_rv==cn) | (share_rv==share_rv[cn])
		
		if len(share_rv)==0:
			continue

		print(cn, len(share_rv))


		print("Handling share_rv=", cn)
		fp = open("output.txt", "w")
		if single_double=="single":
			p = subprocess.Popen(["mpiexec", "-np", "5", "python3", install_path+"single_FitSpectrum_MCMC_multfiles_faster2.py", "ATM", "run_single_one_by_oneMCMC", str(cn)], stdout=fp)
		elif single_double=="double":
			p = subprocess.Popen(["mpiexec", "-np", "5", "python3", install_path+"double_FitSpectrum_MCMC_multfiles_faster.py", "ATM", "run_double_one_by_oneMCMC", str(cn)], stdout=fp)
		p.wait()
		fp.close()
		
		done_wanted_index.append(wanted_index)
		
		finished_index.append(cn)




if len(done_wanted_index) < 1:
	for cn, _ in enumerate(share_rv[share_rv==-1]):
		print("Handling spectrum=", cn)
		fp = open("output.txt", "w")
		if single_double=="single":
			p = subprocess.Popen(["mpiexec", "-np", "5", "python3", install_path+"single_FitSpectrum_MCMC_multfiles_faster2.py", "ATM", "run_single_one_by_oneMCMC", str(cn)], stdout=fp)
		elif single_double=="double":
			p = subprocess.Popen(["mpiexec", "-np", "5", "python3", install_path+"double_FitSpectrum_MCMC_multfiles_faster.py", "ATM", "run_double_one_by_oneMCMC", str(cn)], stdout=fp)
		p.wait()
		fp.close()
		







all_lines_RV1, all_lines_RV2 = [], []
for aresult in natsorted(os.listdir(os.getcwd()+"/out")):
	if ("result" in aresult and ".out" in aresult) and (not aresult=="result.out") and (not aresult=="result_one_by_one.out"):
		with open("out/"+aresult) as resultfile:
			result_lines = resultfile.readlines()
			for iii, lll in enumerate(result_lines):
				if "RV1" in lll:   all_lines_RV1.append(result_lines[iii+1])
				if "RV2" in lll:   all_lines_RV2.append(result_lines[iii+1])




with open("out/result_one_by_one.out", 'w') as output_file:
	output_file.write("RV1:\n")
	for aline in all_lines_RV1:
		output_file.write(aline)
	
	if len(all_lines_RV2)>0:
		output_file.write("RV2:\n")
		for aline in all_lines_RV2:
			output_file.write(aline)


for aresult in os.listdir(os.getcwd()+"/out"):
	if ("result" in aresult and ".out" in aresult) and (not aresult=="result.out") and (not "result_one_by_one" in aresult):
		os.remove(os.getcwd()+"/out/"+aresult)

	
