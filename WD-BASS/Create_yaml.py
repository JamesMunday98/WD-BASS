import os, sys
sys.path.append(os.environ['WD_BASS_INSTALL_DIR']+"/scripts")
from miscAstro import miscAstro
from natsort import natsorted
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.io import fits

##  I assume that the time/HJD/BJD values desired are commented on the first line of the file. Otherwise, I will not place the HJDs into the config file



model_wavelengths_min, model_wavelengths_max = -90, 90
normalise_wavelengths_min, normalise_wavelengths_max = -91, 91
cut_wavelengths_min, cut_wavelengths_max = -100, 100
resolution_all = 8878
resolution_all_blue = 2428
reference_wavelength = 6562.81
reference_wavelength_blue = 4861.35
RVboundaries_1_min, RVboundaries_1_max = -200, 200
RVboundaries_2_min, RVboundaries_2_max = -200, 200


if len(sys.argv) > 1:
    single_or_double = sys.argv[1]
    if not (single_or_double=="single" or single_or_double=="double"): raise ValueError("Use is - 'python Create_dbl_yaml.py double'  /  'python Create_dbl_yaml.py single'")
else:
    raise ValueError("Use is - 'python Create_dbl_yaml.py double'  /  'python Create_dbl_yaml.py single'")

print(os.getcwd())

#science_spectra_directory = "/home/james/Desktop/FitDA_properly/grid_3D/dwd_fit_package/test_data_single/"
science_spectra_directory = os.getcwd()
if science_spectra_directory.endswith("/"): None
else: science_spectra_directory+="/"


string_to_write="spectraHa: "

num_files=0
fnames=[]
string_to_write+="["

vis_uvb=[]
for fi in natsorted(os.listdir(science_spectra_directory)):
	check=True
	if fi.endswith(".dat") and not "blue" in fi:
		try:
			ff=open(science_spectra_directory+fi)
			line1 = ff.readline()
			#print(line1)
			if "OrderedDict" in line1:
				HJD = line1.split("HJD")[1].split(")")[0][3:];    RAdeg = line1.split("RA")[1].split(")")[0][3:];    decdeg = line1.split("Dec")[1].split(")")[0][3:]
			elif line1.startswith("#"):
				line1=line1.split(" ")
				vals = [x for x in line1 if x]
				RAdeg=vals[1];    decdeg=vals[2];    HJD=vals[3][:-1]
			else:
				print("I can not read the HJD of file:  " + str(fi))
				print("Line 1 is:");    print(line1);    print()
				check=False
				
			if check==True:
				num_files+=1
				fnames.append(science_spectra_directory+str(fi))
				string_to_write+='"' + str(fi)+ '", '
		except: None
	elif fi.endswith(".fit") or fi.endswith(".fits"):
		hdul = fits.open(fi)
		header = hdul[0].header
		if header["HIERARCH ESO SEQ ARM"]=="VIS":
			num_files+=1
			fnames.append(science_spectra_directory+str(fi))
			vis_uvb.append("VIS")
			string_to_write+='"' + str(fi)+ '", '
		elif header["HIERARCH ESO SEQ ARM"]=="UVB":
			num_files+=1
			fnames.append(science_spectra_directory+str(fi))
			vis_uvb.append("UVB")
			string_to_write+='"' + str(fi)+ '", '


num_blue=0
for fi in natsorted(os.listdir(science_spectra_directory)):
	check=True
	if fi.endswith(".dat") and "blue" in fi:
		try:
			ff=open(science_spectra_directory+fi)
			line1 = ff.readline()
			#print(line1)
			if "OrderedDict" in line1:
				HJD = line1.split("HJD")[1].split(")")[0][3:];    RAdeg = line1.split("RA")[1].split(")")[0][3:];    decdeg = line1.split("Dec")[1].split(")")[0][3:]
			elif line1.startswith("#"):
				line1=line1.split(" ")
				vals = [x for x in line1 if x]
				RAdeg=vals[1];    decdeg=vals[2];    HJD=vals[3][:-1]
			else:
				print("I can not read the HJD of file:  " + str(fi))
				print("Line 1 is:");    print(line1);    print()
				check=False
				
			if check==True:
				num_files+=1
				fnames.append(science_spectra_directory+str(fi))
				string_to_write+='"' + str(fi)+ '", '
			num_blue+=1
		except: None
	elif fi.endswith(".fit") or fi.endswith(".fits"):
		hdul = fits.open(fi)
		header = hdul[0].header
		if header["HIERARCH ESO SEQ ARM"]=="VIS":
			num_files+=1
			fnames.append(science_spectra_directory+str(fi))
			vis_uvb.append("VIS")
			string_to_write+='"' + str(fi)+ '", '
		elif header["HIERARCH ESO SEQ ARM"]=="UVB":
			num_files+=1
			fnames.append(science_spectra_directory+str(fi))
			vis_uvb.append("UVB")
			string_to_write+='"' + str(fi)+ '", '
			num_blue+=1

for j in range(4):
	for i in fnames:
		if "blue" in i:
			string_to_write+='"' + str(i.split("/")[-1])+ '", '
	

string_to_write=string_to_write[:-2]
string_to_write+="] # full path desired\n"



string_to_write+="spectra_source_type: 'wl_flux_fluxerr' # need to INSERT your source type here. Options are 'wl_flux_fluxerr', 'Xshooter'. Default is 'wl_flux_fluxerr' \n"


string_to_write+="modelHa: ["
for i in range(num_files):
	string_to_write+="["+str(model_wavelengths_min)+","+str(model_wavelengths_max)+"], "
for i in range(4):
	for j in range(num_blue):
		if i==0: string_to_write+="["+str(-60)+","+str(60)+"], "
		elif i==1: string_to_write+="["+str(-50)+","+str(50)+"], "
		elif i==2: string_to_write+="["+str(-35)+","+str(35)+"], "
		elif i==3: string_to_write+="["+str(-15)+","+str(15)+"], "


string_to_write=string_to_write[:-2]+"]  # in angstroms"



string_to_write+="\nnormaliseHa: ["
for i in range(num_files):
	string_to_write+="["+str(normalise_wavelengths_min)+","+str(normalise_wavelengths_max)+"], "

for i in range(4):
	for j in range(num_blue):
		if i==0: string_to_write+="["+str(-61)+","+str(61)+"], "
		elif i==1: string_to_write+="["+str(-51)+","+str(51)+"], "
		elif i==2: string_to_write+="["+str(-36)+","+str(36)+"], "
		elif i==3: string_to_write+="["+str(-16)+","+str(16)+"], "

string_to_write=string_to_write[:-2]+"]  # in angstroms"


string_to_write+="\ncut_Ha_max: ["
for i in range(num_files):
	string_to_write+="["+str(cut_wavelengths_min)+","+str(cut_wavelengths_max)+"], "

for i in range(4):
	for j in range(num_blue):
		if i==0: string_to_write+="["+str(-70)+","+str(70)+"], "
		elif i==1: string_to_write+="["+str(-60)+","+str(60)+"], "
		elif i==2: string_to_write+="["+str(-41)+","+str(41)+"], "
		elif i==3: string_to_write+="["+str(-20)+","+str(20)+"], "

string_to_write=string_to_write[:-2]+"]  # in angstroms"



string_to_write+="\nresolutions: ["
if vis_uvb==[]:
	for i in range(num_files):
		if i<num_files-num_blue:  string_to_write+=str(resolution_all)+", "
		else:  string_to_write+=str(resolution_all_blue)+", "
	for i in range(4):
		for j in range(num_blue):
			if i==0: string_to_write+=str(2168)+", "
			elif i==1: string_to_write+=str(2049)+", "
			elif i==2: string_to_write+=str(1983)+", "
			elif i==3: string_to_write+=str(1942)+", "
else:
	for i, v_uv in zip(fnames, vis_uvb):
		if v_uv=="VIS":
			string_to_write+=str(10000)+", "
		elif v_uv=="UVB":
			string_to_write+=str(5000)+", "



string_to_write=string_to_write[:-2]+"] # R=lambda/dlambda"




string_to_write+="\nrefwave: ["
if vis_uvb==[]:
	for i in range(num_files):
		if i<num_files-num_blue:   string_to_write+=str(reference_wavelength)+", "
		else:   string_to_write+=str(reference_wavelength_blue)+", "
else:
	for i in vis_uvb:
		if i=="VIS":
			string_to_write+=str(6562.81)+", "
		elif i=="UVB":
			string_to_write+=str(4861.35)+", "

for i in range(4):
	for j in range(num_blue):
		if i==0: string_to_write+=str(4340.472)+", "
		elif i==1: string_to_write+=str(4101.734)+", "
		elif i==2: string_to_write+=str(3970.075)+", "
		elif i==3: string_to_write+=str(3889.064)+", "

string_to_write=string_to_write[:-2]+"]"


string_to_write+="\nshare_rv: ["
for i in range(num_files):
	string_to_write+="-1, "
for i in range(4):
	for j in range(num_blue):
		string_to_write+=str(-1)+", "
string_to_write=string_to_write[:-2]+"]  # 0 indexed!  needs all shared rvs to be at the end of the line, otherwise won't work properly!!!!!!"



if single_or_double=="double":
	string_to_write+="\nstarType: ["
	for i in range(2):
		string_to_write+='"DA", '
	string_to_write=string_to_write[:-2]+"]  # 'DA' or 'DBA'"
else:
	string_to_write+="\nstarType: ['DA'] # 'DA' or 'DBA"



string_to_write+="\nRV_boundaries1: ["
for i in range(num_files):
	if single_or_double=="single":
		string_to_write+="["+str(int(RVboundaries_1_min/2))+","+str(int(RVboundaries_1_max/2))+"], "
	else:
		string_to_write+="["+str(RVboundaries_1_min)+","+str(RVboundaries_1_max)+"], "
for i in range(4):
	for j in range(num_blue):
		string_to_write+="["+str(int(RVboundaries_1_min/2))+","+str(int(RVboundaries_1_max/2))+"], "
string_to_write=string_to_write[:-2]+"] # in kms-1"


if single_or_double=="double":
	string_to_write+="\nRV_boundaries2: ["
	for i in range(num_files):
		string_to_write+="["+str(RVboundaries_2_min)+","+str(RVboundaries_2_max)+"], "
	string_to_write=string_to_write[:-2]+"] # in kms-1"


string_to_write+="\n# please note that, where possible, the performance is about 1.4x when identical refwaves have the same normalisation/model/cut."



string_to_write+="\nHJD_Values: ["
for i in fnames:
	if i.endswith(".dat"):
		ff=open(i)
		line1 = ff.readline()
		#print(line1)
		if "OrderedDict" in line1:
			HJD = line1.split("HJD")[1].split(")")[0][3:];    RAdeg = line1.split("RA")[1].split(")")[0][3:];    decdeg = line1.split("Dec")[1].split(")")[0][3:]
		elif line1.startswith("#"):
			line1=line1.split(" ")
			vals = [x for x in line1 if x]
			RAdeg=vals[1];    decdeg=vals[2];    HJD=vals[3][:-1]
		else:
			print("I can not read the HJD of file:  " + str(i))
			print("Line 1 is:");    print(line1);    print()
		string_to_write+= str(HJD) + ", "
	elif i.endswith(".fits"):
		hdul = fits.open(i)
		header = hdul[1].header
		RAdeg = header["RA"];  decdeg = header["DEC"];  MJD=header["TMID"]  # MJD=header["MJD-OBS"] + (header["EXPTIME"]/86400)/2
		print("WARNING WARNING WARNING")
		print("ASSUMED PARANAL LOCATION")
		
		VLT= EarthLocation.from_geodetic(lat=-24.62722, lon=-70.40472, height=2635.43)   
		HJD=miscAstro.jd_corr(MJD, RAdeg, decdeg, VLT, jd_type='bjd')
		string_to_write+= str(HJD.tdb.value) + ", "
		
for i in range(4):
	for j in range(num_blue):
		string_to_write+=str(0) + ", "
		
string_to_write=string_to_write[:-2]+"]"


RAdeg, decdeg = float(RAdeg), float(decdeg)
RAhr, dechr = miscAstro.ra_dec_deg_to_hr(RAdeg,decdeg)

if single_or_double=="double":
	string_to_write+="\nforced_teff: [0, 0]   # number or 0"
	string_to_write+="\np0teff: [[5000,15000], [5000,15000]]   # number"
	string_to_write+="\nforced_logg: [0, 0]   # number or 0"
	string_to_write+="\np0logg: [[6.8,9], [6.8,9]]   # number"
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
for i in range(num_files):
	string_to_write+="4, "
for i in range(4):
	for j in range(num_blue):
		string_to_write+="4, "
string_to_write=string_to_write[:-2]+"]"

string_to_write+='\nparallax: "A"'
string_to_write+='\nparallax_uncertainty: "A"'
string_to_write+='\nphot_min_val: 0'
string_to_write+='\nplot_phot_spectrum: False'
string_to_write+="\n#to-do still = pre-compute the model spectrum when T1+logg1 or T2+logg2 is fixed"
string_to_write+="\n#when multiple lines are desired, precompute the entire model and then trim it/normalise it as is necessary, instead"
string_to_write+="\n#put all np.loadtxt things that are inside of /tmp into .npy format to make it faster https://stackoverflow.com/questions/40868231/fastest-most-optimized-way-to-read-write-ascii-file-in-python"
string_to_write+="\n#incorporate /Music/IstrateCooling/MTR_Istrate.py into my import function (all bits, so it is obvious what has been done in the future, like in def save_tables_output():)"
string_to_write+="\n#Herschel resolutions for Ha Hb Hg He: 8878, 2428, 2168, 2049, 1983, 1942, 1916"
string_to_write+="\n#INT resolution for Ha: 6310"
string_to_write+="\n#6562.81, 4861.35, 4340.472, 4101.734, 3970.075, 3889.064, 3835.397"


if single_or_double=="double":
	with open("example_Config_dbl.yaml", "w+") as f:
		f.write(string_to_write)
else:
	with open("example_Config_sgl.yaml", "w+") as f:
		f.write(string_to_write)


