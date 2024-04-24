import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import griddata

all_massesDB, all_tempDB, all_loggDB, all_radiusDB = [], [], [], []
for fi in os.listdir("/home/james/python_scripts_path/dwd_fit_package/DBmassradius/DB/Z002"):
	if ".trk" in fi:
		logL, logTeff, logTc, logRho_c, logAbundanceHCentre, logAbundanceHeCentre, OuterMassFraction, InnerMassConv, logAge, Mass, LogMdot, modelnum, LogLpp, LogLcno, LogLHe, LogLC, LogL_g, LogLnu, LogMHtot, LogMHeBuffer, MassHfreeCore, MassHefreeCore, Logg, RoverRsun, _, _ = np.loadtxt(fi, comments="#", unpack=True)
		
		mask=(10**logTeff<30000) & (10**logTeff>5000)
		
		#plt.scatter(10**logTeff[mask], RoverRsun[mask])
		#plt.title(str(fi))
		#plt.show()
		
		for mm, tt, lglg, rr in zip(Mass[mask], 10**logTeff[mask], Logg[mask], RoverRsun[mask]):
			all_massesDB.append(mm);   all_tempDB.append(tt);   all_loggDB.append(lglg);   all_radiusDB.append(rr)

all_massesDB, all_tempDB, all_loggDB, all_radiusDB = np.asarray(all_massesDB), np.asarray(all_tempDB), np.asarray(all_loggDB), np.asarray(all_radiusDB)


rad_He = float(griddata(np.array([all_tempDB.ravel(),all_loggDB.ravel()]).T,all_radiusDB.ravel(),np.array([10000, 8]).T, method='linear'))

print(np.amin(all_loggDB), np.amax(all_loggDB))
