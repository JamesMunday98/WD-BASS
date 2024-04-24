import numpy as np
import matplotlib.pyplot as plt

fname, hjd, RV1, RV1err, RV2, RV2err = np.loadtxt("/home/james/Fit_DWDs_with_my_script/antoine_pier_mixed/Systems_hybrid_nogaia/WDJ181058.67+31194/RVfits/RVfit_results.dat", delimiter="\t", dtype=str, unpack=True)
hjd, RV1, RV1err, RV2, RV2err = hjd.astype(float), RV1.astype(float), RV1err.astype(float), RV2.astype(float), RV2err.astype(float)
plt.errorbar(hjd, RV1, yerr=RV1err, fmt='.b', label="RV1")
plt.errorbar(hjd, RV2, yerr=RV2err, fmt='.r', label="RV2")
plt.show()
