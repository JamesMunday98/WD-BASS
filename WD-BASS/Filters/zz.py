import numpy as np
import matplotlib.pyplot as plt

wls, trans = np.loadtxt("GAIA_GAIA3.G.dat", unpack=True)

plt.plot(wls,trans)
plt.show()
