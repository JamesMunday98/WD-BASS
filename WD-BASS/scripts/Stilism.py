#!/usr/bin/env python3

"""stilism_dist.py: Get reddening from stilism web site for stars with l, b and distance"""
__author__ = "Nicolas Leclerc"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Nicolas Leclerc"
__status__ = "Production"

import numpy as np
import pandas, requests, os
from io import StringIO


def complete_with_stilism(input_csv):
    url = "http://stilism.obspm.fr/reddening?frame=galactic&vlong={}&ulong=deg&vlat={}&ulat=deg&distance={}"
    
    df = pandas.read_csv(input_csv)
    df.loc[:, "distance[pc][stilism]"] = np.nan
    df.loc[:, "reddening[mag][stilism]"] = np.nan
    df.loc[:, "distance_uncertainty[pc][stilism]"] = np.nan
    df.loc[:, "reddening_uncertainty_min[mag][stilism]"] = np.nan
    df.loc[:, "reddening_uncertainty_max[mag][stilism]"] = np.nan

    for index, row in df.iterrows():
        print("l:", row["l"], "deg, b:", row["b"], "deg, distance:", row["distance"], "pc")
        res = requests.get(url.format(row["l"], row["b"], row["distance"]), allow_redirects=True)
        if res.ok:
            file = StringIO(res.content.decode("utf-8"))
            dfstilism = pandas.read_csv(file)
            print(dfstilism)
            df.loc[index, "distance[pc][stilism]"] = dfstilism["distance[pc]"][0]
            df.loc[index, "reddening[mag][stilism]"] = dfstilism["reddening[mag]"][0]
            df.loc[index, "distance_uncertainty[pc][stilism]"] = dfstilism["distance_uncertainty[pc]"][0]
            df.loc[index, "reddening_uncertainty_min[mag][stilism]"] = dfstilism["reddening_uncertainty_min[mag]"][0]
            df.loc[index, "reddening_uncertainty_max[mag][stilism]"] = dfstilism["reddening_uncertainty_max[mag]"][0]

    filename, extension = os.path.splitext(input_csv)
    df.to_csv("%s_with_stilism_reddening%s"%(filename, extension), index=False)


    
