#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 20 10:01:23 2022

@author: james
"""

from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u
import numpy as np
from miscAstro import miscAstro

class checkLocalStars(object):
    def find_star_in_gaia_edr3(RAdeg,Decdeg, radius_arcsec, predicted_Gmag=None):
        v = Vizier(columns=['RA_ICRS', 'DE_ICRS', 'Gmag', 'Plx', 'e_Plx', 'BP-RP'])
        obj=coord.SkyCoord(ra=RAdeg, dec=Decdeg, unit=(u.deg, u.deg), frame='icrs')
        
        # grab star information
        result = v.query_region(obj, radius=radius_arcsec * u.arcsec, catalog='I/350/gaiaedr3')[0] # gaia edr3   https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/350&-out.max=50&-out.form=HTML%20Table&-out.add=_r&-out.add=_RAJ,_DEJ&-sort=_r&-oc.form=sexa
        mask1=(result["BP-RP"]<1.5) & (result['Plx']>0)
        result=result[mask1]
        
        if len(result['Gmag'])>1 and predicted_Gmag==None:
        	print("Multiple/No results in Gaia. I am not able to know if there is a nearby and bright contaminant to the photometry. I found results for:")
        	print(result)
        	print("around the coordinates")
        	print(obj)
        
        
        else:
            if len(result['Gmag']) == 1:
                star_mag=result['Gmag']
                star_plx=result['Plx']
                star_e_plx=result['e_Plx']
                RAdeg=result["RA_ICRS"]
                Decdeg=result["DE_ICRS"]
            elif len(result['Gmag'])>1 and predicted_Gmag!=None:
                ## search for most probable target
                minarg = np.argmin(np.abs(result['Gmag']-predicted_Gmag))
                star_mag=result['Gmag'][minarg]
                star_plx=result['Plx'][minarg]
                star_e_plx=result['e_Plx'][minarg]
                RAdeg=result["RA_ICRS"][minarg]
                Decdeg=result["DE_ICRS"][minarg]
        
            
            result_large = v.query_region(obj, radius=7.5 * u.arcsec, catalog='I/350/gaiaedr3')[0] # gaia edr3   https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/350&-out.max=50&-out.form=HTML%20Table&-out.add=_r&-out.add=_RAJ,_DEJ&-sort=_r&-oc.form=sexa
            
            
            for count, row in enumerate(result_large): 
                if row['Gmag'] == star_mag:  result_large.remove_row(count)
            
            warning=False
            
            if len(result_large)>=1 and np.amin(result_large['Gmag'])<15:
                print("ALERT ALERT")
                print("ALERT ALERT")
                print("ALERT ALERT")
                print("The photometry may be compromised. I found a star brighter than mag 15 within 7.5 arcsec of the source")
                warning=True
            
            
            
            
            
            result_very_large = v.query_region(obj, radius=30 * u.arcsec, catalog='I/350/gaiaedr3')[0] # gaia edr3   https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/350&-out.max=50&-out.form=HTML%20Table&-out.add=_r&-out.add=_RAJ,_DEJ&-sort=_r&-oc.form=sexa
        
            for count, row in enumerate(result_very_large): 
                if row['Gmag'] == star_mag:  result_very_large.remove_row(count)
            
            if len(result_very_large)>=1 and np.amin(result_very_large['Gmag'])<11:
                print("ALERT ALERT")
                print("ALERT ALERT")
                print("ALERT ALERT")
                print("The photometry may be compromised. I found a star brighter than mag 11 within 30 arcsec of the source")
                warning=True
            
            
            #return np.asarray(RAdeg)[0], np.asarray(Decdeg)[0], np.asarray(star_mag)[0], np.asarray(star_plx)[0], np.asarray(star_e_plx)[0], warning
            if type(RAdeg)==float:
                return RAdeg, Decdeg, star_mag, star_plx, star_e_plx, warning
            else:
                return np.asarray(RAdeg)[0], np.asarray(Decdeg)[0], np.asarray(star_mag)[0], np.asarray(star_plx)[0], np.asarray(star_e_plx)[0], warning


#checkLocalStars.find_star_in_gaia_edr3(121.59563339554, 15.45861172208, 5, 20)

