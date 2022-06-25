# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 16:05:13 2022

@author: Andrew Ridden-Harper

This script takes a list of planet names as input and queries SIMBAD 
for RA, dec, magnitude, and proper motion. 

It then outputs this information in a target list correctly formatted for
loading into the Gemini Phase 1 Tool (PIT).

This was created to enter hundreds of potential targets at Phase I for our 
"ExoGems" Gemini Large and Long Program (LLP).

"""

import numpy as np
import matplotlib.pyplot as plt 
from astroquery.simbad import Simbad


import pandas as pd


 


#


customSimbad = Simbad()

customSimbad.add_votable_fields('pmdec')
customSimbad.add_votable_fields('pmra')




customSimbad.add_votable_fields('flux_system(B)')
customSimbad.add_votable_fields('flux_system(V)')
customSimbad.add_votable_fields('flux_system(G)')
customSimbad.add_votable_fields('flux_system(J)')
customSimbad.add_votable_fields('flux_system(H)')
customSimbad.add_votable_fields('flux_system(K)')

customSimbad.add_votable_fields('flux(B)')
customSimbad.add_votable_fields('flux(V)')
customSimbad.add_votable_fields('flux(G)')
customSimbad.add_votable_fields('flux(J)')
customSimbad.add_votable_fields('flux(H)')
customSimbad.add_votable_fields('flux(K)')

AlternativeNameDict = {}
AlternativeNameDict['KELT-19 A b'] = 'KELT-19 b'


df = pd.read_csv('2022B_ LLP Transits  - LLP Targets w_ Transits.csv')


PlanetNamesIncludigTOIS = df['Planet Name']
PlanetNamesAsTICs = df['TIC Name Andrew']


for PlanetNameIndex in range(len(PlanetNamesIncludigTOIS)):
    
    TOIName = PlanetNamesIncludigTOIS[PlanetNameIndex]
    TICName = PlanetNamesAsTICs[PlanetNameIndex]    
    
    if 'TOI' in TOIName:
        AlternativeNameDict[TOIName] = TICName
        
    


#customSimbad.add_votable_fields('fluxdata(V)')

ToWriteToFile = 'Name,RAJ2000,DecJ2000,pmRA,pmDec,B,B_sys,V,V_sys,G,G_sys,J,J_sys,H,H_sys,K,K_sys\n'


#Simbad.get_field_description('V')



#df = pd.read_csv('All Possible LLP Targets - All - Observed_NoTOIs.csv')



#for PlanetName in df['Planet']:
    
for PlanetName in df['Planet Name']:

    
#for PlanetName in ['HATS-37 b']:
    
    
    if PlanetName in AlternativeNameDict.keys():
        PlanetName = AlternativeNameDict[PlanetName]
    
    
    if 'Qatar' in PlanetName:
        PlanetName = PlanetName.replace('-',' ')
        
    if 'TOI' in PlanetName:
        pass
        
    
    LookupName = PlanetName[0:-2]
    
    print('Looking up: %s for %s'%(LookupName,PlanetName))

    
    #name = 'CoRoT-11'
    
    t = customSimbad.query_object(LookupName)
    
    
    RAListFormat1 = str(t['RA']).split('\n')
    RAStringToWrite = RAListFormat1[-1].replace(' ',':')
    
    
    DECListFormat1 = str(t['DEC']).split('\n')
    DECStringToWrite = DECListFormat1[-1].replace(' ',':')
    
    
    PMRA = str(t['PMRA']).split('\n')[-1].strip()
    
    PMDEC = str(t['PMDEC']).split('\n')[-1].strip()
    
    
    
    
    
    FLUX_B = str(t['FLUX_B']).split('\n')[-1].strip()
    FLUX_V = str(t['FLUX_V']).split('\n')[-1].strip()
    FLUX_G = str(t['FLUX_G']).split('\n')[-1].strip()
    FLUX_J = str(t['FLUX_J']).split('\n')[-1].strip()
    FLUX_H = str(t['FLUX_H']).split('\n')[-1].strip()
    FLUX_K = str(t['FLUX_K']).split('\n')[-1].strip()
    
    FLUX_SYSTEM_B = str(t['FLUX_SYSTEM_B']).split('\n')[-1].strip()
    FLUX_SYSTEM_V = str(t['FLUX_SYSTEM_V']).split('\n')[-1].strip()
    FLUX_SYSTEM_G = str(t['FLUX_SYSTEM_G']).split('\n')[-1].strip()
    FLUX_SYSTEM_J = str(t['FLUX_SYSTEM_J']).split('\n')[-1].strip()
    FLUX_SYSTEM_H = str(t['FLUX_SYSTEM_H']).split('\n')[-1].strip()
    FLUX_SYSTEM_K = str(t['FLUX_SYSTEM_K']).split('\n')[-1].strip()
    
    
    
    if FLUX_B == '--': 
        FLUX_B = ''
        
    if FLUX_V == '--': 
        FLUX_V = ''        

    if FLUX_G == '--': 
        FLUX_G = ''
        
    if FLUX_J == '--': 
        FLUX_J = ''
        
    if FLUX_H == '--': 
        FLUX_H = ''
        
    if FLUX_K == '--': 
        FLUX_K = ''
        


    
    
    
    LineToAdd = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'%(PlanetName,
                                                                   RAStringToWrite,
                                                                   DECStringToWrite,
                                                                   PMRA,
                                                                   PMDEC,
                                                                   FLUX_B,
                                                                   FLUX_SYSTEM_B,
                                                                   FLUX_V,
                                                                   FLUX_SYSTEM_V,
                                                                   FLUX_G,
                                                                   FLUX_SYSTEM_G,
                                                                   FLUX_J,
                                                                   FLUX_SYSTEM_J,
                                                                   FLUX_H,
                                                                   FLUX_SYSTEM_H,
                                                                   FLUX_K,
                                                                   FLUX_SYSTEM_K)
    
    
    
    ToWriteToFile += LineToAdd
    
                                                                  
    
    f = open('LLP_S2022B_PIT_targets.csv','w')
    f.write(ToWriteToFile)
    f.close()
    
    #t = result_table