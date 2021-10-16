#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:14:43 2020

@author: ariddenharper
"""

import numpy as np
import pandas as pd 
from astropy import units as u
from astropy import constants as const
import matplotlib.pyplot as plt 


def CalcPlanetSignalStrength(NumberOfScaleHeights):
    SignalStrength = (((NumberOfScaleHeights*ScaleHeight_km.to(u.m) + pl_rad_m)**2 - pl_rad_m**2)/(((df['st_rad'].values*u.Rsun).to(u.m))**2)).value
    
    return SignalStrength


def CalcTotalDataSNR(PassedSNR_from_ETC,NumLinesOfSpecies=1,FuncNumPixPerLine=3):
    
    TotalDataSNR = PassedSNR_from_ETC*np.sqrt(NumExpInTrans)*np.sqrt(FuncNumPixPerLine)*np.sqrt(NumTransitsObserved)*np.sqrt(NumLinesOfSpecies)
    
    return TotalDataSNR


#ResolvingPower = 60e3
ResolvingPower = 67.5e3

#MeanMolecularWeightValue = 28.964
MeanMolecularWeightValue = 2.3
#MeanMolecularWeightValue = 10

MeanMoleculeMass = (MeanMolecularWeightValue*u.u).to(u.kg) ## 2.3u for H/He dominated

#### For H-alpha
NumberOfHScaleHeights = 37  #### measured from the top of Rp (ie will be Rp + NH)
NumHSpectralLines = 1 

#NumberOfNaScaleHeights = 20
NumberOfNaScaleHeights = 25
NumNaSpectralLines = 2


### original numbers 
NumberOfH2OScaleHeights = 3#0.5
#NumberOfH2OScaleHeights = 0.5
NumH2OSpectralLines = 1500

NaFWHM_A = 0.52

SpecResEl_A = 5889.950/ResolvingPower

#NumNaPix = NaFWHM_A/SpecResEl_A
NumNaPix = 6 ## rounded up from above 

NumPixPerLine = 2

NumTransitsObserved = 1


# NumberOfH2OScaleHeights = 0#0.01
# NumH2OSpectralLines = 500

#### For He
#NumberOfHScaleHeights = 80 ### As in Nortmann et al 2018


### For Subaru 
#ReadoutTime_s = 60*u.s
### For CARMEMES
#ReadoutTime_s = 40*u.s


#ReadoutTime_s = 35*u.s
ReadoutTime_s = 25*u.s ### From the output in the GRACES ITC 



### For Subaru:
#deltaRadVelPerPix = 2.6*u.km/u.s


deltaRadVelPerPix = (2.998e5/ResolvingPower)*u.km/u.s



### CARMENES https://carmenes.caha.es/ext/publications/CARMENES_201807_SPIE_Quirrenbachetal.pdf

### For calculating number of pixels per line in CARMENES
# =============================================================================
# CARMENES_vis_res = 94600
# KELT9bHaFWHM_kms = 38.3*(u.km/u.s)
# 
# c = const.c.to(u.km/u.s)
# 
# RefWave_A = 6562.8
# res_el_A = RefWave_A/CARMENES_vis_res
# res_el_kms = c/CARMENES_vis_res
# 
# pix_per_res_el = 2.8
# 
# ### For CARMENES 
# deltaRadVelPerPix = res_el_kms/pix_per_res_el
# 
# NumResElInHalphaFWHM = KELT9bHaFWHM_kms/res_el_kms
# 
# NumPixInHalphaFWHM = NumResElInHalphaFWHM*pix_per_res_el
# =============================================================================



### For Subaru 
#NumPixPerLine = 10
## For CARMENES 
#NumPixPerLine = 33  ## calculated by above code 



#FileToLoadName = 'AllNorthernAStarsNoHeaderForCARMENES_Halpha'
#FileToLoadName = 'AllNorthernAStarsNoHeaderForCARMENES_HalphaTimes0-8'

#FileToLoadName = 'AllNorthernAStarsNoHeaderForCARMENES_He'
#FileToLoadName = 'AllNorthernAStarsNoHeaderForCARMENES_HeTimes0-8'

#FileToLoadName = 'AllNorthernAStarsNoHeader_ForSubaru2020BRedo'
#FileToLoadName = 'Gemini_LLP_2020_potential_target_list - TransitingDecGrt30S_10March2020_mag13_0.1%_10RE'
#FileToLoadName = 'Gemini_LLP_2020_potential_target_list - TransitingDecGrt30S_10March2020_mag13_0.1%_10RE_7h'
#FileToLoadName = 'Gemini_LLP_2020_potential_target_list - DecGrt30S_mag130.0pc_10RE_PastedCalcs'

#FileToLoadName = 'Gemini_LLP_2020_potential_target_list - TransitingDecGrt30S_10March2020_mag13_0.1%_SMALL'

#FileToLoadName = 'Gemini_LLP_2020_potential_target_list - FINAL'
#FileToLoadName = 'Gemini_LLP_2020_potential_target_list - Selection-NewNewUpdatedCalcs'
#FileToLoadName = 'Gemini_LLP_2020_potential_target_list - Selection-NewNewUpdatedCalcs_K2-226_correct_transit_duration'


##FileToLoadName = 'Gemini_LLP_2020_potential_target_list - Selection-NewNew'
#FileToLoadName = 'Gemini_LLP_2020_potential_target_list - FinalUpdatedCalculations2Transits'


#FileToLoadName = 'Gemini_LLP_2020_potential_target_list - Selection-NewNewUpdatedCalcs_K2-226_correct_transit_duration_new_SNR_GRACES_ITC_v2p0'

#FileToLoadName = 'Gemini_LLP_2020_potential_target_list - Selection-NewCubed_July15_GRACES_ITC_v2p0'

#FileToLoadName = 'Gemini_LLP_2020_potential_target_list - Selection-NewCubed_July15_GRACES_ITC_v2p0_GRACES_ITC_v2p0'
#FileToLoadName = 'Gemini_LLP_2020_potential_target_list - Selection-NewCubed_July15_GRACES_ITC_v2p0_GRACES_ITC_v2p0'


#FileToLoadName = 'Gemini_LLP_2020_potential_target_list - Selection-NewCubed_July15_GRACES_ITC_v2p0_GRACES_ITC_v2p0_GRACES_ITC_v2p0'

#FileToLoadName = 'Gemini_LLP_2020_potential_target_list - Selection-NewCubed_July15_SeveralGracesITCs'


#FileToLoadName = 'Gemini_LLP_2020_potential_target_list - Selection-NewCubed_July15_SeveralGracesITCs_GRACES_ITC_Incl_v2p1_GRACES_ITC_Incl_v2p1'


#FileToLoadName = 'Gemini_LLP_2020_potential_target_list - Selection-NewCubed_July15_SeveralGracesITCs_Incl_v2p1'

#FileToLoadName = 'Gemini_LLP_2020_potential_target_list_IncGliese486b-Selection-NewCubed_July15_SeveralGracesITCs_New_ITC_From_v2p1_GaiaMagAsR'


#FileToLoadName = 'Gemini_LLP_RefiningMarch2021_ExpTimes'

#FileToLoadName = 'LLP_Target_List_Reloaded-transit depth_0.1pc_FilledValuesV2_WithSNR'
#FileToLoadName = 'LLP_Target_List_Reloaded-transit depth_0.1pc_FilledValuesV4_WithSNR'

#FileToLoadName = 'LLP_Target_List_Reloaded-TOIs_WithSNR'

#FileToLoadName = 'Planets_in_S2021B_but_missing_from_Target_List_Reloaded_WithSNR'

FileToLoadName = 'Confirmed TOIs (10_15_2021) - Reduced_WithSNR'

#FileToLoadName = 'Confirmed TOIs (10_15_2021) - Reduced_LeftOvers_WithSNR'



#TransitDepthThreshold_pc = 0.03
#VMagThreshold = 13
#DecThreshold = 0

df = pd.read_csv('%s.csv'%(FileToLoadName))

#df2 = df[(df['pl_trandep']>TransitDepthThreshold_pc)&(df['st_optmag']<VMagThreshold)&(df['dec']>DecThreshold)]

#df2 = df[df['pl_tranflag']==1]

#df3 = df2[df['dec']>(-30)]

#df2.to_csv('%s_OnlyTransiting.csv'%(FileToLoadName),header=True)
#df3.to_csv('%s_OnlyTransitingDecGtr30S.csv'%(FileToLoadName),header=True)




#df = pd.read_csv('AllNorthernAStarsNoHeaderOldMASCARA1bValues.csv')

# ExposureTime_s = df['exp_time'].values*u.s
# ExposureTime_s = 300*u.s



#SNR_from_ETC = df['SNR_from_ETC']
#SNR_from_ETC = df['SNR_pixel_100']
#SNR_from_ETC = df['GRACES_ITC_v2p0_300s']
#SNR_from_ETC = df['GRACES_ITC_v2p0_100s_IQ85']

#SNR_from_ETC = df['GRACES_ITC_v2p0_100s_IQ85_SB80_CC70']
#SNR_from_ETC = df['GRACES_ITC_v2p0_300s_IQ85']
#SNR_from_ETC = df['GRACES_ITC_v2p0_60s_IQ85']
#SNR_from_ETC = df['GRACES_ITC_v2p0_60s_IQ100']
#SNR_from_ETC = df['GRACES_ITC_v2p0_300s_IQ100']
#SNR_from_ETC = df['GRACES_ITC_v2p1_60s_IQ85']
#SNR_from_ETC = df['GRACES_ITC_v2p1_60s_IQ100']
#SNR_from_ETC = df['GRACES_ITC_v2p1_300s_IQ100']

#SNR_from_ETC = df['GRACES_ITC_v2p1_100s_IQ100']
#SNR_from_ETC = df['GRACES_ITC_v2p1_150s_IQ100']

#SNR_from_ETC = df['GRACES_ITC_v2p1_60s_IQ100']

#SNR_from_ETC = df['GRACES_ITC_v2p1_60s_IQ100']

#SNR_from_ETC = df['GRACES_ITC_v2p1_60s_IQ100_GaiaMagAsR']

#SNR_from_ETC = df['GRACES_ITC_v2p1_60s_IQ100_UsingActualVMag_NaRegion']

#SNR_from_ETC = df['GRACES_ITC_v2p1_60s_IQ100_UsingActualVMag_Overall']

#SNR_from_ETC = df['GRACES_ITC_v2p1_300s_IQ100_UsingActualVMag_NaRegion']

#SNR_from_ETC = df['GRACES_ITC_v2p1_300s_IQ100_UsingActualVMag_Overall']

#SNR_from_ETC = df['GRACES_ITC_v2p1_VmagAsVmag_60s_NaRegion']
#SNR_from_ETC = df['GRACES_ITC_v2p1_VmagAsVmag_60s_OverallAverage']
#SNR_from_ETC = df['GRACES_ITC_v2p1_VmagAsVmag_120s_NaRegion']
#SNR_from_ETC = df['GRACES_ITC_v2p1_VmagAsVmag_120s_OverallAverage']
#SNR_from_ETC = df['GRACES_ITC_v2p1_VmagAsVmag_200s_NaRegion']
#SNR_from_ETC = df['GRACES_ITC_v2p1_VmagAsVmag_300s_NaRegion']

# SNR_from_ETC = df['GRACES_ITC_v2p1_60s_IQ100_OverallSNR']
# SNR_from_ETCNaRegion = df['GRACES_ITC_v2p1_60s_IQ100_NaRegionSNR']

SNR_from_ETC = df['GRACES_ITC_v2p1_300s_IQ100_OverallSNR']
SNR_from_ETCNaRegion = df['GRACES_ITC_v2p1_300s_IQ100_NaRegionSNR']

#ExposureTime_s = 60*u.s 
#ExposureTime_s = 120*u.s 
#ExposureTime_s = 200*u.s 
ExposureTime_s = 300*u.s 


pl_rad_m = (df['pl_radj'].values*u.Rjup).to(u.m)

#pl_dur_day = df['pl_trandur'].values*u.day
pl_dur_day = (df['Updated_trandur'].values*u.hour).to(u.day)



pl_orbper_day = df['pl_orbper'].values*u.day
pl_orbper_s = pl_orbper_day.to(u.s)

NumExpInTrans = (pl_dur_day/(ReadoutTime_s + ExposureTime_s)).decompose().value

df['NumExpInTrans'] = NumExpInTrans

PlanetGravity = const.G*(df['pl_massj'].values*u.Mjup).to(u.kg)/((pl_rad_m)**2)
df['PlanetGravity'] = PlanetGravity.value

#=W3*((Y3*0.00465)/(2*E3))^(1/2)


#=(((6.67408E-11)*(X3*1.988435E+30)*((C3*24*3600)^(2)))/(4*(PI())^2))^(1/3)

CalculatedPlanetSemiMajorAxis = (const.G*((df['st_mass'].values*u.Msun).to(u.kg))*(pl_orbper_s**2)/(4*(np.pi**2)))**(1/3)
CalculatedPlanetSemiMajorAxis_au = CalculatedPlanetSemiMajorAxis.to(u.au)

df['calc_pl_orbsmax'] = CalculatedPlanetSemiMajorAxis_au

CalculatedPlanetEffectiveTemp = (df['st_teff'].values*u.K)*(((df['st_rad'].values*u.Rsun).to(u.km))/(2*(CalculatedPlanetSemiMajorAxis.to(u.km))))**(1/2)
df['calc_pl_eqt'] = CalculatedPlanetEffectiveTemp


#ScaleHeight_km = ((const.k_B.decompose()*df['pl_eqt'].values*u.K)/(MeanMoleculeMass*PlanetGravity)).to(u.km)
ScaleHeight_km = ((const.k_B.decompose()*CalculatedPlanetEffectiveTemp)/(MeanMoleculeMass*PlanetGravity)).to(u.km)

df['calc_ScaleHeight_km'] = ScaleHeight_km.value

HSignalStrength = CalcPlanetSignalStrength(NumberOfHScaleHeights)
NaSignalStrength = CalcPlanetSignalStrength(NumberOfNaScaleHeights)
H2OSignalStrength = CalcPlanetSignalStrength(NumberOfH2OScaleHeights)

df['HSignalStrength'] = HSignalStrength
df['NaSignalStrength'] = NaSignalStrength
df['H2OSignalStrength'] = H2OSignalStrength


#orbvel = 2*np.pi*((df['pl_orbsmax'].values*u.AU).to(u.km))/pl_orbper_day.to(u.s)
orbvel = 2*np.pi*(CalculatedPlanetSemiMajorAxis_au.to(u.km))/pl_orbper_day.to(u.s)



df['orbvel'] = orbvel.value

RadVelShiftInTrans = 2*orbvel*np.sin((df['pl_orbincl_Fake90s'].values*u.deg).to(u.rad))*np.sin((2*np.pi*((pl_dur_day/2)/pl_orbper_day))*u.rad)
df['RadVelShiftInTrans'] = RadVelShiftInTrans.value

RadVelShiftInTrans_pix = RadVelShiftInTrans/deltaRadVelPerPix
df['RadVelShiftInTrans_pix'] = RadVelShiftInTrans_pix.value

RadVelShiftInExp = 2*orbvel*np.sin((df['pl_orbincl_Fake90s'].values*u.deg).to(u.rad))*np.sin((2*np.pi*((ExposureTime_s.to(u.day)/2)/pl_orbper_day))*u.rad)

RadVelShiftInExp2 = 2*orbvel*np.sin((df['pl_orbincl_Fake90s'].values*u.deg).to(u.rad))*np.sin((2*np.pi*((ExposureTime_s.to(u.day)/2)/pl_orbper_day))*u.rad)






RadVelAtStartOfTransit = orbvel*np.sin((df['pl_orbincl_Fake90s'].values*u.deg).to(u.rad))*np.sin((2*np.pi*((pl_dur_day/2)/pl_orbper_day))*u.rad)

RadVelAfterFirstExp = orbvel*np.sin((df['pl_orbincl_Fake90s'].values*u.deg).to(u.rad))*np.sin((2*np.pi*(((pl_dur_day-ExposureTime_s.to(u.day))/2)/pl_orbper_day))*u.rad)

RadVelShiftFirstExp = RadVelAtStartOfTransit - RadVelAfterFirstExp

orbinc_r = (df['pl_orbincl_Fake90s'].values*u.deg).to(u.rad).value
pl_orbper_day_val = pl_orbper_day.value
orbvel_val = orbvel.value
ExposureTime_day = ExposureTime_s.value/(24*3600)

RadVelShiftInExpNoUnits = 2*orbvel_val*np.sin(orbinc_r)*np.sin(np.pi*ExposureTime_day/pl_orbper_day_val)

DesiredDeltaRV = 1 
ExpTimeFor1kms_day = np.arcsin((DesiredDeltaRV/(2*orbvel_val*np.sin(orbinc_r))))*(pl_orbper_day_val/np.pi)
ExpTimeFor1kms_s = ExpTimeFor1kms_day*3600*24
df['ExpTimeFor1kms_s'] = ExpTimeFor1kms_s

#ExpTimeFor1kms = np.arcsin(1/(2*orbvel*np.sin((df['pl_orbincl'].values*u.deg).to(u.rad))))


df['RadVelShiftInExp'] = RadVelShiftInExp.value
#df['RadVelShiftInExp_pix'] = RadVelShiftInExp.value/deltaRadVelPerPix.value

HTotalDataSNR = CalcTotalDataSNR(NumHSpectralLines,NumHSpectralLines)
NaTotalDataSNR = CalcTotalDataSNR(SNR_from_ETCNaRegion,NumNaSpectralLines,NumNaPix)
H2OTotalDataSNR = CalcTotalDataSNR(SNR_from_ETC,NumH2OSpectralLines,NumPixPerLine)


#df['HTotalDataSNR'] = HTotalDataSNR
df['NaTotalDataSNR'] = NaTotalDataSNR
df['H2OTotalDataSNR'] = H2OTotalDataSNR


HTotalDataNoise = 1/HTotalDataSNR
NaTotalDataNoise = 1/NaTotalDataSNR
H2OTotalDataNoise = 1/H2OTotalDataSNR


#df['HTotalDataNoise'] = HTotalDataNoise
df['NaTotalDataSNR'] = NaTotalDataSNR
df['H2OTotalDataSNR'] = H2OTotalDataSNR


HDetectionSNR = HSignalStrength/HTotalDataNoise
NaDetectionSNR = NaSignalStrength/NaTotalDataNoise
H2ODetectionSNR = H2OSignalStrength/H2OTotalDataNoise


#df['HDetectionSNR'] = HDetectionSNR
df['NaDetectionSigma'] = NaDetectionSNR
df['H2ODetectionSigma'] = H2ODetectionSNR

xindexvect = np.arange(len(CalculatedPlanetEffectiveTemp))

plt.scatter(xindexvect,CalculatedPlanetEffectiveTemp,label='calculated value',s=5)
plt.scatter(xindexvect,df['pl_eqt'].values,s=80,facecolors='none',edgecolors='red',label='from NASA exoplanet archive')
plt.xlabel('planet index')
plt.ylabel('equilibrium temp')
plt.legend()

plt.figure()

plt.scatter(xindexvect,CalculatedPlanetSemiMajorAxis_au,label='calculated value',s=5)
plt.scatter(xindexvect,df['pl_orbsmax'].values,s=80,facecolors='none',edgecolors='red',label='from NASA exoplanet archive')
plt.xlabel('planet index')
plt.ylabel('semi-major axis')
plt.legend()


df.to_csv('%s_Calcs_MMW%.2f_ExpTime%.2f.csv'%(FileToLoadName,MeanMolecularWeightValue,ExposureTime_s.value),header=True)


#c = 299792.458 
#
#R = 30e3
#
#
#
#TransitStartPhase = (pl_dur_day/2)/pl_orbper_day  
#ExpTimeInPhase = (ExposureTime_s.to(u.day))/pl_orbper_day

#PhaseVect = np.arange(1-TransitStartPhase,1,ExpTimeInPhase)
#PhaseVect2 = np.arange(0,TransitStartPhase,ExpTimeInPhase)

#FullPhaseVect = []
#for i in PhaseVect1:
#    FullPhaseVect.append(i)
#for i in PhaseVect2:
#    FullPhaseVect.append(i)
#FullPhaseVect = np.array(FullPhaseVect)

#rv_at_each_phase = orbvel_val*np.sin(orbinc_r)*np.sin(2*np.pi*PhaseVect2)

#ChangeInRVAfterEachObservation = np.diff(rv_at_each_phase)






#plt.figure()
#plt.plot(rv_at_each_phase)
#
#plt.plot(ChangeInRVAfterEachObservation)
#
#MeanChangeInRVEachOb = np.mean(ChangeInRVAfterEachObservation)


#dv = c/(2*R) ### Since there seems to be a factor of 2 to make it match the paper