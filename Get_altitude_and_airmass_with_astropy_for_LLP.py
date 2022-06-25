# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 16:51:34 2022

@author: Andrew Ridden-Harper

This script loads transit ingress and egress times such as those generated 
by the NASA Exoplanet Archive Transit and Ephemeris service and calculates 
observation intervals by expanding the transit duration by an out-of-transit 
baseline and an additional time interval for observation scheduling flexibility.

It also calculates the airmass at calculated times, to see if the airmass at
the time of observation is within acceptable limits. 

It also writes a note for the observer, indicating the observing windows 
for our time-critical transit observations.

Note that the airmass calculations are based on this tutorial with astropy:
https://docs.astropy.org/en/stable/generated/examples/coordinates/plot_obs-planning.html

All unusual variable names come from the tutorial above.
The variable names have been kept the same but their values were changed for 
our LLP (e.g. setting the observatory location to that of Gemini-North)



"""

import astropy.units as u
from astropy.time import Time,TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd

def GiveAstropyTimeFromNASATransitFinderTime(NASAtimeStr):

    #calendar = calendarList[PlanetIndex]
    
    calendar = NASAtimeStr
    calendarSplit1 = calendar.split()
    
    datestr1 = calendarSplit1[0]
    datestr1Split = datestr1.split('/')
    
    timestr1 = calendarSplit1[1]
    
    ReformattedTimeStr = '%s-%s-%sT%s:00'%(datestr1Split[2],datestr1Split[0],datestr1Split[1],timestr1)
    
    astropytimeObj = Time(ReformattedTimeStr)
    
    return astropytimeObj

def ConvertDeltaTimeObjectToHourMinuteString(DeltaTimeObj):
    
    TimeDiff_hr = DeltaTimeObj.to(u.hour).value    

    IntNumWholeHours = int(TimeDiff_hr)

    IntNumMins = (TimeDiff_hr - IntNumWholeHours)*60
    
    return '%d hours %.0f mins'%(IntNumWholeHours,IntNumMins) 

def GetAltAndAirmass(LookupName,TimeObj):

    TargetSkyCoordinates = SkyCoord.from_name(LookupName)
        
    frame = AltAz(obstime=TimeObj,
                              location=bear_mountain)
    
    
    target_azs = TargetSkyCoordinates.transform_to(frame)
    
    airmasss = target_azs.secz
    
    altitude = target_azs.alt
    
    return airmasss.value, altitude

#df = pd.read_csv('LLP_2022B_updated_transit_timings - Sheet1.csv')

#df = pd.read_csv('LLP_2022B_updated_transit_timings - Sheet1(1).csv')
#df = pd.read_csv('LLP_2022B_updated_transit_timings - Sheet1.csv')

df = pd.read_csv('Good_TOI-849b_LLP_2022B_updated_transit_timings.csv')







# NameList = df['name'].values
# ingresscalendarList = df['ingresscalendar'].values
# egresscalendarList = df['egresscalendar'].values

NameList = df['Planet'].values
ingresscalendarList = df['Ingress (UT)'].values
egresscalendarList = df['Egress (UT)'].values


### bear_mountain used in the example that was modified for this
bear_mountain = EarthLocation(lat=19.82396*u.deg, lon=-155.46984*u.deg, height=4213*u.m) ## Location of Gemini-N/Mauna Kea
utcoffset = 0*u.hour  # Hawaii Standard time  ### For using UTC, set this to 0

AllNotesString = ''

#PlanetIndex = 0


for PlanetIndex in range(len(NameList)):
#for PlanetIndex in [0,1]:
    

    PlanetName = NameList[PlanetIndex]
    LookupName = PlanetName[0:-2]
    
    print('Doing PlanetIndex %d of %d (%s)'%(PlanetIndex,len(NameList),PlanetName))

    
    if PlanetName == 'Qatar-7 b':
        LookupName = '2MASS J23540364+3701185'
        
    if PlanetName == 'Qatar-5 b':
     	LookupName = 'UCAC3 265-004681'
         
    if PlanetName == 'KELT-19 A b':
        LookupName = 'KELT-19'
    
    
    IngressTimeObj = GiveAstropyTimeFromNASATransitFinderTime(ingresscalendarList[PlanetIndex])
    EgressTimeObj = GiveAstropyTimeFromNASATransitFinderTime(egresscalendarList[PlanetIndex])
    
    Delta45mins = TimeDelta(45*u.min)
    
    Delta30mins = TimeDelta(30*u.min)
    
    Delta1hour = TimeDelta(1*u.hour)
    
    Delta10hour = TimeDelta(10*u.hour)
    
    
    IdealStartTime = IngressTimeObj - Delta45mins
    IdealEndTime = EgressTimeObj + Delta45mins
    
    ProposedStartTime = IdealStartTime - Delta30mins
    ProposedEndTime = IdealEndTime + Delta30mins
    
    TotalProgramParterTime = IdealEndTime - IdealStartTime
    
    #TotalProgramParterTime_hr = TotalProgramParterTime.to(u.hour).value
    
    TotalWindowDuration = ProposedEndTime - ProposedStartTime
    #TotalWindowDuration_hr = TotalWindowDuration.to(u.hour).value
    
    
    TotalProgramPartnerTimeStr = ConvertDeltaTimeObjectToHourMinuteString(TotalProgramParterTime)
    TotalWindowDurationStr = ConvertDeltaTimeObjectToHourMinuteString(TotalWindowDuration)
    
    if str(ProposedStartTime)[0:10] == str(ProposedEndTime)[0:10]:
        IdealTimingWindowStr = '%s %s - %s'%(str(IdealStartTime)[0:10],str(IdealStartTime)[11:16],str(IdealEndTime)[11:16])
        ProposedTimingWindowStr = '%s %s - %s'%(str(ProposedStartTime)[0:10],str(ProposedStartTime)[11:16],str(ProposedEndTime)[11:16])
        IngressTimingStr = '%s %s - %s'%(str(IngressTimeObj)[0:10],str(IngressTimeObj)[11:16],str(EgressTimeObj)[11:16])
        
    if str(ProposedStartTime)[0:10] != str(ProposedEndTime)[0:10]:
        IdealTimingWindowStr = '%s %s - %s %s'%(str(IdealStartTime)[0:10],str(IdealStartTime)[11:16],str(IdealEndTime)[0:10],str(IdealEndTime)[11:16])
        ProposedTimingWindowStr = '%s %s - %s %s'%(str(ProposedStartTime)[0:10],str(ProposedStartTime)[11:16],str(ProposedEndTime)[0:10],str(ProposedEndTime)[11:16])
        IngressTimingStr = '%s %s - %s %s'%(str(IngressTimeObj)[0:10],str(IngressTimeObj)[11:16],str(EgressTimeObj)[0:10],str(EgressTimeObj)[11:16])
    
        
    airmass_ingress, alt_ingress = GetAltAndAirmass(LookupName,IngressTimeObj)
    airmass_egress, alt_egress = GetAltAndAirmass(LookupName,EgressTimeObj)
    
    
    airmass_IdealStart, alt_idealstart = GetAltAndAirmass(LookupName,IdealStartTime)
    airmass_IdealEnd, alt_egress = GetAltAndAirmass(LookupName,IdealEndTime)
    
    airmass_ProposedStart, alt_proposedstart = GetAltAndAirmass(LookupName,ProposedStartTime)
    airmass_ProposedEnd, alt_proposedstart = GetAltAndAirmass(LookupName,ProposedEndTime)
    
    airmass_lim = 2.5
    
    if airmass_ingress >= airmass_lim:
        print('%s Ingress at airmass %.1f'%(PlanetName,airmass_ingress))
        # raise Exception
        
    if airmass_egress >= airmass_lim:
        print('%s Egress at airmass %.1f'%(PlanetName,airmass_egress))
        # raise Exception
    
    
    
    if airmass_IdealStart >= airmass_lim:
        print('%s IdealStart at airmass %.1f'%(PlanetName,airmass_IdealStart))            
        AllNotesString += '!!!! Ideal start is at airmass %.1f \n'%(airmass_IdealStart)

        # raise Exception
    
        
    if airmass_IdealEnd >= airmass_lim:
        print('%s Ideal End at airmass %.1f'%(PlanetName,airmass_IdealEnd))    
        AllNotesString += '!!!! Ideal end is at airmass %.1f \n'%(airmass_IdealEnd)


        # raise Exception
    
    
    if airmass_ProposedStart >= airmass_lim:
        print('%s Proposed Start at airmass %.1f'%(PlanetName,airmass_ProposedStart))    
        # raise Exception
        AllNotesString += '!!!! Default proposed start is at airmass %.1f \n'%(airmass_ProposedStart)
    
    if airmass_ProposedEnd >= airmass_lim:
        print('%s Proposed End at airmass %.1f'%(PlanetName,airmass_ProposedEnd))    
        AllNotesString += '!!!! Default proposed end is at airmass %.1f \n'%(airmass_ProposedEnd)

        # raise Exception
        
    
        
    NoteString = '''%s
    - Observe target once
    - All timing windows include 1.5 hours of total baseline to be disturbed between pre- and post- transit baselines
    - Observations can begin 30 minutes before or after the start of the ideal timing window.
    - The science observations should *NEVER* be started more than 1 hour after the indicated start time in the proposed timing window.
    - Ideal timing windows 45 mins of baseline before and after the transit.
    
    ---------------------------------
    - Total Program Partner Time: %s
    - Total Window Duration: %s
    - Proposed Timing Observing windows (UT): 
      %s 
    - Ideal Observing Window (45 mins baseline before and after transit) (UT)
      %s
    - Transit Ingress (UT), Egress (UT)
      %s'''%(PlanetName,
            TotalProgramPartnerTimeStr,
            TotalWindowDurationStr,
            ProposedTimingWindowStr,
            IdealTimingWindowStr,
            IngressTimingStr) 
    
    AllNotesString += NoteString
    AllNotesString += ' \n '
    AllNotesString += ' \n '
    AllNotesString += '*** End of timing note *** \n'
    AllNotesString += ' \n '



#TimingNoteWriteFile = open('Correct_TOI-849b_LLP_20222B_PhaseII_timing_notes.txt','w')
TimingNoteWriteFile = open('TEST.txt','w')
TimingNoteWriteFile.write(AllNotesString)
TimingNoteWriteFile.close()




#m33 = SkyCoord.from_name('WASP-25')
#m33 = SkyCoord.from_name(LookupName) ## TOI-2109
#time = Time('2022-6-15 00:00:00') - utcoffset

# m33altaz = m33.transform_to(AltAz(obstime=time,location=bear_mountain))
# print(f"M33's Altitude = {m33altaz.alt:.2}")

#midnight = Time('2022-3-22 10:00:00') - utcoffset
# midnight = Time('2022-6-17 10:00:00') - utcoffset ### This should be the middle of the local night in UTC 

# ## 

# #delta_midnight = np.linspace(-2, 10, 100)*u.hour
# #delta_midnight = np.linspace(5, 16, 1000)*u.hour
# delta_midnight = np.linspace(-6, 6, 1000)*u.hour



# frame_July13night = AltAz(obstime=midnight+delta_midnight,
#                           location=bear_mountain)

# m33altazs_July13night = m33.transform_to(frame_July13night)

# airmasss = m33altazs_July13night.secz

# altitude = m33altazs_July13night.alt

# plt.figure()
# plt.plot(delta_midnight, altitude)
# #plt.xlim(-2, 10)
# #plt.ylim(1, 4)
# plt.xlabel('Hours from Midnight')
# plt.ylabel('altiude')


# plt.figure()

# plt.plot(delta_midnight, airmasss)
# #plt.xlim(-2, 10)
# plt.ylim(1, 4)
# plt.xlabel('Hours from Midnight')
# plt.ylabel('Airmass [Sec(z)]')