#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 14:42:01 2022

@author: Andrew Ridden-Harper

This script was created to make it easier to reduce our large initial sample 
of potential target planets (submitted in Phase I) to a much smaller sample
that are visible during the periods of time when GRACES is on the telescope

It also determines which planets can be observed on each night, where night
dates are for the start of the night.  This required careful treatment 
of transit timings that occured after midnight (local time). While this script 
loads and processes transit times in local time (Hawaii standard time), times 
in UTC could also be used by doing the simple conversion (e.g. UTC -10)

"""

from astropy.time import Time,TimeDelta
import astropy.units as u
import itertools


### The basic time constraints (midnight on each date)
# t1 = Time('2022-09-29', format='isot', scale='utc')
# t2 = Time('2022-10-09', format='isot', scale='utc')

# t3 = Time('2023-01-20', format='isot', scale='utc')
# t4 = Time('2023-01-31', format='isot', scale='utc')



### I think more accurate time constraints, considering that the dates on the schedule are start of the night (starting at 16:00 and going to 07:00 on the date after)
# t1 = Time('2022-09-29T16:00:00', format='isot', scale='utc')
# t2 = Time('2022-10-10T07:00:00', format='isot', scale='utc')

# t3 = Time('2023-01-20T16:00:00', format='isot', scale='utc')
# t4 = Time('2023-02-01T07:00:00', format='isot', scale='utc')


#### For finding the additional transits in the adjusted GRACES windows 
t1 = Time('2022-09-22T16:00:00', format='isot', scale='utc')
t2 = Time('2022-09-28T07:00:00', format='isot', scale='utc')

# t3 = Time('2023-01-20T16:00:00', format='isot', scale='utc')
# t4 = Time('2023-02-01T07:00:00', format='isot', scale='utc')

InPeriod2 = True


OutputStr = ''

PlanetPerDayDict = {}
PlanetPerDayDict_FormattedForCopyingText = {}


# IngressEgressStr = 'Ingress (UT)      Egress (UT)\n'


with open('Aug-1-2022-Jan-31-2023-LLP_Transits.txt') as f:
#with open('Only_HAT-P-20b.txt') as f:
    lines = f.readlines()
    
    

    
WrittenPlanetNameList = []
    
for i in range(len(lines)):    
   
    line = lines[i] 
    
    if 'b' in line:        
        LastPlanetName = line
        
        
    elif 'c' in line:        
        LastPlanetName = line

        
    elif 'd' in line:        
        LastPlanetName = line
        

        
    elif 'TIC' in line:        
        LastPlanetName = line
        
        
    elif 'K2-' in line:        
        LastPlanetName = line       

        
    elif line == '\n':
        #OutputStr += line        
        pass
        
    elif line == 'Ingress (UT)      Egress (UT)\n':        
        pass
      
        
    else:
        
        #OutputStr += IngressEgressStr
        
        SplitDateLine = line.split()
        
        IngressDateStr = SplitDateLine[0]        
        SplitIngressDateStr = IngressDateStr.split('/')
        
        NeededFormatForDate = '20%s-%s-%s'%(SplitIngressDateStr[2],SplitIngressDateStr[0],SplitIngressDateStr[1])
        

        IngressTimeStr = SplitDateLine[1]
        
        SplitIngressTimeStr = IngressTimeStr.split(':')
        
        NeededFormatForTime = '%s:%s:00'%(SplitIngressTimeStr[0],SplitIngressTimeStr[1])
        
        NeededFormatForDateTime = '%sT%s'%(NeededFormatForDate,NeededFormatForTime)
       
        t = Time(NeededFormatForDateTime, format='isot', scale='utc')
        t_midnight_on_date = Time('%sT00:00:00'%(NeededFormatForDate), format='isot', scale='utc')
        
        TimeDeltaBack1Day = TimeDelta(-1*u.day)  
        
        MorningCutoffHr = 9.0
        
        ### If transit happens in the first 9 hours of the night (i.e., before 9 am) subtract 1 day from the date at midnight for the dictionary with keys of "night of" or "night starting"
        if (t.mjd - t_midnight_on_date.mjd) <= MorningCutoffHr/24.0:
            
            t_for_night_of_date_for_key = t_midnight_on_date + TimeDeltaBack1Day
            
            PlanetPerDayDictKey = str(t_for_night_of_date_for_key).split('T')[0]
            
        else:
            PlanetPerDayDictKey = str(t_midnight_on_date).split('T')[0]
        
        
        
        
        
        InPeriod1 = ((t.mjd >= t1.mjd) & (t.mjd<=t2.mjd))
        #InPeriod2 = ((t.mjd >= t3.mjd) & (t.mjd<=t4.mjd))
        
        # if LastPlanetName == 'HAT-P-20 b\n':
        #     raise Exception        
        
        
        if InPeriod1 or InPeriod2:
            
            if LastPlanetName not in WrittenPlanetNameList:                
                OutputStr += '\n'
                OutputStr += LastPlanetName
                OutputStr += 'Ingress (HST)    Egress (HST)\n' 
                WrittenPlanetNameList.append(LastPlanetName)

                
                
            OutputStr += line
            
            if PlanetPerDayDictKey not in PlanetPerDayDict.keys():
                PlanetPerDayDict[PlanetPerDayDictKey] = []
                
            PlanetPerDayDict[PlanetPerDayDictKey].append(LastPlanetName.strip('\n'))
            
            
for key in PlanetPerDayDict.keys():
    PlanetListString = ''
    for PlanetName in PlanetPerDayDict[key]:
        #PlanetListString += PlanetName+', '
        PlanetListString += PlanetName+'\n'
        
    #PlanetListString = PlanetListString[0:-2]        
    PlanetPerDayDict_FormattedForCopyingText[key] = PlanetListString
    
PlanetPerDayWriteFile = open('LLP_2022B_planets_per_day.txt','w')

for key in sorted(list(PlanetPerDayDict_FormattedForCopyingText.keys())):
    print()
    print(key+':')
    print(PlanetPerDayDict_FormattedForCopyingText[key])
    
    PlanetPerDayWriteFile.write('\n') 
    PlanetPerDayWriteFile.write(key+':'+'\n')
    PlanetPerDayWriteFile.write(PlanetPerDayDict_FormattedForCopyingText[key]+':'+'\n')
  
PlanetPerDayWriteFile.close()
    

        
#WriteFile = open('LLP_2022B_transits_in_GRACES_periods_BetterLastAndFirstDayV3.txt','w')
#WriteFile = open('test_write.txt','w')
WriteFile = open('LLP_2022B_transits_Sept22-28.txt','w')


WriteFile.write(OutputStr)
WriteFile.close()

PlanetNameListStr = ''
for i in WrittenPlanetNameList:
    PlanetNameListStr += i

with open('PlanetNameList_Sept22-28.txt','w') as PlanetNameListFile:
    PlanetNameListFile.write(PlanetNameListStr)

### Check that all planet names in WrittenPlanetNameList are also in PlanetPerDayDict_FormattedForCopyingText

for PlanetName in WrittenPlanetNameList:
    
    InDict = False    
    
    for key in sorted(list(PlanetPerDayDict.keys())):
        #print(PlanetPerDayDict[key])
        
        if PlanetName[0:-1] in PlanetPerDayDict[key]:
            
            InDict = True 
        
    if not InDict:
        raise Exception
        
### Check the opposite of above 
ListOfDictValues = list(itertools.chain.from_iterable(PlanetPerDayDict.values()))

for PlanetName in ListOfDictValues:
    if PlanetName+'\n' not in WrittenPlanetNameList:
        raise Exception


