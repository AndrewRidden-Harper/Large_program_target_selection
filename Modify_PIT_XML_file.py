# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 11:47:01 2022

@author: Andrew Ridden-Harper

This script opens an XML file created by the Gemini Phase 1 Tool (PIT) and 
modifies it to add targets and observation times.

This script was created to assist with submitting hundreds of possible targets
in our Gemini Long and Large Program (LLP). Note that as GRACES is only 
on the telescope for a few weeks each semester, the large number of potential 
targets submitted at Phase I are significantly reduced to only those that are 
visible during the periods when GRACES is on the telescope

"""

import copy
import pandas as pd

df = pd.read_csv('2022B_ LLP Transits  - LLP Targets w_ Transits.csv')

XML_file_name = 'Andrew_version_LLP2022B_No_Times.xml'

f = open(XML_file_name,'r')
#s = f.read()

# TargetNameString = 'observation band="Band 1/2" enabled="true" target='

# TargetNameIndex = s.find(TargetNameString)


# TargetNumbertartIndex = TargetNameIndex+len(TargetNameString) + 1

# TargetNumberEndIndex = s.find('"',TargetNumbertartIndex)

# print(s[TargetNumbertartIndex:TargetNumberEndIndex])

# #TargetName = s[]

TargetNumberToNameDict = {}

l = f.readlines()

#LastTargetNumberStringList = ['test']

for LineIndex in range(len(l)):
    
    line = l[LineIndex]
    
    if ('sidereal epoch' in line) & ('id=' in line):         
        
        SplitLine = line.split('=')
        TargetNumberStr = SplitLine[-1].strip('\'>\n')
        
        LastTargetNumberString = copy.copy(TargetNumberStr)
        
        NameLine = l[LineIndex+1]
        
        NameLine2 = NameLine.replace('<','>')
        
        SplitNameLine = NameLine2.split('>')
        
        Name = SplitNameLine[2]
        
        TargetNumberToNameDict[TargetNumberStr] = Name        
        
        #raise Exception
        


f.close()

f2 = open(XML_file_name,'r')
s = f2.read()
        
index1 = s.find(LastTargetNumberString) + 100  ### Adding 100 to get it well past the last target info target-216 since it was just finding that instead of the next part in the obserations section

TargetNumberList = list(TargetNumberToNameDict.keys())

#TargetNumberList = TargetNumberList[0:1]
#TargetNumberList = ['target-186']

AlternativeNamesDict = {}

AlternativeNamesDict['KELT-19 b'] = 'KELT-19 A b'


for TargetNumberIndex in range(len(TargetNumberList)):
#for TargetNumberIndex in [216]:
    
    TargetNumberString = TargetNumberList[TargetNumberIndex]
    
    TargetName = TargetNumberToNameDict[TargetNumberString] 
    
    if 'Qatar' in TargetName:
        TargetName = TargetName.replace(' ','-',1)
        
    
    if TargetName in AlternativeNamesDict.keys():
        TargetName = AlternativeNamesDict[TargetName]
        
    if 'TIC' in TargetName:
        
        targetdf = df.loc[df['TIC Name Andrew'] == TargetName]        
        ObTime = targetdf['Observation duration (h)'].values
        
    else:           
        targetdf = df.loc[df['Planet Name'] == TargetName]        
        ObTime = targetdf['Observation duration (h)'].values
        
    # raise Exception
    
    # ObTime = 4.7
    
    
    InsertString = '            <progTime units="hr">%.1f</progTime>\n            <partTime units="hr">0.0</partTime>\n            <time units="hr">%.1f</time>\n'%(ObTime,ObTime)    
    
    
    
    index2 = s.find(TargetNumberString,index1)
    
    index3 = s.find('>',index2) + 2
    
    sBefore = s[0:index3]
    sAfter = s[index3:]
    
    #raise Exception
    
    s = sBefore + InsertString + sAfter
    

outf=open('LLP2022B_with_ObTimes.xml','w')
###outf=open('test.xml','w')
outf.write(s)    
outf.close()
    


        
        
        