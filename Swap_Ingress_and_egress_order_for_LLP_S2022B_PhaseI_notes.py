# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 16:45:08 2022

@author: Andrew Ridden-Harper

A simple script to swap the order of ingress and egress
columns in a .txt file of ingress and egress times 
with each transit being on a new line

"""

f = open('LLP_S2022B_PhaseI_notes_to_modify.txt')
l = f.readlines()
f.close()

headerline = 'Ingress (UT)      Egress (UT)\n'

for LineIndex in range(len(l)):
    
    line = l[LineIndex] 
    
    if 'Egress' in line:
        l[LineIndex] = headerline
        
    else:
        
        l[LineIndex] = l[LineIndex].strip(' \n\t') 
        l[LineIndex] = l[LineIndex].replace('\t',' ')
        l[LineIndex] += '\n'
        
StringToWrite = ''

for line in l:
    StringToWrite += line
    
outf = open('Modified_LLP_S2022B_PhaseI_notes.txt','w')
outf.write(StringToWrite)
outf.close()
