#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 17:46:55 2020

@author: ariddenharper
"""

import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
import pandas as pd 
import matplotlib.cm as cm
import matplotlib.ticker as ticker

plt.rc('font', family='serif', size=20)



HasBeenObserved = ['WASP-52 b', 
'HAT-P-32 b',
'WASP-177 b',
'WASP-76 b',
'GJ 3470 b (OLD)',
'WASP-140 b',
'WASP-49 b',
'GJ 436 b (OLD)',
'WASP-183 b',
'KELT-18 b',
'KELT-23 A b',
'TOI 851.01',
'WASP-70 A b']


ToCircle = ['TOI-1201 b', 'TOI-1728 b', 'TOI-1259 A b', 'TOI-1601 b',\
'TOI-1333 b', 'TOI-954 b', 'HD 202772 A b', 'TOI-4329b', 'TOI-1431 b', 'TOI-1518b',\
'TOI-2109b_Mass1','TOI 1641.01_Mass1']

PlotLabel1 = True
PlotLabel2 = True 

df_prop = pd.read_csv('Gemini LLP Target List - Target List_UpdatedWithPreviousSheets.csv')
df2 = pd.read_csv('Gemini LLP Target List - Target List_UpdatedWithPreviousSheets.csv')




pl_mass_prop = df_prop['Planet Mass (MJup'].values
pl_rad_prop = df_prop['Planet Radius (RJup'].values
pl_teq_prop = df_prop['Planet Temp (K'].values
n_prop = df_prop['pl_name']

pl_trandur = df_prop['Planet Duration (hours']

names = list(df_prop['pl_name'])




pl_mass_prop2 = df2['Planet Mass (MJup'].values
pl_rad_prop2 = df2['Planet Radius (RJup'].values
pl_teq_prop2 = df2['Planet Temp (K'].values
n_prop2 = df2['pl_name']

names2 = df2['pl_name']


strtype = type('')

ObTimeRunningSum = 0
for i in range(len(n_prop)):
    if type(n_prop[i]) == strtype:
        ObTimeRunningSum += pl_trandur[i]*24+(90/60)
        


##################

plt.figure(figsize=(10,5))
#plt.scatter(pl_rad, pl_mass, c=pl_teq,vmin=np.min(pl_teq), vmax=3000, cmap=cm.rainbow,s=50)
plt.scatter(pl_mass_prop, pl_teq_prop, c=pl_rad_prop,s=30)
plt.colorbar(label=r'Planet radius (R$_{Jup}$)',pad=0.01)



#ListOfAllPlanetsToBeNamed = ToCircle+HasBeenObserved
ListOfAllPlanetsToBeNamed = []


namelist3 = []
nv3 = n_prop.values
        

for i in nv3:    
    if i in ListOfAllPlanetsToBeNamed:   
        namelist3.append(i)
        
    else:
        namelist3.append('')




ax = plt.gca()

ax.set_xscale('log')

ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=10))
ax.tick_params(which='major', length=10)
ax.tick_params(which='minor', length=5)



plt.ylabel(r'Planetary temperature (K)')
plt.xlabel(r'Planet mass (M$_{Jup}$)')

# ### To do the labels 
plt.rc('font', family='serif', size=11)


plt.tight_layout()

axes = plt.gca()

ylims = axes.get_ylim()

Mjup = 1
Msat = 0.2994
Mnep = 0.05397 
Mearth = 0.003146

plt.plot([Mjup,Mjup],[ylims[0],2800],color='black',linestyle='--')
plt.plot([Mjup,Mjup],[3200,ylims[1]],color='black',linestyle='--')

plt.plot([Msat,Msat],[ylims[0],2800],color='black',linestyle='--')
plt.plot([Msat,Msat],[3200,ylims[1]],color='black',linestyle='--')

plt.plot([Mnep,Mnep],[ylims[0],2800],color='black',linestyle='--')
plt.plot([Mnep,Mnep],[3200,ylims[1]],color='black',linestyle='--')




plt.ylim(ylims)

plt.annotate('Jupiter',(1,3000),color='black',fontsize=15,horizontalalignment='center',verticalalignment='center')
plt.annotate('Saturn',(Msat,3000),color='black',fontsize=15,horizontalalignment='center',verticalalignment='center')
plt.annotate('Neptune',(Mnep,3000),color='black',fontsize=15,horizontalalignment='center',verticalalignment='center')


CircledPlanetCount = 0
for i in range(len(names)):
    if names[i] in ToCircle:
        CircledPlanetCount += 1
        if PlotLabel1:
            plt.scatter(pl_mass_prop[i], pl_teq_prop[i],edgecolor='grey',facecolors='none',s=300,label='Already observed')
            # PlotLabel1 = False 
        if not PlotLabel1:
            plt.scatter(pl_mass_prop[i], pl_teq_prop[i],edgecolor='grey',facecolors='none',s=300)



for i in range(len(names)):
    if names[i] in HasBeenObserved:
        plt.scatter(pl_mass_prop[i], pl_teq_prop[i],edgecolor='grey',facecolors='none',marker='D',s=300,label='Already observed')

plt.scatter(1.5, 700,edgecolor='grey',facecolors='none',s=200,marker='D')
plt.annotate('Already observed',(1.7,650),color='black',fontsize=10,horizontalalignment='left',verticalalignment='center')

plt.scatter(1.5, 1000,edgecolor='grey',facecolors='none',s=200)
plt.annotate('New TESS planet',(1.7,1000),color='black',fontsize=10,horizontalalignment='left',verticalalignment='center')

plt.savefig('GeminiLLP_15_oct_2021.pdf')



#########################################################
########################################################
#### Attempt at using broken axis 

# If we were to simply plot pts, we'd lose most of the interesting
# details due to the outliers. So let's 'break' or 'cut-out' the y-axis
# into two portions - use the top (ax) for the outliers, and the bottom
# (ax2) for the details of the majority of our data

plt.figure(figsize=(20,10))
BottomRatioValue = 5
f, (ax, ax2) = plt.subplots(2, 1, sharex=True,gridspec_kw={'height_ratios': [1, BottomRatioValue]})


# plot the same data on both axes
ts = ax.scatter(pl_mass_prop, pl_teq_prop, c=pl_rad_prop,s=10)



bs = ax2.scatter(pl_mass_prop, pl_teq_prop, c=pl_rad_prop,s=10)

# zoom-in / limit the view to different portions of the data

#ax1ylims = (4000,4100)
ax1ylims = (3900,4600)
ax2ylims = (400, 2900)

ax.set_ylim(ax1ylims)  # outliers only
ax2.set_ylim(ax2ylims)  # most of the data

# hide the spines between ax and ax2
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()

# This looks pretty good, and was fairly painless, but you can get that
# cut-out diagonal lines look with just a bit more work. The important
# thing to know here is that in axes coordinates, which are always
# between 0-1, spine endpoints are at these locations (0,0), (0,1),
# (1,0), and (1,1).  Thus, we just need to put the diagonals in the
# appropriate corners of each of our axes, and so long as we use the
# right transform and disable clipping.

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d*BottomRatioValue, +d*BottomRatioValue), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d*BottomRatioValue, +d*BottomRatioValue), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

ax.set_xscale('log')

# What's cool about this is that now if we vary the distance between
# ax and ax2 via f.subplots_adjust(hspace=...) or plt.subplot_tool(),
# the diagonal lines will move accordingly, and stay right at the tips
# of the spines they are 'breaking'

cax = f.add_axes([0.85, 0.12, 0.03, 0.76])
f.colorbar(bs, cax=cax, orientation='vertical',label=r'Planet radius (R$_{Jup}$)')

#ax2.set_ylabel(r'Planetary temperature (K)')
# ax2.text(0.015, 1200, r'Planetary temperature (K)',rotation=90)
# ax2.set_xlabel(r'Planet mass (M$_{Jup}$)')

CircledPlanetCount = 0 

for i in range(len(names)):
    if names[i] in ToCircle:
        CircledPlanetCount += 1
        if PlotLabel1:
            ax2.scatter(pl_mass_prop[i], pl_teq_prop[i],edgecolor='grey',facecolors='none',s=150,label='Already observed')
            PlotLabel1 = False 
        if not PlotLabel1:
            ax2.scatter(pl_mass_prop[i], pl_teq_prop[i],edgecolor='grey',facecolors='none',s=150)
            ax.scatter(pl_mass_prop[i], pl_teq_prop[i],edgecolor='grey',facecolors='none',s=150)

            

for i in range(len(names)):
    if names[i] in HasBeenObserved:
        if PlotLabel2:
            ax2.scatter(pl_mass_prop[i], pl_teq_prop[i],edgecolor='grey',facecolors='none',marker='D',s=150,label='Already observed')
            PlotLabel2 = False
        if not PlotLabel2:
            ax2.scatter(pl_mass_prop[i], pl_teq_prop[i],edgecolor='grey',facecolors='none',marker='D',s=150)
            


#### Ax 1 

ax.plot([Mjup,Mjup],[ax1ylims[0],ax1ylims[1]],color='black',linestyle='--',zorder=0)
ax2.plot([Mjup,Mjup],[ax2ylims[0],ax2ylims[1]],color='black',linestyle='--',zorder=0)


ax.plot([Msat,Msat],[ax1ylims[0],ax1ylims[1]],color='black',linestyle='--',zorder=0)
ax2.plot([Msat,Msat],[ax2ylims[0],ax2ylims[1]],color='black',linestyle='--',zorder=0)

ax.plot([Mnep,Mnep],[ax1ylims[0],ax1ylims[1]],color='black',linestyle='--',zorder=0)
ax2.plot([Mnep,Mnep],[ax2ylims[0],ax2ylims[1]],color='black',linestyle='--',zorder=0)




f.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.ylabel(r'Planetary temperature (K)',labelpad=15)
plt.xlabel(r'Planet mass (M$_{Jup}$)',labelpad=5)


plt.annotate('Jupiter',(0.705,0.805),color='black',fontsize=8,horizontalalignment='center',verticalalignment='center')
plt.annotate('Saturn',(0.46,0.805),color='black',fontsize=8,horizontalalignment='center',verticalalignment='center')
plt.annotate('Neptune',(0.13,0.805),color='black',fontsize=8,horizontalalignment='center',verticalalignment='center')


ax2.scatter(1.25, 550,edgecolor='grey',facecolors='none',s=100,marker='D')
ax2.annotate('Already observed',(1.50,550),color='black',fontsize=6,horizontalalignment='left',verticalalignment='center')



plt.gcf().subplots_adjust(right=0.8)

