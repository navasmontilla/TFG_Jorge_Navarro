#!/usr/bin/env python
# coding: utf-8

# ## WFM results
# 

# Libraries are imported and program folders and files are defined:

# In[121]:


import math                    
import numpy as np             
import matplotlib.pyplot as plt 
from sklearn.linear_model import LinearRegression

folder_base="output-files/"
folder_sol=""
folder_out=folder_base

type_model=1 #1: (T,Y), 2: (H,Y)


# Data is read:

# In[122]:


from glob import glob
import os
import re

files = sorted(glob(folder_out + "/listFire0*.out"), key=lambda x: int(re.findall(r'\d+', os.path.basename(x))[0]))
#files = glob(folder_out+"/listFire0*")
lf=len(files)
print("The number of files is", lf)

#print(files)

file1 = open(files[0], 'r')
line = file1.readline()
line = file1.readline()
cells=re.findall(r'\b\d+\b',line)
xcells=int(cells[0])
ycells=int(cells[1])

print("Mesh dimensions are:",xcells,"x",ycells)
print("\nList of files read:\n")

xc=np.zeros((xcells,ycells))
yc=np.zeros((xcells,ycells))
T=np.zeros((xcells,ycells,lf))
Y=np.zeros((xcells,ycells,lf))
H=np.zeros((xcells,ycells,lf))
Rf=np.zeros((xcells,ycells,lf))
Z=np.zeros((xcells,ycells))


j=0
for fname in files:
    print(fname+" file read")
    data     = np.loadtxt(fname, skiprows=2)
    k=0
    for l in range(0,xcells):   
        for m in range(0,ycells):
            #k = l + m*xcells + n*xcells*ycells
            xc[l,m]=data[k,0]
            yc[l,m]=data[k,1]
            T[l,m,j]=data[k,2]
            Y[l,m,j]=data[k,3]
            Rf[l,m,j]=data[k,4]
            Z[l,m]=100.0*np.exp(-((xc[l,m]-200)**2+(yc[l,m]-200)**2)/20000) #topography
            if type_model==2:
                H[l,m,j]=data[k,5]
            k+=1       
    j=j+1   
    
times  = np.loadtxt(folder_out+"/times_fire.out")


# **2D PLOT:**
# 
# This creates a 2D plot, to select the time use `j=...`, the index for the file that takes values from `0` to `lf-1`:

# In[123]:


j=-1 #file number (time)

print("Plot at time:",times[j], "s")

#filename = 

xp = xc[:,0]     #puntos en x
yp = yc[0,:]      #puntos en y
Xp, Yp = np.meshgrid(xp, yp)    #matriz de puntos
ST=np.transpose(T[:,:,j])
SY=np.transpose(Y[:,:,j])
SRf=np.transpose(Rf[:,:,j])
ZT=np.transpose(Z[:,:])

fig, ax = plt.subplots(1,2,figsize=(9,4.5))      #genera el objeto "figura"

levels = np.linspace(np.min(T), np.max(T), 12)
#print(levels)
plot1=ax[0].contour(Xp, Yp, ST, 6,levels=levels,colors="k",linewidths=0.75)  
plot1=ax[0].contourf(Xp, Yp, ST, 200, cmap='YlOrRd')
ax[0].set_title('T(K)')
ax[0].set_xlabel("x(m)") 
ax[0].set_ylabel("y(m)") 
ax[0].set_aspect('equal', 'box')
#ax[0].plot([100,100],[1,199],'k--')
#ax[0].plot([1,199],[100,100],'k--')
#plot1.set_clim( np.min(T), np.max(T) )
plot1.set_clim( 300.0, 1400.0 )

plot2=ax[1].contourf(Xp, Yp, SY, 200, cmap='Greens' )
plot1=ax[1].contour(Xp, Yp, SY, 6,colors="k",linewidths=0.75) 
#ax[1].plot([100,100],[1,199],'k--')
#ax[1].plot([1,199],[100,100],'k--')


ax[1].set_title('Y')
ax[1].set_xlabel("x(m)") 
ax[1].set_ylabel("y(m)") 
ax[1].set_aspect('equal', 'box')

#plot2=ax[2].contourf(Xp, Yp, SRf, 200, cmap='viridis' )
#plot1=ax[2].contour(Xp, Yp, SRf, 12,colors="k",linewidths=0.75) 
#ax[2].set_title('Rf(x,y)')
#ax[2].set_xlabel("x") 
#ax[2].set_ylabel("y") 
#ax[2].set_aspect('equal', 'box')

dx=xc[1,0]-xc[0,0]
burned_area = dx*dx*np.sum(SY < 0.99)

print("Dx is: ",dx)
print("Burned area: ",burned_area)

print("Max T is: ",np.max(ST))
print("Min T is: ",np.min(ST))

print("Max Y is: ",np.max(SY))
print("Min Y is: ",np.min(SY))

print("Max Rf is: ",np.max(SRf))
print("Min Rf is: ",np.min(SRf))





SRf=np.transpose(Rf[:,:,0])
for j in range(1,lf):
    SY=np.transpose(Y[:,:,j])
    plot3=ax[1].contour(Xp, Yp, SY, [0.99],colors="r",linewidths=0.75) 

fig.tight_layout()

fig.savefig("fig_2D_TY.png",dpi=500)  
  
  
  
  
  
  
  
j=-1 #file number (time)

#filename = 

xp = xc[:,0]     #puntos en x
yp = yc[0,:]      #puntos en y
Xp, Yp = np.meshgrid(xp, yp)    #matriz de puntos
ST=np.transpose(T[:,:,j])
SY=np.transpose(Y[:,:,j])
SRf=np.transpose(Rf[:,:,j])
  

fig, ax = plt.subplots(1,2,figsize=(10.5, 4.5))      #genera el objeto "figura"


plot2=ax[0].contourf(Xp, Yp, SY, 200, cmap='Greens' )
plot3=ax[0].contour(Xp, Yp, SY, 6,colors="k",linewidths=0.75) 
plot2.set_clim( 0.0, 1.0 )

ax[0].set_title('Y')
ax[0].set_xlabel("x(m)") 
ax[0].set_ylabel("y(m)") 
ax[0].set_aspect('equal', 'box')


plot2=ax[1].contourf(Xp, Yp, SRf*400.0, 200, cmap='BrBG' )
plot3=ax[1].contour(Xp, Yp, SRf*400.0, 6,colors="k",linewidths=0.75) 
plot2.set_clim( 0.0, 4.00 )

ax[1].set_title('W')
ax[1].set_xlabel("x(m)") 
ax[1].set_ylabel("y(m)") 
ax[1].set_aspect('equal', 'box')

cbar=fig.colorbar(plot2, ax=ax[1])
cbar.set_label('W(kg/m$^2$)')


fig.text(0.15, 0.85, "SC=1.0", fontsize=11)
fig.text(0.35, 0.85, "SC=0.4", fontsize=11)


SRf=np.transpose(Rf[:,:,0])
for j in range(1,lf):
    SY=np.transpose(Y[:,:,j])
    plot4=ax[0].contour(Xp, Yp, SY, [0.99],colors="r",linewidths=0.75) 
    plot5=ax[1].contour(Xp, Yp, SY, [0.99],colors="r",linewidths=0.75) 

fig.tight_layout()

fig.savefig("fig_2D_YW.png",dpi=500)  
  
  
j=-1 #file number (time)

#filename = 

xp = xc[:,0]     #puntos en x
yp = yc[0,:]      #puntos en y
Xp, Yp = np.meshgrid(xp, yp)    #matriz de puntos
ST=np.transpose(T[:,:,j])
SY=np.transpose(Y[:,:,j])
SRf=np.transpose(Rf[:,:,j])
  

fig, ax = plt.subplots(1,2,figsize=(10.5, 4.5))      #genera el objeto "figura"


plot1=ax[0].contour(Xp, Yp, ST, 6,levels=levels,colors="k",linewidths=0.75)  
plot1=ax[0].contourf(Xp, Yp, ST, 200, cmap='YlOrRd') 
#plot1=ax[0].contour(Xp, Yp, ZT, 12,colors="k",linestyles="--",linewidths=0.3)       
ax[0].set_title('T(K)')
ax[0].set_xlabel("x(m)") 
ax[0].set_ylabel("y(m)") 
ax[0].set_aspect('equal', 'box')
plot1.set_clim( 300.0, 1400.0 )

ax[0].set_title('T')
ax[0].set_xlabel("x(m)") 
ax[0].set_ylabel("y(m)") 
ax[0].set_aspect('equal', 'box')

fig.text(0.11, 0.88, "T_max="+str(round(np.max(ST),1))+" K", fontsize=9)
fig.text(0.11, 0.83, "T_min="+str(round(np.min(ST),1))+" K", fontsize=9)


plot2=ax[1].contourf(Xp, Yp, SRf*400.0, 200, cmap='BrBG' )
plot3=ax[1].contour(Xp, Yp, SRf*400.0, 6,colors="k",linewidths=0.75) 
plot2.set_clim( 0.0, 4.00 )

ax[1].set_title('W')
ax[1].set_xlabel("x(m)") 
ax[1].set_ylabel("y(m)") 
ax[1].set_aspect('equal', 'box')

cbar=fig.colorbar(plot2, ax=ax[1])
cbar.set_label('W(kg/m$^2$)')



SRf=np.transpose(Rf[:,:,0])
for j in range(1,lf):
    SY=np.transpose(Y[:,:,j])
    plot4=ax[0].contour(Xp, Yp, SY, [0.99],colors="r",linewidths=0.75) 
    plot5=ax[1].contour(Xp, Yp, SY, [0.99],colors="r",linewidths=0.75) 

fig.tight_layout()

fig.savefig("fig_2D_TW.png",dpi=500)  

    
  
  
  
  
  
