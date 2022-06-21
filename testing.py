# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 14:33:33 2022

@author: Pedro Miranda
"""

from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
from annealD import annealD
global xE,yE,zE     #location of measurements
global gE           #measurements
global G,roV,roB,VS #constants

def plotcost(cost,vmin,vmax,V,vnames,fn,tit=''):
    '''
    Plot cross-sections of the cost function in the vicinity of the "best" solution
    '''
    sh=np.shape(vmin)
    ndim=sh[0]
    nd=range(ndim)
    plans=combinations(nd,2)

    k=0        
    plt.figure(figsize=(12,8))
    kP=0

    nlines=3
    for plan in plans:
        print(plan)
        k=k+1
        xP=plan[0];yP=plan[1];
        nx=201;ny=201;
        xx=np.zeros((nx,ny));yy=np.zeros((nx,ny));
        CCC=np.zeros((nx,ny))
        xmin=vmin[xP];xmax=vmax[xP];
        ymin=vmin[yP];ymax=vmax[yP];
        for ix in range(nx):
            for iy in range(ny):
                xx[ix,iy]=xmin+ix*(xmax-xmin)/nx
                yy[ix,iy]=ymin+iy*(ymax-ymin)/ny
                Vp=np.copy(V)
                Vp[xP]=xx[ix,iy];
                Vp[yP]=yy[ix,iy];
                CCC[ix,iy]=cost(Vp);
        kP=kP+1
        print(kP)
        plt.subplot(nlines,3,kP)
        map=plt.contourf(xx,yy,np.log(CCC)/np.log(10),cmap='jet')
        
        plt.colorbar(map,label=r'$log_10(J)$')
        plt.scatter(V[xP],V[yP],marker='x',color='white')
        plt.xlabel(vnames[xP])
        plt.ylabel(vnames[yP]);

    plt.suptitle(tit)    
    plt.tight_layout()
    return 

def inipol(iseed=10,nE=1000): 
    '''
    Synthetic data (problem dependent)
    Fourth order polynomial 
    '''
    a=5;b=2;c=3;d=2 #solution
    xmin=0;xmax=10
    np.random.seed(iseed)
    xE=xmin+(xmax-xmin)*np.random.sample(nE)
    yO=a*xE**3+b*xE**2+c*xE+d
    return xE,yO

def cost(V): 
    '''
    Cost function (depends on the problem)
    '''
    a,b,c,d=V
    custo=np.sum(((a*xE**3+b*xE**2+c*xE+d)-yO)**2)
    return custo

# For a new problem change inipol and cost, CHECK parameters
xE,yO=inipol()                     # Define synthetic data OR read observations (yO) at locations (xE)
V=np.zeros((4))                    # First guess of the multidimensional solution
vnames=np.array(['a','b','c','d']) # Names of parameters (for output)
vmin=np.zeros(V.shape);vmax=10*np.ones(V.shape) # Domain of search

Jmin=-10.                   # Stop when Jmin (if Jmin>0)
minvstep=(vmax-vmin)/100000 # Cool until range<minvstep

kappa=np.float64(0.1) # Parameters for annealling 
T=10.                 # Parameters for annealling 
outITER=1             # Parameters for annealling (output)

maxITER=1000    # Number of cooling cycles
maxPERT=100000  # Number of trials per cycle
COOL=0.9        # Cooling rate

# Compute best solution
n,V,iTER,path,J=annealD(V,cost,vmin,vmax,Jmin,minvstep,maxITER,maxPERT)
# Solution
print(V)  
# Plot cross-sections of the cost funtion in the vicinity of the solution
plotcost(cost,vmin,vmax,V,vnames,10,tit=r'$ITER=%3i(%4i),PERT=%5i,min\Delta x=%6.5f,a=%4.3f,b=%4.3f,c=%4.3f,d=%4.3f$'\
         % (iTER,maxITER,maxPERT,minvstep[0],V[0],V[1],V[2],V[3]))
