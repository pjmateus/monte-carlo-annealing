# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 14:40:49 2022

@author: Pedro Miranda
"""
import numpy as np
from time import process_time

def annealD(V,cost,vmin,vmax,Jmin,minvstep,maxITER,maxPERT,\
            COOL=0.9,kappa=0.1,T=10,outITER=1,iseed=None,FG=True):
    """
    Monte-Carlo search for the minimum of the multidimensional "cost" function
    Miranda and Mateus (2022), Water vapor tomography with optimized GNSS mapping functions
    Geophysical Research Letters, submitted.

    Input Variables:
    V        : 1D array with the values of the free parameters
    cost     : cost function; cost(V)
    vmin,vmax: 1D arrays defining the domain of the search for each free parameter
    Jmin     : if positive and J<Jmin, stop iteration
    minvstep : stop "cooling" when vstep[0]<minvstep
    maxITER  : maxmimum number of iterations (in the outer cycle) or cooling cycles
    maxPERT  : number of random trials in each iteration
    COOl     : cooling rate, controls the change in "Temperature" (T) between iterations
    T        : initial "Temperature"
    kappa    : parameter for the Boltzmann function
    outITER  : debug iteration interval between prints, 1: prints all iterations
    if FG=True generate random first guess, otherwise use input V 

    Output:
    ndim     : number of free parameters
    V        : "optimal" solution
    iTER     : number of iterations done
    path     : evolution of the solution
    J        : "optimal" (minimum) value of the cost function

    This algorithm starts from a first guess of the set of free parameters (input V),
    and computes the cost function for that solution. Then it tries maxPERT alternative 
    solutions contained in the multidimensional domain defined by [vmin,vmax], accepting, 
    in turn, solutions that lead to a smaller cost function. Once this iteration 
    is concluded the method proceeds by reducing the range in the search domain to a region 
    around the current best solution, using a Boltzmann formula taken from the
    concept of "simulated annealing" (e.g., Press et al, 1984, "Numerical Recipes")
    """    

    t0=process_time()
    path = []
    sh = np.shape(vmin)
    ndim = sh[0]
    
    fmt='%3i,%2i,['
    for idim in range(ndim):
        fmt+='%3.2e,'
    fmt+='],%4.3f,%5.0fs'
    
    vstep = 2*(vmax-vmin)
    minxstep = minvstep[0]
    xstep = vstep[0]
    np.random.seed(iseed)
    rrr = np.random.sample(sh)
    if not FG:
        print('random FG')
        V = vmin+(vmax-vmin)*rrr
    VI = np.zeros(sh)
    J = cost(V)
    print('start:',V,1000*np.sqrt(J),vmin,vmax)

    iTER = 0
    nHIT = 1

    print(maxITER, Jmin, minxstep,maxPERT,ndim)
    while (iTER < maxITER and J > Jmin and xstep > minxstep):
        nHIT = 0
        for iP in range(maxPERT):
            for idim in range(ndim):
                rr = np.random.rand()-0.5
                VI[idim] = V[idim]+vstep[idim]*rr
                while VI[idim] < vmin[idim] or VI[idim] > vmax[idim]:
                    rr = np.random.rand()-0.5
                    VI[idim] = V[idim]+vstep[idim]*rr
            JI = cost(VI)
            if JI < J:
                J = np.copy(JI)
                V = np.copy(VI)
                path.append(V)
                nHIT = nHIT+1

            iP = iP+1
        iTER = iTER+1
        T = T*COOL
        vstep = np.maximum(minvstep, vstep*np.exp(-kappa/T))
        xstep = vstep[0]
        t1=process_time()
        if iTER % outITER == 0:
            # change the number of V[] depending on dimension (only for output)
            print(fmt %\
                      (iTER,nHIT,V[0],V[1],V[2],V[3],1000*np.sqrt(J),t1-t0))
    return ndim, V, iTER, path, J
