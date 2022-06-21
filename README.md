# monte-carlo-annealing
**Monte-Carlo search for the minimum of the multidimensional "cost" function
Miranda and Mateus (2022), Water vapor tomography with optimized GNSS mapping functions,
Geophysical Research Letters, submitted.**

<img src="https://github.com/pjmateus/monte-carlo-annealing/blob/80329255f53d957f4e4ecd28d85fe946c23f048b/logos.png" width="350">

The algorithm was developed and written at the Dom Luiz Institute (IDL), Faculty of Sciences of the University of Lisbon (FCUL), by Pedro Miranda (pmmiranda@fc.ul.pt).
The following scripts (annealD.py and testing.py) contains guidelines to run the code.

This algorithm starts from a first guess of the set of free parameters (input V), and computes the cost function for that solution. Then it tries maxPERT alternative 
solutions contained in the multidimensional domain defined by [vmin, vmax], accepting, in turn, solutions that lead to a smaller cost function. Once this iteration 
is concluded the method proceeds by reducing the range in the search domain to a region around the current best solution, using a Boltzmann formula taken from the
concept of "simulated annealing" (e.g., Press et al, 1984, "Numerical Recipes").

The testing.py script can be used to test the **annealD** algorithm. 

The **inipol** funtion is used here to generate synthetic data (can be replaced by observational data).
In this example we generate synthetic data for a fourth order polynomial equation.
```Python
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
```

The **cost** funtion defines the "cost funtion" (that depend on the problem) 
```Python
def cost(V): 
    '''
    Cost function (depends on the problem)
    '''
    a,b,c,d=V
    custo=np.sum(((a*xE**3+b*xE**2+c*xE+d)-yO)**2)
    return custo
```

```Python
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
# check the size of V (only for output)
plotcost(cost,vmin,vmax,V,vnames,10,tit=r'$ITER=%3i(%4i),PERT=%5i,min\Delta x=%6.5f,a=%4.3f,b=%4.3f,c=%4.3f,d=%4.3f$'\
         % (iTER,maxITER,maxPERT,minvstep[0],V[0],V[1],V[2],V[3]))
```

**If you have any questions do not hesitate to contact me by email pmmiranda@fc.ul.pt**
