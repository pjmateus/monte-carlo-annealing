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



