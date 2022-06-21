# monte-carlo-annealing
**Monte-Carlo search for the minimum of the multidimensional "cost" function
Miranda and Mateus (2022), Water vapor tomography with optimized GNSS mapping functions,
Geophysical Research Letters, submitted.**

This algorithm starts from a first guess of the set of free parameters (input V),
and computes the cost function for that solution. Then it tries maxPERT alternative 
solutions contained in the multidimensional domain defined by [vmin, vmax], accepting, 
in turn, solutions that lead to a smaller cost function. Once this iteration 
is concluded the method proceeds by reducing the range in the search domain to a region 
around the current best solution, using a Boltzmann formula taken from the
concept of "simulated annealing" (e.g., Press et al, 1984, "Numerical Recipes").


