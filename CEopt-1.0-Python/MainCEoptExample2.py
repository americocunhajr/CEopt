import time
import numpy as np
from CEopt import CEopt  # Assumes CEopt.py is in the same directory

print(' ------------------- ')
print(' MainCEoptExample2.py ')
print(' ------------------- ')

# -----------------------------------------------------------------
# Objective function
# -----------------------------------------------------------------
def PeaksFunc(x):
    x = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]

    nSamples = x.shape[0]
    F        = np.zeros((nSamples, 1))

    for i in range(nSamples):
        F[i] = (3*(1 - x1[i])**2 * np.exp(-x1[i]**2 - (x2[i] + 1)**2) -
             10*(x1[i]/5 - x1[i]**3 - x2[i]**5) * np.exp(-x1[i]**2 - x2[i]**2) -
             (1/3)*np.exp(-(x1[i]+1)**2 - x2[i]**2))

    return F
# -----------------------------------------------------------------

# Define the objective function handle.
F = PeaksFunc
    
# Bound for design variables (2D problem)
lb = np.array([-3, -3])
ub = np.array([ 3,  3])
    
# Initialize mean and standard deviation
mu0    = lb + (ub - lb) * np.random.rand(2)
sigma0 = 5 * (ub - lb)
    
# Define parameters for the CE optimizer.
CEstr = {}
CEstr['isVectorized'] = True       # Vectorized function (True)
CEstr['EliteFactor']  = 0.1         # Elite samples percentage
CEstr['Nsamp']        = 50          # Number of samples
CEstr['MaxIter']      = 80          # Maximum number of iterations
CEstr['TolAbs']       = 1.0e-2      # Absolute tolerance
CEstr['TolRel']       = 1.0e-2      # Relative tolerance
CEstr['alpha']        = 0.7         # Smoothing parameter for the mean update
CEstr['beta']         = 0.8         # Smoothing parameter for the std. dev. update
CEstr['q']            = 10          # Exponent for dynamic smoothing parameter
    
# Run the CE optimizer
start_time = time.time()
Xopt, Fopt, ExitFlag, CEstr = CEopt(F, mu0, sigma0, lb, ub, None, CEstr)
elapsed_time = time.time() - start_time
    
print("\nCE optimizer results:")
print("Xopt =", Xopt)
print("Fopt =", Fopt)
print("ExitFlag =", ExitFlag)
print("Elapsed time: {:.4f} seconds".format(elapsed_time))
