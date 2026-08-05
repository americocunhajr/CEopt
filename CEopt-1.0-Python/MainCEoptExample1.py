import time
import numpy as np
from CEopt import CEopt  # Assumes CEopt.py is in the same directory

print(' ------------------- ')
print(' MainCEoptExample1.py ')
print(' ------------------- ')

# Objective function
# -----------------------------------------------------------------
def F(x):
    # x can be a scalar or an array; operations are elementwise.
    return -0.8 * np.exp(-(x - 2)**2) - 0.5 * np.exp(-(x + 2)**2) + 1
# -----------------------------------------------------------------

# Bound for design variables
lb = -5
ub =  5

# Initialize mean and standard deviation vectors.
mu0    = (ub + lb) / 2.0
sigma0 = (ub - lb) / 6.0

# Run the CE optimizer
start_time = time.time()
Xopt, Fopt, ExitFlag, CEobj = CEopt(F, mu0, sigma0, lb, ub)
elapsed_time = time.time() - start_time

print("\nXopt =", Xopt)
print("Fopt =", Fopt)
print("ExitFlag =", ExitFlag)
#print("CEobj =", CEobj)
print("Elapsed time: {:.4f} seconds".format(elapsed_time))
