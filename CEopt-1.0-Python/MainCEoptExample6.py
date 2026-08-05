import time
import numpy as np
from CEopt import CEopt  # Ensure CEopt.py is in the same directory

print(' ------------------- ')
print(' MainCEoptExample6.py ')
print(' ------------------- ')

# Objetive function
# -----------------------------------------------------------------
def PatternSearchFunc(x):
    x        = np.atleast_2d(x)
    x1       = x[:, 0]
    x2       = x[:, 1]
    nSamples = x.shape[0]
    F        = np.zeros((nSamples, 1))
    for i in range(nSamples):
        if x1[i] < -5:
            F[i] = (x1[i] + 5)**2 + np.abs(x2[i])
        elif x1[i] < -3:
            F[i] = -2 * np.sin(x1[i]) + np.abs(x2[i])
        elif x1[i] < 0:
            F[i] = 0.5 * x1[i] + 2 + np.abs(x2[i])
        elif x1[i] >= 0:
            F[i] = 0.3 * np.sqrt(x1[i]) + 5/2 + np.abs(x2[i])

    return F
# -----------------------------------------------------------------

# Constraint function
# -----------------------------------------------------------------
def ConicConstraints(x):
    x  = np.atleast_2d(x)
    x1 = x[:, 0]
    x2 = x[:, 1]
    nSamples = x.shape[0]
    G        = np.zeros((nSamples, 1))
    H        = np.zeros((nSamples, 1))
    for i in range(nSamples):
        G[i]  = 2 * x1[i]**2 + x2[i]**2 - 3
        H[i]  = (x1[i] + 1)**2 - (x2[i] / 2)**4

    return G, H
# -----------------------------------------------------------------

# Objective function and constraints
F       = PatternSearchFunc
nonlcon = ConicConstraints

# Bound for design variables and initial mean
lb  = np.array([-6, -4])
ub  = np.array([ 2,  4])
mu0 = np.array([-4, 2])

# Cross-entropy optimizer struct
CEstr = {}
CEstr['isVectorized'] = True   # vectorized function
CEstr['TolCon'] = 1.0e-6       # relative tolerance for constraints

# Run the CE optimizer
start_time = time.time()
Xopt, Fopt, ExitFlag, CEstr = CEopt(F, None, None, lb, ub, nonlcon, CEstr)
elapsedTime = time.time() - start_time

print("\nCE optimizer results:")
print("Xopt =", Xopt)
print("Fopt =", Fopt)
print("ExitFlag =", ExitFlag)
print("Elapsed time: {:.4f} seconds".format(elapsedTime))
