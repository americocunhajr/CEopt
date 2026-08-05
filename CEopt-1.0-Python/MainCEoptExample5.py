import time
import numpy as np
from CEopt import CEopt  # Ensure CEopt.py is in the same directory

print(' ------------------- ')
print(' MainCEoptExample5.py ')
print(' ------------------- ')

# Objective function
# -----------------------------------------------------------------
def MyObjFunc(x, p):
    xE     = p[0, :]
    yE     = p[1, :]
    thetaA = p[2, :]
    a0     = 1.0
    b0     = x[0]
    c0     = x[0]
    d0     = x[1]
    gamma  = x[2]
    thetaD = x[3]
    theta0 = x[4]
    thetaE = (180 - gamma) / 2.0

    # Use degree-based cosine and sine:
    e0 = 2 * b0 * np.cos(np.deg2rad(thetaE))
    m  = np.sqrt(a0**2 + d0**2 - 2 * a0 * d0 * np.cos(np.deg2rad(thetaA + theta0)))

    # Compute beta in degrees using arcsin (convert from radians)
    beta   = np.rad2deg(np.arcsin((a0 * np.sin(np.deg2rad(thetaA + theta0))) / m))
    thetaB = np.rad2deg(np.arccos((b0**2 + m**2 - c0**2) / (2 * b0 * m))) - beta
    xEc0   = a0 * np.cos(np.deg2rad(thetaA + theta0 + thetaD)) + e0 * np.cos(np.deg2rad(thetaB + thetaE + thetaD))
    yEc0   = a0 * np.sin(np.deg2rad(thetaA + theta0 + thetaD)) + e0 * np.sin(np.deg2rad(thetaB + thetaE + thetaD))
    
    # Normalize data: note that Matlab normalizes yE using the range of xE.
    xE1  = (xE - np.mean(xE)) / (np.max(xE) - np.min(xE))
    yE1  = (yE - np.mean(yE)) / (np.max(xE) - np.min(xE))
    xEc1 = (xEc0 - np.mean(xEc0)) / (np.max(xEc0) - np.min(xEc0))
    yEc1 = (yEc0 - np.mean(yEc0)) / (np.max(xEc0) - np.min(xEc0))
    
    F1 = np.sum((xE1 - np.real(xEc1))**2 + (yE1 - np.real(yEc1))**2)
    F2 = np.sum(np.abs(np.imag(xEc1))**2 + np.abs(np.imag(yEc1))**2 + np.abs(np.imag(thetaB))**2)
    
    return F1 + F2
# -----------------------------------------------------------------

# Read data from CSV file
Xdata = np.loadtxt("Example5data.csv", delimiter=",")
xE    = Xdata[:, 0]
yE    = Xdata[:, 1]

# Generate angle vector (thetaA) from 0 to 360 with length equal to length(xE)
thetaA = np.linspace(0, 360, len(xE))

# Build parameter matrix p (with three rows: xE, yE, thetaA)
p = np.vstack((xE, yE, thetaA))

# Objective function
F = lambda x: MyObjFunc(x, p)

# Bound for design variables
lb = np.array([1.3, 2.0, 120.0, -15.0, -50.0])
ub = np.array([2.5, 3.5, 150.0, -5.0, -35.0])

# Run the CE optimizer
start_time = time.time()
Xopt, Fopt, ExitFlag, CEstr = CEopt(F, None, None, lb, ub, None)
elapsedTime = time.time() - start_time

print("\nCE optimizer results:")
print("Xopt =", Xopt)
print("Fopt =", Fopt)
print("ExitFlag =", ExitFlag)
print("Elapsed time: {:.4f} seconds".format(elapsedTime))
