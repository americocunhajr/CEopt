# -----------------------------------------------------------------
#  MainCEoptExample5Ext.py
# -----------------------------------------------------------------
#  programmer: Jose Geraldo Telles Ribeiro
#              jose.gt.ribeiro@gmail.com
# 
#  Originally programmed in: Mar 28, 2025
#           Last updated in: Mar 28, 2025
# -----------------------------------------------------------------
#  Example 5: Design of a complex mechanism
# -----------------------------------------------------------------

import time
import numpy as np
import matplotlib.pyplot as plt
from CEopt import CEopt  # Ensure CEopt.py is in the same directory

# --- Set random seed for reproducibility ---
np.random.seed(30081984)

# Close all figures
plt.close('all')

print(' ------------------- ')
print(' MainCEoptExample5Ext.py ')
print(' ------------------- ')

# -----------------------------------------------------------------
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

# Create te angle vector (thetaA) from 0 to 360 with length equal to length(xE)
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

# -----------------------------------------------------------------
# Plot optimal mechanism
# -----------------------------------------------------------------
a0     = 1.0
b0     = Xopt[0,0] * a0
c0     = b0
d0     = Xopt[0,1] * a0
gamma  = Xopt[0,2]
thetaD = Xopt[0,3]
theta0 = Xopt[0,4]

# Create a fine vector for thetaA_1 from 0 to 360 in steps of 0.5 degrees
thetaA_1 = np.arange(0, 360.5, 0.5)

# Compute thetaE in degrees
thetaE = (180 - gamma) / 2.0

# cosd and sind: use np.cos(np.deg2rad(...)) and np.sin(np.deg2rad(...))
e0 = 2 * b0 * np.cos(np.deg2rad(thetaE))
m  = np.sqrt(a0**2 + d0**2 - 2*a0*d0*np.cos(np.deg2rad(thetaA_1 + theta0)))

# Compute beta (in degrees) using arcsin (converted from radians)
beta   = np.rad2deg(np.arcsin((a0 * np.sin(np.deg2rad(thetaA_1 + theta0))) / m))
thetaB = np.rad2deg(np.arccos((b0**2 + m**2 - c0**2) / (2 * b0 * m))) - beta

# Compute xEc0 and yEc0 using degree-based trigonometric functions
xEc0 = a0 * np.cos(np.deg2rad(thetaA_1 + theta0 + thetaD)) + e0 * np.cos(np.deg2rad(thetaB + thetaE + thetaD))

# Scale factors: normalize the mechanism dimensions based on the data
a = a0 * (np.max(xE) - np.min(xE)) / (np.max(xEc0) - np.min(xEc0))
b = b0 * (np.max(xE) - np.min(xE)) / (np.max(xEc0) - np.min(xEc0))
c = c0 * (np.max(xE) - np.min(xE)) / (np.max(xEc0) - np.min(xEc0))
d = d0 * (np.max(xE) - np.min(xE)) / (np.max(xEc0) - np.min(xEc0))
e = e0 * (np.max(xE) - np.min(xE)) / (np.max(xEc0) - np.min(xEc0))

# Recompute the designed trajectory
xEc = a * np.cos(np.deg2rad(thetaA_1 + theta0 + thetaD)) + e * np.cos(np.deg2rad(thetaB + thetaE + thetaD))
yEc = a * np.sin(np.deg2rad(thetaA_1 + theta0 + thetaD)) + e * np.sin(np.deg2rad(thetaB + thetaE + thetaD))

# Translate the mechanism to match the data's mean values
x_A = np.mean(xE) - np.mean(xEc)
y_A = np.mean(yE) - np.mean(yEc)
xEc = xEc + x_A
yEc = yEc + y_A

# Set up the figure for mechanism plot
plt.figure(figsize=(8,6))
plt.clf()
plt.plot(xE, yE, 'bo', linewidth=2.0, label='Target Trajectory')
plt.plot(xEc, yEc, 'r-', linewidth=2.0, label='Designed Trajectory')
plt.xlim([0, 50])
plt.ylim([0, 40])
plt.xlabel('x', fontsize=20, fontname='Arial')
plt.ylabel('y', fontsize=20, fontname='Arial')
plt.legend(loc='best', fontsize=12)
plt.box(True)
plt.tight_layout()
plt.savefig('CEoptExample5.png', format='png', dpi=300)
plt.show()
# -----------------------------------------------------------------
