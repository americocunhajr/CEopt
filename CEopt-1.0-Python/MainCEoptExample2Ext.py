# -----------------------------------------------------------------
#  MainCEoptExample2Ext.py
# -----------------------------------------------------------------
#  programmer: Americo Cunha
#              americo.cunhajr@gmail.com
# 
#  Originally programmed in: Mar 27, 2025
#           Last updated in: Mar 28, 2025
# -----------------------------------------------------------------
#  Example 2: The Peaks function in 2D
# -----------------------------------------------------------------

import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
from CEopt import CEopt  # Ensure CEopt.py is in the same directory

# --- Set random seed for reproducibility ---
np.random.seed(30081984)

# Close all figures
plt.close('all')

print(' ------------------- ')
print(' MainCEoptExample2Ext.py ')
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

# Fix the random seed for reproducibility (similar to RandStream in Matlab)
np.random.seed(30081984)

# Define objective function handle
F = PeaksFunc

# Bound for design variables (2D)
lb = np.array([-3, -3])
ub = np.array([ 3,  3])

# Initialize mean and std. dev. vectors
# mu0: randomly initialized within the bounds (as a 2-element array)
mu0    = lb + (ub - lb) * np.random.rand(2)
sigma0 = 5 * (ub - lb)

# Define parameters for the CE optimizer in a dictionary
CEstr = {}
CEstr['isVectorized'] = True        # Vectorized function
CEstr['EliteFactor' ] = 0.1         # Elite samples percentage
CEstr['Nsamp'       ] = 50          # Number of samples
CEstr['MaxIter'     ] = 80          # Maximum number of iterations
CEstr['TolAbs'      ] = 1.0e-2      # Absolute tolerance
CEstr['TolRel'      ] = 1.0e-2      # Relative tolerance
CEstr['alpha'       ] = 0.7         # Smoothing parameter for mean update
CEstr['beta'        ] = 0.8         # Smoothing parameter for std. dev. update
CEstr['q'           ] = 10          # Exponent for dynamic smoothing parameter

# Run the CE optimizer and time the execution
start_time = time.time()
Xopt, Fopt, ExitFlag, CEstr = CEopt(F , mu0, sigma0, lb, ub, None, CEstr)
elapsed_time = time.time() - start_time

print("\nCE optimizer results:")
print("Xopt =", Xopt)
print("Fopt =", Fopt)
print("ExitFlag =", ExitFlag)
print("Elapsed time: {:.4f} seconds".format(elapsed_time))

# Meshgrid for visualization
Ngrid  = 100
xRange = np.linspace(lb[0], ub[0], Ngrid)
yRange = np.linspace(lb[1], ub[1], Ngrid)
X, Y   = np.meshgrid(xRange, yRange)
pts    = np.column_stack((X.ravel(), Y.ravel()))
Z      = F(pts)
Z      = Z.reshape(X.shape)

# Custom colors
MyYellow = [0.9290, 0.6940, 0.1250]
MyBlue   = [0.0000, 0.4470, 0.7410]
MyRed    = [0.6350, 0.0780, 0.1840]

# Custom colormap (using 'viridis' with 50 discrete colors)
Ncolors = 50
MyColorMap1 = colormaps['viridis'].resampled(Ncolors)

# Set up the figure for visualization
plt.figure()
plt.clf()

# Plot F(x) as a filled contour plot
cont = plt.contourf(X, Y, Z, Ncolors, cmap=MyColorMap1)
plt.colorbar()

# Plot the analytical solution
plt.plot(0.22, -1.62, 'o', color=MyYellow, markersize=15,
         markerfacecolor=MyYellow, label='Analytical Solution')

# Plot the initial approximation
plt.plot(mu0[0], mu0[1], 'd', color=MyRed, linewidth=1,
         markersize=7, markerfacecolor=MyRed, label='Initial Approximation')

# Plot CE iterations
for j in range(CEstr['iter']):
    if j == 0:
        plt.plot(CEstr['xmean'][j, 0], CEstr['xmean'][j, 1], 'd',
                 color=MyRed, linewidth=1, markersize=7, label='CE Iterations')
    else:
        plt.plot(CEstr['xmean'][j, 0], CEstr['xmean'][j, 1], 'd',
                 color=MyRed, linewidth=1, markersize=7, label='_nolegend_')

# Plot the CE solution
plt.plot(Xopt[0, 0], Xopt[0, 1], 'x', color=MyBlue, linewidth=3,
         markersize=25, label='CE Optimum')

# Labeling
plt.xlabel('x1', fontsize=20)
plt.ylabel('x2', fontsize=20)

# Legend (placed in the upper right, equivalent to Matlab's 'NorthEast')
plt.legend(loc='upper right', fontsize=12)

# Set additional font and box properties
ax = plt.gca()
for label in ax.get_xticklabels():
    label.set_fontname('Arial')
for label in ax.get_yticklabels():
    label.set_fontname('Arial')
ax.tick_params(labelsize=18)
# Box on
ax.set_frame_on(True)

# Save the figure as PNG
plt.savefig('CEoptExample2.png', format='png', dpi=300)

# Show the figure
plt.show()
